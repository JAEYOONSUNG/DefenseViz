"""Find native DNMBsuite sample roots and read source-bound display metadata.

This lightweight inventory never imports calls, executes R or infers protein
function. An accepted root is a boundary: its modules and plasmids are not
additional samples. Explicit manifests remain the way to select repeated runs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

GENBANK_SUFFIXES = frozenset({".gbff", ".gbk", ".gb"})
ASSEMBLY = re.compile(r"(?<![A-Za-z0-9_])GC[AF]_\d+\.\d+(?![A-Za-z0-9_.])")


def _source_metadata(path):
    """Read GenBank DBLINK and source qualifiers, without loading sequences."""
    before = path.stat()
    digest = hashlib.sha256()
    records = []
    record = None
    field = None
    feature = None
    qualifier = None
    pieces = []

    def flush_qualifier():
        nonlocal qualifier, pieces
        if qualifier is not None:
            value = " ".join(" ".join(pieces).split()).strip('"')
            if value:
                record["source"].setdefault(qualifier, []).append(value)
        qualifier, pieces = None, []

    def flush_record():
        nonlocal record
        flush_qualifier()
        if record is not None:
            records.append(record)
        record = None

    with path.open("rb") as stream:
        for raw in stream:
            digest.update(raw)
            line = raw.decode("utf-8-sig").rstrip("\r\n")
            if line.startswith("LOCUS "):
                flush_record()
                record = dict(assemblies=[], organism=None, source={})
                field, feature = None, None
            if record is None:
                continue
            if line.startswith("//"):
                flush_record()
                continue
            key = line[:12].strip()
            if line.startswith("FEATURES "):
                field, feature = "FEATURES", None
                continue
            if line.startswith(("ORIGIN", "CONTIG")):
                flush_qualifier()
                field, feature = "SEQUENCE", None
                continue
            if field == "FEATURES":
                name = line[5:21].strip() if len(line) >= 21 else ""
                if name:
                    flush_qualifier()
                    feature = name
                if feature == "source":
                    content = line[21:].strip()
                    if content.startswith("/"):
                        flush_qualifier()
                        match = re.match(r'/([^=]+)=(.*)', content)
                        if match and match[1] in {"organism", "strain", "isolate"}:
                            qualifier, pieces = match[1], [match[2]]
                    elif qualifier is not None:
                        pieces.append(content)
                continue
            if field == "SEQUENCE":
                continue
            if key:
                field = key
            if field == "DBLINK":
                # A reference paper or annotation comment mentioning an
                # accession is not the assembly identity of this record.
                match = re.search(r"\bAssembly:\s*(.*)", line[12:])
                if match:
                    record["assemblies"].extend(ASSEMBLY.findall(match[1]))
            if key == "ORGANISM":
                record["organism"] = " ".join(line[12:].split()) or None
    flush_record()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"{path}: GenBank metadata source changed while reading")
    if not records:
        raise ValueError(f"{path}: no GenBank records for automatic sample metadata; use a valid source or explicit manifest")
    return dict(path=str(path), sha256=digest.hexdigest(), records=records)


def _display_metadata(sources, directory):
    records = [record for source in sources for record in source["records"]]
    assemblies = sorted({value for record in records for value in record["assemblies"]})
    if len(assemblies) > 1:
        raise ValueError(f"{directory}: conflicting GenBank assembly identities {assemblies}; use an explicit manifest to select a sample")
    labels = set()
    all_strains = list(dict.fromkeys(alias.strip() for record in records
                                    for value in record["source"].get("strain", [])
                                    for alias in value.split(";") if alias.strip()))
    for record in records:
        organisms = set(record["source"].get("organism", []))
        if len(organisms) > 1:
            raise ValueError(f"{directory}: conflicting source organism qualifiers")
        organism = next(iter(organisms), record["organism"])
        if not organism:
            continue
        strains = list(dict.fromkeys(alias.strip() for value in record["source"].get("strain", [])
                                     for alias in value.split(";") if alias.strip()))
        if not strains:
            # A plasmid/contig may omit the strain qualifier that is present
            # on the chromosome. Its organism still must agree after removal
            # of an exact, declared strain suffix.
            strains = [s for s in all_strains if organism.endswith(" " + s)]
        if not strains:
            strains = list(dict.fromkeys(record["source"].get("isolate", [])))
        # Prefer a declared alias already present at the end of the organism;
        # otherwise retain the first declared strain. Never parse the folder.
        matched = [s for s in strains if organism.endswith(" " + s)]
        strain = max(matched, key=len) if matched else (strains[0] if strains else "")
        taxon = organism[:-len(strain)].rstrip() if matched else organism
        if matched and taxon.endswith(" strain"):
            taxon = taxon[:-7].rstrip()
        labels.add((taxon, strain))
    # In one assembly, plasmids/contigs can declare only the organism. Fill a
    # missing strain only from one unambiguous explicit strain of that exact
    # taxon; retain all original source qualifiers in the inventory. Neither
    # a different taxon nor conflicting explicit strains can be reconciled.
    complete = {(taxon, strain) for taxon, strain in labels if strain}
    normalized = set()
    for taxon, strain in labels:
        matches = {pair for pair in complete if pair[0] == taxon}
        normalized.add(next(iter(matches)) if not strain and len(matches) == 1 else (taxon, strain))
    labels = normalized
    if len(labels) > 1:
        raise ValueError(f"{directory}: conflicting GenBank organism/strain metadata; use an explicit manifest")
    result = dict(sample_id=assemblies[0] if assemblies else directory.name,
                  path=str(directory), label=directory.name)
    if labels:
        taxon, strain = next(iter(labels))
        result.update(label=" ".join(x for x in (taxon, strain) if x),
                      display=dict(taxon_label=taxon, strain_label=strain))
        parts = taxon.split(maxsplit=1)
        if len(parts) == 2:
            result.update(genus=parts[0], species=parts[1])
    return result, ("genbank_dblink_assembly" if assemblies else "directory_name")


def _native_root(path):
    return any(p.name.startswith("dnmb_module_") and p.is_dir() for p in path.iterdir())


def _generated_root(path):
    """Names such as output/reviews alone never prove that a tree is disposable."""
    if not (path / "analysis.json").is_file() or not (path / "receipt.json").is_file():
        return False
    try:
        receipt = json.loads((path / "receipt.json").read_text(encoding="utf-8"))
    except (OSError, ValueError, UnicodeError):
        return False
    return isinstance(receipt, dict) and receipt.get("tool") == "DefenseViz"


def _symlink_boundary(path, root):
    target = path.resolve()
    if target == root or not target.is_relative_to(root):
        return None
    relative = target.relative_to(root)
    if any(part.startswith(".") for part in relative.parts):
        return "hidden_target"
    for parent in target.parents:
        if parent == root:
            break
        if _native_root(parent):
            return "native_root_descendant"
        if _generated_root(parent):
            return "generated_output_descendant"
    return None


def discover_suite(root):
    """Return deterministic sample specs and their exact metadata-source hashes.

    Symbolic links to the same directory are deduplicated, cycles are skipped,
    and different native directories yielding the same sample ID are rejected.
    Hidden trees and identified DefenseViz output roots are excluded. Ordinary
    folder names, including output/outputs/reviews, remain discoverable.
    Directory names remain valid identities when no assembly metadata exists.
    No cohort/target/transformability label is invented.
    """
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"{root}: suite root must be an existing directory")
    samples, skipped, seen, identities = [], [], set(), {}
    pending = [root]
    while pending:
        directory = pending.pop()
        stat = directory.stat()
        identity = (stat.st_dev, stat.st_ino)
        if identity in seen:
            skipped.append(dict(path=str(directory), reason="duplicate_or_cycle"))
            continue
        seen.add(identity)
        children = sorted(directory.iterdir(), key=lambda p: (p.name.casefold(), p.name))
        subdirs = [p for p in children if p.is_dir()]
        if any(p.name.startswith("dnmb_module_") for p in subdirs):
            directory = directory.resolve()
            genomes = sorted({p.resolve() for p in children if p.is_file() and p.suffix.lower() in GENBANK_SUFFIXES})
            sources = [_source_metadata(p) for p in genomes]
            spec, origin = _display_metadata(sources, directory)
            sid = spec["sample_id"]
            if sid in identities:
                raise ValueError(f"duplicate suite sample_id {sid!r}: {identities[sid]} and {directory}; select the intended run with an explicit manifest")
            identities[sid] = str(directory)
            samples.append(dict(spec=spec, directory=str(directory), identity_origin=origin,
                                metadata_sources=sources))
            continue
        if _generated_root(directory):
            skipped.append(dict(path=str(directory), reason="generated_output"))
            continue
        for child in reversed(subdirs):
            if child.name.startswith(".") or child.name.startswith("dnmb_module_"):
                skipped.append(dict(path=str(child), reason="excluded_directory"))
                continue
            boundary = _symlink_boundary(child, root) if child.is_symlink() else None
            if boundary:
                skipped.append(dict(path=str(child), reason=boundary))
                continue
            pending.append(child)
    if not samples:
        raise ValueError(f"{root}: no native suite folders discovered")
    return dict(schema_version=1, root=str(root), samples=samples, skipped=skipped)

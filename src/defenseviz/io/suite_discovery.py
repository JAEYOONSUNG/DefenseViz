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


def _declared_strains(record):
    return list(dict.fromkeys(alias.strip() for value in record["source"].get("strain", [])
                             for alias in re.split(r"[;=]", value) if alias.strip()))


def _organism_strain_suffix(organism, strains):
    """Separate only declared suffixes or anchored explicit alias chains."""
    parts = re.split(r"\s*[=;]\s*", organism)
    collection = r"[A-Z]{2,}[A-Z0-9_-]*\s*\d[A-Za-z0-9._/-]*"
    if len(parts) > 1:
        # '=' is source-declared equivalence, but an unknown free-text prefix
        # must not become a guessed strain or a truncated taxon. The last
        # alias must exactly match a source strain qualifier.
        infer_collection = "=" in organism and ";" not in organism
        if parts[-1] not in strains or any(
                part not in strains and not (infer_collection and re.fullmatch(collection, part))
                for part in parts[1:]):
            return organism, []
    else:
        infer_collection = False
    matched = [strain for strain in strains if parts[0].endswith(" " + strain)]
    first = max(matched, key=len) if matched else None
    if first is None and infer_collection:
        match = re.search(r"(?:^|\s)(" + collection + r")$", parts[0])
        first = match[1] if match else None
    if first is None:
        return organism, []
    taxon = parts[0][:-len(first)].rstrip()
    if taxon.endswith(" strain"):
        taxon = taxon[:-7].rstrip()
    return (taxon, [first, *parts[1:]]) if taxon else (organism, [])


def _display_metadata(sources, directory):
    records = [record for source in sources for record in source["records"]]
    assemblies = sorted({value for record in records for value in record["assemblies"]})
    if len(assemblies) > 1:
        raise ValueError(f"{directory}: conflicting GenBank assembly identities {assemblies}; use an explicit manifest to select a sample")
    labels = set()
    equivalences, preferences = [], []
    all_strains = list(dict.fromkeys(alias for record in records for alias in _declared_strains(record)))
    source_alias_groups = []
    for record in records:
        declared = _declared_strains(record)
        organism = next(iter(record["source"].get("organism", [])), record["organism"])
        _, matched = _organism_strain_suffix(organism, declared) if organism else (None, [])
        source_alias_groups.append(set([*declared, *matched]))
        all_strains.extend(alias for alias in matched if alias not in all_strains)
    for record in records:
        organisms = set(record["source"].get("organism", []))
        if len(organisms) > 1:
            raise ValueError(f"{directory}: conflicting source organism qualifiers")
        organism = next(iter(organisms), record["organism"])
        if not organism:
            continue
        declared = _declared_strains(record)
        strains = declared
        if declared:
            aliases = set(declared)
            while True:
                expanded = aliases | set().union(*(group for group in source_alias_groups if group & aliases))
                if expanded == aliases:
                    break
                aliases = expanded
            strains = [alias for alias in all_strains if alias in aliases]
        else:
            # A plasmid/contig may omit the strain qualifier that is present
            # on the chromosome. Its organism still must agree after removal
            # of an exact, declared strain suffix.
            _, inherited = _organism_strain_suffix(organism, all_strains)
            strains = inherited
        if not strains:
            strains = list(dict.fromkeys(record["source"].get("isolate", [])))
        taxon, matched = _organism_strain_suffix(organism, strains)
        if not matched and re.split(r"\s*[=;]\s*", organism)[-1] in strains and re.search(r"[=;]", organism):
            raise ValueError(f"{directory}: ambiguous GenBank organism/strain alias chain; use an explicit manifest")
        strain = matched[0] if matched else (strains[0] if strains else "")
        labels.add((taxon, strain))
        if strain:
            equivalences.append((taxon, set([*declared, *matched, strain])))
            preferences.append((taxon, matched or [strain]))
    # Equivalent chromosome/plasmid aliases get one consistent display name.
    # Prefer the first alias in an explicit organism chain, then an exact
    # organism suffix, without merging unrelated declared strains.
    reconciled = set()
    for taxon, strain in labels:
        aliases = {strain} if strain else set()
        while True:
            expanded = aliases | set().union(*(group for name, group in equivalences
                                               if name == taxon and group & aliases))
            if expanded == aliases:
                break
            aliases = expanded
        ordered = sorted(preferences, key=lambda item: len(item[1]) < 2)
        preferred = next((alias for name, group in ordered if name == taxon
                          for alias in group if alias in aliases), strain)
        reconciled.add((taxon, preferred))
    labels = reconciled
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


def _project_root(path):
    # Its presence declares a separate project even if its JSON is malformed.
    # Discovery must not silently consume that project's samples as primary.
    return (path / "defenseviz_project.json").is_file()


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
    if target == root:
        return None
    # A link to a project's descendant must not bypass its declared boundary,
    # including when the linked project lives outside the discovery root.
    for parent in (target, *target.parents):
        if parent == root:
            break
        if _project_root(parent):
            return "nested_project" if parent == target else "nested_project_descendant"
    if not target.is_relative_to(root):
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
    Hidden trees, nested projects and identified DefenseViz output roots are
    excluded. A project at the discovery root remains eligible. Ordinary
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
        if directory.resolve() != root and _project_root(directory):
            skipped.append(dict(path=str(directory), reason="nested_project"))
            continue
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

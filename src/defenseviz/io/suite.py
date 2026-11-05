"""Import saved DNMBsuite folders, retaining native identities and decisions.

No detector, model scoring or recognition transfer is invoked here. Native
completion metadata is separate from the consistency of the selected live files.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import time
from collections import defaultdict
from importlib.resources import files
from pathlib import Path

from .project import Reader, _canonical, _family_maps, _string

MODULES = {"defensefinder": "DefenseFinder", "padloc": "PADLOC",
           "defensepredictor": "DefensePredictor", "rebasefinder": "REBASE"}
NATIVE_DEFAULTS = {
    "defensefinder": {"calls": "defensefinder_best_solution_genes.tsv", "systems": "defensefinder_systems.tsv",
                      "id_map": "defensefinder_id_map.tsv", "profiles": "defensefinder_hmmer.tsv"},
    "padloc": {"calls": "padloc_query_proteins_padloc.csv", "id_map": "padloc_query_proteins_id_map.tsv",
               "fasta": "padloc_query_proteins.faa", "profiles": "padloc_query_proteins.domtblout"},
    "defensepredictor": {"calls": "defense_predictor_output.csv", "id_map": "defense_predictor_query_id_map.tsv",
                        "fasta": "defense_predictor_query_protein.faa"},
    "rebasefinder": {"calls": "DNMB_REBASEfinder_augmented_hits.tsv", "workbook": "R-M_REBASE_analysis.xlsx"},
}
MISSING = {"", "NA", "NaN", "nan", "null", "None", "-"}
RM = {f"RM_Type_{t}": f"Type {t}" for t in ("I", "II", "III", "IV")}
RM.update({f"RM_type_{t}": f"Type {t}" for t in ("I", "II", "III", "IV")})
RM.update({f"Type {t}": f"Type {t}" for t in ("I", "II", "III", "IV")})
RM.update({name: "Type II" for name in ("RM_Type_IIG", "RM_type_IIG", "Type IIG")})
# Exact system taxonomy metadata, not protein function or catalytic-role rules.
# This legacy model is defined in the official tagged DefenseFinder 2.0.2 tree.
DF_RM_ALIASES = {"RM_Type_IV_1": "Type IV"}
DF_RM_ALIAS_REFERENCE = "https://github.com/mdmparis/defense-finder-models/blob/2.0.2/definitions/RM/RM/RM_Type_IV_1.xml"
# Named native systems whose taxonomy is not spelled out by the parent label.
# This is model taxonomy only: it confers neither a CDS catalytic role nor a
# retained DNMB annotation. McrBC/GmrSD XML explicitly use Type IV profiles;
# EcoKMcrA is a modification-dependent restriction system (linked primary work).
DF_NAMED_TYPE_IV_MODELS = {
    "McrBC": ("defense-finder-models/RM/McrBC/McrBC", "definitions/RM/McrBC/McrBC.xml"),
    "GmrSD_RM_Type_IV": ("defense-finder-models/RM/GmrSD_RM_Type_IV/GmrSD_RM_Type_IV",
                        "definitions/RM/GmrSD_RM_Type_IV/GmrSD_RM_Type_IV.xml"),
    "EcoKMcrA": ("defense-finder-models/DefenseFinder/EcoKMcrA/EcoKMcrA",
                "definitions/DefenseFinder/EcoKMcrA/EcoKMcrA.xml"),
}
DF_NAMED_TYPE_IV_SOURCE = "https://github.com/mdmparis/defense-finder-models/blob/afb0e5a8b466be53586b13266f5d38d98c3ac268/"
DF_ECOKMCRA_REFERENCE = "https://doi.org/10.1093/nar/gky731"
REFERENCE_PAIRS = (
    ("original_raw", "raw_hit_label", "raw_blast_identity", ("raw_family_id", "family_id", "rm_type", "final_family")),
    ("legacy_blast", "blast_match", "blast_identity", ("family_id", "rm_type", "final_family")),
    ("native_blast", "hit_label", "blast_identity", ("family_id", "rm_type", "final_family")),
    ("curated_blast", "curated_blast_match", "curated_blast_identity", ("curated_blast_family", "family_id", "rm_type", "final_family")),
    # Supplemental matches may disagree with the context candidate's family.
    # Only a subject-specific exported family can classify this observation.
    ("supplemental_blast", "supplemental_blast_match", "supplemental_blast_identity", ("supplemental_blast_family",)),
)


def _value(v):
    return None if v is None or (isinstance(v, str) and v.strip() in MISSING) else v


def _bool(v, where):
    if _value(v) is None:
        return None
    if isinstance(v, bool):
        return v
    if str(v).lower() in {"true", "1"}:
        return True
    if str(v).lower() in {"false", "0"}:
        return False
    raise ValueError(f"{where}: expected explicit true/false")


def _number(v, where):
    if _value(v) is None:
        return None
    try:
        result = float(v)
    except (ValueError, TypeError):
        raise ValueError(f"{where}: expected finite numeric value") from None
    if not math.isfinite(result):
        raise ValueError(f"{where}: expected finite numeric value")
    return result


def _protein(sequence, where):
    sequence = "".join(str(sequence).split()).upper()
    if not sequence or any(c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ*" for c in sequence):
        raise ValueError(f"{where}: invalid protein sequence")
    return hashlib.sha256(sequence.encode("ascii")).hexdigest(), len(sequence)


def _native_context_label(row, label):
    """Recognize DNMB's reserved non-BLAST candidate labels, not free text."""
    if row.get("source") != "rebasefinder" or not isinstance(label, str):
        return False
    if row.get("evidence_mode") in {"operon_context", "operon_context_rebase_blast"}:
        parts = label.split(":", 2)
        return (len(parts) == 3 and parts[0] in {"typeI_context", "typeII_context", "typeIII_context"}
                and parts[1] in {"R", "M", "S"} and bool(parts[2]))
    return (row.get("evidence_mode") == "annotation_motif_candidate"
            and label == "typeIV_context:Mrr_or_modified_DNA_REase")


def _reference_pair(row, where):
    """Use one coherent exported subject/measurement pair, without borrowing."""
    for origin, subject_key, identity_key, family_keys in REFERENCE_PAIRS:
        subject, identity = _value(row.get(subject_key)), _value(row.get(identity_key))
        if subject is not None and identity is not None and not _native_context_label(row, subject):
            return origin, subject, _number(identity, f"{where} {identity_key}"), family_keys
    return None


def _reference_observation_state(row, pair):
    """Missing BLAST measurements are unresolved observations, not no-hit rows."""
    if pair is not None:
        return "coherent"
    explicit_subjects = tuple(subject for _, subject, _, _ in REFERENCE_PAIRS if subject != "hit_label")
    explicit_identities = tuple(identity for _, _, identity, _ in REFERENCE_PAIRS)
    if (any(_value(row.get(k)) is not None and not _native_context_label(row, row[k]) for k in explicit_subjects)
            or any(_value(row.get(k)) is not None for k in explicit_identities)):
        return "incomplete"
    # Generic hit_label is also used by explicitly HMM-only native exports.
    hmm_only = (str(row.get("source", "")).lower() in {"hmm", "hmmer", "hmm_only"}
                or row.get("evidence_kind") == "hmm_hit")
    if _value(row.get("hit_label")) is not None and not hmm_only and not _native_context_label(row, row["hit_label"]):
        return "incomplete"
    return "absent"


def _legacy_df_taxonomy(row, where):
    """Recover missing export columns from the exact native model hierarchy.

    This only projects taxonomy already present in a recognized FQN namespace.
    It supplies neither a protein role nor a producer/model version.
    """
    fqn = _value(row.get("model_fqn"))
    parts = fqn.split("/") if isinstance(fqn, str) else []
    if (len(parts) == 4 and parts[0] == "defense-finder-models"
            and parts[1] in {"DefenseFinder", "RM", "AntiDefenseFinder"}
            and all(part and part.strip() == part for part in parts)):
        parent = "RM" if parts[1] == "RM" else parts[2]
    elif len(parts) == 2 and parts[0] == "CasFinder" and parts[1] and parts[1].strip() == parts[1]:
        parent = "CAS"
    else:
        raise ValueError(f"{where}: missing DefenseFinder type requires an exact recognized model_fqn hierarchy")
    subtype = _value(row.get("subtype"))
    if subtype is not None and subtype != parts[-1]:
        raise ValueError(f"{where}: DefenseFinder subtype conflicts with exact model_fqn leaf")
    return dict(rule="native_model_fqn_taxonomy_projection_v1", source_field="model_fqn",
                model_fqn=fqn, type=parent, subtype=parts[-1])


def _native_family_display(provider, family, native, evidence_kind):
    """Bounded producer taxonomy for display, never catalytic-role inference."""
    if evidence_kind != "final_member":
        return None
    display = category = rule = None
    reference = None
    if provider == "DefenseFinder":
        parent, fqn = _value(native.get("type")), _value(native.get("model_fqn"))
        lineage = fqn.split("/") if isinstance(fqn, str) else []
        cas_namespace = len(lineage) == 2 and lineage[0] == "CasFinder" and lineage[-1] == family
        cas_legacy = (fqn is None and parent in {None, "CAS", "CRISPR-Cas", family}
                      and (family == "CAS_Cluster" or re.fullmatch(
                          r"CAS_(?:Class1-Subtype-(?:I|III|IV)|Class2-Subtype-(?:II|V|VI))(?:-[A-Z][0-9]*)?", family)))
        if cas_namespace or cas_legacy or (fqn is None and parent in {"CAS", "CRISPR-Cas"}):
            display, category, rule = "CRISPR-Cas", "CRISPR-Cas", "defensefinder_cas_family_display_v1"
            reference = "https://github.com/macsy-models/CasFinder/tree/3.1.0"
        cbass_namespace = lineage == ["defense-finder-models", "DefenseFinder", "CBASS", family]
        cbass_name = re.fullmatch(r"CBASS_(?:I|II|III|IV)(?:_(?:Cap4|Endonuclease|HNH|Other|PNP|Phospholipase|TIR|TM))?", family)
        if cbass_name and ((parent == "CBASS" and fqn is None) or cbass_namespace):
            if parent not in {None, "CBASS", family}:
                raise ValueError("DefenseFinder CBASS namespace contradicts native parent")
            display, category, rule = "CBASS", "Abortive infection", "defensefinder_cbass_family_display_v1"
            reference = "https://github.com/mdmparis/defense-finder-models/tree/2.0.2/definitions/DefenseFinder/CBASS"
        pago = {"pAgo_LongA": ("pAgo LongA", "Nuclease"),
                "pAgo_LongB_SIRim": ("pAgo LongB-SIRim", "Abortive infection")}
        if family in pago and ((fqn is None and parent in {"pAgo", family}) or
                               lineage == ["defense-finder-models", "DefenseFinder", "pAgo", family]):
            if parent not in {None, "pAgo", family}:
                raise ValueError("DefenseFinder pAgo namespace contradicts native parent")
            display, category = pago[family]
            rule = "defensefinder_exact_pago_family_display_v1"
            reference = DF_NAMED_TYPE_IV_SOURCE + "definitions/DefenseFinder/pAgo/" + family + ".xml"
    elif provider == "PADLOC":
        if re.fullmatch(r"cas_type_(?:(?:I|II|III|IV|V|VI)(?:-[A-Z][0-9]*)?|arrays|other)", family):
            display, category, rule = "CRISPR-Cas", "CRISPR-Cas", "padloc_cas_family_display_v1"
        elif re.fullmatch(r"cbass_(?:type_(?:I|II|III|IV|IIs)|other)", family):
            display, category, rule = "CBASS", "Abortive infection", "padloc_cbass_family_display_v1"
        elif match := re.fullmatch(r"PDC-([MS])[0-9]{2}", family):
            display, category, rule = f"PDC-{match[1]} (pooled)", "Unresolved", "padloc_pdc_family_display_v1"
        elif family == "argonaute_solo":
            # This native model accepts several pAgo/PIWI classes. It is not
            # equivalent to DefenseFinder's specific LongA model.
            display, category, rule = "Argonaute", "Unresolved", "padloc_untyped_argonaute_display_v1"
        if rule:
            reference = "https://github.com/padlocbio/padloc-db/blob/v2.0.0/sys/" + family + ".yaml"
    if rule is None:
        return None
    return dict(display_family=display, category=category, rule=rule, reference=reference,
                inputs={k: native[k] for k in ("type", "subtype", "model_fqn", "system") if k in native})


def _sheets(path):
    import zipfile
    import xml.etree.ElementTree as ET
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("xl/workbook.xml"))
    return [x.attrib["name"] for x in root.findall("{*}sheets/{*}sheet")]


class SuiteReader:
    def __init__(self, root, manifest, rscript, strict):
        self.root = root
        self.reader = Reader(manifest or root / "suite-project.json")
        self.manifest = manifest
        self.rscript = rscript
        self.strict = strict
        self.data = {k: [] for k in ("samples", "calls", "genes", "coverage", "saved_decisions", "annotations",
                                     "effector_evidence", "groups", "complexes", "source_metadata", "evidence_records")}
        self.raw = {k: [] for k in ("sources", "identity_maps", "coverage", "features", "contexts", "motifs",
                                    "structures", "profiles", "scores", "annotations", "reference_matches", "conflicts",
                                    "capabilities", "warnings", "commands", "unmapped_records", "unresolved_contexts", "software_resources", "table_transforms")}
        self.raw.update(schema_version=1, source_mode="saved_native_import")
        self.data["suite"] = self.raw
        self.genes = {}
        self.source_index = {}
        self.call_index = {}
        self.saved_index = {}
        self.mapping_records = {}
        self.data["input_roles"] = {}

    def pin(self, path):
        path = Path(path).resolve()
        before = path.stat()
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"{path}: input changed while hashing")
        digest = digest.hexdigest()
        old = self.reader.hashes.get(str(path))
        if old and old != digest:
            raise ValueError(f"{path}: input changed while loading suite")
        self.reader.hashes[str(path)] = digest
        if path not in self.reader.inputs:
            self.reader.inputs.append(path)
        return {"path": self.reader.label(path), "sha256": digest}

    def source(self, path, sample, module, role, *, sheet=None):
        ref = self.pin(path)
        if sheet:
            ref["sheet"] = sheet
        resolved = str(Path(path).resolve())
        native_result = role in {"native_calls", "saved_candidates", "membership_crosscheck", "hmm_hits", "hmm_domains", "domain_sites", "domain_hits", "functional_annotations"} or role.startswith("DNMB_REBASEfinder_") or role == "foldseek_results.tsv"
        input_role = dict(sample_id=sample, module=module, kind=role, native_result=native_result)
        roles = self.data["input_roles"].setdefault(resolved, [])
        if input_role not in roles:
            roles.append(input_role)
        key = (resolved, sheet)
        if key not in self.source_index:
            sid = "source-" + str(len(self.raw["sources"]) + 1)
            self.source_index[key] = sid
            self.raw["sources"].append(dict(source_id=sid, sample_id=sample, module=module, role=role,
                                            selected=True, resolved_path=resolved, native_result=native_result, **ref))
        return dict(ref)

    def table(self, path, sample, module, role, *, sheet=None):
        ref = self.source(path, sample, module, role, sheet=sheet)
        previous_limit = csv.field_size_limit()
        try:
            csv.field_size_limit(max(previous_limit, 64 * 1024 * 1024))
            if Path(path).suffix.lower() == ".tsv":
                # DNMB writes native TSV with quote=FALSE. Literal quote characters
                # in annotation text are data, not CSV escape delimiters.
                with Path(path).open(encoding="utf-8-sig", newline="") as stream:
                    lines = iter(stream)
                    first = next(lines, "")
                    while first.startswith("# ") or (first and not first.strip()):
                        first = next(lines, "")
                    import itertools
                    parser = csv.DictReader(itertools.chain([first], lines), delimiter="\t", quoting=csv.QUOTE_NONE)
                    self.reader.headers(parser.fieldnames, path)
                    headers = set(parser.fieldnames)
                    rows = list(parser)
                    if any(None in row or any(v is None for v in row.values()) for row in rows):
                        raise ValueError(f"{path}: row/header column count differs")
                self.pin(path)
            else:
                options = {"literal_hyperlinks": True} if Path(path).suffix.lower() == ".xlsx" and role == "gene_annotations" else {}
                rows, headers = self.reader.table(Path(path), sheet=sheet, **options)
                for transform in getattr(self.reader, "table_transforms", []):
                    if transform not in self.raw["table_transforms"]:
                        self.raw["table_transforms"].append(dict(transform))
        finally:
            csv.field_size_limit(previous_limit)
        for row in rows:
            if _value(row.get("sample_id")) is not None and row["sample_id"] != sample:
                raise ValueError(f"{path}: row sample_id conflicts with the declared sample unit")
        # Reader rehashes the exact bytes it parses, enforcing the pin above.
        item = next(x for x in self.raw["sources"] if x["path"] == ref["path"] and x.get("sheet") == sheet)
        item.update(row_count=len(rows), columns=sorted(headers))
        return rows, headers, ref

    def warn(self, message):
        self.raw["warnings"].append(message)

    def conflict(self, sample, module, code, observed, expected, ref=None):
        action = {"malformed_gene_tsv_fallback": "valid_native_workbook_selected", "unbound_native_candidate": "unbound_record_quarantined",
                  "unbound_eggnog_query": "unbound_record_quarantined", "conflicting_eggnog_query": "all_raw_observations_preserved",
                  "saved_dp_flag_score_rule_mismatch": "native_saved_flag_preserved"}.get(code, "selected_live_artifact_preserved")
        c = dict(sample_id=sample, module=module, code=code, severity="warning", observed=observed,
                 expected=expected, action=action, source_ref=ref)
        self.raw["conflicts"].append(c)
        self.warn(f"{sample}/{module}: {code}; {action.replace(chr(95), chr(32))}")
        if self.strict:
            raise ValueError(f"{sample}/{module}: {code}; strict input consistency failed")

    def gene(self, sample, gene_id, native=None, ref=None):
        gene_id = _string(gene_id, f"{sample} gene ID")
        key = sample, gene_id
        if key not in self.genes:
            self.genes[key] = dict(sample_id=sample, gene_id=gene_id)
        result = self.genes[key]
        native = native or {}
        fields = {"contig": native.get("contig"), "strand": native.get("direction", native.get("strand")),
                  "start": native.get("start"), "end": native.get("end")}
        for field, value in fields.items():
            if _value(value) is None and not (field == "strand" and value == "-"):
                continue
            if field in {"start", "end"}:
                n = _number(value, f"{sample}/{gene_id} {field}")
                if n < 1 or n != int(n):
                    raise ValueError(f"{sample}/{gene_id}: coordinates must be positive integers")
                value = int(n)
            if field == "strand" and value not in {"+", "-"}:
                raise ValueError(f"{sample}/{gene_id}: strand must be + or -")
            if field in result and result[field] != value:
                raise ValueError(f"{sample}/{gene_id}: conflicting {field} in exact identity mapping")
            result[field] = value
        if result.get("start", 0) > result.get("end", float("inf")):
            raise ValueError(f"{sample}/{gene_id}: start exceeds end")
        if _value(native.get("translation")) is not None:
            digest, length = _protein(native["translation"], f"{sample}/{gene_id}")
            if "protein_sha256" in result and result["protein_sha256"] != digest:
                raise ValueError(f"{sample}/{gene_id}: conflicting protein sequence")
            normalized = "".join(str(native["translation"]).split()).upper()
            result.update(protein_sha256=digest, sequence_sha256=digest, protein_length=length,
                          protein_md5=hashlib.md5(normalized.encode("ascii")).hexdigest())
        if ref:
            refs = result.setdefault("source_refs", [])
            if ref not in refs:
                refs.append(ref)
        return result

    def resolve(self, sample, value, mapping, where):
        value = _string(value, where)
        if mapping is not None:
            if value not in mapping:
                raise ValueError(f"{where}: ID map has no entry for {value!r}")
            return mapping[value]
        if (sample, value) not in self.genes:
            raise ValueError(f"{where}: no exact gene binding; supply the native ID map or an explicit manifest map")
        return value

    def metadata(self, directory, sample):
        paths = [directory / ".dnmb_module_stage_cache.rds"]
        for mod in MODULES:
            md = directory / f"dnmb_module_{mod}"
            paths.extend([md / f".dnmb_{mod}_complete", md / ".input_signature.rds"])
        existing = [p for p in paths if p.is_file()]
        for p in existing:
            self.source(p, sample, "suite", "completion_metadata")
        if not existing:
            return {}
        executable = self.rscript or shutil.which("Rscript")
        if not executable:
            self.warn(f"{sample}: RDS completion metadata unavailable without Rscript; live files remain importable")
            return {}
        helper = Path(str(files("defenseviz").joinpath("data/read_suite_metadata.R")))
        self.raw["software_resources"].append(dict(resource="package:defenseviz/data/read_suite_metadata.R", sha256=hashlib.sha256(helper.read_bytes()).hexdigest()))
        argv = [str(executable), "--vanilla", str(helper), str(directory)]
        started = time.monotonic()
        try:
            proc = subprocess.run(argv, capture_output=True, text=True, timeout=60, check=False)
        except (OSError, subprocess.TimeoutExpired) as exc:
            if self.rscript:
                raise ValueError(f"RDS metadata helper failed: {exc}") from exc
            self.warn(f"{sample}: RDS metadata helper unavailable: {exc}")
            return {}
        self.raw["commands"].append(dict(argv=argv, exit_code=proc.returncode,
                                         elapsed_s=round(time.monotonic() - started, 6), stderr=proc.stderr))
        if proc.returncode:
            raise ValueError(f"{sample}: RDS metadata helper failed: {proc.stderr.strip()}")
        try:
            metadata = json.loads(proc.stdout)
        except ValueError as exc:
            raise ValueError(f"{sample}: RDS metadata helper did not return JSON") from exc
        if metadata.get("reader_schema_version") != 1:
            raise ValueError(f"{sample}: unsupported RDS helper schema")
        for p in existing:
            self.pin(p)
        self.data["source_metadata"].append(dict(sample_id=sample, provider="DNMB", source_kind="native_suite_metadata",
                                                 producer_commit=None, metadata=metadata))
        return metadata

    def module_path(self, directory, spec, module, kind):
        override = spec.get("modules", {}).get(module, {})
        if not isinstance(override, dict):
            raise ValueError(f"{module}: module override must be an object")
        if kind in override:
            raw = _string(override[kind], f"{module}.{kind}")
            path = Path(raw).expanduser()
            return (path if path.is_absolute() else self.reader.base / path).resolve()
        default = NATIVE_DEFAULTS[module].get(kind)
        return directory / f"dnmb_module_{module}" / default if default else None

    def mapping(self, directory, spec, module, sample):
        path = self.module_path(directory, spec, module, "id_map")
        if path is None or not path.is_file():
            if "id_map" in spec.get("modules", {}).get(module, {}):
                raise ValueError(f"{path}: explicit ID map is missing")
            return None
        rows, headers, ref = self.table(path, sample, module, "identity_map")
        source = next((c for c in ("source_id", "defensefinder_id", "query_id") if c in headers), None)
        dest = next((c for c in ("gene_id", "locus_tag") if c in headers), None)
        if not source or not dest:
            raise ValueError(f"{path}: ID map requires source_id/defensefinder_id/query_id and gene_id/locus_tag")
        mapping, destinations = {}, set()
        for number, row in enumerate(rows, 2):
            a = _string(row[source], f"{path}:{number} source ID")
            b = _string(row[dest], f"{path}:{number} gene ID")
            if a in mapping or b in destinations:
                raise ValueError(f"{path}:{number}: ID map must be one-to-one within sample")
            mapping[a] = b
            self.mapping_records[(sample, module, a)] = row
            destinations.add(b)
            self.gene(sample, b, row, {**ref, "row": number})
            self.raw["identity_maps"].append(dict(sample_id=sample, namespace=module, source_id=a, gene_id=b,
                                                    source_ref={**ref, "row": number}, native_fields=row))
        return mapping

    def fasta(self, path, sample, mapping=None):
        if not path or not path.is_file():
            return
        ref = self.source(path, sample, "suite", "protein_sequences")
        name, parts, seen = None, [], set()
        def flush():
            if name is None:
                return
            if name in seen:
                raise ValueError(f"{path}: duplicate FASTA identifier {name!r}")
            seen.add(name)
            gene = self.resolve(sample, name, mapping, str(path))
            self.gene(sample, gene, {"translation": "".join(parts)}, {**ref, "record_id": name})
        with path.open(encoding="utf-8") as stream:
            for line in stream:
                if line.startswith(">"):
                    flush()
                    name = _string(line[1:].strip().split()[0] if line[1:].strip() else "", str(path))
                    parts = []
                elif line.strip():
                    if name is None:
                        raise ValueError(f"{path}: sequence before FASTA header")
                    parts.append(line.strip())
            flush()
        self.pin(path)

    def gene_annotations(self, directory, spec, sample):
        selected = spec.get("genes")
        sheet = None
        if selected:
            path = self.reader.path(selected)
        else:
            path = directory / "dnmb_module_rebasefinder/rebasefinder_input.tsv"
            if not path.is_file():
                choices = sorted(directory.glob("*_total.xlsx"))
                if len(choices) > 1:
                    raise ValueError(f"{directory}: multiple total workbooks; specify manifest genes")
                if not choices:
                    return
                path = choices[0]
                sheet = "1.GenBank_table"
        try:
            rows, headers, ref = self.table(path, sample, "suite", "gene_annotations", sheet=sheet)
        except ValueError as exc:
            # Historical R quote=FALSE TSVs can contain literal line breaks in
            # annotation cells. Do not guess where to stitch broken records.
            choices = sorted(directory.glob("*_total.xlsx")) if not selected else []
            if "row/header column count differs" not in str(exc) or len(choices) != 1:
                raise
            self.conflict(sample, "suite", "malformed_gene_tsv_fallback", str(exc), "valid gene worksheet", self.pin(path))
            for source in self.raw["sources"]:
                if source["resolved_path"] == str(path.resolve()):
                    source.update(selected=False, selection_reason="malformed TSV; using full gene workbook")
            path, sheet = choices[0], "1.GenBank_table"
            rows, headers, ref = self.table(path, sample, "suite", "gene_annotations", sheet=sheet)
        idcol = "locus_tag" if "locus_tag" in headers else "gene_id"
        if idcol not in headers:
            raise ValueError(f"{path}: gene table needs locus_tag or gene_id")
        seen = set()
        for number, row in enumerate(rows, 2):
            gene = _string(row[idcol], f"{path}:{number} gene ID")
            if gene in seen:
                raise ValueError(f"{path}:{number}: duplicate gene annotation identity")
            seen.add(gene)
            self.gene(sample, gene, row, {**ref, "row": number})
            self.raw["annotations"].append(dict(sample_id=sample, gene_id=gene, source_ref={**ref, "row": number}, native_fields=row))

    def coverage(self, sample, module, spec, metadata, present, rows, selected_ref):
        provider = MODULES[module]
        cache = metadata.get("cache") or {}
        cached = (cache.get("modules") or {}).get("REBASEfinder" if module == "rebasefinder" else provider, {})
        sentinel = (metadata.get("sentinels") or {}).get("REBASEfinder" if module == "rebasefinder" else provider) or {}
        override = spec.get("modules", {}).get(module, {})
        declared = override.get("status")
        if declared is not None and declared not in {"complete", "partial", "missing", "not_run"}:
            raise ValueError(f"{sample}/{module}: invalid manifest status")
        valid_sentinel = sentinel.get("sentinel_version") == 1 and sentinel.get("alias", "").lower() == ("REBASEfinder" if module == "rebasefinder" else provider).lower()
        ok = valid_sentinel or cached.get("ok") is True
        counts = cached.get("table_counts") or {}
        completed_zero = (valid_sentinel and sentinel.get("n_hits") == 0) or (cached.get("ok") is True and counts.get("hits") == 0)
        if module == "rebasefinder":
            # Modern n_hits counts retained calls, not all candidate reference
            # observations. Retention zero alone cannot prove reference zero.
            completed_zero = cached.get("ok") is True and counts.get("all_hits") == 0
        completion = "sentinel" if valid_sentinel else "stage_cache" if cached.get("ok") is True else "live_artifact" if present else "unavailable"
        status = "complete" if present or completed_zero else "partial" if ok or cached.get("ok") is False else "missing"
        if declared:
            if present and declared in {"missing", "not_run"}:
                raise ValueError(f"{sample}/{module}: declared {declared} conflicts with supplied live artifact")
            status = declared
            completion = "manifest"
        if sentinel and not valid_sentinel:
            self.warn(f"{sample}/{module}: unrecognized or unreadable completion sentinel")
        if present and cached.get("ok") is False:
            status = "partial"
            self.conflict(sample, module, "cached_failure_with_live_artifact", True, False, selected_ref)
        if present and valid_sentinel and module != "rebasefinder" and sentinel.get("n_hits") is not None and sentinel["n_hits"] != rows:
            self.conflict(sample, module, "sentinel_live_row_count_mismatch", rows, sentinel["n_hits"], selected_ref)
        self.data["coverage"].append(dict(sample_id=sample, caller=provider, status=status,
                                             endpoint="reference_matches" if module == "rebasefinder" else "native_calls"))
        detailed = dict(sample_id=sample, module=module, run_status="complete" if ok else "failed" if cached.get("ok") is False else "unknown",
                        completion_source=completion, artifact_status="present" if present else "completed_zero" if completed_zero else "missing",
                        selected_source=selected_ref, cached_table_counts=cached.get("table_counts", {}), sentinel=sentinel or None)
        self.raw["coverage"].append(detailed)
        return detailed

    def call(self, sample, gene, provider, family, system, ref, *, native=None, evidence_kind="final_member", model=None, mapping=None,
             reference_family=None):
        native = native or {}
        typ = reference_family if reference_family is not None else native.get("type", family)
        rm_type = RM.get(typ) or RM.get(family)
        category = "R-M" if rm_type else None
        mapping_rule = "native_exact_rm_type_v1" if rm_type else None
        mapping_reference = None
        mapping_references = None
        if provider == "DefenseFinder":
            parent = _value(native.get("type"))
            fqn = _value(native.get("model_fqn"))
            lineage = fqn.split("/") if isinstance(fqn, str) else []
            rm_namespace = len(lineage) >= 3 and lineage[:2] == ["defense-finder-models", "RM"]
            cas_namespace = len(lineage) == 2 and lineage[0] == "CasFinder"
            if (rm_namespace or cas_namespace) and _value(native.get("subtype")) is not None and native["subtype"] != lineage[-1]:
                raise ValueError(f"{ref['path']}: DefenseFinder subtype conflicts with exact model_fqn leaf")
            if parent == "RM" or rm_namespace:
                category = "R-M"
                hierarchy_type = RM.get(lineage[2]) if rm_namespace else None
                alias_type = DF_RM_ALIASES.get(family)
                if rm_type and hierarchy_type and rm_type != hierarchy_type:
                    raise ValueError(f"{ref['path']}: contradictory DefenseFinder R-M taxonomy")
                rm_type = rm_type or hierarchy_type or alias_type
                mapping_rule = "defensefinder_rm_parent_v1"
                if hierarchy_type:
                    mapping_rule = "defensefinder_exact_rm_hierarchy_v1"
                if alias_type:
                    if rm_type != alias_type:
                        raise ValueError(f"{ref['path']}: contradictory DefenseFinder R-M subtype alias")
                    mapping_rule = "defensefinder_2.0.2_exact_rm_subtype_v1"
                    mapping_reference = DF_RM_ALIAS_REFERENCE
            named = DF_NAMED_TYPE_IV_MODELS.get(family) if evidence_kind == "final_member" else None
            if named is not None and (fqn is None or fqn == named[0]):
                if parent not in {None, "RM", family, "Type IV", "RM_Type_IV", "RM_type_IV"} or (rm_type and rm_type != "Type IV"):
                    raise ValueError(f"{ref['path']}: named DefenseFinder Type IV model contradicts native taxonomy")
                category, rm_type = "R-M", "Type IV"
                mapping_rule = "defensefinder_exact_named_type_iv_model_v1"
                mapping_reference = DF_NAMED_TYPE_IV_SOURCE + named[1]
                mapping_references = [mapping_reference]
                if family == "EcoKMcrA":
                    mapping_references.append(DF_ECOKMCRA_REFERENCE)
            if parent in {"CAS", "CRISPR-Cas"} or cas_namespace:
                if category == "R-M":
                    raise ValueError(f"{ref['path']}: contradictory DefenseFinder RM/CAS taxonomy")
                category = "CRISPR-Cas"
                mapping_rule = "defensefinder_casfinder_namespace_v1" if cas_namespace else "defensefinder_cas_parent_v1"
                mapping_reference = "https://github.com/macsy-models/CasFinder/tree/3.1.0"
        raw = dict(sample_id=sample, gene_id=gene, caller=provider, family=family, system_id=system,
                   role=None, rm_type=rm_type, category=category,
                   evidence_kind=evidence_kind, native_family=family, native_type=typ, native_model=model,
                   source=ref["path"], source_ref=ref, run_id=ref["sha256"], run_id_kind="selected_source_sha256", selected=evidence_kind != "hmm_hit")
        display = _native_family_display(provider, family, native, evidence_kind)
        effective_mapping = mapping or {}
        explicit_display = any(m.get("source") == "family_map" for key in ((provider, family), (None, family))
                               if (m := effective_mapping.get(key)) is not None)
        if display and not explicit_display:
            effective_mapping = {**effective_mapping, (provider, family): dict(
                category=display["category"], display_family=display["display_family"], rm_type=None,
                priority=0.5, source="native_taxonomy", caller=provider)}
        elif display and explicit_display and rm_type is None:
            # The optional user classification takes precedence over display
            # taxonomy defaults. Native fields and source provenance stay intact.
            raw["category"] = None
        raw = _canonical(raw, effective_mapping, ref["path"])
        if display and not explicit_display:
            raw["family_mapping_rule"] = display["rule"]
            raw["family_mapping_inputs"] = display["inputs"]
            raw["family_mapping_reference"] = display["reference"]
        if mapping_rule:
            raw["mapping_rule"] = mapping_rule
            raw["mapping_inputs"] = {key: native[key] for key in ("type", "subtype", "model_fqn") if key in native}
            if mapping_reference:
                raw["mapping_reference"] = mapping_reference
            if mapping_references:
                raw["mapping_references"] = mapping_references
            if raw["category_source"] in {"explicit", "native_type"}:
                raw["category_source"] = "native_system_taxonomy"
        raw["native_fields"] = native
        if provider == "DefenseFinder" and _value(native.get("replicon")) is not None:
            raw["replicon_id"] = native["replicon"]
        if self.genes[(sample, gene)].get("protein_sha256"):
            raw["protein_sha256"] = self.genes[(sample, gene)]["protein_sha256"]
        identity = tuple(raw.get(k) for k in ("sample_id", "gene_id", "caller", "run_id", "replicon_id", "system_id", "native_family", "native_model", "evidence_kind"))
        previous = self.call_index.get(identity)
        if previous is not None:
            previous.setdefault("source_refs", [previous["source_ref"]]).append(ref)
            previous.setdefault("native_observations", [previous["native_fields"]]).append(native)
            return previous
        self.call_index[identity] = raw
        self.data["calls"].append(raw)
        return raw

    def dp_decision_rule(self, override, sample):
        """Resolve a score rule without replacing a saved classification.

        Older DNMB wrappers used >= 4; other producers use > 4. The legacy
        fallback is an explicit assumption, not verified producer metadata.
        """
        source_rule = {}
        metadata_ref = None
        if "decision_metadata" in override:
            path = self.reader.path(override["decision_metadata"])
            metadata_ref = self.source(path, sample, "defensepredictor", "decision_metadata")
            native = self.reader.json(path)
            if not isinstance(native, dict) or not {"threshold", "operator"} <= native.keys():
                raise ValueError(f"{path}: DefensePredictor decision_metadata requires threshold and operator")
            source_rule = {"score_threshold": native["threshold"], "score_operator": native["operator"]}
        for source, candidate in (("manifest", override), ("decision_metadata", source_rule)):
            if "score_threshold" in candidate:
                threshold = candidate["score_threshold"]
                if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(threshold):
                    raise ValueError(f"{sample}/{source}: score_threshold must be finite numeric")
            if "score_operator" in candidate and candidate["score_operator"] not in (">", ">="):
                raise ValueError(f"{sample}/{source}: score_operator must be '>' or '>='")
        for key in ("score_threshold", "score_operator"):
            if key in override and key in source_rule and override[key] != source_rule[key]:
                raise ValueError(f"{sample}: manifest {key} contradicts DefensePredictor decision_metadata")
        rule = {}
        for key, output_key, default, assumed_source in (
            ("score_threshold", "threshold", 4.0, "legacy_default_assumption_4.0"),
            ("score_operator", "score_operator", ">=", "legacy_inclusive_assumption"),
        ):
            rule[output_key] = override.get(key, source_rule.get(key, default))
            rule[output_key + "_source"] = ("manifest" if key in override else
                                             "decision_metadata" if key in source_rule else assumed_source)
        rule["score_rule_explicit"] = all(key in override or key in source_rule
                                          for key in ("score_threshold", "score_operator"))
        if metadata_ref is not None:
            rule["decision_metadata_source_ref"] = metadata_ref
        return rule

    def comparator(self, directory, spec, module, sample, metadata, families):
        provider = MODULES[module]
        mapping = self.mapping(directory, spec, module, sample)
        override = spec.get("modules", {}).get(module, {})
        dp_rule = self.dp_decision_rule(override, sample) if module == "defensepredictor" else None
        selected = override.get("calls")
        multipart = isinstance(selected, list)
        if multipart:
            if module != "defensefinder" or not selected:
                raise ValueError("Only DefenseFinder calls accepts a nonempty list of artifact paths")
            paths = [self.reader.path(_string(p, "DefenseFinder calls path")) for p in selected]
            if len(set(paths)) != len(paths):
                raise ValueError("DefenseFinder calls paths must be unique")
        else:
            paths = [self.module_path(directory, spec, module, "calls")]
        for path in paths:
            if not path.is_file() and "calls" in override:
                raise ValueError(f"{path}: explicit call file is missing")
        systems_path = self.module_path(directory, spec, module, "systems")
        if module == "defensefinder" and not multipart and not paths[0].is_file() and systems_path and systems_path.is_file():
            paths = [systems_path]
        artifacts = [(path, *self.table(path, sample, module, "native_calls")) for path in paths if path.is_file()]
        refs = [ref for path, rows, headers, ref in artifacts]
        selected_ref = refs[0] if len(refs) == 1 else None
        detailed = self.coverage(sample, module, spec, metadata, bool(refs),
                                 sum(len(rows) for path, rows, headers, ref in artifacts), selected_ref)
        if multipart:
            detailed.update(selected_sources=refs, artifact_count=len(refs),
                            source_scope="explicit_native_artifact_list")
        self.fasta(self.module_path(directory, spec, module, "fasta"), sample, mapping)
        groups = defaultdict(list)
        group_family = {}
        group_refs = defaultdict(list)
        for path, rows, headers, ref in artifacts:
            required = {"defensefinder": {"sys_id"}, "padloc": {"system.number", "seqid", "system", "target.name"},
                        "defensepredictor": {"locus_tag", "mean_log_odds"}}[module]
            if module == "defensefinder" and not ({"type", "model_fqn"} & headers):
                raise ValueError(f"{path}: DF requires type or model_fqn")
            if module == "defensefinder" and not ({"hit_id", "protein_in_syst"} & headers):
                raise ValueError(f"{path}: DF requires hit_id or protein_in_syst")
            if not required <= headers:
                raise ValueError(f"{path}: missing required columns {sorted(required - headers)}")
        predicted = set()
        if dp_rule is not None and not dp_rule["score_rule_explicit"] and any(
                rows and "predicted_general_defense" not in headers for path, rows, headers, ref in artifacts):
            self.warn(f"{sample}/defensepredictor: no saved classification or fully explicit score rule; "
                      "legacy threshold/operator assumptions recorded in decision_rule")
        entries = ((path, headers, ref, number, row) for path, rows, headers, ref in artifacts
                   for number, row in enumerate(rows, 2))
        for path, headers, ref, number, row in entries:
            rr = {**ref, "row": number}
            taxonomy_projection = None
            if module == "defensefinder":
                if "protein_in_syst" in headers:
                    ids = [x.strip() for x in row["protein_in_syst"].split(",")]
                    if not all(ids) or len(set(ids)) != len(ids):
                        raise ValueError(f"{path}:{number}: invalid system member list")
                    if _value(row.get("genes_count")) is not None and _number(row["genes_count"], str(path)) != len(ids):
                        raise ValueError(f"{path}:{number}: genes_count disagrees with members")
                elif "hit_id" in headers:
                    ids = [row["hit_id"]]
                else:
                    raise ValueError(f"{path}: DF requires hit_id or protein_in_syst")
                if "type" not in headers or (_value(row.get("type")) is None and _value(row.get("subtype")) is None):
                    taxonomy_projection = _legacy_df_taxonomy(row, f"{path}:{number}")
                    family = taxonomy_projection["subtype"]
                else:
                    family = _value(row.get("subtype")) or row["type"]
                system = _string(row["sys_id"], str(path))
                model = _value(row.get("gene_name"))
            elif module == "padloc":
                ids = [row["target.name"]]
                family = row["system"]
                system = json.dumps([row["seqid"], row["system.number"]], separators=(",", ":"))
                model = _value(row.get("hmm.name"))
            else:
                gene = self.resolve(sample, row["locus_tag"], None, str(path))
                if gene in predicted:
                    raise ValueError(f"{path}: duplicate DefensePredictor gene prediction")
                predicted.add(gene)
                score = _number(row["mean_log_odds"], str(path))
                if score is None:
                    raise ValueError(f"{path}: missing DefensePredictor score")
                threshold = dp_rule["threshold"]
                score_selected = score > threshold if dp_rule["score_operator"] == ">" else score >= threshold
                if "predicted_general_defense" in headers:
                    selected = _bool(row["predicted_general_defense"], f"{path}:{number}: predicted_general_defense")
                    if selected is None:
                        raise ValueError(f"{path}:{number}: predicted_general_defense requires explicit true/false; missing is unassessed")
                    selection_source = "native_predicted_general_defense"
                    if dp_rule["score_rule_explicit"] and selected != score_selected:
                        self.conflict(sample, module, "saved_dp_flag_score_rule_mismatch", selected,
                                      {"score": score, "operator": dp_rule["score_operator"],
                                       "threshold": threshold, "selected": score_selected}, rr)
                else:
                    selected = score_selected
                    selection_source = "explicit_score_rule" if dp_rule["score_rule_explicit"] else "legacy_score_rule_assumption"
                self.raw["scores"].append(dict(sample_id=sample, gene_id=gene, score_name="mean_log_odds", score=score,
                    sd=_number(row.get("sd_log_odds"), str(path)), min=_number(row.get("min_log_odds"), str(path)),
                    max=_number(row.get("max_log_odds"), str(path)),
                    **{k: v for k, v in dp_rule.items() if k != "decision_metadata_source_ref"},
                    score_rule_selected=score_selected, selection_source=selection_source,
                    selected=selected, source_ref=rr, native_fields=row))
                if selected:
                    call = self.call(sample, gene, provider, "General defense", None, rr, native=row, mapping=families)
                    call["selection_source"] = selection_source
                continue
            for source_id in ids:
                gene = self.resolve(sample, source_id, mapping, str(path))
                binding = "unassessed"
                replicon = _value(row.get("replicon")) if module == "defensefinder" else None
                if module == "defensefinder" and replicon is not None:
                    mapped = self.mapping_records.get((sample, module, source_id), {})
                    expected_contig = _value(mapped.get("contig_id"))
                    if expected_contig is not None and replicon != expected_contig:
                        raise ValueError(f"{path}:{number}: DefenseFinder replicon contradicts the exact gene/contig map")
                    if expected_contig == replicon or self.genes[(sample, gene)].get("contig") == replicon:
                        binding = "verified"
                if module == "padloc":
                    self.gene(sample, gene, {"start": row.get("start"), "end": row.get("end"), "strand": row.get("strand")}, rr)
                    mapped = self.mapping_records.get((sample, module, source_id), {})
                    expected_contig = _value(mapped.get("contig_id"))
                    known_contigs = {g.get("contig") for (sid, _), g in self.genes.items() if sid == sample and g.get("contig")}
                    if expected_contig is None and row["seqid"] in known_contigs:
                        expected_contig = self.genes[(sample, gene)].get("contig")
                    if expected_contig is not None:
                        if row["seqid"] != expected_contig:
                            raise ValueError(f"{path}:{number}: PADLOC replicon contradicts the exact gene/contig map")
                        binding = "verified"
                call = self.call(sample, gene, provider, family, system, rr, native=row, model=model, mapping=families)
                if taxonomy_projection is not None:
                    call["native_taxonomy_projection"] = taxonomy_projection
                    call["native_type"] = row.get("type")
                call["replicon_binding"] = binding
                if replicon is not None:
                    call["replicon_id"] = replicon
                group_key = ref["sha256"], system, replicon
                groups[group_key].append(gene)
                if group_key in group_family and group_family[group_key] != family:
                    raise ValueError(f"{path}:{number}: contradictory family for one native system")
                group_family[group_key] = family
                if ref not in group_refs[group_key]:
                    group_refs[group_key].append(ref)
        if module == "defensepredictor":
            mapped = set(mapping.values()) if mapping is not None else {g for s, g in self.genes if s == sample}
            detailed.update(mapped_gene_count=len(mapped), scored_gene_count=len(predicted), missing_gene_count=len(mapped - predicted))
        for group_key, members in groups.items():
            run_id, system, replicon = group_key
            system_refs = group_refs[group_key]
            unique = list(dict.fromkeys(members))
            contigs = {self.genes[(sample, g)].get("contig") for g in unique}
            bound_calls = [c for c in self.data["calls"] if c["sample_id"] == sample and c["caller"] == provider and c["system_id"] == system and c["run_id"] == run_id and c["evidence_kind"] == "final_member" and (replicon is None or c.get("replicon_id") == replicon)]
            verified = bool(bound_calls) and all(c.get("replicon_binding") == "verified" for c in bound_calls) if module == "padloc" or replicon is not None else None not in contigs and len(contigs) == 1
            group_id = [sample, provider, run_id, system]
            if replicon is not None:
                group_id.insert(-1, replicon)
            self.data["groups"].append(dict(group_id=json.dumps(group_id), sample_id=sample,
                replicon_binding="verified" if verified else "unassessed", replicon_id=replicon,
                contigs=sorted(c for c in contigs if c is not None), provider=provider, run_id=run_id, native_system_id=system, system_id=system, membership_kind="native_system", family=group_family[group_key],
                members=[self.member(sample, g) for g in unique], source_refs=system_refs))
        # Check redundant native DF outputs without adding duplicate memberships.
        if module == "defensefinder" and refs and systems_path not in paths and systems_path and systems_path.is_file():
            sr, sh, sf = self.table(systems_path, sample, module, "membership_crosscheck")
            if not {"sys_id", "protein_in_syst"} <= sh:
                raise ValueError(f"{systems_path}: missing DF system membership columns")
            expected = {(r["sys_id"], self.resolve(sample, i.strip(), mapping, str(systems_path))) for r in sr for i in r["protein_in_syst"].split(",")}
            observed = {(group_key[1], g) for group_key, gs in groups.items() for g in gs}
            if expected != observed:
                self.conflict(sample, module, "genes_systems_membership_mismatch", len(observed), len(expected), selected_ref)
        self.native_profiles(directory, spec, module, sample, mapping, families)
        versions = (((metadata.get("cache") or {}).get("signature") or {}).get("db_state") or {}).get(provider, {})
        source_metadata = dict(sample_id=sample, provider=provider, source_kind="native_module",
            versions=versions, source_ref=selected_ref, source_refs=refs, producer_commit=None)
        if dp_rule is not None:
            source_metadata["decision_rule"] = dict(dp_rule, selection_precedence="saved_flag_then_score_rule",
                saved_flag_present=any("predicted_general_defense" in headers for path, rows, headers, ref in artifacts))
        self.data["source_metadata"].append(source_metadata)
        for call in self.data["calls"]:
            if call["sample_id"] == sample and call["caller"] == provider:
                call["source_version"] = versions
                native = call.get("native_fields", {})
                for field in ("model_version", "database_version"):
                    observations = {_value(r.get(field)) for r in call.get("native_observations", [native])
                                    if _value(r.get(field)) is not None}
                    if len(observations) > 1:
                        raise ValueError(f"{sample}/{provider}: contradictory duplicate native {field}")
                row_model = _value(native.get("model_version"))
                row_database = _value(native.get("database_version"))
                for observed, expected, field in ((row_model, versions.get("models_version"), "model_version"),
                                                  (row_database, versions.get("db_version"), "database_version")):
                    if observed is not None and expected is not None and str(observed) != str(expected):
                        raise ValueError(f"{sample}/{provider}: native {field} conflicts with saved module metadata")
                call["model_version"] = row_model or versions.get("models_version") or row_database or versions.get("db_version")
                call["model_version_source"] = ("native_row.model_version" if row_model is not None else
                    "saved_module_metadata.models_version" if versions.get("models_version") is not None else
                    "native_row.database_version" if row_database is not None else
                    "saved_module_metadata.db_version" if versions.get("db_version") is not None else "unavailable")
                database = row_database or versions.get("db_version")
                if database is not None:
                    call["database_version"] = database
                    call["database_version_source"] = ("native_row.database_version" if row_database is not None
                                                       else "saved_module_metadata.db_version")

    def member(self, sample, gene):
        result = {"gene_id": gene}
        if self.genes[(sample, gene)].get("protein_sha256"):
            result["protein_sha256"] = self.genes[(sample, gene)]["protein_sha256"]
        return result

    def native_profiles(self, directory, spec, module, sample, mapping, families):
        path = self.module_path(directory, spec, module, "profiles")
        if module in {"defensefinder", "padloc"}:
            self.data["coverage"].append(dict(sample_id=sample, caller=MODULES[module],
                status="complete" if path and path.is_file() else "missing", endpoint="hmm_hits"))
        if not path or not path.is_file():
            return
        if module == "defensefinder":
            rows, headers, ref = self.table(path, sample, module, "hmm_hits")
            if not {"hit_id", "gene_name"} <= headers:
                raise ValueError(f"{path}: DF HMM table requires hit_id and gene_name")
            for number, row in enumerate(rows, 2):
                gene = self.resolve(sample, row["hit_id"], mapping, str(path))
                rr = {**ref, "row": number}
                self.raw["profiles"].append(dict(sample_id=sample, gene_id=gene, namespace="DefenseFinder", source_ref=rr, native_fields=row))
                self.call(sample, gene, "DefenseFinder", row["gene_name"], None, rr, native=row,
                          evidence_kind="hmm_hit", model=row["gene_name"], mapping=families)
        else:
            ref = self.source(path, sample, module, "hmm_domains")
            with path.open() as stream:
                for number, line in enumerate(stream, 1):
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip().split(maxsplit=22)
                    if len(fields) < 22:
                        raise ValueError(f"{path}:{number}: malformed HMMER domain record")
                    # Native PADLOC releases use both HMMER orientations. Resolve only
                    # against the explicit map; ambiguous bindings are an error.
                    bound = [x for x in (fields[0], fields[3]) if (mapping is not None and x in mapping) or (mapping is None and (sample, x) in self.genes)]
                    if len(bound) != 1:
                        raise ValueError(f"{path}:{number}: HMMER protein identity is ambiguous or unmapped")
                    gene = self.resolve(sample, bound[0], mapping, str(path))
                    self.raw["profiles"].append(dict(sample_id=sample, gene_id=gene, namespace="PADLOC", source_ref={**ref, "row": number},
                        native_fields=dict(zip(("target_name", "target_accession", "target_length", "query_name", "query_accession", "query_length",
                        "evalue", "score", "bias", "domain_number", "domain_count", "conditional_evalue", "independent_evalue", "domain_score",
                        "domain_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "accuracy", "description"), fields))))
            self.pin(path)

    def rebase(self, directory, spec, sample, metadata):
        module = "rebasefinder"
        path = self.module_path(directory, spec, module, "calls")
        if "calls" in spec.get("modules", {}).get(module, {}) and not path.is_file():
            raise ValueError(f"{path}: explicit REBASE candidate file missing")
        sheet = None
        if not path.is_file():
            path = self.module_path(directory, spec, module, "workbook")
            if path.is_file():
                self.pin(path)
                sheets = _sheets(path)
                sheet = "All_Augmented_Evidence" if "All_Augmented_Evidence" in sheets else "RM_Comprehensive"
        rows, headers, ref = self.table(path, sample, module, "saved_candidates", sheet=sheet) if path.is_file() else ([], set(), None)
        detail = self.coverage(sample, module, spec, metadata, ref is not None, len(rows), ref)
        reference_coverage = next(c for c in self.data["coverage"]
                                  if c["sample_id"] == sample and c["caller"] == "REBASE" and c["endpoint"] == "reference_matches")
        reference_schema = any({subject, identity} <= headers for _, subject, identity, _ in REFERENCE_PAIRS)
        if ref is not None and not reference_schema:
            # A curation/HMM-only export cannot establish a reference-match zero.
            reference_coverage["status"] = "partial"
        unbound_count = unbound_reference_count = incomplete_reference_count = 0
        retained_count = 0
        for number, row in enumerate(rows, 2):
            value = row.get("query", row.get("locus_tag", row.get("gene_id")))
            rr = {**ref, "row": number}
            # Establish the reference endpoint before checking candidate binding.
            # A no-hit candidate can be unresolved without losing a reference hit.
            pair = _reference_pair(row, f"{path}:{number}")
            reference_state = _reference_observation_state(row, pair)
            if reference_state == "incomplete":
                incomplete_reference_count += 1
                reference_coverage["status"] = "partial"
            if (sample, value) not in self.genes:
                unbound_count += 1
                unbound_reference_count += reference_state != "absent"
                self.raw["unmapped_records"].append(dict(sample_id=sample, module=module, source_ref=rr, native_fields=row,
                    assessment_status="unassessed_unbound_identity", has_coherent_reference_pair=pair is not None,
                    reference_observation_state=reference_state))
                self.conflict(sample, module, "unbound_native_candidate", value, "exact canonical gene identity", rr)
                if pair is not None:
                    reference_coverage["status"] = "partial"
                continue
            gene = self.resolve(sample, value, None, str(path))
            self.raw["features"].append(dict(sample_id=sample, gene_id=gene, source_ref=rr, native_fields=row))
            keep = _bool(row.get("curation_keep"), f"{path}:{number} curation_keep")
            state = "unassessed" if keep is None else "retained" if keep else "not_retained"
            retained_count += keep is True
            family = _value(row.get("final_family")) if keep is not None else None
            native_role = _value(row.get("final_role")) if keep is not None else None
            role = native_role if native_role in {"R", "M", "S", "RM", "Aux", "A", "C", "F", "unknown", "Unknown"} else None
            record = dict(sample_id=sample, gene_id=gene, state=state, family=family, role=role,
                rm_type=family if family in {"Type I", "Type II", "Type III", "Type IV"} else None,
                source_ref=rr, score_total=_number(row.get("curation_score"), str(path)),
                decision_reason=_value(row.get("curation_reason")), final_component=_value(row.get("final_component")),
                curation_tier=_value(row.get("curation_tier")), native_role=native_role)
            if self.genes[(sample, gene)].get("protein_sha256"):
                record["protein_sha256"] = self.genes[(sample, gene)]["protein_sha256"]
            context_payload = {k: _value(v) for k, v in row.items()
                               if k.startswith(("typei_", "typeii_", "typeiii_")) or k in {"context_state", "orphan"}}
            record["native_context"] = context_payload
            previous = self.saved_index.get((sample, gene))
            if previous is not None:
                if previous.get("native_context", {}) != context_payload:
                    raise ValueError(f"{sample}/{gene}: contradictory saved context payloads; select one active candidate record")
                semantic = ("state", "family", "role", "rm_type", "score_total", "final_component", "curation_tier", "native_role")
                if any(previous.get(k) != record.get(k) for k in semantic):
                    raise ValueError(f"{sample}/{gene}: contradictory saved decision rows; select one decision source")
                previous.setdefault("source_refs", [previous["source_ref"]]).append(rr)
            else:
                self.saved_index[(sample, gene)] = record
                self.data["saved_decisions"].append(record)
            # A matched reference is metadata, never a retained annotation by itself.
            # Keep each reference subject paired with its own measurements.
            # Raw/original REBASE observations take precedence over a later
            # curated alternative; an incomplete pair never borrows identity.
            if pair is not None:
                origin, subject, identity, family_keys = pair
                native_family = next((_value(row.get(k)) for k in family_keys if _value(row.get(k)) is not None), "Unresolved")
                self.raw["reference_matches"].append(dict(sample_id=sample, gene_id=gene, reference_namespace="REBASE",
                    reference_id=subject, identity=identity, match_origin=origin, scope="exported_reference_observations",
                    source_ref=rr, native_fields=row))
                reference_call = self.call(sample, gene, "REBASE", native_family, None, rr, native=row, evidence_kind="reference_hit", model=subject,
                                           reference_family=native_family if origin == "supplemental_blast" else None)
                reference_call.update(reference_id=subject, match_identity=identity, match_origin=origin,
                                      reference_scope="exported_reference_observations")
                if origin == "supplemental_blast":
                    reference_call["reference_family_source"] = "supplemental_blast_family" if _value(row.get("supplemental_blast_family")) is not None else "unassessed"
            for prefix, typ in (("typei", "Type I"), ("typeii", "Type II"), ("typeiii", "Type III")):
                context = {k[len(prefix)+1:]: v for k, v in row.items() if k.startswith(prefix + "_")}
                if not context or all(_value(v) is None for v in context.values()):
                    continue
                self.raw["contexts"].append(dict(sample_id=sample, gene_id=gene, family=typ, source_ref=rr, native_fields=context))
                if family == typ and _value(context.get("context_status")):
                    record["context_state"] = context["context_status"]
                partners = _value(context.get("context_partners"))
                if not partners:
                    continue
                member_ids = [gene]
                unresolved = False
                known = [g for s, g in self.genes if s == sample]
                for partner in str(partners).split(" | "):
                    matches = [g for g in known if partner == g or partner.startswith(g + ":")]
                    if len(matches) != 1:
                        unresolved = True
                        break
                    member_ids.append(matches[0])
                if unresolved:
                    self.raw["unresolved_contexts"].append(dict(sample_id=sample, gene_id=gene,
                        family=typ, source_ref=rr, partner_ids=str(partners).split(" | "),
                        resolved_member_ids=list(dict.fromkeys(member_ids)),
                        failed_partner_id=partner, matching_gene_ids=matches,
                        reason="context_partner_requires_one_exact_canonical_identity"))
                    self.warn(f"{sample}/{gene}: saved {typ} context partner identity unresolved; group unavailable")
                    continue
                group_id = json.dumps([sample, "DNMB", ref["sha256"], number, prefix])
                self.data["groups"].append(dict(group_id=group_id, sample_id=sample, provider="DNMB", run_id=ref["sha256"],
                    membership_kind="saved_context", family=typ, context_state=_value(context.get("context_status")),
                    members=[self.member(sample, g) for g in dict.fromkeys(member_ids)], source_refs=[rr]))
        if ref:
            cached = detail["cached_table_counts"]
            # all_hits is the comparable all-candidate table in newer source versions.
            expected = cached.get("all_hits", cached.get("hits"))
            if expected is not None and len(rows) != expected:
                self.conflict(sample, module, "cached_live_candidate_count_mismatch", len(rows), expected, ref)
            sentinel = detail.get("sentinel") or {}
            if "curation_keep" in headers and sentinel.get("n_hits") is not None and retained_count != sentinel["n_hits"]:
                self.conflict(sample, module, "sentinel_live_retained_count_mismatch", retained_count, sentinel["n_hits"], ref)
        binding_status = "partial" if unbound_count else "complete" if ref is not None else "unassessed"
        detail.update(reference_pair_schema=reference_schema, candidate_binding_status=binding_status,
            candidate_row_count=len(rows), bound_candidate_row_count=len(rows)-unbound_count,
            unbound_candidate_row_count=unbound_count, unbound_reference_row_count=unbound_reference_count,
            incomplete_reference_row_count=incomplete_reference_count)
        self.raw["capabilities"].append(dict(sample_id=sample, saved_retention="curation_keep" in headers,
            native_candidates=ref is not None, candidate_count=len(rows), candidate_binding_status=binding_status,
            unbound_candidate_row_count=unbound_count, reference_scope="exported_candidates",
            retained_row_count=retained_count if "curation_keep" in headers else None))
        for name, destination in (("DNMB_REBASEfinder_motif_hits.tsv", "motifs"), ("DNMB_REBASEfinder_motif_hits_raw.tsv", "motifs"),
            ("DNMB_REBASEfinder_structure_coverage.tsv", "structures"), ("DNMB_REBASEfinder_motif_geometry.tsv", "structures"),
            ("DNMB_REBASEfinder_homology_templates.tsv", "structures"), ("foldseek_results.tsv", "structures")):
            extra = directory / "dnmb_module_rebasefinder" / name
            if not extra.is_file():
                continue
            er, eh, ef = self.table(extra, sample, module, name)
            idcol = next((c for c in ("locus_tag", "query", "query_id") if c in eh), None)
            if not idcol:
                self.warn(f"{sample}: {name} lacks supported exact ID column; preserved source only")
                continue
            for number, row in enumerate(er, 2):
                # Structure query labels can include filename suffixes; no guessing.
                value = row[idcol]
                if (sample, value) not in self.genes:
                    self.raw["unmapped_records"].append(dict(sample_id=sample, module=module, source_ref={**ef, "row": number}, native_fields=row, assessment_status="unassessed_unbound_identity"))
                    self.warn(f"{sample}: {name} row {number} has unbound structure identifier; raw observation quarantined")
                    continue
                self.raw[destination].append(dict(sample_id=sample, gene_id=value, evidence_type=name,
                    source_ref={**ef, "row": number}, native_fields=row))

    def eggnog(self, directory, sample):
        """Preserve native orthology annotations without inferring gene roles."""
        path = directory / "dnmb_module_eggnog" / "eggnog_out.emapper.annotations"
        if not path.is_file():
            self.raw["capabilities"].append(dict(sample_id=sample, raw_eggnog_status="missing"))
            return
        ref = self.source(path, sample, "eggnog", "functional_annotations")
        headers, comments, seen = None, [], {}
        row_count = bound_count = unbound_count = conflict_count = 0
        with path.open(encoding="utf-8-sig", newline="") as stream:
            for number, line in enumerate(stream, 1):
                text = line.rstrip("\r\n")
                if not text.strip():
                    continue
                cells = text.split("\t")
                if cells[0] in {"#query", "query"}:
                    if headers is not None:
                        raise ValueError(f"{path}:{number}: repeated EggNOG header")
                    self.reader.headers(cells, path)
                    if "seed_ortholog" not in cells:
                        raise ValueError(f"{path}:{number}: EggNOG header lacks seed_ortholog")
                    headers = cells
                    continue
                if text.startswith("#"):
                    comments.append(dict(row=number, text=text))
                    continue
                if headers is None:
                    raise ValueError(f"{path}:{number}: EggNOG data precedes a query header")
                if len(cells) != len(headers):
                    raise ValueError(f"{path}:{number}: EggNOG row/header column count differs")
                row = dict(zip(headers, cells))
                gene = cells[0]
                if _value(row.get("sample_id")) is not None and row["sample_id"] != sample:
                    raise ValueError(f"{path}:{number}: row sample_id conflicts with the declared sample unit")
                row_count += 1
                rr = {**ref, "row": number}
                if (sample, gene) not in self.genes:
                    unbound_count += 1
                    self.raw["unmapped_records"].append(dict(sample_id=sample, module="eggnog", namespace="EggNOG",
                        source_ref=rr, native_fields=row, assessment_status="unassessed_unbound_identity"))
                    self.conflict(sample, "eggnog", "unbound_eggnog_query", gene, "exact existing gene ID", rr)
                    continue
                if gene in seen and seen[gene] != row:
                    conflict_count += 1
                    self.conflict(sample, "eggnog", "conflicting_eggnog_query", gene, "one coherent annotation row per query", rr)
                seen.setdefault(gene, row)
                bound_count += 1
                self.raw["annotations"].append(dict(sample_id=sample, gene_id=gene, namespace="EggNOG",
                    source_ref=rr, native_fields=row))
        if headers is None:
            raise ValueError(f"{path}: EggNOG query header is missing")
        self.pin(path)
        item = next(x for x in self.raw["sources"] if x["path"] == ref["path"])
        item.update(row_count=row_count, columns=headers, comments=comments)
        self.raw["capabilities"].append(dict(sample_id=sample, raw_eggnog_scope="all_exact_bound_queries",
            raw_eggnog_status="partial" if unbound_count or conflict_count else "complete", raw_eggnog_rows=row_count,
            raw_eggnog_bound_rows=bound_count, raw_eggnog_unbound_rows=unbound_count, raw_eggnog_conflicting_rows=conflict_count))

    def domains(self, directory, sample):
        self.eggnog(directory, sample)
        candidates = _system_review_gene_scope(self.data, sample)
        found = False
        for name in INTERPRO_FIELDS:
            choices = [directory / folder / name for folder in ("dnmb_module_interproscan", "dnmb_interproscan")]
            existing = [p for p in choices if p.is_file()]
            if len(existing) > 1:
                raise ValueError(f"{directory}: ambiguous InterPro exports in two native layouts")
            if not existing:
                continue
            found = True
            path = existing[0]
            ref = self.source(path, sample, "interproscan", "domain_sites" if name.endswith("sites") else "domain_hits")
            _attach_interpro_file(self.raw, self.genes, sample, candidates, path, ref, strict=self.strict)
            self.pin(path)
        capability = dict(sample_id=sample,
            raw_domain_scope="all_observed_calls_DNMB_candidates_and_explicit_context_members")
        if not found:
            capability["raw_interpro_status"] = "missing"
        self.raw["capabilities"].append(capability)


INTERPRO_FIELDS = {
    "interproscan_results.tsv": ("protein_id", "sequence_md5", "protein_length", "analysis", "signature_accession", "signature_description",
        "start", "end", "score", "status", "date", "interpro_accession", "interpro_description", "go_annotations", "pathway_annotations"),
    "interproscan_results.tsv.sites": ("protein_id", "sequence_md5", "protein_length", "analysis", "signature_accession", "start", "end",
        "site_count", "residue", "site_start", "site_end", "site_description"),
}


def _system_review_gene_scope(data, sample):
    """Observed identities only; membership is not a functional assignment."""
    result = {r["gene_id"] for table in ("calls", "saved_decisions", "annotations", "effector_evidence")
              for r in data.get(table, []) if r.get("sample_id") == sample and r.get("gene_id")}
    result.update(m["gene_id"] for g in data.get("groups", []) if g.get("sample_id") == sample
                  for m in g.get("members", []) if m.get("gene_id"))
    return result


def _file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _attach_interpro_file(raw, genes, sample, scope, path, ref, *, strict=False):
    """Read a pinned observation source; quarantine contradictions without role inference."""
    path = Path(path)
    if _file_sha256(path) != ref.get("sha256"):
        raise ValueError(f"{path}: InterPro source changed from its pinned hash")
    fields = INTERPRO_FIELDS[path.name]
    profiles, quarantine, conflicts = [], [], []
    parsed = excluded = bound = sequence_unassessed = 0
    scoped_genes = set()
    with path.open(encoding="utf-8-sig", newline="") as stream:
        for number, row in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            parsed += 1
            rr = {**ref, "row": number}
            if len(row) != len(fields):
                raise ValueError(f"{path}:{number}: unsupported InterPro column count {len(row)}")
            native = dict(zip(fields, row))
            gene = genes.get((sample, row[0]))
            reason = None
            state = "unassessed"
            if gene is None:
                reason = "unbound_interpro_query"
            else:
                bound += 1
                try:
                    length = int(native["protein_length"])
                    start, end = int(native["start"]), int(native["end"])
                    valid = length > 0 and 1 <= start <= end <= length
                    if path.name.endswith(".sites"):
                        valid = valid and int(native["site_count"]) >= 0 and 1 <= int(native["site_start"]) <= int(native["site_end"]) <= length
                    if not valid:
                        reason = "invalid_interpro_coordinates"
                except (ValueError, TypeError):
                    reason = "invalid_interpro_coordinates"
                reported_md5 = _value(native["sequence_md5"])
                if reason is None and gene.get("protein_length") is not None and length != gene["protein_length"]:
                    reason = "interpro_sequence_length_mismatch"
                if reason is None and reported_md5 is not None:
                    if not re.fullmatch(r"[0-9a-fA-F]{32}", str(reported_md5)):
                        reason = "invalid_interpro_sequence_md5"
                    elif gene.get("protein_md5"):
                        if str(reported_md5).lower() != gene["protein_md5"]:
                            reason = "interpro_sequence_md5_mismatch"
                        else:
                            state = "verified"
            observation = dict(sample_id=sample, gene_id=row[0], namespace="InterProScan", source_ref=rr,
                native_fields=native, evidence_type="site_annotation" if path.name.endswith(".sites") else "domain_hit")
            if gene and gene.get("protein_sha256"):
                observation["protein_sha256"] = gene["protein_sha256"]
            if reason:
                observation.update(module="interproscan", assessment_status="unassessed_quarantined",
                                   sequence_binding="conflicting" if gene else "unbound", reason=reason)
                quarantine.append(observation)
                conflicts.append(dict(sample_id=sample, module="interproscan", code=reason, severity="warning",
                    observed=dict(gene_id=row[0], sequence_md5=native["sequence_md5"], protein_length=native["protein_length"]),
                    expected="exact current gene, sequence and valid protein coordinates",
                    action="raw_observation_quarantined", source_ref=rr))
            elif row[0] not in scope:
                excluded += 1
            else:
                observation.update(sequence_binding=state,
                    coordinate_binding="verified" if gene.get("protein_length") is not None else "reported_sequence_only",
                    qualification="unassessed", activity="unassessed")
                sequence_unassessed += state != "verified"
                profiles.append(observation)
                scoped_genes.add(row[0])
    if _file_sha256(path) != ref["sha256"]:
        raise ValueError(f"{path}: InterPro source changed while reading")
    if strict and conflicts:
        raise ValueError(f"{sample}/interproscan: {conflicts[0]['code']}; strict input consistency failed")
    def same_source(r):
        source = r.get("source_ref") or {}
        return r.get("sample_id") == sample and source.get("sha256") == ref["sha256"] and source.get("path") == ref["path"]
    for table, new in (("profiles", profiles), ("unmapped_records", quarantine), ("conflicts", conflicts)):
        raw[table] = [r for r in raw.get(table, []) if not (same_source(r) and
            (r.get("namespace") == "InterProScan" or r.get("module") == "interproscan"))] + new
    capability = dict(sample_id=sample, source_ref=ref,
        raw_interpro_scope="all_observed_calls_DNMB_candidates_and_explicit_context_members",
        raw_interpro_status="partial" if quarantine or sequence_unassessed else "complete",
        raw_interpro_parsed_rows=parsed, raw_interpro_bound_rows=bound,
        raw_interpro_retained_rows=len(profiles), raw_interpro_retained_genes=len(scoped_genes),
        raw_interpro_excluded_scope_rows=excluded, raw_interpro_excluded_reason="outside_declared_observed_identity_scope",
        raw_interpro_quarantined_rows=len(quarantine), raw_interpro_sequence_unassessed_rows=sequence_unassessed)
    assert parsed == len(profiles) + excluded + len(quarantine)
    raw["capabilities"] = [c for c in raw.get("capabilities", []) if not ("raw_interpro_scope" in c and same_source(c))] + [capability]
    warning_prefix = f"{sample}/interproscan [{ref['path']}]: "
    raw["warnings"] = [w for w in raw.get("warnings", []) if not str(w).startswith(warning_prefix)]
    if quarantine:
        raw["warnings"].append(warning_prefix + f"{len(quarantine)} raw observations quarantined; exact binding or coordinates conflict")
    for source in raw.get("sources", []):
        if source.get("sample_id") == sample and source.get("path") == ref["path"] and source.get("sha256") == ref["sha256"]:
            source.update(row_count=parsed, columns=list(fields), evidence_schema="InterProScan native TSV",
                          retained_row_count=len(profiles), excluded_scope_row_count=excluded, quarantined_row_count=len(quarantine))


def attach_system_review_annotations(data, root, *, strict_inputs=False):
    """Expand pinned InterPro observations after supplements; mutate and return data.

    Only already imported InterPro files are reopened. Current exact gene IDs,
    source hashes and protein coordinates are checked; no role is inferred.
    ``root`` resolves relative paths in serialized source inventories.
    """
    root = Path(root).expanduser().resolve()
    raw = data.get("suite", {})
    genes = {(g["sample_id"], g["gene_id"]): g for g in data.get("genes", [])}
    samples = {s["sample_id"] for s in data.get("samples", [])}
    for source in raw.get("sources", []):
        if source.get("module") != "interproscan" or source.get("role") not in {"domain_hits", "domain_sites"}:
            continue
        sample = source["sample_id"]
        if sample not in samples:
            raise ValueError("InterPro source belongs to an unknown sample")
        path = Path(source.get("resolved_path") or source["path"])
        path = path if path.is_absolute() else root / path
        ref = {k:source[k] for k in ("path", "sha256")}
        _attach_interpro_file(raw, genes, sample, _system_review_gene_scope(data, sample), path, ref,
                             strict=strict_inputs)
    return data


def load_suite(root, *, manifest=None, rscript=None, strict_inputs=False):
    """Load one native folder or an immediate-child cohort without preprocessing."""
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"{root}: suite root must be an existing directory")
    manifest = Path(manifest).expanduser().resolve() if manifest is not None else None
    r = SuiteReader(root, manifest, rscript, strict_inputs)
    config = r.reader.json(manifest) if manifest else {}
    if manifest:
        if not isinstance(config, dict) or type(config.get("schema_version")) is not int or config["schema_version"] != 1:
            raise ValueError("suite manifest schema_version must be integer 1")
        unknown = set(config) - {"schema_version", "samples", "family_map", "options"}
        if unknown:
            raise ValueError(f"unknown suite manifest keys: {sorted(unknown)}")
        specs = config.get("samples")
        if not isinstance(specs, list) or not specs:
            raise ValueError("suite manifest samples must be a nonempty list")
    else:
        def native(p):
            return any(x.is_dir() and x.name.startswith("dnmb_module_") for x in p.iterdir())
        excluded = {"staging", "receipts", "output", "outputs", "reviews", "artifacts"}
        dirs = [root] if native(root) else sorted(p for p in root.iterdir() if p.is_dir() and not p.name.startswith(".") and p.name.lower() not in excluded and native(p))
        if not dirs:
            raise ValueError(f"{root}: no native suite folders discovered")
        specs = [{"sample_id": p.name, "path": str(p), "label": p.name} for p in dirs]
    families = _family_maps(r.reader, config)
    catalog_path = Path(str(files("defenseviz").joinpath("data/families.json"))).resolve()
    if str(catalog_path) in r.reader.hashes:
        r.raw["software_resources"].append(dict(resource="package:defenseviz/data/families.json", sha256=r.reader.hashes.pop(str(catalog_path))))
        r.reader.inputs.remove(catalog_path)
    seen = set()
    directories = set()
    for spec in specs:
        if not isinstance(spec, dict) or set(spec) - {"sample_id", "path", "label", "highlight", "genus", "species", "cohort", "display", "modules", "genes", "fasta"}:
            raise ValueError("invalid or unknown suite sample fields")
        sample = _string(spec.get("sample_id"), "suite sample_id")
        if sample in seen:
            raise ValueError(f"duplicate suite sample_id {sample!r}")
        seen.add(sample)
        directory = r.reader.path(spec.get("path"))
        if not directory.is_dir():
            raise ValueError(f"{directory}: sample directory missing")
        if directory in directories:
            raise ValueError(f"{directory}: sample directory bound more than once")
        directories.add(directory)
        metadata_sample = {k: spec[k] for k in ("genus", "species", "cohort", "display") if k in spec}
        if "display" in spec and not isinstance(spec["display"],dict):
            raise ValueError("sample display must be a metadata object")
        r.data["samples"].append(dict(sample_id=sample, label=_string(spec.get("label") or sample, "sample label"),
            highlight=_bool(spec.get("highlight", False), "sample highlight") or False, **metadata_sample))
        modules = spec.get("modules", {})
        if not isinstance(modules, dict) or set(modules) - set(MODULES):
            raise ValueError("suite modules must use defensefinder,padloc,defensepredictor,rebasefinder")
        for name, override in modules.items():
            if not isinstance(override, dict) or set(override) - {"calls", "systems", "id_map", "profiles", "fasta", "workbook", "status", "score_threshold", "score_operator", "decision_metadata"}:
                raise ValueError(f"{name}: invalid module override")
            if name != "defensepredictor" and set(override) & {"score_operator", "decision_metadata"}:
                raise ValueError(f"{name}: score_operator and decision_metadata are DefensePredictor-only")
        r.gene_annotations(directory, spec, sample)
        metadata = r.metadata(directory, sample)
        for module in ("defensefinder", "padloc", "defensepredictor"):
            r.comparator(directory, spec, module, sample, metadata, families)
        # Canonical locus-labelled suite FAA supplies hashes even without gene table.
        if spec.get("fasta"):
            fasta = r.reader.path(spec["fasta"])
        else:
            candidates = sorted(directory.glob("*.faa"))
            if len(candidates) > 1:
                raise ValueError(f"{directory}: ambiguous root protein FASTA; supply manifest fasta")
            fasta = candidates[0] if candidates else None
        r.fasta(fasta, sample)
        r.rebase(directory, spec, sample, metadata)
        r.domains(directory, sample)
        r.data["source_metadata"].append(dict(sample_id=sample, provider="DNMB", source_kind="saved_native_import",
            source_directory=r.reader.label(directory), protein_hash_method="sha256 of uppercase sequence with whitespace removed; terminal stop preserved"))
    # Attach hashes loaded later from root FASTA to already imported memberships.
    for collection in (r.data["calls"], r.data["saved_decisions"]):
        for row in collection:
            digest = r.genes[(row["sample_id"], row["gene_id"])].get("protein_sha256")
            if digest:
                row["protein_sha256"] = digest
                row["protein_binding_origin"] = "importer_attached"
    for group in r.data["groups"]:
        for member in group["members"]:
            digest = r.genes[(group["sample_id"], member["gene_id"])].get("protein_sha256")
            if digest:
                member["protein_sha256"] = digest
    r.data.update(genes=list(r.genes.values()), inputs=r.reader.inputs, input_sha256=r.reader.hashes,
                  options=config.get("options", {}), evidence_extension=None)
    if not isinstance(r.data["options"], dict):
        raise ValueError("suite options must be an object")
    if manifest:
        declaration_ref = r.pin(manifest)
        declarations = config["options"] if "options" in config else {}
        global_scope = declarations.get("assessment_coverage", {})
        sample_scopes = declarations.get("assessment_coverage_by_sample", {})
        for sample in r.data["samples"]:
            sid = sample["sample_id"]
            for scope, axes, option_path in (("global", global_scope, "options.assessment_coverage"),
                    ("sample", sample_scopes.get(sid, {}) if isinstance(sample_scopes, dict) else {},
                     f"options.assessment_coverage_by_sample[{sid!r}]")):
                if isinstance(axes, dict) and axes:
                    r.data["source_metadata"].append(dict(sample_id=sid, provider="DNMB",
                        source_kind="assessment_coverage_declaration", declaration_scope=scope,
                        axes=dict(axes), option_path=option_path, source_ref=declaration_ref))
    return r.data

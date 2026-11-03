"""Read saved calls without guessing sample identities or biological activity.

The adapters accept concrete table schemas, preserve native system membership,
and keep ID translation separate from feature interpretation.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
from importlib.resources import files
from pathlib import Path
from typing import Any

CATEGORIES = frozenset({"R-M", "Modification", "Nuclease", "CRISPR-Cas",
                        "Abortive infection", "Other mechanisms", "Unresolved"})
KINDS = frozenset({"final_member", "hmm_hit", "reference_hit"})
ROLES = frozenset({"R", "M", "S", "RM", "C", "Aux"})
STATUSES = frozenset({"complete", "missing", "not_run"})
# Exact native type labels carry classification; substring/keyword inference is forbidden.
RM_TYPES = {**{f"RM_Type_{t}": f"Type {t}" for t in ("I", "II", "III", "IV")},
            **{f"RM_type_{t}": f"Type {t}" for t in ("I", "II", "III", "IV")}}
DEFAULT_CALLERS = {"defensefinder": "DefenseFinder", "padloc": "PADLOC",
                   "defensepredictor": "DefensePredictor"}


def _error(where: str, message: str) -> ValueError:
    return ValueError(f"{where}: {message}")


def _string(value: Any, where: str, *, optional: bool = False) -> str | None:
    if value is None or value == "":
        if optional:
            return None
        raise _error(where, "required nonempty text is missing")
    if not isinstance(value, str) or not value.strip():
        raise _error(where, "expected nonempty text (identifiers must be strings)")
    if value != value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise _error(where, "text must not contain surrounding whitespace, NUL, or line/tab separators")
    return value


def _boolean(value: Any, where: str) -> bool:
    if isinstance(value, bool):
        return value
    if value in (None, ""):
        return False
    if isinstance(value, str) and value.lower() in {"true", "false", "1", "0"}:
        return value.lower() in {"true", "1"}
    raise _error(where, "expected true/false or 1/0")


def _integer(value: Any, where: str) -> int:
    if isinstance(value, bool):
        raise _error(where, "expected a positive integer")
    try:
        parsed = int(value)
    except (ValueError, TypeError):
        raise _error(where, "expected a positive integer") from None
    if str(parsed) != str(value) or parsed < 1:
        raise _error(where, "expected a positive integer")
    return parsed


def _sequence_hash(value: Any, where: str) -> str:
    value = _string(value, where)
    if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise _error(where, "sequence_sha256 must be 64 lowercase hexadecimal characters")
    return value


def _unique_object(pairs):
    value = {}
    for key, val in pairs:
        if key in value:
            raise ValueError(f"duplicate JSON key {key!r}")
        value[key] = val
    return value


def _nonfinite(token):
    raise ValueError(f"non-finite JSON number {token!r} is not permitted")


class Reader:
    def __init__(self, config: Path):
        self.base = config.parent
        self.inputs: list[Path] = []
        self.hashes: dict[str, str] = {}

    def path(self, value: Any) -> Path:
        raw = _string(value, "file path")
        p = Path(raw).expanduser()
        return (p if p.is_absolute() else self.base / p).resolve()

    def label(self, p: Path) -> str:
        """Portable provenance reference; resolved paths stay in the receipt only."""
        return Path(os.path.relpath(p, self.base)).as_posix()

    def read(self, p: Path) -> bytes:
        if not p.is_file():
            raise _error(str(p), "input file does not exist or is not a regular file")
        raw = p.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if str(p) in self.hashes and self.hashes[str(p)] != digest:
            raise _error(str(p), "input changed while loading project")
        self.hashes[str(p)] = digest
        if p not in self.inputs:
            self.inputs.append(p)
        return raw

    def json(self, p: Path):
        raw = self.read(p)
        try:
            return json.loads(raw.decode("utf-8-sig"),
                              object_pairs_hook=_unique_object, parse_constant=_nonfinite)
        except (UnicodeError, ValueError) as exc:
            raise _error(str(p), f"invalid JSON: {exc}") from exc

    def table(self, p: Path, *, sheet: str | None = None):
        raw = self.read(p)
        if p.suffix.lower() == ".json":
            values = self.json(p)
            if not isinstance(values, list) or any(not isinstance(x, dict) for x in values):
                raise _error(str(p), "table JSON must be a list of row objects")
            headers = set().union(*(x.keys() for x in values)) if values else set()
            return [{k: row.get(k) for k in headers} for row in values], headers
        if p.suffix.lower() == ".xlsx":
            try:
                import openpyxl
            except ImportError as exc:
                raise _error(str(p), "XLSX support requires pip install 'defenseviz[excel]'") from exc
            workbook = openpyxl.load_workbook(io.BytesIO(raw), read_only=True, data_only=False)
            try:
                if sheet is None and len(workbook.sheetnames) != 1:
                    raise _error(str(p), "multiple worksheets: specify sheet explicitly")
                if sheet is not None and sheet not in workbook.sheetnames:
                    raise _error(str(p), f"worksheet {sheet!r} does not exist")
                ws = workbook[sheet or workbook.sheetnames[0]]
                data = []
                for row in ws.iter_rows():
                    if any(c.data_type == "f" for c in row):
                        raise _error(str(p), "formulas are not imported; export literal values first")
                    data.append([c.value for c in row])
            finally:
                workbook.close()
            if not data:
                raise _error(str(p), "empty workbook lacks a header")
            headers = data[0]
            self.headers(headers, p)
            rows = []
            for i, row in enumerate(data[1:], 2):
                if all(x is None for x in row):
                    continue
                if len(row) != len(headers):
                    raise _error(f"{p}:{i}", "row/header column count differs")
                rows.append(dict(zip(headers, ("" if x is None else x for x in row))))
            return rows, set(headers)
        if p.suffix.lower() not in {".tsv", ".csv", ".txt"}:
            raise _error(str(p), "supported tables are .tsv, .csv, .txt (tab-separated), .json, or .xlsx")
        try:
            content = raw.decode("utf-8-sig")
        except UnicodeError as exc:
            raise _error(str(p), "table must use UTF-8 encoding") from exc
        if "\x00" in content:
            raise _error(str(p), "NUL byte in table")
        lines = content.splitlines(keepends=True)
        # DefenseFinder's systems files carry a leading # metadata preamble.
        while lines and (not lines[0].strip() or lines[0].startswith("# ")):
            lines.pop(0)
        if not lines:
            raise _error(str(p), "empty file lacks a header; use a header-only file for a completed zero-call run")
        delimiter = "," if p.suffix.lower() == ".csv" else "\t"
        try:
            reader = csv.DictReader(io.StringIO("".join(lines)), delimiter=delimiter, strict=True)
            headers = reader.fieldnames
            self.headers(headers, p)
            rows = []
            for i, row in enumerate(reader, 2):
                if None in row or any(v is None for v in row.values()):
                    raise _error(f"{p}:{i}", "row/header column count differs")
                rows.append(row)
            return rows, set(headers)
        except csv.Error as exc:
            raise _error(str(p), f"malformed delimited table: {exc}") from exc

    @staticmethod
    def headers(headers, p):
        if not headers or any(not isinstance(h, str) or not h.strip() for h in headers):
            raise _error(str(p), "all column headers must be nonempty strings")
        if len(headers) != len(set(headers)):
            raise _error(str(p), "duplicate column header")
        if any(h != h.strip() for h in headers):
            raise _error(str(p), "column headers must not have surrounding whitespace")


def _require(headers, names, where):
    missing = set(names) - headers
    if missing:
        raise _error(str(where), "missing required columns: " + ", ".join(sorted(missing)))


def _descriptor(value, where):
    if isinstance(value, str):
        return {"path": value}
    if isinstance(value, dict) and "path" in value:
        return value
    raise _error(where, "expected a filename or an object containing path")


def _bound_sample(row, spec, known, where):
    supplied = spec.get("sample_id")
    observed = row.get("sample_id")
    if supplied and observed and supplied != observed:
        raise _error(where, "row sample_id conflicts with source sample_id")
    sample = _string(supplied or observed, where + " sample_id")
    if sample not in known:
        raise _error(where, f"unknown sample_id {sample!r}")
    return sample


def _maps(reader, spec, known):
    if not spec.get("id_map"):
        return None
    path = reader.path(spec["id_map"])
    rows, headers = reader.table(path)
    # Exported maps may use query_id/locus_tag, but never parse a protein label.
    source = "source_id" if "source_id" in headers else "query_id"
    dest = "gene_id" if "gene_id" in headers else "locus_tag"
    _require(headers, [source, dest], path)
    mappings, destinations = {}, set()
    for i, row in enumerate(rows, 2):
        where = f"{path}:{i}"
        sample = _bound_sample(row, spec, known, where)
        a = _string(row[source], where + " source_id")
        b = _string(row[dest], where + " gene_id")
        key = (sample, a)
        if key in mappings or (sample, b) in destinations:
            raise _error(where, "ID map must be one-to-one within each sample; duplicate source or destination")
        mappings[key] = b
        destinations.add((sample, b))
    return mappings


def _translate(gene, sample, mapping, where):
    gene = _string(gene, where + " gene_id")
    if mapping is None:
        return gene
    if (sample, gene) not in mapping:
        raise _error(where, f"ID map has no entry for sample {sample!r}, source ID {gene!r}")
    return mapping[sample, gene]


def _family_maps(reader, config):
    mappings = {}
    # The optional shipped catalog is declarative classification metadata, never
    # a biological-role or scoring rule. Explicit user mappings take precedence.
    resource = files("defenseviz").joinpath("data/families.json")
    if resource.is_file():
        p = Path(str(resource)).resolve()
        defaults = reader.json(p)
        if not isinstance(defaults, list) or any(not isinstance(r, dict) for r in defaults):
            raise _error(str(p), "default family catalog must be a list of mapping objects")
        groups = [(p, defaults)]
    else:
        groups = []
    if config.get("family_map"):
        desc = _descriptor(config["family_map"], "family_map")
        p = reader.path(desc["path"])
        rows, headers = reader.table(p, sheet=desc.get("sheet"))
        _require(headers, ["family", "category"], p)
        groups.append((p, rows))
    for priority, (p, rows) in enumerate(groups):
        mapping_source = "builtin_catalog" if p.name == "families.json" and p == Path(str(resource)).resolve() else "family_map"
        seen = set()
        for i, row in enumerate(rows, 2):
            where = f"{p}:{i}"
            family = _string(row.get("family"), where + " family")
            caller = _string(row.get("caller"), where + " caller", optional=True)
            category = _string(row.get("category"), where + " category")
            if category not in CATEGORIES:
                raise _error(where, f"unknown category {category!r}")
            key = (caller, family)
            if key in seen:
                raise _error(where, "duplicate family/caller mapping")
            seen.add(key)
            mappings[key] = dict(category=category, priority=priority, source=mapping_source, caller=caller,
                                 display_family=_string(row.get("display_family"), where, optional=True),
                                 rm_type=_string(row.get("rm_type"), where, optional=True))
    return mappings


def _canonical(raw, family_map, where):
    for field in ("sample_id", "gene_id", "caller", "family"):
        raw[field] = _string(raw.get(field), where + " " + field)
    raw.setdefault("native_family", raw["family"])
    candidates = [family_map[key] for key in ((raw["caller"], raw["family"]), (None, raw["family"])) if key in family_map]
    mapping = max(candidates, key=lambda m: (m["priority"], m["caller"] is not None)) if candidates else None
    explicit_category = raw.pop("_category_explicit", False)
    if explicit_category and mapping and mapping["source"] == "builtin_catalog":
        mapping = None  # User-supplied normalized classification outranks shipped defaults.
    raw["category_source"] = "native_type" if raw.get("native_type") and raw.get("category") == "R-M" else "explicit" if raw.get("category") else "unresolved"
    category = _string(raw.get("category"), where + " category", optional=True) or "Unresolved"
    if mapping:
        explicit = raw.get("category")
        if explicit and explicit != "Unresolved" and explicit != mapping["category"]:
            raise _error(where, "explicit category conflicts with family_map")
        category = mapping["category"]
        raw["category_source"] = mapping["source"]
        if mapping["display_family"]:
            raw["native_family"] = raw["family"]
            raw["family"] = mapping["display_family"]
        if mapping["rm_type"]:
            if raw.get("rm_type") and raw["rm_type"] != mapping["rm_type"]:
                raise _error(where, "native rm_type conflicts with family_map")
            raw["rm_type"] = mapping["rm_type"]
    if category not in CATEGORIES:
        raise _error(where, f"unknown category {category!r}; supply an exact family_map")
    raw["category"] = category
    raw["evidence_kind"] = _string(raw.get("evidence_kind"), where + " evidence_kind", optional=True) or "final_member"
    if raw["evidence_kind"] not in KINDS:
        raise _error(where, "evidence_kind must be final_member, hmm_hit, or reference_hit")
    for field in ("system_id", "role", "rm_type"):
        raw[field] = _string(raw.get(field), where + " " + field, optional=True)
    if raw["role"] is not None and raw["role"] not in ROLES:
        raise _error(where, f"unsupported role {raw['role']!r}")
    if raw["rm_type"] is not None and raw["rm_type"] not in {"Type I", "Type II", "Type III", "Type IV"}:
        raise _error(where, f"unsupported rm_type {raw['rm_type']!r}")
    if raw["rm_type"] and category != "R-M":
        raise _error(where, "rm_type requires category R-M")
    return raw


def _calls(reader, spec, known, family_map):
    if not isinstance(spec, dict):
        raise _error("calls", "each source must be an object")
    accepted_keys = {"format", "path", "sample_id", "caller", "id_map", "sheet", "score_threshold", "score_column", "id_column"}
    if set(spec) - accepted_keys:
        raise _error("calls", "unknown source keys: " + ", ".join(sorted(set(spec) - accepted_keys)))
    for key in ("sample_id", "caller", "id_map", "sheet", "score_column", "id_column"):
        if key in spec:
            _string(spec[key], "calls source " + key)
    fmt = spec.get("format")
    if fmt not in {"normalized", "defensefinder", "padloc", "defensepredictor"}:
        raise _error("calls", f"unsupported format {fmt!r}")
    if fmt != "defensepredictor" and set(spec) & {"score_threshold", "score_column", "id_column"}:
        raise _error("calls", "score_threshold, score_column and id_column are DefensePredictor-only source options")
    p = reader.path(spec.get("path"))
    rows, headers = reader.table(p, sheet=spec.get("sheet"))
    mapping = _maps(reader, spec, known)
    fixed_caller = spec.get("caller") or DEFAULT_CALLERS.get(fmt)
    if fixed_caller is not None:
        fixed_caller = _string(fixed_caller, str(p) + " caller")
    if spec.get("sample_id"):
        _bound_sample({}, spec, known, str(p))
    if fmt == "normalized":
        _require(headers, ["gene_id", "family"], p)
        if not spec.get("sample_id"):
            _require(headers, ["sample_id"], p)
        if not fixed_caller:
            _require(headers, ["caller"], p)
    elif fmt == "defensefinder":
        is_systems = "protein_in_syst" in headers
        _require(headers, ["sys_id", "type", "protein_in_syst" if is_systems else "hit_id"], p)
    elif fmt == "padloc":
        _require(headers, ["system.number", "seqid", "system", "target.name"], p)
    else:
        score_column = spec.get("score_column", "mean_log_odds")
        id_column = spec.get("id_column", "locus_tag")
        _require(headers, [score_column, id_column], p)
        threshold = spec.get("score_threshold")
        if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(threshold):
            raise _error(str(p), "DefensePredictor requires an explicit finite numeric score_threshold")
    calls, observed_callers, prediction_ids = [], set(), set()
    for i, row in enumerate(rows, 2):
        where = f"{p}:{i}"
        sample = _bound_sample(row, spec, known, where)
        if row.get("caller") and fixed_caller and row["caller"] != fixed_caller:
            raise _error(where, "row caller conflicts with source caller")
        caller = fixed_caller or _string(row.get("caller"), where + " caller")
        observed_callers.add(caller)
        base = {"sample_id": sample, "caller": caller, "source": reader.label(p)}
        if fmt == "normalized":
            base.update({k: row.get(k) for k in ("family", "category", "system_id", "role", "rm_type", "evidence_kind")})
            base["_category_explicit"] = row.get("category") not in (None, "")
            ids = [row["gene_id"]]
        elif fmt == "defensefinder":
            native_type = _string(row["type"], where + " type")
            base.update(family=row.get("subtype") or native_type, system_id=_string(row["sys_id"], where + " sys_id"),
                        native_type=native_type, native_model=row.get("gene_name"),
                        rm_type=RM_TYPES.get(native_type), category="R-M" if native_type in RM_TYPES else "Unresolved")
            if is_systems:
                value = _string(row["protein_in_syst"], where + " protein_in_syst")
                ids = [x.strip() for x in value.split(",")]
                if not all(ids) or len(ids) != len(set(ids)):
                    raise _error(where, "system member list contains empty or duplicate protein IDs")
                if row.get("genes_count") and _integer(row["genes_count"], where + " genes_count") != len(ids):
                    raise _error(where, "genes_count disagrees with explicit member list")
            else:
                ids = [row["hit_id"]]
        elif fmt == "padloc":
            family = _string(row["system"], where + " system")
            # Tuple encoding avoids collisions between local system numbers on different contigs.
            system = json.dumps([_string(row["seqid"], where + " seqid"),
                                 _string(row["system.number"], where + " system.number")], separators=(",", ":"))
            base.update(family=family, system_id=system, native_model=row.get("hmm.name"),
                        rm_type=RM_TYPES.get(family), category="R-M" if family in RM_TYPES else "Unresolved")
            ids = [row["target.name"]]
        else:
            try:
                score = float(row[score_column])
            except (ValueError, TypeError):
                raise _error(where, f"{score_column} must be a finite number") from None
            if not math.isfinite(score):
                raise _error(where, f"{score_column} must be a finite number")
            # Validate even unselected identifiers; a broken row is never silently skipped.
            canonical_id = _translate(row[id_column], sample, mapping, where)
            if (sample, canonical_id) in prediction_ids:
                raise _error(where, "duplicate DefensePredictor sample/gene prediction")
            prediction_ids.add((sample, canonical_id))
            if score < threshold:
                continue
            base.update(family="General defense", category="Unresolved", system_id=None,
                        native_score=score, score_threshold=threshold, score_column=score_column)
            ids = [row[id_column]]
        for gene in ids:
            raw = dict(base, gene_id=_translate(gene, sample, mapping, where))
            calls.append(_canonical(raw, family_map, where))
    inferred = []
    if spec.get("sample_id"):
        # Caller must still be identifiable for a header-only completed run.
        if fmt == "normalized" and rows:
            callers = {c["caller"] for c in calls if c["evidence_kind"] == "final_member"}
        else:
            callers = {fixed_caller} if fixed_caller else observed_callers
        if not callers and not rows:
            raise _error(str(p), "empty normalized source needs caller in its descriptor to identify completed coverage")
        inferred = [{"sample_id": spec["sample_id"], "caller": c, "status": "complete"} for c in sorted(callers)]
    return calls, inferred


def _evidence(reader, value, known):
    if value is None:
        return []
    spec = _descriptor(value, "evidence")
    p = reader.path(spec["path"])
    if spec.get("format") == "dnmb_features":
        return _feature_evidence(reader, spec, known)
    rows = reader.json(p)
    if not isinstance(rows, list) or any(not isinstance(r, dict) for r in rows):
        raise _error(str(p), "evidence JSON must be a list of record objects")
    fmt = spec.get("format", "normalized")
    bindings = None
    if fmt == "dnmb":
        if not spec.get("id_map"):
            raise _error(str(p), "DNMB target_id export requires id_map with target_id,sample_id,gene_id")
        mp = reader.path(spec["id_map"])
        mr, mh = reader.table(mp)
        _require(mh, ["target_id", "sample_id", "gene_id"], mp)
        bindings, genes = {}, set()
        for i, r in enumerate(mr, 2):
            where = f"{mp}:{i}"
            sample = _bound_sample(r, {}, known, where)
            target = _string(r["target_id"], where + " target_id")
            gene = _string(r["gene_id"], where + " gene_id")
            if target in bindings or (sample, gene) in genes:
                raise _error(where, "DNMB evidence ID map must be one-to-one")
            supplied_hash = r.get("sequence_sha256")
            if supplied_hash not in (None, ""):
                supplied_hash = _sequence_hash(supplied_hash, where + " sequence_sha256")
            else:
                supplied_hash = None
            bindings[target] = (sample, gene, supplied_hash)
            genes.add((sample, gene))
    elif fmt != "normalized":
        raise _error(str(p), f"unsupported evidence format {fmt!r}")
    result, seen = [], set()
    for i, row in enumerate(rows, 1):
        where = f"{p}:record {i}"
        if bindings is not None:
            target = _string(row.get("target_id"), where + " target_id")
            if target not in bindings:
                raise _error(where, f"no explicit ID binding for target {target!r}")
            sample, gene, mapped_hash = bindings[target]
            row_hash = row.get("sequence_sha256")
            if row_hash not in (None, ""):
                row_hash = _sequence_hash(row_hash, where + " sequence_sha256")
            else:
                row_hash = None
            if row_hash and mapped_hash and row_hash != mapped_hash:
                raise _error(where, "record sequence_sha256 conflicts with evidence ID map")
            provenance = {"source": reader.label(p), "target_id": target,
                          "upstream_adjudication": row.get("evidence_adjudication", {})}
            if row_hash or mapped_hash:
                provenance["sequence_sha256"] = row_hash or mapped_hash
        else:
            sample = _bound_sample(row, {}, known, where)
            gene = _string(row.get("gene_id"), where + " gene_id")
            provenance = row.get("provenance", {})
        if (sample, gene) in seen:
            raise _error(where, "duplicate sample/gene evidence record")
        seen.add((sample, gene))
        if "evidence" not in row:
            raise _error(where, "missing evidence key; use explicit null for an unassessed candidate")
        evidence = row.get("evidence")
        if evidence is not None and not isinstance(evidence, dict):
            raise _error(where, "evidence must be an object; old decision fields are never used for rescoring")
        if not isinstance(provenance, dict):
            raise _error(where, "provenance must be an object")
        result.append(dict(sample_id=sample, gene_id=gene, evidence=evidence, provenance=provenance))
    return result


def load_project(config_path: Path) -> dict:
    """Load a versioned project. All references resolve relative to its config.

    Native exports need an explicit per-source sample_id or a sample_id column.
    Omitted caller/sample combinations are marked missing, never a zero result.
    Scientific policy validation happens in :func:`defenseviz.analysis.integrate`.
    """
    config_path = Path(config_path).expanduser().resolve()
    reader = Reader(config_path)
    config = reader.json(config_path)
    if not isinstance(config, dict) or type(config.get("schema_version")) is not int or config["schema_version"] != 1:
        raise _error(str(config_path), "schema_version must be integer 1")
    allowed = {"schema_version", "samples", "calls", "evidence", "genes", "coverage", "family_map", "options", "evidence_extension"}
    if set(config) - allowed:
        raise _error(str(config_path), "unknown config keys: " + ", ".join(sorted(set(config) - allowed)))
    desc = _descriptor(config.get("samples"), "samples")
    sp = reader.path(desc["path"])
    rows, headers = reader.table(sp, sheet=desc.get("sheet"))
    _require(headers, ["sample_id"], sp)
    samples, known = [], set()
    for i, row in enumerate(rows, 2):
        where = f"{sp}:{i}"
        sid = _string(row.get("sample_id"), where + " sample_id")
        if sid in known:
            raise _error(where, f"duplicate sample_id {sid!r}")
        known.add(sid)
        sample = {"sample_id": sid, "label": _string(row.get("label") or sid, where + " label"),
                  "highlight": _boolean(row.get("highlight"), where + " highlight")}
        for key in ("genus", "species"):
            value = _string(row.get(key), where + " " + key, optional=True)
            if value is not None:
                sample[key] = value
        samples.append(sample)
    if not samples:
        raise _error(str(sp), "at least one sample is required")
    family_map = _family_maps(reader, config)
    sources = config.get("calls", [])
    if not isinstance(sources, list):
        raise _error("calls", "expected a list of source descriptors")
    calls, coverage, seen = [], {}, set()
    for spec in sources:
        imported, inferred = _calls(reader, spec, known, family_map)
        for call in imported:
            key = tuple(call[k] for k in ("sample_id", "gene_id", "caller", "family", "system_id", "evidence_kind"))
            if key in seen:
                raise _error(call["source"], f"duplicate call identity {key!r}; do not import both systems and genes for the same caller run")
            seen.add(key)
            calls.append(call)
        for row in inferred:
            coverage[(row["sample_id"], row["caller"])] = row
    if config.get("coverage"):
        desc = _descriptor(config["coverage"], "coverage")
        cp = reader.path(desc["path"])
        cr, ch = reader.table(cp, sheet=desc.get("sheet"))
        _require(ch, ["sample_id", "caller", "status"], cp)
        explicit = set()
        for i, row in enumerate(cr, 2):
            where = f"{cp}:{i}"
            sid = _bound_sample(row, {}, known, where)
            caller = _string(row["caller"], where + " caller")
            status = row["status"]
            if status not in STATUSES:
                raise _error(where, "status must be complete, missing, or not_run")
            key = sid, caller
            if key in explicit:
                raise _error(where, "duplicate sample/caller coverage")
            if key in coverage and status != "complete":
                raise _error(where, "coverage conflicts with an explicitly supplied completed per-sample source")
            explicit.add(key)
            coverage[key] = dict(sample_id=sid, caller=caller, status=status)
    callers = list(dict.fromkeys([c["caller"] for c in calls if c["evidence_kind"] != "reference_hit"] + [k[1] for k in coverage]))
    for sample in samples:
        for caller in callers:
            coverage.setdefault((sample["sample_id"], caller), dict(sample_id=sample["sample_id"], caller=caller, status="missing"))
    for call in calls:
        key = call["sample_id"], call["caller"]
        if call["evidence_kind"] == "final_member" and coverage[key]["status"] != "complete":
            raise _error(call["source"], f"calls for {key!r} need explicit complete coverage; multi-sample tables cannot imply omitted-sample completeness")
    evidence = _evidence(reader, config.get("evidence"), known)
    genes, seen = [], set()
    if config.get("genes"):
        desc = _descriptor(config["genes"], "genes")
        gp = reader.path(desc["path"])
        gr, gh = reader.table(gp, sheet=desc.get("sheet"))
        _require(gh, ["sample_id", "gene_id", "contig", "start", "end", "strand"], gp)
        for i, row in enumerate(gr, 2):
            where = f"{gp}:{i}"
            sid = _bound_sample(row, {}, known, where)
            gene = _string(row["gene_id"], where + " gene_id")
            if (sid, gene) in seen:
                raise _error(where, "duplicate sample/gene coordinate record")
            seen.add((sid, gene))
            start, end = _integer(row["start"], where + " start"), _integer(row["end"], where + " end")
            if start > end or row["strand"] not in {"+", "-"}:
                raise _error(where, "coordinates require start <= end and strand + or -")
            gene_row = dict(sample_id=sid, gene_id=gene, contig=_string(row["contig"], where + " contig"), start=start, end=end, strand=row["strand"])
            if row.get("sequence_sha256"):
                digest = _string(row["sequence_sha256"], where + " sequence_sha256")
                if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                    raise _error(where, "sequence_sha256 must be 64 lowercase hexadecimal characters")
                gene_row["sequence_sha256"] = digest
            genes.append(gene_row)
    options = config.get("options", {})
    if not isinstance(options, dict):
        raise _error("options", "expected an object")
    extension = None
    if config.get("evidence_extension") is not None:
        spec = config["evidence_extension"]
        if not isinstance(spec, dict) or set(spec) - {"registry", "candidates", "admission_rules", "source_root"}:
            raise _error("evidence_extension", "expected registry, candidates, admission_rules and optional source_root")
        for key in ("registry", "candidates", "admission_rules"):
            if key not in spec:
                raise _error("evidence_extension", f"missing {key}")
        registry_path = reader.path(spec["registry"])
        extension = {"registry": reader.json(registry_path),
                     "candidates": reader.json(reader.path(spec["candidates"])),
                     "admission_rules": reader.json(reader.path(spec["admission_rules"])),
                     "registry_dir": registry_path.parent,
                     "source_root": reader.path(spec["source_root"]) if spec.get("source_root") else reader.base}
        if not extension["source_root"].is_dir():
            raise _error("evidence_extension.source_root", "source root must be an existing directory")
    return dict(samples=samples, calls=calls, evidence_records=evidence, genes=genes,
                coverage=list(coverage.values()), inputs=reader.inputs, input_sha256=reader.hashes,
                options=options, evidence_extension=extension)


FEATURE_NUMERIC = (
    "blast_identity", "blast_query_coverage", "blast_reference_coverage", "blast_evalue",
    "curated_blast_identity", "curated_blast_query_coverage", "curated_blast_reference_coverage", "curated_blast_evalue",
)
FEATURE_FLAGS = (
    "blast_role_compatible", "hmm_role_ambiguous", "rease_nuclease_verified",
    "rease_motor_verified", "mtase_motif_verified", "numeric_structure_admissible_structural_evidence",
    "curated_blast_role_compatible", "curated_blast_family_compatible",
    *("factor_" + key for key in ("strong", "moderate", "profile", "annotation", "motif", "context", "alternative_annotation", "conflict", "partial")),
)
FEATURE_REQUIRED = ("family_id", "enzyme_role", "sequence_length_aa", "blast_family_compatible",
                    "rease_motif_hits_raw", "profiles", "divergent_context", *FEATURE_NUMERIC, *FEATURE_FLAGS)


def _feature_evidence(reader, spec, known):
    """Adapt declared native factors, preserving incomplete rows as unassessed.

    Individual profile records must be supplied; independent per-column maxima
    cannot reconstruct a coherent HMM hit. Recognition fields are never read.
    """
    from dataclasses import asdict
    from defenseviz.analysis.policy import from_feature, Profile, truth
    from defenseviz.analysis.evidence import parse_evidence
    p = reader.path(spec["path"])
    rows, headers = reader.table(p, sheet=spec.get("sheet"))
    mapping = _maps(reader, spec, known)
    id_column = spec.get("id_column") or next((k for k in ("gene_id", "locus_tag", "query_id", "target_id") if k in headers), None)
    if not id_column:
        raise _error(str(p), "DNMB feature table requires gene_id/locus_tag/query_id or an explicit id_column")
    _require(headers, [id_column], p)
    if not spec.get("sample_id"):
        _require(headers, ["sample_id"], p)
    results, seen = [], set()
    for i, native in enumerate(rows, 2):
        where = f"{p}:{i}"
        sid = _bound_sample(native, spec, known, where)
        gene = _translate(native[id_column], sid, mapping, where)
        if (sid, gene) in seen:
            raise _error(where, "duplicate sample/gene native feature record")
        seen.add((sid, gene))
        row = dict(native)
        if "family_id" not in row and "rm_type" in row:
            row["family_id"] = row["rm_type"]
        if "enzyme_role" not in row and "subunit" in row:
            row["enzyme_role"] = row["subunit"]
        missing = [k for k in FEATURE_REQUIRED if k not in row]
        # A missing boolean is not evidence that the predicate was evaluated false.
        missing.extend(k for k in (*FEATURE_FLAGS, "divergent_context") if k in row and (row[k] is None or str(row[k]).strip().lower() in {"", "na", "nan", "none", "null", "-"}))
        if row.get("profiles") in (None, "") and "profiles" not in missing:
            missing.append("profiles")
        provenance = {"source": reader.label(p), "format": "dnmb_features", "source_row": i,
                      "family": row.get("family_id"), "role": row.get("enzyme_role"),
                      "missing_fields": sorted(set(missing)),
                      "imported_native_fields": row}
        if row.get("sequence_sha256") not in (None, ""):
            provenance["sequence_sha256"] = _sequence_hash(row["sequence_sha256"], where + " sequence_sha256")
        if missing:
            results.append(dict(sample_id=sid, gene_id=gene, evidence=None, status="unassessed", provenance=provenance))
            continue
        try:
            for key in FEATURE_NUMERIC:
                observed = row[key]
                if observed is None or str(observed).strip().lower() in {"", "na", "nan", "none", "null", "-"}:
                    continue
                if isinstance(observed, bool):
                    raise ValueError(f"{key} must be numeric, not boolean")
                number = float(observed)
                if not math.isfinite(number):
                    raise ValueError(f"{key} must be finite")
            profiles = row["profiles"]
            if isinstance(profiles, str):
                profiles = json.loads(profiles, object_pairs_hook=_unique_object, parse_constant=_nonfinite)
            if not isinstance(profiles, list) or any(not isinstance(x, dict) for x in profiles):
                raise ValueError("profiles must be a JSON list of coherent Profile objects")
            profile_values = [Profile(**x) for x in profiles]
            value = asdict(from_feature(row, profile_values, divergent=truth(row["divergent_context"])))
            value["profiles"] = list(value["profiles"])
            value["raw_nuclease_classes"] = sorted(value["raw_nuclease_classes"])
            value["reference_nuclease_anchor_profile_indices"] = list(value["reference_nuclease_anchor_profile_indices"])
            parse_evidence(value)
        except (ValueError, TypeError, KeyError) as exc:
            raise _error(where, f"invalid declared DNMB features: {exc}") from exc
        results.append(dict(sample_id=sid, gene_id=gene, evidence=value, provenance=provenance))
    return results

"""Load explicit reviews bound to the saved inputs and current protein sequences.

This module checks provenance and merges records. Biological interpretation is
owned by analysis.saved; no identifier, label or review wording is interpreted.
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from .io.project import Reader
from .pipeline import digest


TABLES = frozenset({"annotations", "effector_evidence", "groups", "complexes"})


def _sha(value, where):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"{where}: expected a lowercase SHA-256")
    return value


def apply_supplement(data: dict, path: str | Path) -> dict:
    """Return a copy with additive source-bound records; originals stay unchanged.

    The supplement must bind at least one exact imported native artifact. Every
    gene-specific record must additionally match the current whole-protein hash.
    Group/complex membership binds each member through explicit member hashes.
    Source documents are hashed without assuming they constitute valid biology.
    """
    path = Path(path).resolve()
    reader = Reader(path)
    doc = reader.json(path)
    if not isinstance(doc, dict) or type(doc.get("schema_version")) is not int or doc["schema_version"] != 1:
        raise ValueError("Supplement requires schema_version: 1")
    unknown = set(doc) - {"schema_version", "input_sha256", "description", *TABLES}
    if unknown:
        raise ValueError(f"Unknown supplement fields: {sorted(unknown)}")
    bindings = doc.get("input_sha256")
    if not isinstance(bindings, dict) or not bindings:
        raise ValueError("Supplement requires nonempty input_sha256 bindings to imported native artifacts")
    imported = {str(Path(p).resolve()): sha for p, sha in data.get("input_sha256", {}).items()}
    if not imported:
        raise ValueError("Native input hashes are required before applying a supplement")
    bound_samples = set()
    input_roles = data.get("input_roles", {})
    for raw, expected in bindings.items():
        p = reader.path(raw)
        _sha(expected, f"input_sha256[{raw!r}]")
        if imported.get(str(p)) != expected:
            raise ValueError(f"Supplement input binding does not match imported native artifact: {raw}")
        reader.read(p)
        if reader.hashes[str(p)] != expected:
            raise ValueError(f"Supplement input changed: {raw}")
        roles = input_roles.get(str(p), [])
        bound_samples.update(r["sample_id"] for r in roles
                             if r.get("native_result") is True and r.get("sample_id"))
    if not bound_samples:
        raise ValueError("Supplement bindings must include native result artifacts, not only manifests or maps")

    known = {}
    for gene in data.get("genes", []):
        key = (gene["sample_id"], gene["gene_id"])
        sha = gene.get("protein_sha256", gene.get("sequence_sha256"))
        if key in known and known[key] != sha:
            raise ValueError(f"Conflicting native protein hash: {key}")
        known[key] = sha

    def check_gene(sid, gid, sha, where):
        key = (sid, gid)
        _sha(sha, where + ".protein_sha256")
        if key not in known or known[key] is None or known[key] != sha:
            raise ValueError(f"{where}: current protein hash does not match native gene {key}")

    def check_ref(ref, where):
        if not isinstance(ref, dict) or not {"path", "sha256"} <= set(ref):
            raise ValueError(f"{where}: source reference requires path and sha256")
        p = reader.path(ref["path"])
        sha = _sha(ref["sha256"], where + ".sha256")
        reader.read(p)
        if reader.hashes[str(p)] != sha:
            raise ValueError(f"{where}: source hash mismatch: {ref['path']}")
        return {**ref, "path": str(p)}

    result = deepcopy(data)
    native_group_keys = {(g["sample_id"], g["group_id"]) for g in data.get("groups", [])}
    added_group_keys = set()
    for table in ("annotations", "effector_evidence", "groups", "complexes"):
        records = doc.get(table, [])
        if not isinstance(records, list) or any(not isinstance(r, dict) for r in records):
            raise ValueError(f"Supplement {table} must be a list of objects")
        for i, record in enumerate(records):
            row = deepcopy(record)
            where = f"Supplement {table}[{i}]"
            sid = row.get("sample_id")
            if not isinstance(sid, str) or not sid:
                raise ValueError(f"{where}: sample_id is required")
            if sid not in bound_samples:
                raise ValueError(f"{where}: no bound native result artifact for sample {sid!r}")
            if "origin" in row and row["origin"] != "supplement":
                raise ValueError(f"{where}: supplement cannot claim a native origin")
            row["origin"] = "supplement"
            if table in {"annotations", "effector_evidence"}:
                check_gene(sid, row.get("gene_id"), row.get("protein_sha256"), where)
            elif table == "groups":
                if row.get("membership_kind") != "saved_context":
                    raise ValueError(f"{where}: supplement groups must use saved_context, never native_system")
                if row.get("provider", "Review") != "Review" or "run_id" in row:
                    raise ValueError(f"{where}: supplement groups cannot claim a native provider/run")
                row["provider"] = "Review"
                row["run_id"] = reader.hashes[str(path)]
                gid = row.get("group_id")
                if not isinstance(gid, str) or not gid:
                    raise ValueError(f"{where}: group_id is required")
                key = (sid, gid)
                if key in native_group_keys or key in added_group_keys:
                    raise ValueError(f"{where}: duplicate group; choose an explicit new review group ID")
                added_group_keys.add(key)
                members = row.get("members")
                if not isinstance(members, list) or not members:
                    raise ValueError(f"{where}: nonempty members are required")
                for j, member in enumerate(members):
                    if not isinstance(member, dict):
                        raise ValueError(f"{where}: member must be an object")
                    if "retained" in member:
                        raise ValueError(f"{where}: review members cannot claim native retention")
                    check_gene(sid, member.get("gene_id"), member.get("protein_sha256"), f"{where}.members[{j}]")
            else:
                core = row.get("core_gene_ids")
                hashes = row.pop("core_protein_sha256", None)
                if not isinstance(core, list) or not core or not isinstance(hashes, dict) or set(hashes) != set(core):
                    raise ValueError(f"{where}: core_protein_sha256 must bind every core_gene_id")
                for gid in core:
                    check_gene(sid, gid, hashes[gid], where)
                if (sid, row.get("group_id")) not in native_group_keys | added_group_keys:
                    raise ValueError(f"{where}: complex group is not present")
                # Retain the explicit bindings for reproducibility.
                row["core_protein_sha256"] = hashes
            if "source_ref" in row:
                row["source_ref"] = check_ref(row["source_ref"], where + ".source_ref")
            refs = row.get("source_refs")
            if refs is not None:
                if not isinstance(refs, list) or not refs:
                    raise ValueError(f"{where}: source_refs must be a nonempty list")
                row["source_refs"] = [check_ref(ref, where + ".source_refs") for ref in refs]
            if "source_ref" not in row and not refs:
                raise ValueError(f"{where}: a hashed source reference is required")
            row["supplement_source"] = {"path": str(path), "sha256": reader.hashes[str(path)], "table": table, "index": i}
            result.setdefault(table, []).append(row)
    result["inputs"] = sorted({Path(p).resolve() for p in data.get("inputs", [])} | set(reader.inputs))
    result["input_sha256"] = {**imported, **reader.hashes}
    # Recheck bound files after all reads, including shared references.
    for p, sha in reader.hashes.items():
        if digest(Path(p)) != sha:
            raise ValueError(f"Supplement source changed during import: {p}")
    return result

"""Load optional exact review rules and evaluate source-bound domain evidence.

This adapter supplies evidence to the additive system review. It never changes
native detections, accepted annotations, or prioritization membership.
"""
from __future__ import annotations

from pathlib import Path

from .io.project import Reader
from .pipeline import digest


def apply_system_review_inputs(data, *, rules_path=None, evidence_path=None):
    """Attach pinned review inputs to an imported, mutable suite bundle."""
    readers = []
    additions = {}
    if rules_path is not None:
        from .analysis.system_review import validate_system_review_rules

        path = Path(rules_path).resolve()
        reader = Reader(path)
        rules = reader.json(path)
        validate_system_review_rules(rules)
        for rule in rules["rules"]:
            contract = rule.get("architecture_contract")
            if contract is not None:
                ref = contract["source_ref"]
                ref_path = reader.path(ref["path"])
                if digest(ref_path) != ref["sha256"]:
                    raise ValueError("System architecture contract source hash mismatch")
                reader.hashes[str(ref_path)] = ref["sha256"]
                reader.inputs.append(ref_path)
        additions["system_review_rules"] = rules
        additions["system_review_rule_source"] = {
            "path": str(path), "sha256": reader.hashes[str(path)]}
        additions["system_review_rules_base"] = str(path.parent)
        readers.append(reader)

    if evidence_path is not None:
        from .analysis.registry import assess_registry

        path = Path(evidence_path).resolve()
        reader = Reader(path)
        pack = reader.json(path)
        allowed = {"schema_version", "registry", "candidates", "source_root"}
        if (not isinstance(pack, dict) or type(pack.get("schema_version")) is not int
                or pack["schema_version"] != 1 or set(pack) - allowed
                or not {"registry", "candidates"} <= set(pack)):
            raise ValueError("System review evidence requires schema_version 1, registry and candidates paths")
        registry_path = reader.path(pack["registry"])
        candidates_path = reader.path(pack["candidates"])
        registry = reader.json(registry_path)
        candidates = reader.json(candidates_path)
        if not isinstance(candidates, list):
            raise ValueError("System review candidates must be a list")
        genes = {(r["sample_id"], r["gene_id"]): r for r in data.get("genes", [])}
        for row in candidates:
            if not isinstance(row, dict):
                raise ValueError("System review candidate must be an object")
            sid, gid = row.get("sample_id"), row.get("query_id")
            if not isinstance(sid, str) or not isinstance(gid, str):
                raise ValueError("System review candidate requires exact sample_id and query_id")
            gene = genes.get((sid, gid))
            if gene is None or not gene.get("protein_sha256"):
                raise ValueError("System review candidate lacks a current suite protein identity")
            if row.get("sequence_sha256") != gene["protein_sha256"]:
                raise ValueError("System review candidate protein differs from current suite gene")
        source_root = reader.path(pack["source_root"]) if "source_root" in pack else registry_path.parent
        # Pin external reference sequences before evaluation as well as after it.
        if isinstance(registry, dict):
            for family in registry.get("families", []):
                for reference in family.get("references", []):
                    if reference.get("sequence_file"):
                        ref = Path(reference["sequence_file"])
                        ref = (registry_path.parent / ref).resolve() if not ref.is_absolute() else ref.resolve()
                        reader.hashes[str(ref)] = digest(ref)
                        reader.inputs.append(ref)
        assessments = assess_registry(registry, candidates, registry_dir=registry_path.parent,
                                      source_root=source_root, verify_sources=True)
        for assessment in assessments:
            for name, sha in assessment.get("source_sha256", {}).items():
                if name in reader.hashes and reader.hashes[name] != sha:
                    raise ValueError(f"System review input changed during evaluation: {name}")
                reader.hashes[name] = sha
                ref = Path(name)
                if ref not in reader.inputs:
                    reader.inputs.append(ref)
        additions["system_review_evidence"] = {
            "schema_version": 1, "registry": registry, "candidates": candidates,
            "assessments": assessments, "input_sha256": dict(reader.hashes),
            "registry_sha256": reader.hashes[str(registry_path)],
            "candidates_sha256": reader.hashes[str(candidates_path)],
            "registry_path": str(registry_path), "candidates_path": str(candidates_path),
            "registry_dir": str(registry_path.parent), "source_root": str(source_root),
            "source_files_verified": True,
            "operation": "evaluation_of_supplied_evidence_not_detector_run"}
        readers.append(reader)

    pinned = dict(data.get("input_sha256", {}))
    for reader in readers:
        for name, sha in reader.hashes.items():
            if name in pinned and pinned[name] != sha:
                raise ValueError(f"System review input conflicts with imported artifact: {name}")
            if digest(Path(name)) != sha:
                raise ValueError(f"System review input changed during loading: {name}")
            pinned[name] = sha
    data.update(additions)
    data["input_sha256"] = pinned
    inputs = data.setdefault("inputs", [])
    for reader in readers:
        for path in reader.inputs:
            if path not in inputs:
                inputs.append(path)
    return data

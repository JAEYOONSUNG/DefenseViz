"""Source-bound re-evaluation of explicit saved R-M Evidence for suite imports.

This bridge runs the packaged, frozen numerical policy. It does not reconstruct
missing axes, search sequences, infer partners, or alter original tool calls.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, fields
from pathlib import Path
import hashlib
import json
from math import isfinite
import time

from .analysis.evidence import evidence_dict, parse_evidence, parse_policy, policy_metadata
from .analysis.policy import Evidence, Policy, evaluate
from .io.project import Reader
from .pipeline import digest

RM_TYPES = {"Type I", "Type II", "Type III", "Type IV"}
SCORABLE_TYPES = RM_TYPES | {"NA", "Unknown", "Unresolved"}
SCORABLE_ROLES = {"R", "M", "S", "RM"}
CANONICAL_ROLES = SCORABLE_ROLES | {"Aux", "A", "C", "F", "unknown", "Unknown"}
CONTEXT_STATES = {"complete", "met", "met_in_saved_context", "missing",
                  "missing_retained_partner_in_saved_context", "unassessed"}


def _sha(value, where):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"{where}: expected lowercase SHA-256")
    return value


def _key(row, samples, where):
    sid, gid = row.get("sample_id"), row.get("gene_id")
    if not isinstance(sid, str) or sid not in samples or not isinstance(gid, str) or not gid:
        raise ValueError(f"{where}: exact known sample_id and nonempty gene_id required")
    return sid, gid


def _candidate_integrity_issues(data, sid, native_decisions):
    """Keep candidate-universe defects distinct from unrelated raw hit gaps."""
    candidate_paths = set()
    candidate_modules = set()
    for raw, roles in data.get("input_roles", {}).items():
        for role in roles:
            if role.get("sample_id") == sid and role.get("kind") == "saved_candidates":
                candidate_paths.add(str(Path(raw).resolve()))
                if role.get("module"):
                    candidate_modules.add(role["module"])
    for key, row in native_decisions.items():
        if key[0] == sid:
            for ref in [row.get("source_ref"), *row.get("source_refs", [])]:
                if isinstance(ref, dict) and isinstance(ref.get("path"), str):
                    candidate_paths.add(str(Path(ref["path"]).resolve()))

    def candidate_source(row):
        ref = row.get("source_ref") or {}
        if ref.get("path"):
            return str(Path(ref["path"]).resolve()) in candidate_paths
        return row.get("module") in candidate_modules

    suite = data.get("suite") or {}
    issues = []
    explicit_codes = {"unbound_native_candidate", "cached_live_candidate_count_mismatch",
                      "sentinel_live_retained_count_mismatch"}
    for row in suite.get("conflicts", []):
        if row.get("sample_id") == sid and (row.get("code") in explicit_codes or
                (row.get("code") in {"sentinel_live_row_count_mismatch", "cached_failure_with_live_artifact"}
                 and candidate_source(row))):
            issues.append(dict(kind="native_candidate_conflict", record=deepcopy(row)))
    for row in suite.get("unmapped_records", []):
        if row.get("sample_id") == sid and candidate_source(row):
            issues.append(dict(kind="unmapped_native_candidate", record=deepcopy(row)))
    return issues


def apply_rm_evidence(data: dict, path: str | Path) -> dict:
    """Return a copied suite bundle with evaluated decisions and full provenance.

    Schema v1 has input_sha256, the full selected_policy parameter object,
    records and candidate_scope. Every record binds an imported whole protein
    and hashed source_refs. Optional context_state/orphan require separate
    context_source_refs; old_context never supplies either statement.

    Default evaluate_unassessed mode rejects contradictory assessed decisions.
    Explicit replace_scoped mode supersedes all old decisions in scoped samples:
    omitted candidates remain present as unassessed, never stale retained calls.
    Optional expected_decision verifies a hashed historical keep/score record;
    new evaluations need not have a historical result.
    """
    if not isinstance(data, dict):
        raise ValueError("R-M evaluation needs a canonical suite bundle")
    if data.get("rm_evaluation") is not None:
        raise ValueError("R-M evaluation already applied; start from the native bundle")
    started = time.monotonic()
    path = Path(path).resolve()
    reader = Reader(path)
    doc = reader.json(path)
    allowed = {"schema_version", "input_sha256", "selected_policy", "records", "candidate_scope", "description", "decision_mode"}
    if not isinstance(doc, dict) or type(doc.get("schema_version")) is not int or doc["schema_version"] != 1:
        raise ValueError("R-M evidence requires schema_version: 1")
    if set(doc) - allowed:
        raise ValueError(f"Unknown R-M evidence fields: {sorted(set(doc)-allowed)}")
    decision_mode = doc.get("decision_mode", "evaluate_unassessed")
    if not isinstance(decision_mode, str) or decision_mode not in {"evaluate_unassessed", "replace_scoped"}:
        raise ValueError("decision_mode must be evaluate_unassessed or replace_scoped")
    parameters = doc.get("selected_policy")
    if not isinstance(parameters, dict) or set(parameters) != {f.name for f in fields(Policy)}:
        raise ValueError("selected_policy must supply the full policy parameter object")
    policy = parse_policy(parameters)
    if asdict(policy) != asdict(Policy()):
        raise ValueError("Suite R-M evaluation requires the unchanged selected policy")
    policy_sha = hashlib.sha256(json.dumps(parameters, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    samples = {s["sample_id"] for s in data.get("samples", [])}
    if not samples:
        raise ValueError("R-M evaluation requires imported samples")
    bindings = doc.get("input_sha256")
    imported = {str(Path(p).resolve()): sha for p, sha in data.get("input_sha256", {}).items()}
    if not isinstance(bindings, dict) or not bindings or not imported:
        raise ValueError("R-M evidence needs nonempty native artifact input_sha256 bindings")
    bound_samples = set()
    for raw, expected in bindings.items():
        p = reader.path(raw)
        _sha(expected, f"input_sha256[{raw!r}]")
        if imported.get(str(p)) != expected:
            raise ValueError(f"R-M native artifact binding mismatch: {raw}")
        reader.read(p)
        if reader.hashes[str(p)] != expected:
            raise ValueError(f"R-M native artifact changed: {raw}")
        bound_samples.update(r.get("sample_id") for r in data.get("input_roles", {}).get(str(p), [])
                             if r.get("native_result") is True and r.get("sample_id") in samples)
    if not bound_samples:
        raise ValueError("R-M bindings must include a native result artifact, not only a manifest/map")

    def refs(values, where):
        if not isinstance(values, list) or not values:
            raise ValueError(f"{where}: nonempty hashed source_refs required")
        out = []
        for ref in values:
            if not isinstance(ref, dict) or not {"path", "sha256"} <= set(ref):
                raise ValueError(f"{where}: source reference needs path and sha256")
            p = reader.path(ref["path"])
            expected = _sha(ref["sha256"], where)
            reader.read(p)
            if reader.hashes[str(p)] != expected:
                raise ValueError(f"{where}: source hash mismatch: {ref['path']}")
            out.append({**deepcopy(ref), "path": str(p)})
        return out

    genes = {}
    for g in data.get("genes", []):
        key = _key(g, samples, "native gene")
        if key in genes:
            raise ValueError(f"Duplicate native gene identity: {key}")
        genes[key] = g.get("protein_sha256", g.get("sequence_sha256"))
    native_decisions = {}
    for r in data.get("saved_decisions", []):
        key = _key(r, samples, "saved decision")
        if key in native_decisions:
            raise ValueError(f"Duplicate saved decision; select one active source: {key}")
        if key not in genes:
            raise ValueError(f"Saved decision identity is absent from native genes: {key}")
        if r.get("protein_sha256") is not None and r["protein_sha256"] != genes[key]:
            raise ValueError(f"Saved decision protein hash contradicts native gene: {key}")
        native_decisions[key] = r
    records = doc.get("records")
    scopes = doc.get("candidate_scope")
    if not isinstance(records, list) or any(not isinstance(r, dict) for r in records):
        raise ValueError("R-M records must be a list of complete Evidence records")
    if not isinstance(scopes, list) or not scopes or any(not isinstance(s, dict) for s in scopes):
        raise ValueError("R-M candidate_scope must explicitly describe every evaluated sample")
    scope_by_sample = {}
    for i, scope in enumerate(scopes):
        if set(scope) - {"sample_id", "status", "gene_ids", "source_refs", "description"}:
            raise ValueError("Unknown candidate_scope fields")
        sid = scope.get("sample_id")
        ids = scope.get("gene_ids")
        if not isinstance(sid, str) or sid not in bound_samples or sid in scope_by_sample:
            raise ValueError("candidate_scope needs one native-bound record per sample")
        if scope.get("status") not in {"complete", "partial"}:
            raise ValueError("candidate_scope status must be complete or partial")
        if (not isinstance(ids, list) or any(not isinstance(x, str) or not x for x in ids)
                or len(set(ids)) != len(ids) or any((sid,x) not in genes for x in ids)):
            raise ValueError("candidate_scope gene_ids require distinct known native identities")
        scope_by_sample[sid] = {**deepcopy(scope), "source_refs": refs(scope.get("source_refs"), f"candidate_scope[{i}]")}
    result = deepcopy(data)
    decisions = deepcopy(native_decisions)
    evaluated = []
    seen = set()
    evaluation_id = reader.hashes[str(path)]
    if decision_mode == "replace_scoped":
        for key, row in decisions.items():
            if key[0] not in scope_by_sample:
                continue
            row.update(state="unassessed", score_total=None,
                       decision_reason="native_decision_superseded_by_explicit_scoped_evaluation",
                       origin="stored_evidence_evaluation", evaluation_id=evaluation_id,
                       native_source_ref=deepcopy(native_decisions[key].get("source_ref")),
                       native_source_refs=deepcopy(native_decisions[key].get("source_refs", [])))
            for name in ("context_state", "orphan", "context_source_refs", "integrated_R_component_support"):
                row.pop(name, None)
    for i, r in enumerate(records):
        where = f"R-M records[{i}]"
        permitted = {"sample_id", "gene_id", "protein_sha256", "evidence", "source_refs", "context_state", "orphan", "context_source_refs", "provenance", "expected_decision"}
        if set(r) - permitted:
            raise ValueError(f"{where}: unknown fields {sorted(set(r)-permitted)}")
        key = _key(r, samples, where)
        if key in seen:
            raise ValueError(f"{where}: duplicate evaluated identity {key}")
        seen.add(key)
        sid, gid = key
        if sid not in scope_by_sample or gid not in scope_by_sample[sid]["gene_ids"]:
            raise ValueError(f"{where}: identity outside its declared candidate_scope")
        protein = _sha(r.get("protein_sha256"), where + ".protein_sha256")
        if genes.get(key) != protein:
            raise ValueError(f"{where}: current native protein hash mismatch")
        sources = refs(r.get("source_refs"), where)
        raw_evidence = r.get("evidence")
        if not isinstance(raw_evidence, dict) or set(raw_evidence) != {f.name for f in fields(Evidence)}:
            raise ValueError(f"{where}: full Evidence fieldset required; absent axes cannot become false")
        evidence = parse_evidence(raw_evidence)
        context = {}
        if "context_state" in r or "orphan" in r:
            context["context_source_refs"] = refs(r.get("context_source_refs"), where + ".context")
        elif "context_source_refs" in r:
            raise ValueError(f"{where}: context_source_refs without an explicit assessment")
        if "context_state" in r:
            if r["context_state"] not in CONTEXT_STATES:
                raise ValueError(f"{where}: unsupported explicit context_state")
            context["context_state"] = r["context_state"]
        if "orphan" in r:
            if type(r["orphan"]) is not bool:
                raise ValueError(f"{where}: orphan must be an explicit boolean")
            context["orphan"] = r["orphan"]
        if "provenance" in r and not isinstance(r["provenance"], dict):
            raise ValueError(f"{where}: provenance must be an object")
        decision = None
        if evidence.family not in SCORABLE_TYPES or evidence.role not in SCORABLE_ROLES:
            state, reason = "out_of_scope", "selected_policy_type_or_role_out_of_scope"
        elif not evidence.valid_sequence:
            state, reason = "unassessed", "invalid_or_missing_sequence_cannot_supply_policy_support"
        else:
            decision = evaluate(evidence, policy)
            # Tuples from the historical predicate become ordinary JSON arrays.
            decision = json.loads(json.dumps(decision, allow_nan=False))
            state = "retained" if decision["keep"] else "not_retained"
            reason = "evaluated_unchanged_selected_policy"
        expected_decision = None
        if "expected_decision" in r:
            expected = r["expected_decision"]
            if not isinstance(expected, dict) or set(expected) != {"keep", "score_total", "source_refs"}:
                raise ValueError(f"{where}: expected_decision needs exactly keep, score_total and source_refs")
            score = expected["score_total"]
            if type(expected["keep"]) is not bool or isinstance(score, bool) or not isinstance(score, (int, float)) or not isfinite(score):
                raise ValueError(f"{where}: expected_decision needs a boolean keep and finite numerical score_total")
            expected_decision = {**deepcopy(expected), "source_refs": refs(expected["source_refs"], where + ".expected_decision")}
            if decision is None:
                raise ValueError(f"{where}: cannot verify expected_decision for an unscored record")
            if expected["keep"] != decision["keep"] or score != decision["score_total"]:
                raise ValueError(f"{where}: expected_decision parity mismatch")
        previous = native_decisions.get(key)
        canonical_role = evidence.role if evidence.role in CANONICAL_ROLES else None
        if decision_mode == "evaluate_unassessed" and previous is not None and previous.get("state") != "unassessed":
            if previous.get("state") != state or any(previous.get(k) not in {None, value} for k,value in
                    [("family", evidence.family), ("role", canonical_role)]) or (
                    previous.get("rm_type") in RM_TYPES and previous["rm_type"] != evidence.family):
                raise ValueError(f"{where}: conflicting assessed saved decision; select one current source")
        generated_ref = {"path": str(path), "sha256": evaluation_id, "record_index": i}
        row = {**deepcopy(previous or {}), "sample_id": sid, "gene_id": gid,
               "protein_sha256": protein, "state": state, "family": evidence.family,
               "rm_type": evidence.family if evidence.family in RM_TYPES else None, "role": canonical_role,
               "native_evidence_role": evidence.role,
               "score_total": decision["score_total"] if decision is not None else None,
               "decision_reason": reason, "origin": "stored_evidence_evaluation",
               "evaluation_id": evaluation_id, "policy_sha256": policy_sha,
               "source_ref": generated_ref, "source_refs": sources, **context}
        row["native_source_ref"] = deepcopy((previous or {}).get("source_ref"))
        row["native_source_refs"] = deepcopy((previous or {}).get("source_refs", []))
        if decision is not None:
            row["integrated_R_component_support"] = decision["integrated_R_component_support"]
        else:
            row.pop("integrated_R_component_support", None)
        # A pre-retention context cannot silently become a current assessment.
        for name in ("context_state", "orphan", "context_source_refs"):
            if name not in context:
                row.pop(name, None)
        decisions[key] = row
        evaluated.append(dict(sample_id=sid, gene_id=gid, protein_sha256=protein, evidence=evidence_dict(evidence),
            decision=decision, state=state, reason=reason, source_refs=sources, evaluation_source=generated_ref,
            expected_decision=expected_decision, parity_status="verified" if expected_decision is not None else "not_supplied",
            previous_saved_decision=deepcopy(previous), provenance=deepcopy(r.get("provenance", {})), **context))
    scope_results = []
    for sid, scope in scope_by_sample.items():
        active = {k for k in decisions if k[0] == sid}
        declared = {(sid,gid) for gid in scope["gene_ids"]}
        native_candidates = {k for k in native_decisions if k[0] == sid}
        integrity_issues = _candidate_integrity_issues(data, sid, native_decisions)
        if scope["status"] == "complete" and integrity_issues:
            raise ValueError(f"{sid}: complete candidate_scope conflicts with unresolved native candidate identity/count integrity")
        if scope["status"] == "complete" and declared - native_candidates:
            raise ValueError(f"{sid}: complete candidate_scope adds genes outside the original native candidate set; use partial scope")
        if scope["status"] == "complete" and (active != declared or any(decisions[k]["state"] == "unassessed" for k in active)):
            raise ValueError(f"{sid}: complete candidate_scope omits or leaves unassessed native candidates")
        status = scope["status"]
        scope_results.append({**scope, "coverage_status": status, "active_candidate_count": len(active),
            "native_candidate_count": len(native_candidates), "native_candidate_integrity_issues": integrity_issues,
            "outside_native_candidate_gene_ids": sorted(k[1] for k in active-native_candidates),
            "evaluated_count": sum(k[0] == sid for k in seen),
            "remaining_unassessed_gene_ids": sorted(k[1] for k in active if decisions[k]["state"] == "unassessed"),
            "outside_declared_scope_gene_ids": sorted(k[1] for k in active-declared)})
    result["saved_decisions"] = list(decisions.values())
    native_saved_coverage = [deepcopy(c) for c in data.get("coverage", [])
                             if c.get("endpoint") in {"saved_retention", "saved_decisions"}]
    result["coverage"] = [c for c in result.get("coverage", []) if not
        (c.get("sample_id") in scope_by_sample and c.get("endpoint") in {"saved_retention", "saved_decisions"})]
    result["coverage"].extend(dict(sample_id=s["sample_id"],caller="DNMB",endpoint="saved_retention",
        status=s["coverage_status"],origin="stored_evidence_evaluation",source_ref={"path":str(path),"sha256":evaluation_id}) for s in scope_results)
    policy_evaluations = sum(r["decision"] is not None for r in evaluated)
    transitions = [dict(sample_id=k[0], gene_id=k[1], evaluated=k in seen,
                        previous_state=native_decisions.get(k, {}).get("state"), current_state=r["state"],
                        previous_score_total=native_decisions.get(k, {}).get("score_total"),
                        current_score_total=r.get("score_total"))
                   for k,r in decisions.items() if k[0] in scope_by_sample]
    result["rm_evaluation"] = dict(schema_version=1, operation="fresh_evaluation_of_stored_RM_evidence",
        decision_mode=decision_mode, scoped_samples=sorted(scope_by_sample), decision_transitions=transitions,
        scoring_runs=int(policy_evaluations > 0), policy_evaluations=policy_evaluations, detector_runs=0,
        structure_prediction_runs=0, thresholds_changed=False, selected_policy=asdict(policy),
        policy=policy_metadata(policy), policy_sha256=policy_sha, evaluation_id=evaluation_id,
        records=evaluated, candidate_scope=scope_results, native_saved_decisions=deepcopy(data.get("saved_decisions", [])),
        native_saved_retention_coverage=native_saved_coverage,
        expected_decisions_verified=sum(r["expected_decision"] is not None for r in evaluated),
        elapsed_seconds=time.monotonic()-started, activity_inferred=False)
    result["inputs"] = sorted({Path(p).resolve() for p in data.get("inputs", [])} | set(reader.inputs))
    result["input_sha256"] = {**imported, **reader.hashes}
    for p, sha in reader.hashes.items():
        if digest(Path(p)) != sha:
            raise ValueError(f"R-M evidence source changed during evaluation: {p}")
    return result

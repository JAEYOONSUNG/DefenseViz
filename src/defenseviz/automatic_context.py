"""Reassess retained R-M composition within exact saved DNMB boundaries.

No proximity search, system boundary inference, activity inference, or change to
the R-M decision policy occurs here. Original groups are retained for audit;
new Review groups describe only current retained members of those boundaries.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import hashlib
import json
from pathlib import Path


_ROLES = {"R", "M", "S", "RM", "Aux", "A", "C", "F"}
_TYPES = {"Type I", "Type II", "Type III", "Type IV"}


def _key(row):
    return row["sample_id"], row["gene_id"]


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _refs(row):
    out = []
    for ref in [row.get("source_ref"), *row.get("source_refs", [])]:
        if isinstance(ref, dict) and ref not in out:
            out.append(deepcopy(ref))
    return out


def _merge_refs(rows):
    out = []
    for row in rows:
        for ref in _refs(row):
            if ref not in out:
                out.append(ref)
    return out


def _bound(refs, hashes):
    return bool(refs) and all(isinstance(r.get("path"), str)
                              and isinstance(r.get("sha256"), str) and len(r["sha256"]) == 64
                              and hashes.get(r["path"]) == r["sha256"] for r in refs)


def _roles(role):
    return {"R", "M"} if role == "RM" else {"A"} if role == "Aux" else {role} if role in _ROLES else set()


def _requirements(group, typ):
    explicit = group.get("required_roles")
    if explicit is not None:
        if not isinstance(explicit, (list, tuple)) or not explicit or any(r not in _ROLES for r in explicit):
            raise ValueError("Current-context required_roles must name explicit known component roles")
        return set().union(*(_roles(role) for role in explicit))
    return {"Type I": {"R", "M", "S"}, "Type III": {"R", "M"}}.get(typ)


def apply_current_context_review(data):
    """Return a copied bundle with source-bound current retained contexts.

    Untouched large raw annotation tables are shared; groups, saved decisions,
    and evaluation provenance use copy-on-write. The input is never mutated.
    Existing explicit context/orphan assessments are not overwritten.
    """
    out = dict(data)
    evaluation = data.get("rm_evaluation")
    if not evaluation:
        return out
    original = [g for g in data.get("groups", []) if g.get("provider") == "DNMB"
                and g.get("membership_kind") == "saved_context"]
    decisions = {_key(r): r for r in data.get("saved_decisions", [])}
    genes = {_key(g): g for g in data.get("genes", [])}
    if len(decisions) != len(data.get("saved_decisions", [])):
        raise ValueError("Duplicate R-M decision identity in current-context review")
    signature = _digest({"evaluation_id": evaluation.get("evaluation_id"), "groups": original,
                         "derived_artifact_sha256": data.get("derived_artifact_sha256", {}),
                         "decisions": [{k: r.get(k) for k in ("sample_id", "gene_id", "protein_sha256", "state", "family", "rm_type", "role", "source_ref", "source_refs")}
                                       for _, r in sorted(decisions.items())]})
    prior = evaluation.get("current_context_review")
    hashes = dict(data.get("input_sha256", {}))
    # Generated, source-bound alignments belong to output provenance rather
    # than the original-input list. The caller explicitly registers those
    # artifacts; their bytes must still match before they can bind a context.
    for name, digest in data.get("derived_artifact_sha256", {}).items():
        if not isinstance(name, str) or not isinstance(digest, str):
            raise ValueError("Invalid derived current-context artifact registry")
        path = Path(name)
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError("Derived current-context artifact hash mismatch")
        if name in hashes and hashes[name] != digest:
            raise ValueError("Derived current-context artifact conflicts with original input hash")
        hashes[name] = digest
    # load_suite records display-relative source paths alongside their exact
    # resolved identities. Bind only that explicit path+SHA association; never
    # search by basename or infer a root from the process working directory.
    for source in (data.get("suite") or {}).get("sources", []):
        path, resolved, digest = source.get("path"), source.get("resolved_path"), source.get("sha256")
        if not all(isinstance(value, str) for value in (path, resolved, digest)):
            continue
        if hashes.get(resolved) != digest:
            continue
        if path in hashes and hashes[path] != digest:
            raise ValueError("Saved context source path has conflicting exact hash bindings")
        hashes[path] = digest
    if prior:
        if prior.get("input_signature") != signature:
            raise ValueError("Current-context review is stale after changed boundaries or decisions")
        return out
    evaluated = {_key(r) for r in evaluation.get("records", [])}
    merged = {}
    for group in original:
        sid = group["sample_id"]
        typ = group.get("rm_type") or group.get("family")
        keys = {(sid, m["gene_id"]) for m in group.get("members", [])}
        if not keys & evaluated:
            continue
        ident = (sid, group.get("run_id"), typ, tuple(sorted(keys)))
        required = _requirements(group, typ)
        if ident not in merged:
            merged[ident] = dict(group=group, group_ids=[], source_refs=[], required_roles=required)
        elif merged[ident]["required_roles"] != required:
            raise ValueError("Reciprocal saved boundaries have contradictory required roles")
        target = merged[ident]
        target["group_ids"].append(group["group_id"])
        for ref in _refs(group):
            if ref not in target["source_refs"]:
                target["source_refs"].append(ref)

    records = []
    for (sid, run, typ, keys), source in sorted(merged.items(), key=lambda item: repr(item[0])):
        group, required = source["group"], source["required_roles"]
        reasons, dispositions, current_roles, retained = [], [], set(), []
        refs = deepcopy(source["source_refs"])
        boundary_bound = bool(run and keys and typ in _TYPES and _bound(refs, hashes))
        if not boundary_bound:
            reasons.append("saved_boundary_identity_or_source_unbound")
        original_members = {m["gene_id"]: m for m in group.get("members", [])}
        if len(original_members) != len(group.get("members", [])):
            boundary_bound = False
            reasons.append("duplicate_original_boundary_member")
        roles_known, dispositions_known, original_R = True, True, False
        for key in keys:
            decision, gene = decisions.get(key, {}), genes.get(key, {})
            role, family = decision.get("role"), decision.get("rm_type") or decision.get("family")
            state = decision.get("state", "unassessed")
            protein = gene.get("protein_sha256")
            original_protein = original_members.get(key[1], {}).get("protein_sha256")
            identity_bound = bool(protein and original_protein == protein
                                  and decision.get("protein_sha256") == protein)
            decision_refs = _refs(decision)
            decision_bound = identity_bound and _bound(decision_refs, hashes)
            if not identity_bound:
                boundary_bound = False
                reasons.append("boundary_member_protein_identity_unbound")
            known = state in {"retained", "not_retained", "out_of_scope"} and decision_bound
            dispositions_known &= known
            role_set = _roles(role) if family == typ else set()
            roles_known &= bool(role_set)
            original_R |= "R" in role_set
            if state == "retained" and identity_bound and decision_bound:
                retained.append(dict(gene_id=key[1], protein_sha256=protein))
                current_roles |= role_set
            dispositions.append(dict(sample_id=key[0], gene_id=key[1], state=state,
                                     family=family, role=role, protein_sha256=protein,
                                     disposition_source_bound=bool(known), current_retained=state == "retained"))
            if _bound(decision_refs, hashes):
                for ref in decision_refs:
                    if ref not in refs:
                        refs.append(ref)
        if not dispositions_known:
            reasons.append("original_boundary_disposition_unassessed")
        if not roles_known:
            reasons.append("original_boundary_role_or_type_unresolved")
        if required is None:
            state = "unassessed"
            reasons.append("no_explicit_required_role_contract" if typ == "Type II" else
                           "type_iv_has_no_inferred_methylation_partner_requirement")
        elif not (boundary_bound and dispositions_known and roles_known):
            state = "unassessed"
        elif required <= current_roles:
            state = "met"
        else:
            state = "missing_retained_partner_in_saved_context"
        identity = [sid, run, typ, keys]
        records.append(dict(sample_id=sid, group_id="automatic_current_context:" + _digest(identity),
                            original_group_ids=sorted(source["group_ids"]), run_id=run, family=typ,
                            original_member_keys=[list(k) for k in keys], original_dispositions=dispositions,
                            retained_members=retained, retained_roles=sorted(current_roles),
                            required_roles=sorted(required) if required is not None else None,
                            context_state=state, evaluated_context_state=state, boundary_source_bound=boundary_bound,
                            dispositions_complete=bool(dispositions_known), roles_complete=bool(roles_known),
                            original_R_present=bool(original_R), reason_codes=sorted(set(reasons)), source_refs=refs))
    # Alternative overlapping boundaries cannot impose incompatible current
    # partner conclusions on the same query. Preserve both contracts in audit.
    by_member = defaultdict(list)
    for record in records:
        for member in record["retained_members"]:
            by_member[(record["sample_id"], member["gene_id"], record["family"])].append(record)
    ambiguous = set()
    for shared in by_member.values():
        states = {r["context_state"] for r in shared} - {"unassessed"}
        if len(states) > 1:
            ambiguous.update(r["group_id"] for r in shared)
    for record in records:
        if record["group_id"] in ambiguous:
            record["context_state"] = "unassessed"
            record["reason_codes"].append("competing_exact_saved_contexts")

    new_groups, gene_contexts = [], defaultdict(list)
    for record in records:
        if not record["boundary_source_bound"] or not record["retained_members"]:
            continue
        derived = dict(group_id=record["group_id"], sample_id=record["sample_id"], provider="Review",
                       run_id=record["run_id"], membership_kind="saved_context", family=record["family"],
                       context_state=record["context_state"], members=deepcopy(record["retained_members"]),
                       source_refs=deepcopy(record["source_refs"]), origin="automatic_current_context_review",
                       original_group_ids=record["original_group_ids"],
                       context_scope="current_retained_members_of_exact_saved_boundary", activity="unassessed")
        if record["required_roles"] is not None:
            derived["required_roles"] = record["required_roles"]
        new_groups.append(derived)
        for member in record["retained_members"]:
            gene_contexts[(record["sample_id"], member["gene_id"])].append(record)
    new_decisions = []
    for key, source in decisions.items():
        decision = deepcopy(source)
        contexts = gene_contexts.get(key, [])
        states = {r["context_state"] for r in contexts}
        if contexts and len(states) == 1 and "unassessed" not in states and "context_state" not in decision:
            decision["context_state"] = next(iter(states))
            decision["context_source_refs"] = _merge_refs(contexts)
            decision["context_review_origin"] = "automatic_current_context"
            # Scope orphanhood to a known original R-bearing context whose R
            # was removed. No absent global partner or singleton is invented.
            if "orphan" not in decision and _roles(decision.get("role")) & {"M", "S"}:
                if all("R" in c["retained_roles"] for c in contexts):
                    decision["orphan"] = False
                    decision["orphan_scope"] = "within_exact_retained_saved_context"
                elif all(c["original_R_present"] and "R" not in c["retained_roles"]
                         and c["dispositions_complete"] and c["roles_complete"] for c in contexts):
                    decision["orphan"] = True
                    decision["orphan_scope"] = "within_exact_retained_saved_context_not_genome_wide_absence"
        new_decisions.append(decision)
    out["saved_decisions"] = new_decisions
    out["groups"] = [*data.get("groups", []), *new_groups]
    out["rm_evaluation"] = dict(evaluation, current_context_review=dict(
        schema_version=1, operation="reassessment_of_exact_saved_boundaries_after_retention",
        input_signature=signature, records=records, reviewed_boundary_count=len(records),
        created_review_groups=len(new_groups), detector_runs=0, scoring_runs=0,
        new_boundaries_inferred=0, activity_inferred=False,
        original_groups_preserved=True, policy_decisions_changed=False))
    return out

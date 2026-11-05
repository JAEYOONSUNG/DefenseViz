"""Account for the automatic functional review of every observed defense CDS.

This stage consumes the exact component/complex evaluations already performed
by ``analyze_saved``. It does not run a detector, infer function from a family
name, or turn a successful review execution into complete functional evidence.
Positive membership is preserved; only explicit role exclusions can resolve a
negative effector endpoint. Raw InterPro/eggNOG annotations remain context until
a source-bound functional rule qualifies them.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy


_POSITIVE = {"supported", "accepted", "supported_annotation",
             "supported_putative_effector_annotation"}
_ENDPOINTS = {"effector": "effector", "associated": "associated",
              "motif": "motif_site_supported"}


def _key(row):
    return row["sample_id"], row["gene_id"]


def _endpoint(state, reason, evidence_ids=(), rule_ids=()):
    return {"status": state,
            "coverage": "complete" if state in {"supported", "not_applicable"} else "unassessed",
            "reason_codes": [reason], "evidence_ids": sorted(set(evidence_ids)),
            "rule_ids": sorted(set(rule_ids))}


def _explicit_exclusions(ledger, review, evaluations):
    """Exclude only an explicit query assignment, never a family-level fact."""
    rules = set()
    for row in ledger.get("annotations", []):
        if row.get("status") in _POSITIVE and row.get("assignment_kind") == "noncatalytic_component":
            rules.add(row["rule_id"])
    # Cas4/Cas6-type processing assignments are distinct from immune effectors.
    for row in ledger.get("effector_evidence", []):
        if row.get("status") in _POSITIVE and row.get("assignment_kind") == "nuclease_other_role":
            rules.add(row["rule_id"])
    ids = []
    if review.get("axes", {}).get("component_role", {}).get("status") != "conflicting":
        for row in evaluations:
            conclusion = row.get("conclusion", {})
            if (row.get("status") == "supported" and row.get("endpoint") == "component_role"
                    and (conclusion.get("role") == "noncatalytic"
                         or "nuclease" in conclusion.get("excludes_roles", [])
                         or "nuclease_effector" in conclusion.get("excludes_roles", []))):
                rules.add(row["rule_id"])
                ids.append(row["evaluation_id"])
    return rules, ids


def summarize_assessments(gene_assessments, member_keys=None, *, inventory_status="unassessed"):
    """Aggregate endpoint uncertainty over this exact CDS set, not a cohort.

    A positive CDS is classified even when its catalytic motif is unassessed.
    Any unresolved member keeps that endpoint partial/unassessed for this set.
    Counts describe the observed, integrated universe, never biological absence.
    """
    scope = None if member_keys is None else {tuple(k) for k in member_keys}
    rows = [r for r in gene_assessments if r["in_integrated_universe"]
            and (scope is None or _key(r) in scope)]
    found = {_key(r) for r in rows}
    if scope is not None and scope - found:
        raise ValueError("Automatic assessment scope contains an unaccounted integrated CDS")
    out = {"member_keys": [list(_key(r)) for r in sorted(rows, key=_key)],
           "CDS_count": len(rows), "metric_status": {}, "unresolved_CDS": {},
           "not_applicable_CDS": {}}
    for axis in _ENDPOINTS:
        states = Counter(r["endpoints"][axis]["status"] for r in rows)
        unresolved = sum(r["endpoints"][axis]["coverage"] != "complete" for r in rows)
        assessed = len(rows) - unresolved
        out[axis + "_CDS"] = states["supported"]
        out["unresolved_CDS"][axis] = unresolved
        out["not_applicable_CDS"][axis] = states["not_applicable"]
        out["metric_status"][axis] = ("complete" if not unresolved else
                                        "partial" if assessed or states["supported"] else "unassessed")
        if not rows and inventory_status != "complete":
            # Empty imported scope cannot establish absence when its original
            # source inventory is missing or partial.
            out["metric_status"][axis] = "unassessed"
    return out


def _membership_proposals(data, ledger, review):
    """Bridge only a rule explicitly establishing an individual immune effector.

    A supported nuclease domain or generic 'nuclease' role is intentionally
    insufficient. The existing review verifies source, protein and version
    binding. This returns proposals; the caller chooses whether to materialize
    them and must rerun normal canonical validation after doing so.
    """
    rules = {r["rule_id"]: r for r in (data.get("system_review_rules") or {}).get("rules", [])}
    source = data.get("system_review_rule_source")
    proposals = {}
    reviews = {_key(r): r for r in review.get("gene_reviews", [])}
    by_gene = defaultdict(list)
    for evaluation in review.get("rule_evaluations", []):
        if evaluation.get("gene_id") is not None:
            by_gene[_key(evaluation)].append(evaluation)
    for row in review.get("rule_evaluations", []):
        if row.get("gene_id") is None:
            continue
        key = _key(row)
        existing = ledger.get(key)
        rule = rules.get(row.get("rule_id"), {})
        if (not existing or existing.get("effector") or not existing.get("protein_sha256")
                or not source or row.get("status") != "supported"
                or row.get("route") != "native_component" or row.get("endpoint") != "component_role"
                or row.get("defense_scope") != "host_defense"
                or row.get("conclusion", {}).get("role") not in {"individual_nuclease_effector", "nuclease_effector"}
                or rule.get("conclusion") != row.get("conclusion")
                or reviews.get(key, {}).get("axes", {}).get("component_role", {}).get("status") == "conflicting"):
            continue
        excluded, _ = _explicit_exclusions(existing, reviews.get(key, {}), by_gene[key])
        if excluded:
            continue
        refs = row.get("version_provenance", {}).get("source_refs", [])
        if not refs:
            continue
        rid = "automatic_system_review:" + row["rule_id"]
        if any(e.get("rule_id") == rid for e in existing.get("effector_evidence", [])):
            continue
        proposal = dict(sample_id=key[0], gene_id=key[1], protein_sha256=existing["protein_sha256"],
                        rule_id=rid, assignment_kind="individual_nuclease_effector", status="supported",
                        evidence_grade="source_bound_exact_effector_component_rule",
                        source_refs=deepcopy([*refs, source]), primary_sources=deepcopy(row["primary_sources"]),
                        substrate=row["conclusion"].get("substrate", "unresolved"),
                        motif_state="unassessed", activity="unassessed", automatic=True,
                        review_evaluation_ids=[row["evaluation_id"]])
        ident = (key, rid)
        if ident in proposals:
            prior = proposals[ident]
            prior["review_evaluation_ids"] = sorted(set(prior["review_evaluation_ids"] + proposal["review_evaluation_ids"]))
            for ref in proposal["source_refs"]:
                if ref not in prior["source_refs"]:
                    prior["source_refs"].append(ref)
        else:
            proposals[ident] = proposal
    return [proposals[k] for k in sorted(proposals)]


def review_systems(data, preliminary_result):
    """Return a lossless per-CDS review and scoped count/uncertainty metadata.

    Inputs are not modified. Existing ``system_review`` is reused so raw
    annotations are not reparsed and exact component assignments are not
    duplicated. The caller can append ``proposals.effector_evidence`` and rerun
    ``analyze_saved`` once to apply explicitly qualified new component rules.
    """
    review = preliminary_result.get("system_review")
    if review is None:
        from .analysis.system_review import review_all_systems
        review = review_all_systems(data, preliminary_result, rules=data.get("system_review_rules"))
    ledger_rows = preliminary_result.get("gene_ledger", [])
    ledger = {_key(r): r for r in ledger_rows}
    gene_reviews = {_key(r): r for r in review.get("gene_reviews", [])}
    if len(ledger) != len(ledger_rows) or len(gene_reviews) != len(review.get("gene_reviews", [])):
        raise ValueError("Duplicate CDS identity in automatic system review")
    observed = {_key(r) for r in data.get("calls", [])}
    if not (observed | set(ledger)) <= set(gene_reviews):
        raise ValueError("Automatic system review lost an observed or integrated CDS")
    evaluations = defaultdict(list)
    for row in review.get("rule_evaluations", []):
        if row.get("gene_id") is not None:
            evaluations[_key(row)].append(row)
    rows = []
    for key, reviewed in sorted(gene_reviews.items()):
        accepted = ledger.get(key, {})
        axes = reviewed.get("axes", {})
        excluded, exclusion_ids = _explicit_exclusions(accepted, reviewed, evaluations[key])
        context = Counter((r.get("namespace") or r.get("table") or "unknown")
                          for r in reviewed.get("context_evidence", []))
        endpoints = {}
        for axis, flag in _ENDPOINTS.items():
            scientific_axis = "motif_site" if axis == "motif" else "component_role"
            state = axes.get(scientific_axis, {})
            ids = state.get("supporting_evidence_ids", [])
            rules = state.get("evaluated_rule_ids", [])
            if accepted.get(flag):
                reason = ("accepted_complete_complex_core" if axis == "effector" and accepted.get("complex_core")
                          and not accepted.get("individual_effector") else
                          "accepted_exact_effector_assignment" if axis == "effector" else
                          "explicit_selected_group_or_effector" if axis == "associated" else
                          "accepted_graded_motif_or_site_evidence")
                endpoints[axis] = _endpoint("supported", reason, ids, rules)
                if state.get("status") == "conflicting":
                    # Preserve accepted membership, while refusing to describe
                    # contradictory new evidence as a complete classification.
                    endpoints[axis]["coverage"] = "partial"
                    endpoints[axis]["reason_codes"].append("accepted_membership_preserved_despite_conflicting_review")
                    endpoints[axis]["evidence_ids"] = sorted(set(ids + state.get("conflicting_evidence_ids", [])))
            elif state.get("status") == "conflicting":
                endpoints[axis] = _endpoint("conflicting", "conflicting_query_functional_evidence",
                                             state.get("conflicting_evidence_ids", []), rules)
            elif excluded and axis == "effector":
                endpoints[axis] = _endpoint("not_applicable", "explicit_non_effector_role_assignment", exclusion_ids, excluded)
            else:
                reason = ("no_qualified_motif_or_site_assignment" if axis == "motif" else
                          "no_qualified_effector_seed_in_explicit_context" if axis == "associated" else
                          "component_identity_requires_complete_complex" if any(
                              e.get("assignment_kind") == "module_component_only" for e in accepted.get("effector_evidence", [])) else
                          "context_requires_qualified_functional_rule" if context else
                          "no_applicable_query_effector_rule")
                endpoints[axis] = _endpoint("unresolved", reason, ids, rules)
        rows.append(dict(sample_id=key[0], gene_id=key[1], protein_sha256=reviewed.get("protein_sha256"),
                         in_integrated_universe=key in ledger, review_execution="completed",
                         family_ids=list(accepted.get("family_ids", [])), unit_ids=list(reviewed.get("unit_ids", [])),
                         endpoints=endpoints, context_evidence_counts=dict(sorted(context.items())),
                         raw_context_inference="not_transferred_without_qualified_rule",
                         activity="unassessed", individual_catalysis=accepted.get("individual_catalysis", "unassessed")))
    families = defaultdict(list)
    samples = defaultdict(list)
    for row in rows:
        if row["in_integrated_universe"]:
            samples[row["sample_id"]].append(list(_key(row)))
            for family in row["family_ids"]:
                families[family].append(list(_key(row)))
    inventory = preliminary_result.get("availability", {}).get("inventory_by_sample", {})
    global_inventory = preliminary_result.get("global_summary", {}).get("metric_status", {}).get("inventory", "unassessed")
    summary = summarize_assessments(rows, inventory_status=global_inventory)
    summary.update(reviewed_CDS=len(rows), integrated_CDS=len(ledger),
                   outside_integrated_CDS=len(set(gene_reviews) - set(ledger)),
                   accounted_detection_records=review.get("accounting", {}).get("accounted_detection_count", len(data.get("calls", []))),
                   review_execution="completed", classification_scope="observed_source_bound_annotations",
                   exhaustive_functional_discovery=False)
    return dict(schema_version=1, analysis_kind="automatic_saved_system_evidence_review",
                gene_assessments=rows, summary=summary,
                counts_by_family=[dict(family_id=f, **summarize_assessments(rows, keys)) for f, keys in sorted(families.items())],
                counts_by_sample=[dict(sample_id=s["sample_id"], **summarize_assessments(rows, samples[s["sample_id"]], inventory_status=inventory.get(s["sample_id"], "unassessed")))
                                  for s in data.get("samples", [])],
                proposals={"effector_evidence": _membership_proposals(data, ledger, review)},
                provenance=dict(detector_runs=0, scoring_runs=0,
                                component_registry_sha256=preliminary_result.get("provenance", {}).get("component_registry_sha256"),
                                system_review_rules_sha256=review.get("provenance", {}).get("rules_sha256")),
                policy="Review execution is complete; unresolved functional evidence remains unresolved. Native calls and accepted membership are unchanged.")

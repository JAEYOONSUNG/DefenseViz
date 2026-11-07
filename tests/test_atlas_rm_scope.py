"""Current observed assessment must not inherit unrelated historical gaps."""
from copy import deepcopy
import hashlib

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.plot.publication.manuscript_model import enrich_manuscript_model

REF = dict(path="synthetic_saved.json", sha256="a" * 64)


def bundle():
    genes = [dict(sample_id="sample", gene_id=g, protein_sha256=hashlib.sha256(g.encode()).hexdigest())
             for g in ["r", "m", "s", "historical", "iv"]]
    records = {g["gene_id"]: g for g in genes}
    data = dict(samples=[dict(sample_id="sample", label="Synthetic")], genes=genes,
                calls=[], coverage=[], saved_decisions=[], annotations=[], groups=[],
                options={"manuscript_atlas": {"group_coverage": {"sample": "complete"}}})
    for caller in ["NativeA", "NativeB", "NativeC"]:
        data["coverage"].append(dict(sample_id="sample", caller=caller, endpoint="native_calls", status="complete"))
        for gid in ["r", "m", "s"]:
            data["calls"].append(dict(**records[gid], caller=caller, evidence_kind="final_member",
                                     family="Type I", rm_type="Type I", category="R-M", source_ref=REF))
    data["coverage"].append(dict(sample_id="sample", caller="Reference", endpoint="reference_matches", status="complete"))
    data["calls"].append(dict(**records["r"], caller="Reference", evidence_kind="reference_hit",
                             family="Type I", rm_type="Type I", category="R-M", source_ref=REF))
    for gid in ["r", "m", "s"]:
        data["saved_decisions"].append(dict(**records[gid], state="retained", family="Type I", role=gid.upper(), source_ref=REF))
    data["saved_decisions"].append(dict(**records["historical"], state="unassessed", family=None, role=None, source_ref=REF))
    data["coverage"].append(dict(sample_id="sample", caller="DNMB", endpoint="saved_retention", status="partial"))
    data["rm_evaluation"] = dict(candidate_scope=[dict(sample_id="sample", gene_ids=["r", "m", "s"], status="partial")])
    data["groups"].append(dict(sample_id="sample", group_id="current", provider="Review", run_id="review-run",
                               membership_kind="saved_context", family="Type I", context_state="complete",
                               members=[dict(gene_id=g, protein_sha256=records[g]["protein_sha256"]) for g in ["r", "m", "s"]], source_refs=[REF]))
    return data


def project(data):
    result = enrich_manuscript_model(data, analyze_saved(data))
    return result, result["manuscript_atlas"]["rows"][0]


def test_complete_current_selection_preserves_broad_partial_audit_and_all_counts():
    data = bundle(); before = deepcopy(data)
    result, row = project(data)
    assert data == before
    assert result["availability"]["saved_decisions_by_sample"][0]["status"] == "partial"
    current = result["availability"]["current_saved_decisions_by_sample"][0]
    assert current["status"] == "complete" and current["record_count"] == 3
    assert current["exhaustive_discovery_claim"] is False
    assert result["sample_summary"][0]["retained_counts_status"] == "partial"
    assert row["retained"]["R_CDS"] == 1 and row["retained"]["total_CDS"] == 3
    assert row["retained"]["R_status"] == row["retained"]["total_status"] == row["retained"]["status"] == "complete"
    assert row["rm_original"]["input_pool_CDS"] == row["rm_original"]["intersection_CDS"] == 3
    assert row["rm_original"]["pool_status"] == row["rm_original"]["intersection_count_status"] == "complete"
    assert all(row["rm_roles"]["I_"+role]["status"] == "complete" for role in ["R", "M", "S"])
    assert row["compositions"]["I_RMS"]["count"] == 1
    assert row["compositions"]["I_RMS"]["status"] == "complete"
    assert result["sample_summary"][0]["rm_pool_broad_assessment_status"] == "partial"
    without_history = deepcopy(data)
    without_history["saved_decisions"].pop()
    clean, _ = project(without_history)
    assert result["gene_ledger"] == clean["gene_ledger"]
    assert result["source_sets"] == clean["source_sets"]
    assert result["global_summary"] == clean["global_summary"]


@pytest.mark.parametrize("defect", ["unassessed", "missing_decision", "unknown_role", "unknown_type"])
def test_unknown_within_current_retention_is_not_promoted(defect):
    data = bundle()
    if defect == "missing_decision": data["saved_decisions"].pop(0)
    elif defect == "unassessed": data["saved_decisions"][0]["state"] = "unassessed"
    elif defect == "unknown_role": data["saved_decisions"][0]["role"] = None
    else: data["saved_decisions"][0]["family"] = None
    result, row = project(data)
    assert row["retained"]["R_status"] == "partial"
    assert row["rm_roles"]["I_R"]["status"] == "partial"
    assert row["rm_original"]["pool_status"] == "complete"  # the source union itself is bound
    if defect == "unknown_role":
        assert row["retained"]["total_status"] == "complete"
    else:
        assert row["retained"]["total_status"] == "partial"


def test_unassessed_sequence_in_current_scope_preserves_known_typed_membership():
    data = bundle(); data["saved_decisions"][0].update(state="unassessed", decision_reason="invalid_or_missing_sequence")
    result, row = project(data)
    assert row["retained"]["total_status"] == "partial"
    assert row["rm_original"]["pool_status"] == "complete"
    assert result["sample_summary"][0]["rm_candidate_membership_status"] == "complete"


@pytest.mark.parametrize("known_out_of_scope", [False, True])
def test_untyped_current_member_is_not_assumed_classified(known_out_of_scope):
    data = bundle(); data["rm_evaluation"]["candidate_scope"][0]["gene_ids"].append("historical")
    if known_out_of_scope: data["saved_decisions"][-1]["state"] = "out_of_scope"
    result, row = project(data)
    assert row["rm_original"]["input_pool_CDS"] == 3
    assert row["rm_original"]["pool_status"] == ("complete" if known_out_of_scope else "partial")
    assert row["retained"]["total_status"] == ("complete" if known_out_of_scope else "partial")
    assert data["saved_decisions"][-1]["family"] is None


@pytest.mark.parametrize("provider,state", [("Reference", "missing"), ("Reference", "partial"), ("NativeC", "missing"), ("NativeC", "partial")])
def test_actual_source_gaps_remain_partial_with_separate_native_numerator(provider, state):
    data = bundle()
    for c in data["coverage"]:
        if c["caller"] == provider: c["status"] = state
    if state == "missing": data["calls"] = [c for c in data["calls"] if c["caller"] != provider]
    _, row = project(data)
    assert row["rm_original"]["pool_status"] == "partial"
    assert row["rm_original"]["intersection_count_status"] == ("complete" if provider == "Reference" else "partial" if state == "partial" else "unassessed")
    assert row["retained"]["total_status"] == "complete"


def test_no_explicit_scope_keeps_broader_retention_state():
    data = bundle(); del data["rm_evaluation"]
    _, row = project(data)
    assert row["retained"]["total_status"] == "partial"
    assert row["rm_roles"]["I_R"]["status"] == "partial"


def test_typeIV_unknown_review_and_unresolved_native_type_pool_are_not_zeroed():
    data = bundle(); gene = next(g for g in data["genes"] if g["gene_id"] == "iv")
    data["calls"].append(dict(**gene, caller="NativeA", evidence_kind="final_member", family="Type IV", rm_type="Type IV", category="R-M", source_ref=REF))
    data["annotations"].append(dict(**gene, status="unassessed", assignment_kind="reviewed_RM_component",
                                    family="Type IV", role="R", include_in_union=False, rule_id="unknown_review", source_refs=[REF]))
    result, row = project(data)
    assert row["rm_roles"]["IV_R"]["count"] == 0
    assert row["rm_roles"]["IV_R"]["status"] == "unassessed"
    typed = next(t for t in result["type_tables"] if t["rm_type"] == "Type IV")
    assert typed["pool_CDS"] == 1 and typed["within_pool_counts_status"] != "complete"


@pytest.mark.parametrize("covered", [False, True])
def test_unresolved_retained_context_requires_current_bound_review(covered):
    data = bundle()
    if not covered: data["groups"][0]["provider"] = "DNMB"
    data["suite"] = {"unresolved_contexts": [dict(sample_id="sample", gene_id="r", family="Type I", source_ref=REF)]}
    _, row = project(data)
    assert row["compositions"]["I_RMS"]["status"] == ("complete" if covered else "partial")
    assert row["compositions"]["III_MR"]["status"] == "complete"


def test_rejected_and_outscope_unresolved_context_cannot_change_current_counts():
    data = bundle(); data["suite"] = {"unresolved_contexts": [dict(sample_id="sample", gene_id="historical", family="Type I", source_ref=REF)]}
    _, row = project(data)
    assert row["compositions"]["I_RMS"]["count"] == 1
    assert row["compositions"]["I_RMS"]["status"] == "complete"


def test_displaced_context_boundaries_are_preserved_without_current_status_poisoning():
    data = bundle(); old = deepcopy(data["groups"][0]); old.update(group_id="historical-boundary", provider="DNMB", run_id="old")
    data["groups"].append(old)
    data["rm_evaluation"].update(decision_mode="replace_scoped")
    result, row = project(data)
    groups = result["manuscript_atlas"]["group_records"]
    assert len(groups) == 2 and sum(g["current_composition"] for g in groups) == 1
    assert row["compositions"]["I_RMS"]["count"] == 1
    assert row["compositions"]["I_RMS"]["status"] == "complete"
    assert next(g for g in groups if g["group_id"] == "historical-boundary")["context_state"] == "unassessed"


@pytest.mark.parametrize("selected", [[], ["r"], ["historical"]])
def test_retained_display_scope_explicitly_includes_preserved_positive_decisions(selected):
    data = bundle(); data["rm_evaluation"]["candidate_scope"][0]["gene_ids"] = selected
    result, row = project(data); summary = result["sample_summary"][0]
    current = result["availability"]["current_saved_decisions_by_sample"][0]
    assert current["record_count"] == len(selected)
    assert current["status"] == ("unassessed" if selected == ["historical"] else "complete")
    assert row["retained"]["total_CDS"] == 3
    assert row["retained"]["total_status"] == ("partial" if selected == ["historical"] else "complete")
    assert row["retained"]["scope"] == "current_policy_selected_or_reviewed_RM_components"
    assert "plus explicitly retained active decisions outside" in summary["current_retention_scope"]
    assert set(map(tuple, summary["current_retention_assessment_member_keys"])) == {
        ("sample", g) for g in set(selected) | {"r", "m", "s"}}
    assert set(map(tuple, summary["current_retained_outside_candidate_member_keys"])) == {
        ("sample", g) for g in {"r", "m", "s"} - set(selected)}


def test_unknown_role_on_preserved_retained_positive_stays_partial():
    data = bundle(); data["rm_evaluation"]["candidate_scope"][0]["gene_ids"] = []
    data["saved_decisions"][0]["role"] = None
    _, row = project(data)
    assert row["retained"]["total_CDS"] == 3 and row["retained"]["total_status"] == "complete"
    assert row["retained"]["R_status"] == row["rm_roles"]["I_R"]["status"] == "partial"


@pytest.mark.parametrize("field", ["role", "family"])
def test_validated_review_resolves_missing_retained_assignment(field):
    data = bundle(); data["saved_decisions"][0][field] = None
    gene = next(g for g in data["genes"] if g["gene_id"] == "r")
    data["annotations"].append(dict(**gene, status="supported", assignment_kind="reviewed_RM_component",
                                    family="Type I", role="R", include_in_union=True, rule_id="resolved_assignment", source_refs=[REF]))
    _, row = project(data)
    assert row["rm_roles"]["I_R"]["count"] == 1
    assert row["retained"]["R_status"] == row["rm_roles"]["I_R"]["status"] == "complete"


def test_result_only_fallback_uses_current_retained_status_fields():
    from defenseviz.plot.publication.manuscript_model import build_manuscript_atlas_model
    result = analyze_saved(bundle())
    assert result["sample_summary"][0]["retained_ratio_status"] == "partial"
    row = build_manuscript_atlas_model(result)["rows"][0]
    assert row["retained"]["R_status"] == row["retained"]["total_status"] == row["retained"]["status"] == "complete"


def test_empty_typeIV_pool_does_not_erase_explicit_unknown_review():
    data = bundle(); gene = next(g for g in data["genes"] if g["gene_id"] == "iv")
    data["annotations"].append(dict(**gene, status="unassessed", assignment_kind="reviewed_RM_component",
                                    family="Type IV", role="R", include_in_union=False, rule_id="unknown_review", source_refs=[REF]))
    other = dict(sample_id="other", gene_id="typed", protein_sha256="b" * 64)
    data["samples"].append(dict(sample_id="other", label="Another synthetic sample"))
    data["genes"].append(other)
    data["calls"].append(dict(**other, caller="NativeA", evidence_kind="final_member", family="Type IV", rm_type="Type IV", category="R-M", source_ref=REF))
    data["coverage"].append(dict(sample_id="other", caller="NativeA", endpoint="native_calls", status="complete"))
    data["saved_decisions"].append(dict(**other, state="not_retained", family="Type IV", role="R", source_ref=REF))
    result, row = project(data)
    typed = next(t for t in result["type_tables"] if t["sample_id"] == "sample" and t["rm_type"] == "Type IV")
    assert typed["pool_CDS"] == 0 and typed["within_pool_counts_status"] == "complete"
    assert row["rm_roles"]["IV_R"]["status"] == "unassessed"


def test_exported_intersections_are_independent_of_unclassified_pool_members():
    data = bundle(); data["rm_evaluation"]["candidate_scope"][0]["gene_ids"].append("historical")
    result, row = project(data)
    summary = result["sample_summary"][0]
    assert summary["rm_pool_CDS"] == 3 and summary["rm_pool_status"] == "partial"
    assert summary["rm_intersection_CDS"] == summary["rm_all_source_intersection_CDS"] == 1
    assert summary["rm_intersection_status"] == summary["rm_all_source_intersection_status"] == "complete"
    assert summary["rm_intersection_sources"] == ["NativeA", "NativeB", "NativeC", "Reference"]
    assert summary["rm_native_intersection_CDS"] == row["rm_original"]["intersection_CDS"] == 3
    assert summary["rm_native_intersection_status"] == "complete"
    assert summary["rm_native_intersection_sources"] == ["NativeA", "NativeB", "NativeC"]
    assert row["rm_original"]["intersection_status"] == "partial"


def test_native_intersection_projects_only_generic_predictions_and_excludes_profiles():
    data = bundle()
    for c in data["calls"]:
        if c["caller"] == "NativeC":
            c.update(family="General defense", category="Unresolved", rm_type=None)
    gene = next(g for g in data["genes"] if g["gene_id"] == "historical")
    data["calls"].append(dict(**gene, caller="NativeC", evidence_kind="final_member", family="General defense", source_ref=REF))
    data["calls"].append(dict(**gene, caller="Profiles", evidence_kind="hmm_hit", family="Type I", rm_type="Type I", source_ref=REF))
    data["coverage"].append(dict(sample_id="sample", caller="Profiles", endpoint="hmm_hits", status="complete"))
    result, row = project(data); summary = result["sample_summary"][0]
    assert summary["rm_native_intersection_CDS"] == row["rm_original"]["intersection_CDS"] == 3
    assert summary["rm_native_intersection_status"] == "complete"
    assert summary["rm_native_intersection_sources"] == ["NativeA", "NativeB", "NativeC"]
    # Compatibility: the old configured-source endpoint still includes Profiles.
    assert summary["rm_intersection_CDS"] == 0 and "Profiles" in summary["rm_intersection_sources"]
    assert summary["rm_intersection_status"] == "unassessed"


@pytest.mark.parametrize("provider,state", [("Reference", "partial"), ("Reference", "missing"),
                                             ("NativeC", "partial"), ("NativeC", "missing")])
def test_exported_native_and_reference_intersections_keep_separate_coverage(provider, state):
    data = bundle()
    for c in data["coverage"]:
        if c["caller"] == provider: c["status"] = state
    if state == "missing": data["calls"] = [c for c in data["calls"] if c["caller"] != provider]
    result, row = project(data); summary = result["sample_summary"][0]
    assert summary["rm_intersection_status"] == ("partial" if state == "partial" else "unassessed")
    expected = "complete" if provider == "Reference" else "partial" if state == "partial" else "unassessed"
    assert summary["rm_native_intersection_status"] == row["rm_original"]["intersection_count_status"] == expected


def test_mixed_provider_reference_does_not_fill_native_intersection_or_coverage():
    data = bundle()
    data["calls"] = [c for c in data["calls"] if not (c["caller"] == "NativeC" and c["gene_id"] == "m")]
    gene = next(g for g in data["genes"] if g["gene_id"] == "m")
    data["calls"].append(dict(**gene, caller="NativeC", evidence_kind="reference_hit", family="Type I", rm_type="Type I", source_ref=REF))
    data["coverage"] = [c for c in data["coverage"] if c["caller"] != "NativeC"]
    data["coverage"].append(dict(sample_id="sample", caller="NativeC", endpoint="reference_matches", status="complete"))
    result, row = project(data); summary = result["sample_summary"][0]
    assert summary["rm_native_intersection_CDS"] == row["rm_original"]["intersection_CDS"] == 2
    assert summary["rm_native_intersection_status"] == row["rm_original"]["intersection_count_status"] == "partial"
    from defenseviz.plot.publication.manuscript_model import build_manuscript_atlas_model
    fallback = build_manuscript_atlas_model(analyze_saved(data))["rows"][0]
    assert fallback["rm_original"]["intersection_CDS"] is None


@pytest.mark.parametrize("field,bad_value", [("rm_native_intersection_CDS", 99),
                                            ("rm_native_intersection_status", "partial")])
def test_atlas_crosschecks_same_native_endpoint_export(field, bad_value):
    data = bundle(); result = analyze_saved(data)
    result["sample_summary"][0][field] = bad_value
    with pytest.raises(ValueError, match="native intersection"):
        enrich_manuscript_model(data, result)
    # A deliberately configured subset is a different endpoint, so compute it
    # from its own native rows without comparing an all-native summary count.
    data["options"]["manuscript_atlas"]["source_order"] = ["NativeA", "NativeB"]
    assert enrich_manuscript_model(data, result)["manuscript_atlas"]["rows"][0]["rm_original"]["intersection_CDS"] == 3


@pytest.mark.parametrize("partial", [False, True])
def test_result_only_analysis_to_fraction_preserves_counts_and_audited_coverage(partial):
    from defenseviz.plot.publication.drawing import Canvas, Fonts
    from defenseviz.plot.publication.manuscript_atlas import _fraction, _fraction_layout
    from defenseviz.plot.publication.manuscript_model import build_manuscript_atlas_model
    data = bundle()
    if partial:
        data["coverage"][2]["status"] = "partial"
        data["saved_decisions"][0]["state"] = "unassessed"
    result = analyze_saved(data)
    row = build_manuscript_atlas_model(result)["rows"][0]
    assert result["sample_summary"][0]["retained_counts_status"] == "partial"
    fonts = Fonts(); canvas = Canvas(80, 30, fonts, "scope-fallback", "synthetic")
    for y, name, record, n, den, ns, ds in [
        (8, "intersection", row["rm_original"], "intersection_CDS", "input_pool_CDS", "intersection_count_status", "pool_status"),
        (16, "retained", row["retained"], "R_CDS", "total_CDS", "R_status", "total_status")]:
        layout = _fraction_layout(fonts, [(record[n], record[den], record[ns], record[ds])], name)
        _fraction(canvas, 5, y, layout["width_mm"], record[n], record[den], record[ns], name, 3, record[ds], layout)
    canvas.validate()
    labels = {name:[t["text"] for t in canvas.ledger["texts"] if (t.get("source_ref") or {}).get("endpoint") == name]
              for name in ["intersection", "retained"]}
    assert labels["intersection"] == ["3", "/", "3"]
    assert labels["retained"] == (["0", "/", "2"] if partial else ["1", "/", "3"])

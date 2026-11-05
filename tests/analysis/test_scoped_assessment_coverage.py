"""Functional review coverage must follow declared sample scope, not cohort labels."""
from copy import deepcopy
import pytest
from defenseviz.analysis.saved import analyze_saved
from tests.analysis.test_saved import REF, gene


def project():
    genes = [gene("g", s) for s in ("reviewed", "new")]
    return dict(samples=[dict(sample_id="reviewed", label="A", cohort="external"),
                         dict(sample_id="new", label="B", cohort="primary")],
        genes=genes, calls=[dict(sample_id=s, gene_id="g", caller="Tool", family="Unknown system",
            category="Unresolved", evidence_kind="final_member", source_ref=REF) for s in ("reviewed", "new")],
        coverage=[dict(sample_id="reviewed", caller="Tool", status="complete"),
                  dict(sample_id="new", caller="Tool", status="partial")],
        effector_evidence=[dict(**genes[0], status="supported", assignment_kind="individual_nuclease_effector",
            rule_id="specific_review", source_refs=[REF])], options={})


def sample_rows(result):
    return {r["sample_id"]: r for r in result["sample_summary"]}


def test_local_defaults_do_not_leak_positive_assessment_or_inventory():
    data = project()
    result = analyze_saved(data)
    rows = sample_rows(result)
    assert rows["reviewed"]["metric_status"] == dict(inventory="complete", effector="partial", associated="partial", motif="unassessed")
    assert rows["new"]["metric_status"] == dict(inventory="partial", effector="unassessed", associated="unassessed", motif="unassessed")
    for cell in result["cells"]:
        assert cell["metric_status"] == rows[cell["sample_id"]]["metric_status"]
        assert cell["counts_status"] == rows[cell["sample_id"]]["metric_status"]["inventory"]
    assert result["global_summary"]["metric_status"]["inventory"] == "partial"


def test_explicit_sample_axes_override_global_fallback_and_preserve_counts():
    data = project()
    baseline = analyze_saved(data)
    data["options"] = dict(assessment_coverage=dict(effector="partial", motif="complete"),
        assessment_coverage_by_sample={"reviewed": dict(effector="complete", associated="complete"),
                                       "new": dict(motif="unassessed")})
    result = analyze_saved(data)
    rows = sample_rows(result)
    assert rows["reviewed"]["metric_status"]["effector"] == "complete"
    assert rows["new"]["metric_status"]["effector"] == "partial"
    assert rows["new"]["metric_status"]["motif"] == "unassessed"
    assert rows["new"]["metric_status"]["associated"] == "unassessed"
    for name in ("gene_ledger", "priority_ledger", "groups_selected"):
        assert result[name] == baseline[name]
    assert result["global_summary"]["effector_CDS"] == baseline["global_summary"]["effector_CDS"]
    proof = result["availability"]["assessment_coverage_provenance_by_sample"]
    assert proof["reviewed"]["effector"]["origin"] == "sample_declaration"
    assert proof["new"]["effector"]["origin"] == "explicit_global_fallback"
    assert proof["new"]["associated"]["origin"] == "local_observed_evidence_default"


def test_declaration_provenance_is_sample_and_axis_bound():
    data = project()
    data["options"]["assessment_coverage_by_sample"] = {"reviewed": {"effector": "complete"}}
    data["source_metadata"] = [dict(sample_id="reviewed", provider="DNMB", source_kind="assessment_coverage_declaration",
        declaration_scope="sample", axes={"effector": "complete"}, source_ref=REF)]
    proof = analyze_saved(data)["availability"]["assessment_coverage_provenance_by_sample"]
    assert proof["reviewed"]["effector"]["source_refs"] == [REF]
    assert proof["reviewed"]["associated"]["source_refs"] == []
    assert proof["new"]["effector"]["source_refs"] == []


@pytest.mark.parametrize("options", [
    {"assessment_coverage_by_sample": {"not_a_sample": {"effector": "complete"}}},
    {"assessment_coverage_by_sample": []},
    {"assessment_coverage_by_sample": {"new": []}},
    {"assessment_coverage_by_sample": {"new": {"activity": "complete"}}},
    {"assessment_coverage_by_sample": {"new": {"effector": "negative"}}},
    {"assessment_coverage_by_sample": {"new": {"effector": []}}},
    {"assessment_coverage": {"typo": "complete"}},
    {"assessment_coverage": {"motif": False}},
])
def test_malformed_coverage_scopes_rejected(options):
    data = project(); data["options"] = options
    with pytest.raises(ValueError, match="assessment_coverage"):
        analyze_saved(data)


def test_cohort_label_or_sample_order_does_not_decide_evidence_status():
    data = project()
    before = sample_rows(analyze_saved(data))
    data["samples"].reverse()
    for sample in data["samples"]:
        sample["cohort"] = "unrelated_label"
    after = sample_rows(analyze_saved(data))
    assert {s:r["metric_status"] for s,r in before.items()} == {s:r["metric_status"] for s,r in after.items()}

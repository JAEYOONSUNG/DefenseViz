"""Integration controls for automatic exact-member functional coverage."""
from copy import deepcopy

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.automatic_systems_output import apply_system_review
from tests.analysis.test_saved import base, call, cov
from tests.analysis.test_scoped_assessment_coverage import project
from tests.test_automatic_systems import rule_data


def counts(value):
    """Snapshot numeric scientific endpoints, ignoring only status provenance."""
    if isinstance(value, dict):
        return {k: counts(v) for k, v in value.items()
                if k not in {"automatic_system_review", "automatic_assessment", "availability", "system_review", "metric_status"}}
    if isinstance(value, list):
        return [counts(v) for v in value]
    return value


def test_exact_family_scope_updates_without_numeric_or_native_changes():
    data = base("r", "unknown")
    data["calls"] = [call("r", family="RloC", model="RloC__RloC"), call("unknown", family="Future")]
    data["coverage"] = cov("DefenseFinder")
    result = analyze_saved(data)
    original_data, baseline = deepcopy(data), counts(result)
    output_data, output = apply_system_review(data, result)
    assert output_data == original_data
    assert counts(output) == baseline
    rloc = next(r for r in output["inventory_rows"]["families"] if r["id"] == "RloC")
    assert rloc["metric_status"]["effector"] == "complete"
    assert rloc["metric_status"]["motif"] == "unassessed"
    assert output["global_summary"]["metric_status"]["effector"] == "partial"
    assert output["automatic_system_review"]["materialization"]["canonical_reanalysis_runs"] == 0


def test_sample_scopes_do_not_inherit_unrelated_unknowns_and_inventory_unchanged():
    data = project()
    _, output = apply_system_review(data, analyze_saved(data))
    sample = {r["sample_id"]: r for r in output["sample_summary"]}
    assert sample["reviewed"]["metric_status"]["effector"] == "complete"
    assert sample["new"]["metric_status"]["effector"] == "unassessed"
    assert sample["new"]["metric_status"]["inventory"] == "partial"
    assert output["availability"]["assessment_coverage_by_sample"]["reviewed"]["effector"] == "complete"
    assert output["availability"]["assessment_coverage_provenance_by_sample"]["reviewed"]["effector"]["origin"] == "automatic_exact_member_review"


@pytest.mark.parametrize("declaration", [
    {"assessment_coverage": {"effector": "complete"}},
    {"assessment_coverage_by_sample": {"reviewed": {"motif": "partial"}}},
])
def test_explicit_coverage_declarations_preserve_original_statuses(declaration):
    data = project()
    data["options"] = declaration
    result = analyze_saved(data)
    before = deepcopy(result)
    _, output = apply_system_review(data, result)
    assert output["automatic_system_review"]["coverage_application"]["status"] == "explicit_declarations_preserved"
    for field in ("global_summary", "sample_summary", "inventory_rows", "cells", "availability"):
        assert output[field] == before[field]


def test_positive_proposals_materialize_once_through_canonical_analysis(tmp_path):
    data, result = rule_data(tmp_path, "nuclease_effector")
    native = deepcopy(data["calls"])
    data, result = apply_system_review(data, result)
    assert data["calls"] == native
    assert result["global_summary"]["effector_CDS"] == 1
    assert result["automatic_system_review"]["materialization"]["canonical_reanalysis_runs"] == 1
    assert result["automatic_system_review"]["materialization"]["new_effector_evidence_records"] == 1
    data, again = apply_system_review(data, result)
    assert len(data["effector_evidence"]) == 1
    assert again["automatic_system_review"]["materialization"]["canonical_reanalysis_runs"] == 0


def test_scope_count_disagreement_fails_instead_of_hiding_bad_join():
    data = base("r")
    data["calls"] = [call("r", family="RloC", model="RloC__RloC")]
    result = analyze_saved(data)
    result["cells"][0]["effector_CDS"] += 1
    with pytest.raises(ValueError, match="disagrees"):
        apply_system_review(data, result)


def test_missing_exact_member_scope_fails_closed():
    data = base("r")
    data["calls"] = [call("r", family="RloC", model="RloC__RloC")]
    result = analyze_saved(data)
    result["cells"][0].pop("member_keys")
    with pytest.raises(ValueError, match="lacks exact"):
        apply_system_review(data, result)


@pytest.mark.parametrize("inventory,expected", [("complete", "complete"), ("partial", "unassessed"), ("missing", "unassessed"), ("not_run", "unassessed")])
def test_empty_scope_requires_complete_original_inventory(inventory, expected):
    data = base()
    data["coverage"] = [dict(sample_id="s", caller="DefenseFinder", status=inventory)]
    data["options"] = {"family_catalog": [dict(id="Empty", label="Empty", category="Unresolved")]}
    _, result = apply_system_review(data, analyze_saved(data))
    for row in [result["global_summary"], *result["sample_summary"], *result["cells"],
                *result["inventory_rows"]["families"], *result["inventory_rows"]["categories"]]:
        assert {row["metric_status"][axis] for axis in ("effector", "associated", "motif")} == {expected}
    assert result["automatic_system_review"]["summary"]["metric_status"]["effector"] == expected

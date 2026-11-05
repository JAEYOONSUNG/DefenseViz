"""Functional and identity controls for the automatic system-review stage."""
from copy import deepcopy
import hashlib
import json

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.analysis.system_review import review_all_systems
from defenseviz.automatic_systems import review_systems, summarize_assessments
from tests.analysis.test_saved import base, call, cov, effect, gene, group, REF
from tests.analysis.test_system_review import saved, component, registry


def assessment(report, gid):
    return next(r for r in report["gene_assessments"] if r["gene_id"] == gid)


def test_jetd_exact_positive_does_not_promote_wadjet_accessory_or_motif():
    data = base("d", "a")
    data["calls"] = [call("d", family="Wadjet", model="Wadjet__JetD_II", system="w", run_id="r"),
                     call("a", family="Wadjet", model="Wadjet__JetA_II", system="w", run_id="r")]
    data["coverage"] = cov("DefenseFinder")
    result = analyze_saved(data)
    before = deepcopy((data, result))
    report = review_systems(data, result)
    assert (data, result) == before
    assert report["summary"]["effector_CDS"] == 1
    assert assessment(report, "d")["endpoints"]["effector"]["status"] == "supported"
    assert assessment(report, "a")["endpoints"]["effector"]["status"] == "unresolved"
    assert assessment(report, "d")["endpoints"]["motif"]["status"] == "unresolved"
    assert report["summary"]["metric_status"]["effector"] == "partial"
    assert report["summary"]["metric_status"]["associated"] == "complete"
    assert report["proposals"]["effector_evidence"] == []


def test_gene_scope_is_independent_of_unreviewed_other_families():
    data = base("d", "future")
    data["calls"] = [call("d", family="Wadjet", model="Wadjet__JetD_II"), call("future", family="Novel")]
    data["coverage"] = cov("DefenseFinder")
    report = review_systems(data, analyze_saved(data))
    exact = summarize_assessments(report["gene_assessments"], [["s", "d"]])
    assert exact["metric_status"]["effector"] == "complete"
    assert exact["metric_status"]["motif"] == "unassessed"
    assert report["summary"]["metric_status"]["effector"] == "partial"
    assert assessment(report, "future")["review_execution"] == "completed"


@pytest.mark.parametrize("model", ["unknown_HNH_nuclease", "Wadjet__JetD_IV", "Wadjet__JetD_II_extra"])
def test_near_model_names_and_family_function_do_not_transfer(model):
    data = base("g")
    data["calls"] = [call("g", family="Wadjet nuclease", model=model)]
    report = review_systems(data, analyze_saved(data))
    assert report["summary"]["effector_CDS"] == 0
    assert assessment(report, "g")["endpoints"]["effector"]["status"] == "unresolved"


def test_hmm_only_observation_is_accounted_without_invented_membership():
    data = base("g")
    data["calls"] = [call("g", model="Wadjet__JetD_II", kind="hmm_hit")]
    report = review_systems(data, analyze_saved(data))
    assert report["summary"]["integrated_CDS"] == 0
    assert report["summary"]["outside_integrated_CDS"] == 1
    assert not assessment(report, "g")["in_integrated_universe"]
    assert report["proposals"]["effector_evidence"] == []


def test_processing_nuclease_is_not_immune_effector():
    data = base("c")
    data["calls"] = [call("c", model="cas4_I_II_III_IV_V_VI_1")]
    report = review_systems(data, analyze_saved(data))
    state = assessment(report, "c")["endpoints"]["effector"]
    assert state["status"] == "not_applicable"
    assert state["rule_ids"] == ["rna_component_3"]
    assert report["summary"]["effector_CDS"] == 0
    assert report["summary"]["metric_status"]["effector"] == "complete"


def test_explicit_noncatalytic_query_review_is_preserved():
    data = base("regulator")
    data["calls"] = [call("regulator")]
    data["annotations"] = [dict(**gene("regulator"), assignment_kind="noncatalytic_component", status="supported",
                                rule_id="specific_query_review", source_refs=[REF])]
    report = review_systems(data, analyze_saved(data))
    assert assessment(report, "regulator")["endpoints"]["effector"]["status"] == "not_applicable"
    # A role exclusion does not establish missing motifs or a missing system.
    assert assessment(report, "regulator")["endpoints"]["motif"]["status"] == "unresolved"


@pytest.mark.parametrize("complete", [True, False])
def test_joint_core_counts_only_complete_exact_system(complete):
    ids = ["f", "g", "h", "accessory"] if complete else ["f", "g", "accessory"]
    data = base(*ids)
    models = {"f": "SspBCDE__SspF", "g": "SspBCDE__SspG", "h": "SspBCDE__SspH"}
    data["calls"] = [call(g, family="Ssp", model=models.get(g), system="ssp", run_id="r") for g in ids]
    data["groups"] = [group("ssp", ids, system_id="ssp")]
    report = review_systems(data, analyze_saved(data))
    assert report["summary"]["effector_CDS"] == (3 if complete else 0)
    assert assessment(report, "accessory")["endpoints"]["effector"]["status"] == "unresolved"
    f = assessment(report, "f")
    assert f["endpoints"]["effector"]["status"] == ("supported" if complete else "unresolved")
    if not complete:
        assert f["endpoints"]["effector"]["reason_codes"] == ["component_identity_requires_complete_complex"]
    else:
        assert f["individual_catalysis"] == "unassigned_complex_core"


def test_explicit_graded_motif_is_separate_endpoint():
    data = base("e")
    data["calls"] = [call("e")]
    data["effector_evidence"] = [effect("e", motif_state="supported", motif_grade="two_bound_reference_sites")]
    report = review_systems(data, analyze_saved(data))
    assert report["summary"]["motif_CDS"] == report["summary"]["effector_CDS"] == 1
    assert assessment(report, "e")["activity"] == "unassessed"


def test_raw_interpro_and_eggnog_are_checked_but_not_sufficient(tmp_path):
    data, result = saved(tmp_path, "a")
    data["suite"] = {"profiles": [dict(sample_id="s", gene_id="a", namespace="InterProScan",
                                      native_fields={"interpro_description": "nuclease domain"})],
                     "annotations": [dict(sample_id="s", gene_id="a", namespace="eggNOG",
                                          native_fields={"Description": "Wadjet JetD nuclease"})]}
    report = review_systems(data, result)
    row = assessment(report, "a")
    assert row["context_evidence_counts"] == {"InterProScan": 1, "eggNOG": 1}
    assert row["endpoints"]["effector"]["reason_codes"] == ["context_requires_qualified_functional_rule"]
    assert report["proposals"]["effector_evidence"] == []


def rule_data(tmp_path, role):
    data, result = saved(tmp_path, "a")
    rules = registry(component(defense_scope="host_defense", conclusion=dict(component="a", role=role, substrate="DNA")))
    source = tmp_path / "rules.json"
    source.write_text(json.dumps(rules))
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    data.update(system_review_rules=rules, system_review_rule_source=dict(path=str(source), sha256=digest))
    data["input_sha256"][str(source)] = digest
    result["system_review"] = review_all_systems(data, result, rules)
    return data, result


def test_explicit_bound_effector_rule_proposal_validates_and_is_idempotent(tmp_path):
    data, result = rule_data(tmp_path, "individual_nuclease_effector")
    report = review_systems(data, result)
    proposals = report["proposals"]["effector_evidence"]
    assert len(proposals) == 1
    assert proposals[0]["motif_state"] == "unassessed"
    data["effector_evidence"] = proposals
    applied = analyze_saved(data)
    assert applied["global_summary"]["effector_CDS"] == 1
    assert review_systems(data, applied)["proposals"]["effector_evidence"] == []


def test_generic_qualified_nuclease_role_is_not_effector_proposal(tmp_path):
    data, result = rule_data(tmp_path, "nuclease")
    assert result["system_review"]["rule_evaluations"][0]["status"] == "supported"
    report = review_systems(data, result)
    assert report["proposals"]["effector_evidence"] == []
    assert report["summary"]["effector_CDS"] == 0


def test_explicit_nuclease_exclusion_blocks_effector_alias_proposal(tmp_path):
    data, result = rule_data(tmp_path, "individual_nuclease_effector")
    # A separate explicit negative contract must not be bypassed merely because
    # 'individual_nuclease_effector' is a more specific string than 'nuclease'.
    negative = deepcopy(result["system_review"]["rule_evaluations"][0])
    negative.update(rule_id="excluded", evaluation_id="negative", conclusion={"role": "noncatalytic"})
    result["system_review"]["rule_evaluations"].append(negative)
    assert review_systems(data, result)["proposals"]["effector_evidence"] == []


def test_version_mismatch_cannot_become_automatic_effector(tmp_path):
    data, _ = rule_data(tmp_path, "individual_nuclease_effector")
    data["calls"][0]["database_version"] = "unseen_version"
    _, result = saved(tmp_path, "a")
    report = review_systems(data, result)
    assert report["proposals"]["effector_evidence"] == []


def test_coverage_accounting_rejects_loss_and_unknown_scope():
    data = base("a")
    data["calls"] = [call("a")]
    result = analyze_saved(data)
    report = review_systems(data, result)
    with pytest.raises(ValueError, match="unaccounted"):
        summarize_assessments(report["gene_assessments"], [["s", "missing"]])
    result["system_review"]["gene_reviews"] = []
    with pytest.raises(ValueError, match="lost"):
        review_systems(data, result)


def test_conflict_preserves_accepted_count_but_not_complete_classification():
    data = base("d")
    data["calls"] = [call("d", model="Wadjet__JetD_II")]
    result = analyze_saved(data)
    result["system_review"]["gene_reviews"][0]["axes"]["component_role"]["status"] = "conflicting"
    report = review_systems(data, result)
    assert report["summary"]["effector_CDS"] == 1
    assert report["summary"]["metric_status"]["effector"] == "partial"
    assert report["summary"]["unresolved_CDS"]["effector"] == 1

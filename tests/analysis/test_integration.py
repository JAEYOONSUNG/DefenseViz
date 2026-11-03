from copy import deepcopy
from dataclasses import asdict
import json
import pytest
from defenseviz.analysis import integrate, evaluate, parse_evidence, parse_policy
from defenseviz.analysis.policy import Evidence, Profile, Policy


def evidence(**updates):
    out = dict(family="Type II", role="R", valid_sequence=True, identity=None,
        query_coverage=None, reference_coverage=None, evalue=None, role_compatible=False,
        family_compatible=None, hmm_role_ambiguous=False, profiles=[], raw_nuclease_classes=[],
        legacy_nuclease=False, legacy_motor=False, legacy_mtase_motif=False,
        old_strong=False, old_moderate=False, old_profile=False, annotation=False,
        old_motif=False, old_context=False, strict_structure=False,
        alternative_annotation=False, conflict=False, partial=False)
    out.update(updates)
    return out


def call(gene="g", caller="alpha", sample="s", **updates):
    out = dict(sample_id=sample, gene_id=gene, caller=caller, family="new-family",
               category="R-M", evidence_kind="final_member", role="R", rm_type="Type II", system_id=None)
    out.update(updates)
    return out


SAMPLES = [dict(sample_id="s", label="Arbitrary α strain", highlight=True)]


def test_native_union_overlap_and_caller_scoped_systems():
    calls = [call(system_id="1"), call(caller="beta", system_id="1"),
             call(gene="g2", system_id="1"), call(gene="g2", family="another", system_id="2")]
    r = integrate(SAMPLES, calls, [])
    assert r["sample_summary"][0]["union_cds"] == 2
    assert r["sample_summary"][0]["unassessed_cds"] == 2
    assert r["cells"][0]["cds_count"] == 2
    assert r["cells"][0]["system_count"] == 2
    assert sum(x["cds_count"] for x in r["cells"]) == 3
    assert {tuple(x["callers"]): x["count"] for x in r["overlap"]} == {("alpha","beta"):1,("alpha",):1}
    assert "complete_systems" not in r
    json.dumps(r, allow_nan=False)


def test_hmm_reference_and_native_are_different_endpoints():
    r = integrate(SAMPLES, [call("h", evidence_kind="hmm_hit"), call("r", evidence_kind="reference_hit"),call("n")], [])
    s = r["sample_summary"][0]
    assert s["union_cds"] == s["native_union_cds"] == 1
    assert s["rm_pool_cds"] == s["unassessed_cds"] == 2
    assert r["cells"][0]["cds_count"] == 1
    assert len(r["evidence_calls"]) == 2
    assert r["rm_components"][0]["count"] == 2
    assert r["rm_components"][0]["supported_count"] == 0


def test_hmm_only_does_not_infer_final_member_zero():
    r=integrate(SAMPLES,[call(evidence_kind="hmm_hit")],[])
    assert r["sample_summary"][0]["caller_counts"] == {"alpha":None}
    assert len(r["evidence_calls"]) == 1
    r=integrate(SAMPLES,[call(evidence_kind="reference_hit")],[],coverage=[dict(sample_id="s",caller="alpha",status="not_run")])
    assert r["sample_summary"][0]["caller_counts"] == {"alpha":None}
    assert r["sample_summary"][0]["rm_pool_cds"] == 1


def test_reference_provider_is_not_an_implicit_native_caller():
    reference=call(evidence_kind="reference_hit",caller="Arbitrary reference source")
    r=integrate(SAMPLES,[reference],[])
    assert r["callers"] == []
    assert r["sample_summary"][0]["caller_counts"] == {}
    assert r["evidence_sources"] == ["Arbitrary reference source"]
    assert r["sample_summary"][0]["rm_pool_cds"] == 1
    assert len(r["evidence_calls"]) == 1
    # Explicit native coverage, unlike a provider name, establishes caller scope.
    r=integrate(SAMPLES,[reference],[],coverage=[dict(sample_id="s",caller="Arbitrary reference source",status="missing")])
    assert r["callers"] == ["Arbitrary reference source"]
    assert r["sample_summary"][0]["caller_counts"] == {"Arbitrary reference source":None}
    # HMM evidence names an actual caller, while its final output remains missing.
    r=integrate(SAMPLES,[call(evidence_kind="hmm_hit",caller="Novel detector")],[])
    assert r["callers"] == ["Novel detector"]
    assert r["evidence_sources"] == ["Novel detector"]
    assert r["sample_summary"][0]["caller_counts"] == {"Novel detector":None}


def test_explicit_control_reference_is_out_of_scope_without_policy_data():
    r=integrate(SAMPLES,[call(evidence_kind="reference_hit",role="C")],[])
    assert r["decisions"][0]["status"] == "out_of_scope"
    assert r["sample_summary"][0]["out_of_scope_cds"] == 1


def test_retained_evidence_only_adds_integrated_not_native_caller():
    ev = [dict(sample_id="s",gene_id="outside",evidence=evidence(old_strong=True,old_context=True))]
    r = integrate(SAMPLES, [call()], ev)
    s = r["sample_summary"][0]
    assert (s["union_cds"],s["native_union_cds"],s["retained_dnmb_only_cds"]) == (2,1,1)
    assert s["caller_counts"] == {"alpha":1}
    assert r["cells"][0]["cds_count"] == 1
    assert r["cell_scope"] == "native_final_members"


def test_missing_coverage_is_na_but_complete_empty_is_zero():
    samples = SAMPLES + [dict(sample_id="s2",label="Other")]
    coverage = [dict(sample_id="s",caller="alpha",status="complete"),
                dict(sample_id="s",caller="beta",status="complete"),
                dict(sample_id="s2",caller="alpha",status="not_run")]
    r = integrate(samples,[call()],[],coverage=coverage)
    assert r["sample_summary"][0]["caller_counts"] == {"alpha":1,"beta":0}
    assert r["sample_summary"][1]["caller_counts"] == {"alpha":None,"beta":None}
    assert r["cells"][1]["count_status"] == "unassessed"
    assert r["overlap"][0]["all_callers_complete"] is True
    with pytest.raises(ValueError,match="coverage"):
        integrate(SAMPLES,[call()],[],coverage=[])


@pytest.mark.parametrize("field,value", [("identity",float("nan")),("identity",2),
    ("evalue",-1),("old_motif","false"),("old_context",None),("valid_sequence",1)])
def test_bad_measurements_never_silently_supply_support(field,value):
    with pytest.raises(ValueError): parse_evidence(evidence(**{field:value}))


def test_missing_boolean_is_error_null_record_is_unassessed():
    ev = evidence();del ev["old_context"]
    with pytest.raises(ValueError,match="missing"):
        parse_evidence(ev)
    r = integrate(SAMPLES,[],[dict(sample_id="s",gene_id="g",evidence=None,provenance={"missing_fields":["old_context"]})])
    assert r["decisions"][0]["status"] == "unassessed"
    assert r["decisions"][0]["provenance"]["missing_fields"] == ["old_context"]


def test_control_and_auxiliary_not_scored_but_in_pool():
    ev = [dict(sample_id="s",gene_id=role,evidence=evidence(role=role,old_strong=True,old_context=True)) for role in ["C","Aux"]]
    r = integrate(SAMPLES,[],ev)
    assert all(x["status"] == "out_of_scope" for x in r["decisions"])
    assert r["sample_summary"][0]["rm_pool_cds"] == 2
    assert r["sample_summary"][0]["scored_cds"] == 0
    assert r["sample_summary"][0]["out_of_scope_cds"] == 2


def test_missing_system_identity_is_unassessed_not_zero_systems():
    r=integrate(SAMPLES,[call()],[])
    assert r["cells"][0]["system_count"] is None
    assert r["native_system_counts"][0]["count"] is None
    assert r["native_system_counts"][0]["count_status"] == "unassessed_missing_system_ids"


def test_invalid_sequence_with_legacy_support_is_unassessed():
    r = integrate(SAMPLES,[],[dict(sample_id="s",gene_id="g",evidence=evidence(valid_sequence=False,old_strong=True,old_context=True))])
    assert r["decisions"][0]["status"] == "unassessed"
    assert r["sample_summary"][0]["retained_cds"] == 0


def test_tuple_identity_no_delimiter_collision_or_same_sequence_merge():
    samples = [dict(sample_id="a__b",label="one"),dict(sample_id="a",label="two")]
    r = integrate(samples,[call(gene="c",sample="a__b"),call(gene="b__c",sample="a")],[])
    assert [x["union_cds"] for x in r["sample_summary"]] == [1,1]
    with pytest.raises(ValueError,match="Duplicate evidence"):
        integrate(SAMPLES,[],[dict(sample_id="s",gene_id="g",evidence=None)]*2)
    with pytest.raises(ValueError,match="duplicate call"):
        integrate(SAMPLES,[call(),call()],[])


def test_native_system_collision_cross_family_rejected_not_merged():
    with pytest.raises(ValueError,match="conflicting families"):
        integrate(SAMPLES,[call(system_id="1"),call(gene="other",system_id="1",family="other")],[])


def test_evidence_component_support_never_transfers_to_different_role():
    r = integrate(SAMPLES,[call(role="M")],[dict(sample_id="s",gene_id="g",evidence=evidence(legacy_nuclease=True))])
    by_role = {x["role"]:x for x in r["rm_components"]}
    assert by_role["M"]["supported_count"] == 0
    assert by_role["R"]["supported_count"] == 1


def test_policy_capped_sequence_axis_and_penalty_no_double_count():
    e = parse_evidence(evidence(old_strong=True,old_moderate=True,old_profile=True,annotation=True,
        old_motif=True,old_context=True,strict_structure=True,alternative_annotation=True,conflict=True,partial=True))
    d = evaluate(e,parse_policy())
    assert d["score_sequence"] == 35
    assert d["score_conflict_penalty"] == 0
    assert d["score_total"] == 75
    assert d["keep"] is True


def test_same_profile_gate_cannot_stitch_coverage_between_hits():
    profiles=[dict(family_compatible=True,role_compatible=True,evalue=1e-40,query_coverage=.8,model_coverage=.1),
              dict(family_compatible=True,role_compatible=True,evalue=1e-40,query_coverage=.1,model_coverage=.8)]
    d = evaluate(parse_evidence(evidence(profiles=profiles)),parse_policy())
    assert not d["profile_support"]


def test_numeric_r_selected_preset_and_independent_context_requirement():
    e=evidence(profiles=[dict(family_compatible=True,role_compatible=True,evalue=1e-20,query_coverage=.9,model_coverage=.9)],
        curated_identity=.6,curated_query_coverage=.9,curated_reference_coverage=.9,curated_evalue=1e-6,
        curated_role_compatible=True,curated_family_compatible=True,raw_nuclease_classes=["PD-ExK"])
    d=evaluate(parse_evidence(e),parse_policy())
    assert d["numeric_R_component_admission"]
    assert not d["integrated_R_component_support"]
    assert not d["new_motif_credit"]
    e["old_context"]=True
    assert evaluate(parse_evidence(e),parse_policy())["new_motif_credit"]
    e["raw_nuclease_classes"]=["PLD-HKD"]
    assert not evaluate(parse_evidence(e),parse_policy())["numeric_R_component_admission"]


def test_parameter_types_and_frozen_weights():
    assert parse_policy().r_identity_min == .6
    for value in [{"numeric_r_enabled":"false"},{"r_identity_min":True},{"motif":99},{"unknown":1},{"retain_cutoff":0}]:
        with pytest.raises(ValueError):parse_policy(value)


def test_empty_inputs_are_explicit_empty_lists():
    r=integrate([],[],[],coverage=[])
    for k in ["samples","callers","families","cells","sample_summary","decisions","overlap","rm_components","native_system_counts","evidence_calls"]:
        assert r[k] == []


def test_geometry_never_assigns_context_from_adjacency():
    g=dict(sample_id="s",gene_id="g",contig="one",start=1,end=20,strand="+",sequence_sha256="0"*64)
    r=integrate(SAMPLES,[call()],[dict(sample_id="s",gene_id="g",evidence=evidence())],genes=[g])
    assert not r["decisions"][0]["decision"]["context_support"]
    bad=deepcopy(g);bad["start"]=True
    with pytest.raises(ValueError,match="coordinates"):
        integrate(SAMPLES,[],[],genes=[bad])


def test_geometry_without_sequence_hash_is_valid_but_cannot_grant_evidence():
    g=dict(sample_id="s",gene_id="g",contig="one",start=1,end=20,strand="+")
    r=integrate(SAMPLES,[call()],[],genes=[g])
    assert r["decisions"][0]["status"] == "unassessed"
    assert r["sample_summary"][0]["retained_cds"] == 0


def test_custom_policy_is_not_mislabeled_selected_or_calibrated():
    default=integrate([],[],[])["policy_metadata"]
    assert default["is_selected_preset"] is True
    custom=integrate([],[],[],policy={"r_identity_min":.8})["policy_metadata"]
    assert custom["is_selected_preset"] is False
    assert custom["preset"] == "custom_parameter_override"
    assert custom["status"] == "custom_parameters_not_independent_validation"
    assert custom["parameter_overrides"] == {"r_identity_min":{"selected":.6,"value":.8}}


def test_conflicting_explicit_geometry_and_policy_sequence_hashes_rejected():
    g=dict(sample_id="s",gene_id="g",contig="one",start=1,end=20,strand="+",sequence_sha256="0"*64)
    r=dict(sample_id="s",gene_id="g",evidence=evidence(old_strong=True,old_context=True),provenance={"sequence_sha256":"1"*64})
    with pytest.raises(ValueError,match="conflicts with gene geometry"):
        integrate(SAMPLES,[],[r],genes=[g])
    r["provenance"]["sequence_sha256"]="malformed"
    with pytest.raises(ValueError,match="SHA-256"):
        integrate(SAMPLES,[],[r])
    r["provenance"]["sequence_sha256"]="0"*64
    assert integrate(SAMPLES,[],[r],genes=[g])["sample_summary"][0]["retained_cds"]==1

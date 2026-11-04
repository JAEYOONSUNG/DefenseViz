"""Current retained contexts must not inherit stale original completeness."""
from copy import deepcopy
import hashlib

import pytest

from defenseviz.analysis.saved import analyze_saved


REF = {"path": "saved-contexts.json", "sha256": "a" * 64}


def bundle(typ="Type I", roles=None, excluded=()):
    roles = roles or {"r": "R", "m": "M", "s": "S"}
    genes = [dict(sample_id="sample", gene_id=g,
                  protein_sha256=hashlib.sha256(g.encode()).hexdigest()) for g in roles]
    return dict(samples=[dict(sample_id="sample", label="Example")], genes=genes,
        calls=[], coverage=[], annotations=[],
        saved_decisions=[dict(**g, state="not_retained" if g["gene_id"] in excluded else "retained",
            family=typ, role=roles[g["gene_id"]], source_ref=REF) for g in genes],
        groups=[dict(sample_id="sample", group_id="context", provider="DNMB", run_id="run",
            family=typ, membership_kind="saved_context", context_state="complete",
            members=[{k: v for k, v in g.items() if k != "sample_id"} for g in genes], source_refs=[REF])])


def type_row(result, typ):
    return next(r for r in result["type_tables"] if r["rm_type"] == typ)


@pytest.mark.parametrize("typ,roles,excluded", [
    ("Type I", {"r": "R", "m": "M", "s": "S"}, ("m", "s")),
    ("Type III", {"r": "R", "m": "M"}, ("m",)),
])
def test_filtered_partner_loss_cannot_keep_original_complete_context(typ, roles, excluded):
    data = bundle(typ, roles, excluded); original = deepcopy(data)
    result = analyze_saved(data); row = type_row(result, typ)
    assert data == original
    assert row["R"] == row["incomplete_R_context"] == 1
    assert row["complete_contexts"] == result["sample_summary"][0]["complete_saved_contexts"] == 0
    assert result["global_summary"]["priority_CDS"] == 1  # retention is not activity
    assessment = result["context_assessments"][0]
    assert assessment["original_context_state"] == "complete"
    assert assessment["context_state"] == "missing_retained_partner_in_saved_context"
    assert len(assessment["removed_member_keys"]) == len(excluded)
    assert row["activity"] == "unassessed"


def test_unknown_current_role_does_not_establish_a_missing_partner():
    data = bundle(roles={"r": "R", "m": "unknown", "s": "S", "extra": "Aux"}, excluded=("extra",))
    result = analyze_saved(data); row = type_row(result, "Type I")
    assert result["context_assessments"][0]["context_state"] == "unassessed"
    assert row["context_state"] == "unassessed"
    assert row["incomplete_R_context"] is None


def test_lost_accessory_does_not_invalidate_preserved_required_roles():
    data = bundle(roles={"r": "R", "m": "M", "s": "S", "extra": "Aux"}, excluded=("extra",))
    result = analyze_saved(data)
    assert result["context_assessments"][0]["context_state"] == "complete"
    assert type_row(result, "Type I")["incomplete_R_context"] == 0


def test_type_iv_has_no_inferred_methylation_partner_requirement():
    data = bundle("Type IV", {"r": "R", "extra": "Aux"}, ("extra",))
    result = analyze_saved(data); row = type_row(result, "Type IV")
    assert result["context_assessments"][0]["context_state"] == "complete"
    assert result["context_assessments"][0]["required_roles"] is None
    assert row["incomplete_R_context"] is None
    assert row["R"] == 1


def test_reciprocal_original_contexts_count_once_in_every_summary():
    data = bundle(); reciprocal = deepcopy(data["groups"][0]); reciprocal["group_id"] = "reciprocal"
    data["groups"].append(reciprocal)
    result = analyze_saved(data)
    assert len(result["context_assessments"]) == 1
    assert result["context_assessments"][0]["group_ids"] == ["context", "reciprocal"]
    assert result["sample_summary"][0]["complete_saved_contexts"] == 1
    assert type_row(result, "Type I")["complete_contexts"] == 1
    assert type_row(result, "Total")["complete_contexts"] == 1
    reciprocal["run_id"] = "other-run"
    assert analyze_saved(data)["sample_summary"][0]["complete_saved_contexts"] == 2


def test_different_original_boundaries_do_not_merge_after_filtering():
    data = bundle(roles={"r": "R", "m": "M", "s": "S", "extra": "Aux"}, excluded=("extra",))
    smaller = deepcopy(data["groups"][0]); smaller["group_id"] = "smaller"
    smaller["members"] = [m for m in smaller["members"] if m["gene_id"] != "extra"]
    data["groups"].append(smaller)
    result = analyze_saved(data)
    assert len(result["context_assessments"]) == 2
    assert result["sample_summary"][0]["complete_saved_contexts"] == 2


def test_conflicting_reciprocal_context_assessments_fail():
    data = bundle(); reciprocal = deepcopy(data["groups"][0]); reciprocal["group_id"] = "reciprocal"
    reciprocal["context_state"] = "missing"
    data["groups"].append(reciprocal)
    with pytest.raises(ValueError, match="Contradictory exact saved-context"):
        analyze_saved(data)


def test_explicit_current_review_context_remains_distinct_from_saved_retention():
    data = bundle("Type III", {"r": "R", "m": "M"}, ("m",))
    data["groups"][0]["provider"] = "Review"
    data["calls"] = [dict(sample_id="sample", gene_id="m", caller="DefenseFinder",
        evidence_kind="final_member", family="RM_Type_III", category="R-M", source_ref=REF)]
    result = analyze_saved(data)
    assessment = result["context_assessments"][0]
    assert assessment["assessment_scope"] == "explicit_reviewed_current_context"
    assert assessment["context_state"] == "complete"
    assert assessment["removed_member_keys"] == []
    assert result["global_summary"]["retained_CDS"] == 1


def test_explicit_current_context_conflicting_with_filtered_group_fails():
    data = bundle("Type III", {"r": "R", "m": "M"}, ("m",))
    data["saved_decisions"][0]["context_state"] = "complete"
    with pytest.raises(ValueError, match="Conflicting context assessments"):
        analyze_saved(data)


def test_pre_retention_context_never_becomes_current_by_retaining_every_member():
    data = bundle()
    data['rm_evaluation'] = {'native_saved_decisions': [dict(r, state='unassessed') for r in data['saved_decisions']]}
    result = analyze_saved(data)
    assert result['context_assessments'][0]['context_state'] == 'unassessed'
    assert result['context_assessments'][0]['original_context_state'] == 'complete'
    assert type_row(result, 'Type I')['context_state'] == 'unassessed'
    assert result['sample_summary'][0]['complete_saved_contexts'] == 0


def test_replaced_native_run_does_not_transfer_original_complete_context():
    data = bundle()
    data['rm_evaluation'] = {'decision_mode': 'replace_scoped', 'candidate_scope': [{'sample_id': 'sample'}],
                            'native_saved_decisions': deepcopy(data['saved_decisions'])}
    result = analyze_saved(data)
    assert result['context_assessments'][0]['assessment_scope'] == 'pre_retention_context_not_reassessed'
    assert result['context_assessments'][0]['context_state'] == 'unassessed'


def test_unmet_policy_R_route_does_not_remove_annotation_based_priority():
    data = bundle('Type IV', {'r': 'R'})
    data['saved_decisions'][0]['integrated_R_component_support'] = False
    result = analyze_saved(data)
    assert result['global_summary']['retained_CDS'] == 1
    assert result['global_summary']['restriction_role_CDS'] == 1
    assert type_row(result, 'Type IV')['R'] == 1
    assert result['global_summary']['priority_CDS'] == 1
    assert result['global_summary']['priority_support_qualified_CDS'] == 0
    assert result['global_summary']['R_only_policy_unmet_CDS'] == 1
    assert result['gene_ledger'][0]['policy_R_component_support'] is False
    assert result['gene_ledger'][0]['activity'] == 'unassessed'


def test_independent_effector_review_can_prioritize_despite_failed_R_gate():
    data = bundle('Type IV', {'r': 'R'})
    data['saved_decisions'][0]['integrated_R_component_support'] = False
    gene = data['genes'][0]
    data['effector_evidence'] = [dict(**gene, status='supported', assignment_kind='individual_nuclease_effector',
                                     rule_id='independent_saved_review', source_refs=[REF])]
    result = analyze_saved(data)
    assert result['global_summary']['priority_CDS'] == 1
    assert result['priority_ledger'][0]['selection_route'] == 'R_and_effector'
    assert result['global_summary']['priority_support_qualified_CDS'] == 1


@pytest.mark.parametrize('mode', ['evaluate_unassessed', 'replace_scoped'])
def test_old_candidate_context_cannot_expand_current_associated_set(mode):
    data = bundle()
    data['rm_evaluation'] = {'decision_mode': mode, 'candidate_scope': [{'sample_id': 'sample'}],
        'native_saved_decisions': [dict(r, state='unassessed') for r in data['saved_decisions']]}
    gene = data['genes'][0]
    data['effector_evidence'] = [dict(**gene, status='supported', assignment_kind='individual_nuclease_effector',
                                     rule_id='independent_saved_review', source_refs=[REF])]
    result = analyze_saved(data)
    assert result['global_summary']['retained_CDS'] == 3
    assert result['global_summary']['associated_CDS'] == 1
    assert result['global_summary']['effector_CDS'] == 1
    assert result['context_assessments'][0]['original_member_keys'] == [['sample', 'm'], ['sample', 'r'], ['sample', 's']]
    assert result['groups_selected'] == []
    review = deepcopy(data['groups'][0]); review.update(group_id='current', provider='Review')
    data['groups'].append(review)
    assert analyze_saved(data)['global_summary']['associated_CDS'] == 3

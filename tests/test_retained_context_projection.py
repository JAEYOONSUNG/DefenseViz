"""The retained context panel consumes current source-bound assessments."""
from copy import deepcopy
import hashlib

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.plot.publication.manuscript_model import enrich_manuscript_model


REF = {"path": "contexts.json", "sha256": "a" * 64}


def bundle(provider="DNMB", excluded=()):
    genes = [dict(sample_id="sample", gene_id=g,
                  protein_sha256=hashlib.sha256(g.encode()).hexdigest()) for g in ["r", "m", "s"]]
    return dict(samples=[dict(sample_id="sample", label="Example")], genes=genes,
        calls=[], coverage=[], annotations=[],
        saved_decisions=[dict(**g, state="not_retained" if g["gene_id"] in excluded else "retained",
            family="Type I", role=g["gene_id"].upper(), source_ref=REF) for g in genes],
        groups=[dict(sample_id="sample", group_id="context", provider=provider, run_id="run",
            family="Type I", membership_kind="saved_context", context_state="complete",
            origin="supplement" if provider == "Review" else "native",
            members=[{k: v for k, v in g.items() if k != "sample_id"} for g in genes], source_refs=[REF])])


def project(data):
    return enrich_manuscript_model(data, analyze_saved(data))["manuscript_atlas"]


def test_projection_uses_current_state_without_overwriting_original_source():
    data = bundle(excluded=("m", "s")); original = deepcopy(data)
    result = project(data); group = result["group_records"][0]
    assert data == original
    assert group["context_state"] == "missing_retained_partner_in_saved_context"
    assert group["original_context_state"] == "complete"
    assert group["context_origin"] == "native_saved_context"
    assert group["rm_R_support"] == group["rm_reference_link"] == "unassessed"


def test_source_bound_review_context_is_visible_without_fabricated_support_badges():
    data = bundle("Review"); result = project(data)
    group = result["group_records"][0]
    assert group["provider"] == "Review" and group["origin"] == "supplement"
    assert group["context_origin"] == "reviewed_context"
    assert group["composition_class"] == "I_RMS"
    composition = result["rows"][0]["compositions"]["I_RMS"]
    assert composition["count"] == 1
    assert composition["R_support_all"] is None and composition["reference_any"] is None
    assert composition["status"] == "partial"


def test_review_context_with_native_only_member_is_not_recast_as_retained_group():
    data = bundle("Review", excluded=("s",))
    data["calls"].append(dict(sample_id="sample", gene_id="s", caller="DefenseFinder",
        family="Type I", category="R-M", evidence_kind="final_member", source_ref=REF))
    result = project(data)
    assert result["group_records"] == []
    assert any("nonretained members" in w for w in result["warnings"])


def test_stale_current_context_binding_is_rejected_before_projection():
    data = bundle(); result = analyze_saved(data)
    result["context_assessments"][0]["member_keys"].pop()
    with pytest.raises(ValueError, match="differs from exact saved group"):
        enrich_manuscript_model(data, result)


def test_unknown_context_does_not_gain_reference_badge():
    data = bundle(); data["groups"][0]["context_state"] = "unassessed"
    result = project(data)
    assert result["group_records"][0]["context_state"] == "unassessed"
    assert result["rows"][0]["compositions"]["I_RMS"]["reference_any"] is None


def test_unresolved_imported_context_caps_declared_complete_composition():
    data = bundle()
    data['options'] = {'manuscript_atlas': {'group_coverage': {'sample': 'complete'}}}
    data['suite'] = {'unresolved_contexts': [{'sample_id': 'sample', 'source_ref': REF,
                                           'reason': 'context_partner_requires_one_exact_canonical_identity'}]}
    atlas = project(data)
    assert atlas['rows'][0]['compositions']['I_RMS']['count'] == 1
    assert atlas['rows'][0]['compositions']['I_RMS']['status'] == 'partial'


def test_prior_context_cannot_duplicate_a_current_reviewed_composition():
    data = bundle()
    data['rm_evaluation'] = {'native_saved_decisions': [dict(r, state='unassessed') for r in data['saved_decisions']]}
    review = deepcopy(data['groups'][0]); review.update(group_id='reviewed', provider='Review', origin='supplement')
    data['groups'].append(review)
    atlas = project(data)
    assert len(atlas['group_records']) == 2  # both sources remain auditable
    assert sum(g['current_composition'] for g in atlas['group_records']) == 1
    assert atlas['rows'][0]['compositions']['I_RMS']['count'] == 1
    assert atlas['rows'][0]['compositions']['I_RMS']['group_ids'] == ['reviewed']


def test_unevaluated_legacy_context_is_preserved_without_current_composition_count():
    data = bundle()
    data['rm_evaluation'] = {'native_saved_decisions': [dict(r, state='unassessed') for r in data['saved_decisions']]}
    atlas = project(data)
    assert len(atlas['group_records']) == 1
    assert atlas['group_records'][0]['composition_class'] == 'I_RMS'
    assert atlas['rows'][0]['compositions']['I_RMS']['count'] == 0
    assert atlas['rows'][0]['compositions']['I_RMS']['status'] == 'unassessed'


def test_current_RM_pool_does_not_absorb_unevaluated_historical_candidates():
    data = bundle()
    data['saved_decisions'][1]['state'] = 'unassessed'
    data['saved_decisions'][2]['state'] = 'unassessed'
    data['groups'] = []
    data['rm_evaluation'] = {'candidate_scope': [{'sample_id': 'sample', 'gene_ids': ['r']}],
                            'native_saved_decisions': deepcopy(data['saved_decisions'])}
    data['calls'] = [dict(sample_id='sample', gene_id='r', caller='DF', evidence_kind='final_member',
                         family='RM_Type_I', category='R-M', source_ref=REF),
                    dict(sample_id='sample', gene_id='m', caller='DP', evidence_kind='final_member',
                         family='General defense', category='Unresolved', source_ref=REF)]
    result = enrich_manuscript_model(data, analyze_saved(data))
    assert result['sample_summary'][0]['rm_pool_CDS'] == 1
    assert result['manuscript_atlas']['rows'][0]['rm_original']['input_pool_CDS'] == 1
    assert result['manuscript_atlas']['rows'][0]['rm_original']['source_counts']['DP'] == 0
    assert result['global_summary']['native_union_CDS'] == 2
    assert len(data['saved_decisions']) == 3
    assert result['rm_candidate_scope']['unselected_historical_member_keys'] == [['sample', 'm'], ['sample', 's']]


def test_untyped_reference_is_counted_without_inventing_RM_taxonomy():
    data = bundle()
    data['groups'] = []
    data['saved_decisions'] = []
    data['calls'] = [dict(sample_id='sample', gene_id='r', caller='Reference', evidence_kind='reference_hit',
                         family='Unresolved', category='Unresolved', reference_id='untyped-subject', source_ref=REF)]
    data['coverage'] = [dict(sample_id='sample', caller='Reference', status='complete', endpoint='reference_matches')]
    result = enrich_manuscript_model(data, analyze_saved(data))
    assert result['sample_summary'][0]['rm_source_counts']['Reference'] == 1
    assert result['manuscript_atlas']['rows'][0]['rm_original']['reference_counts']['Reference'] == 1
    assert result['global_summary']['retained_CDS'] == 0
    assert result['global_summary']['final_RM_CDS'] == 0
    assert data['calls'][0]['category'] == 'Unresolved'

"""Type-table assessment is conditional on its observed assigned-type pool."""
from copy import deepcopy

from defenseviz.analysis.saved import analyze_saved
from defenseviz.plot.publication.atlas import _type_display
from defenseviz.plot.publication.manuscript_model import enrich_manuscript_model
from tests.analysis.test_saved import base, call, cov, decision


def test_outside_unassessed_nominee_does_not_make_a_fully_assessed_type_partial():
    data = base('r', 'other')
    data['calls'] = [dict(call('r'), family='Type I', rm_type='Type I', category='R-M')]
    data['coverage'] = cov('DefenseFinder')
    data['saved_decisions'] = [decision('r', typ='Type I'), decision('other', typ='NA', role='unknown', state='unassessed')]
    result = analyze_saved(data)
    row = next(r for r in result['type_tables'] if r['rm_type'] == 'Type I')
    assert row['counts_status'] == 'partial'  # broader candidate scope remains honest
    assert row['within_pool_counts_status'] == 'complete'
    assert row['within_pool_role_states']['R'] == 'complete'
    assert _type_display(row, 'DNMB')[0] == '1'
    assert result['sample_summary'][0]['retained_counts_status'] == 'partial'
    assert len(result['rm_candidate_scope']['member_keys']) == 2


def test_unknown_member_inside_type_pool_preserves_partial_audit_and_observed_count():
    data = base('r', 'other')
    data['calls'] = [dict(call(g), family='Type I', rm_type='Type I', category='R-M') for g in ['r', 'other']]
    data['coverage'] = cov('DefenseFinder')
    data['saved_decisions'] = [decision('r', typ='Type I'), decision('other', typ='Type I', state='unassessed')]
    row = next(r for r in analyze_saved(data)['type_tables'] if r['rm_type'] == 'Type I')
    assert (row['pool_CDS'], row['unique_CDS']) == (2, 1)
    assert row['within_pool_counts_status'] == 'partial'
    assert _type_display(row, 'DNMB')[0] == '1'


def test_reference_and_typed_candidates_define_rm_input_not_untyped_no_hit_nominees():
    data = base('r', 'untyped', 'ref')
    data['saved_decisions'] = [decision('r', typ='Type I'), decision('untyped', typ='NA', role='unknown', state='not_retained')]
    data['calls'] = [dict(call('r'), family='Type I', category='R-M', rm_type='Type I'),
                     call('untyped', caller='DP', family='General defense'),
                     call('ref', caller='REBASE', family='Unresolved', kind='reference_hit')]
    data['coverage'] = cov('DefenseFinder', 'DP', 'REBASE')
    data['rm_evaluation'] = {'candidate_scope': [{'sample_id': 's', 'gene_ids': ['r', 'untyped']}]}
    original = deepcopy(data)
    result = enrich_manuscript_model(data, analyze_saved(data))
    assert data == original
    assert result['rm_candidate_scope']['member_keys'] == [['s', 'r'], ['s', 'untyped']]
    assert result['rm_candidate_scope']['typed_member_keys'] == [['s', 'r']]
    assert result['rm_candidate_scope']['untyped_member_keys'] == [['s', 'untyped']]
    assert result['sample_summary'][0]['rm_pool_CDS'] == 2
    assert result['manuscript_atlas']['rows'][0]['rm_original']['input_pool_CDS'] == 2
    assert result['manuscript_atlas']['rows'][0]['rm_original']['source_counts']['DP'] == 0
    assert next(r for r in result['source_sets'] if r['source_id'] == 'DP')['count'] == 1
    assert result['global_summary']['retained_CDS'] == 1

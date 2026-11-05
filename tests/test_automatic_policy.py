"""Bounds must preserve unknown measurements and the frozen policy decisions."""
from copy import deepcopy
from dataclasses import asdict
import itertools

import pytest

from defenseviz.analysis.evidence import evidence_dict, parse_evidence
from defenseviz.analysis.policy import Evidence, Policy, evaluate
from defenseviz.automatic import policy_bounds


def facts(**updates):
    d = dict(family='Type II', role='R', valid_sequence=True,
             identity=None, query_coverage=None, reference_coverage=None, evalue=None,
             role_compatible=False, family_compatible=None, hmm_role_ambiguous=False,
             profiles=(), raw_nuclease_classes=frozenset(), legacy_nuclease=False,
             legacy_motor=False, legacy_mtase_motif=False, old_strong=False,
             old_moderate=False, old_profile=False, annotation=False, old_motif=False,
             old_context=False, strict_structure=False, alternative_annotation=False,
             conflict=False, partial=False)
    d.update(updates)
    return evidence_dict(Evidence(**d))


def test_complete_evidence_matches_original_policy():
    for flags in itertools.product((False, True), repeat=7):
        d = facts(**dict(zip(('old_strong', 'old_profile', 'old_motif', 'old_context',
                              'strict_structure', 'alternative_annotation', 'partial'), flags)))
        b = policy_bounds(d, [])
        ref = evaluate(parse_evidence(d), Policy())
        assert b['score_min'] == b['score_max'] == ref['score_total']
        assert b['keep'] == ref['keep']
        assert b['integrated_R_component_support'] == ref['integrated_R_component_support']


def test_missing_penalties_are_enumerated_with_mutual_suppression():
    names = ('alternative_annotation', 'conflict', 'partial')
    d = facts(old_strong=True, old_context=True, **dict.fromkeys(names, None))
    b = policy_bounds(d, names)
    assert (b['score_min'], b['score_max']) == (30, 55)
    assert b['keep'] is None
    assert all(d[n] is None for n in names)


@pytest.mark.parametrize('unknown', [
    ['old_strong', 'old_moderate', 'annotation', 'old_context'],
    ['old_motif', 'strict_structure', 'alternative_annotation', 'conflict', 'partial'],
    ['legacy_nuclease', 'hmm_role_ambiguous', 'old_context', 'legacy_motor'],
])
def test_bounds_cover_every_boolean_completion(unknown):
    d = facts(curated_identity=1., curated_query_coverage=1., curated_reference_coverage=1.,
              curated_evalue=0., curated_role_compatible=True, curated_family_compatible=True)
    d['profiles'] = [dict(family_compatible=True,role_compatible=True,evalue=0.,query_coverage=1.,model_coverage=1.)]
    d['raw_nuclease_classes'] = ['PD-ExK']
    for n in unknown: d[n] = None
    b = policy_bounds(d, unknown)
    decisions=[]
    for flags in itertools.product((False,True), repeat=len(unknown)):
        full=deepcopy(d);full.update(zip(unknown,flags))
        r=evaluate(parse_evidence(full),Policy()); decisions.append(r)
        assert b['score_min'] <= r['score_total'] <= b['score_max']
        assert b['capped_axes_min'] <= r['capped_evidence_axis_count'] <= b['capped_axes_max']
    if b['keep'] is not None: assert all(r['keep']==b['keep'] for r in decisions)


def test_unmeasured_structure_can_be_irrelevant_to_robust_retention():
    d=facts(old_strong=True, old_motif=True, old_context=True, strict_structure=None)
    b=policy_bounds(d,['strict_structure'])
    assert b['keep'] is True and b['score_total'] is None
    assert b['assessment']=='conservative_policy_bounds'
    assert b['unknown_axes']==['strict_structure']
    assert d['strict_structure'] is None


def test_unknown_anchor_can_link_hypothetical_profile():
    d=facts(old_context=True,curated_identity=1.,curated_query_coverage=1.,
            curated_reference_coverage=1.,curated_evalue=0.,curated_role_compatible=True,
            curated_family_compatible=True, reference_nuclease_anchor_support=None)
    unknown=['reference_nuclease_anchor_support','profiles','reference_nuclease_anchor_profile_indices']
    b=policy_bounds(d,unknown)
    assert b['score_min']==20 and b['score_max']==65
    assert b['integrated_R_component_support'] is None


def test_fixed_policy_and_missing_boolean_validation():
    with pytest.raises(ValueError,match='unchanged'): policy_bounds(facts(),[],Policy(retain_cutoff=30))
    with pytest.raises(ValueError,match='Unassessed boolean'): policy_bounds(facts(old_motif=None),[])
    with pytest.raises(ValueError,match='Known family'): policy_bounds(facts(),['family'])


def test_automatic_decisions_preserve_native_memberships_and_unknowns(monkeypatch):
    import defenseviz.automatic_rm as compiler
    from defenseviz.automatic import apply_automatic_rm
    native_ref={'path':'native.tsv','sha256':'a'*64,'row':2}
    genes=[{'sample_id':'s','gene_id':g,'protein_sha256':'b'*64} for g in ('kept','pending','excluded')]
    data={'samples':[{'sample_id':'s'}], 'genes':genes,
          'calls':[{'sample_id':'s','gene_id':'pending','caller':'DF','family':'Type II','source_ref':native_ref}],
          'saved_decisions':[{'sample_id':'s','gene_id':g['gene_id'],'state':'unassessed','source_ref':native_ref} for g in genes],
          'coverage':[], 'suite':{}, 'input_roles':{}}
    records=[]
    for g in genes:
        known={'kept':dict(old_strong=True,old_motif=True,old_context=True), 'pending':dict(old_profile=True), 'excluded':{}}[g['gene_id']]
        unknown=['strict_structure'] if g['gene_id']!='excluded' else []
        e=facts(**known)
        for n in unknown:e[n]=None
        records.append({**g,'assignment_status':'supported','evidence':e,'unknown_axes':unknown,'source_refs':[native_ref]})
    compiled={'records':records,'candidate_scope':[{'sample_id':'s','status':'complete'}]}
    monkeypatch.setattr(compiler,'compile_rm_evidence',lambda d:deepcopy(compiled))
    before=deepcopy(data);out=apply_automatic_rm(data)
    assert data==before and out['calls']==before['calls']
    assert {r['gene_id']:r['state'] for r in out['saved_decisions']}=={'kept':'retained','pending':'unassessed','excluded':'not_retained'}
    assert out['rm_evaluation']['candidate_scope'][0]['coverage_status']=='partial'
    assert out['rm_evaluation']['policy_evaluations']==3
    assert out['saved_decisions'][0]['score_total'] is None


def test_candidate_integrity_prevents_false_complete_scope(monkeypatch):
    import defenseviz.automatic_rm as compiler
    from defenseviz.automatic import apply_automatic_rm
    ref={'path':'native.tsv','sha256':'a'*64}
    gene={'sample_id':'s','gene_id':'g','protein_sha256':'b'*64}
    data={'samples':[{'sample_id':'s'}],'genes':[gene],'calls':[],
          'saved_decisions':[{'sample_id':'s','gene_id':'g','state':'retained','source_ref':ref}],
          'coverage':[], 'suite':{'conflicts':[{'sample_id':'s','code':'unbound_native_candidate'}]}}
    rec={**gene,'assignment_status':'supported','evidence':facts(),'unknown_axes':[],'source_refs':[ref]}
    monkeypatch.setattr(compiler,'compile_rm_evidence',lambda d:{'records':[rec],'candidate_scope':[{'sample_id':'s','status':'complete'}]})
    out=apply_automatic_rm(data)
    assert out['saved_decisions']==data['saved_decisions']
    scope=out['rm_evaluation']['candidate_scope'][0]
    assert scope['coverage_status']=='partial' and scope['source_integrity_issues']
    assert scope['reviewed_count']==1 and scope['evaluated_count']==0 and scope['assessed_count']==1


def test_evaluation_identity_is_relocatable_and_order_stable():
    from defenseviz.automatic import evaluation_identity
    d={'records':[{'sample_id':'s','gene_id':'a','assignment_support':[{'role':'M'},{'role':'R'}],
                   'evidence':facts(), 'source_refs':[{'path':'/old/root/native.tsv','sha256':'a'*64,'row':2}]}],
       'sources':[{'path':'native.tsv','resolved_path':'/old/root/native.tsv','sha256':'a'*64}]}
    copy=deepcopy(d)
    copy['records'][0]['assignment_support'].reverse()
    copy['records'][0]['source_refs'][0]['path']='/new/root/native.tsv'
    copy['sources'][0]['resolved_path']='/new/root/native.tsv'
    assert evaluation_identity(d)==evaluation_identity(copy)
    copy['records'][0]['source_refs'][0]['sha256']='b'*64
    assert evaluation_identity(d)[0]!=evaluation_identity(copy)[0]


@pytest.mark.parametrize('status,expected',[('complete','complete'),('missing','missing'),('not_run','not_run')])
def test_automatic_completed_zero_is_not_an_unrun_module(monkeypatch,status,expected):
    import defenseviz.automatic_rm as compiler
    from defenseviz.automatic import apply_automatic_rm
    from defenseviz.analysis.saved import analyze_saved
    data={'samples':[{'sample_id':'s','label':'s'}],'genes':[],'calls':[],'saved_decisions':[],
          'coverage':[{'sample_id':'s','caller':'DNMB','endpoint':'saved_retention','status':status}],
          'suite':{'capabilities':[{'sample_id':'s','native_candidates':True,'candidate_count':0,'saved_retention':True}] if status=='complete' else []}}
    monkeypatch.setattr(compiler,'compile_rm_evidence',lambda d:{'records':[],'candidate_scope':[{'sample_id':'s','status':'complete'}]})
    out=apply_automatic_rm(data)
    assert out['rm_evaluation']['candidate_scope'][0]['coverage_status']==expected
    assert out['rm_evaluation']['operation']=='automatic_features_with_saved_decision_replay'
    assert out['rm_evaluation']['policy_evaluations']==0
    result=analyze_saved(out)
    assert result['selected_policy'] is None
    assert result['availability']['saved_decisions']==('complete' if status=='complete' else 'unassessed')

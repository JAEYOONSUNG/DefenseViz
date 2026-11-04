from copy import deepcopy
import hashlib
import pytest
from defenseviz.analysis.saved import analyze_saved

REF = {"path":"saved.tsv", "sha256":"a"*64, "row":2}

def gene(g, s="s"):
    return dict(sample_id=s,gene_id=g,protein_sha256=hashlib.sha256((s+g).encode()).hexdigest())

def base(*ids):
    return dict(samples=[dict(sample_id="s",label="Synthetic")],genes=[gene(g) for g in ids],calls=[],coverage=[],saved_decisions=[],annotations=[],effector_evidence=[],groups=[],complexes=[])

def call(g, caller="DefenseFinder", family="Novel", model=None, system=None, kind="final_member", **extra):
    return dict(sample_id="s",gene_id=g,caller=caller,family=family,category="Unresolved",evidence_kind=kind,native_model=model,system_id=system,source_ref=REF,**extra)

def cov(*names):
    return [dict(sample_id="s",caller=n,status="complete")for n in names]

def effect(g, **extra):
    return dict(**gene(g),status="supported",assignment_kind="individual_nuclease_effector",source_refs=[REF],rule_id="review_specific_component",**extra)

def group(gid, members, **extra):
    return dict(group_id=gid,sample_id="s",provider="DefenseFinder",run_id="r",membership_kind="native_system",members=[{k:v for k,v in gene(g).items()if k!='sample_id'}for g in members],source_refs=[REF],**extra)

def decision(g, role="R", typ="Type II", state="retained", **extra):
    return dict(sample_id="s",gene_id=g,state=state,family=typ,role=role,source_ref=REF,**extra)

def test_saved_retention_never_rescores_missing_features(monkeypatch):
    import defenseviz.analysis.policy as p
    monkeypatch.setattr(p,'evaluate',lambda *a,**k:pytest.fail('must not score'))
    d=base('r');d['saved_decisions']=[decision('r')]
    r=analyze_saved(d)
    assert r['global_summary']['retained_CDS']==1
    assert r['global_summary']['priority_CDS']==1
    assert r['global_summary']['effector_CDS']==0
    assert r['provenance']['scoring_runs']==0


def test_exact_union_retained_reviewed_reference_distinct():
    d=base('native','retained','iv','ref')
    d['calls']=[call('native'),call('ref','REBASE',kind='reference_hit')];d['coverage']=cov('DefenseFinder','REBASE')
    d['saved_decisions']=[decision('retained',role='M')]
    d['annotations']=[dict(**gene('iv'),status='supported',assignment_kind='reviewed_RM_component',family='Type IV',role='R_component',include_in_union=True,source_refs=[REF],rule_id='typed_domain_review')]
    r=analyze_saved(d);g=r['global_summary']
    assert (g['integrated_CDS'],g['native_union_CDS'],g['retained_CDS'],g['reviewed_CDS'])==(3,1,1,1)
    assert (g['final_RM_CDS'],g['restriction_role_CDS'],g['effector_CDS'])==(2,1,0)
    assert 'ref' not in {x['gene_id']for x in r['gene_ledger']}


def test_reference_not_native_and_priority_outside_source_mask_preserved():
    d=base('a','r');d['calls']=[call('a','REBASE',kind='reference_hit')];d['coverage']=cov('REBASE')
    d['saved_decisions']=[decision('r')]
    r=analyze_saved(d);o=r['priority_overlap_regions']
    assert o['outside_source_memberships']==1
    assert next(x for x in o['regions']if x['mask']=='0')['count']==1
    assert r['global_summary']['native_union_CDS']==0


def test_no_transitive_group_expansion_and_no_fabricated_singletons():
    d=base('a','b','c','u');d['calls']=[call(x)for x in ['a','b','c','u']];d['coverage']=cov('DefenseFinder')
    d['effector_evidence']=[effect('a'),effect('u')];d['groups']=[group('ab',['a','b']),group('bc',['b','c'])]
    r=analyze_saved(d);g=r['global_summary']
    assert (g['associated_CDS'],g['system_member_CDS'],g['ungrouped_effector_CDS'])==(3,2,1)
    assert [x['group_id']for x in r['groups_selected']]==['ab']
    assert not next(x for x in r['gene_ledger']if x['gene_id']=='c')['associated']


def test_retained_context_uses_only_retained_members():
    d=base('a','b');d['calls']=[call('a'),call('b')];d['coverage']=cov('DefenseFinder')
    d['saved_decisions']=[decision('a'),decision('b',state='not_retained')];d['effector_evidence']=[effect('a')]
    d['groups']=[dict(group('g',['a','b']),membership_kind='saved_context')]
    r=analyze_saved(d)
    assert r['global_summary']['associated_CDS']==1


def test_complete_complex_core_not_individual_catalysis_and_no_accessory_credit():
    d=base('f','g','h','m');d['calls']=[call(x)for x in ['f','g','h','m']];d['coverage']=cov('DefenseFinder')
    d['groups']=[group('mod',['f','g','h','m'])]
    d['complexes']=[dict(complex_id='complete',group_id='mod',sample_id='s',status='complete',core_gene_ids=['f','g','h'],source_refs=[REF],rule_id='review_module')]
    r=analyze_saved(d);g=r['global_summary']
    assert (g['individual_effector_CDS'],g['complex_core_CDS'],g['effector_CDS'],g['associated_CDS'])==(0,3,3,4)
    assert next(x for x in r['gene_ledger']if x['gene_id']=='f')['individual_catalysis']=='unassigned_complex_core'
    assert not next(x for x in r['gene_ledger']if x['gene_id']=='m')['effector']
    d['complexes'][0]['status']='incomplete';r=analyze_saved(d)
    assert r['global_summary']['effector_CDS']==0


def test_complex_cannot_borrow_other_group_or_duplicate_gene():
    d=base('f','g','h');d['calls']=[call(x)for x in ['f','g','h']];d['coverage']=cov('DefenseFinder');d['groups']=[group('fg',['f','g']),group('h',['h'])]
    d['complexes']=[dict(complex_id='x',group_id='fg',sample_id='s',status='complete',core_gene_ids=['f','g','h'],source_refs=[REF],rule_id='module')]
    with pytest.raises(ValueError,match='outside exact group'):analyze_saved(d)
    d['complexes'][0]['core_gene_ids']=['f','f','g']
    with pytest.raises(ValueError,match='unique'):analyze_saved(d)


def test_known_jetd_has_effector_annotation_without_motif_and_abc_not_promoted():
    d=base('d','a');d['calls']=[call('d',model='Wadjet__JetD_II'),call('a',model='Wadjet__JetA_II')];d['coverage']=cov('DefenseFinder')
    r=analyze_saved(d)
    assert (r['global_summary']['effector_CDS'],r['global_summary']['motif_CDS'])==(1,0)
    assert r['availability']['motif']=='unassessed'


def test_generic_name_and_hmm_only_do_not_grant_effector():
    d=base('a','b');d['calls']=[call('a',family='Wadjet nuclease',model='unknown_HNH_nuclease'),call('b',model='Wadjet__JetD_II',kind='hmm_hit')]
    r=analyze_saved(d)
    assert r['global_summary']['effector_CDS']==0
    assert r['global_summary']['native_union_CDS']==1


def test_component_source_required_and_noncatalytic_veto():
    d=base('d');c=call('d',model='Wadjet__JetD_II');c.pop('source_ref');d['calls']=[c];d['coverage']=cov('DefenseFinder')
    assert analyze_saved(d)['global_summary']['effector_CDS']==0
    d['calls'][0]['source_ref']=REF;d['annotations']=[dict(**gene('d'),status='supported',assignment_kind='noncatalytic_component',family='X',role='Aux',source_refs=[REF],rule_id='resolved_motor_only')]
    assert analyze_saved(d)['global_summary']['effector_CDS']==0
    d['effector_evidence']=[effect('d')]
    with pytest.raises(ValueError,match='Conflicting noncatalytic'):analyze_saved(d)


@pytest.mark.parametrize('model',['cas4_I_II_III_IV_V_VI_1','cas5_I-C_11'])
def test_processing_adaptation_not_effector(model):
    d=base('x');d['calls']=[call('x',model=model)];d['coverage']=cov('DefenseFinder')
    assert analyze_saved(d)['global_summary']['effector_CDS']==0


def test_typeiv_not_inactive_without_M_and_fused_unique_total():
    d=base('iv','f');d['saved_decisions']=[decision('iv',typ='Type IV'),decision('f',role='RM')]
    r=analyze_saved(d);rows={x['rm_type']:x for x in r['type_tables']}
    assert rows['Type IV']['incomplete_R_context'] is None
    assert rows['Type IV']['context_state']=='unassessed'
    assert (rows['Type II']['unique_CDS'],rows['Type II']['R'],rows['Type II']['M'])==(1,1,1)
    assert r['global_summary']['restriction_role_CDS']==2


def test_incomplete_context_subset_does_not_remove_retained_R():
    d=base('r','m');d['saved_decisions']=[decision('r',typ='Type III',context_state='missing_retained_partner_in_saved_context'),decision('m',typ='Type III',role='M',orphan=True)]
    r=analyze_saved(d);t=r['type_tables'][0]
    assert (t['R'],t['incomplete_R_context'],t['M'],t['orphan_M'])==(1,1,1,1)
    assert r['global_summary']['priority_CDS']==1


def test_bad_review_binding_missing_source_or_stale_hash():
    d=base('g');d['calls']=[call('g')];d['coverage']=cov('DefenseFinder');e=effect('g');e['protein_sha256']='0'*64;d['effector_evidence']=[e]
    with pytest.raises(ValueError,match='hash mismatch'):analyze_saved(d)
    e['protein_sha256']=gene('g')['protein_sha256'];e['source_refs']=[]
    with pytest.raises(ValueError,match='source_refs'):analyze_saved(d)


def test_duplicate_and_unknown_gene_fail():
    d=base('g');d['saved_decisions']=[decision('g'),decision('g')]
    with pytest.raises(ValueError,match='Duplicate saved decision'):analyze_saved(d)
    d['saved_decisions']=[decision('bad')]
    with pytest.raises(ValueError,match='unknown gene'):analyze_saved(d)


def test_partial_source_scope_not_absence_bits():
    d=base('a');d['calls']=[call('a')];d['coverage']=cov('DefenseFinder')+[dict(sample_id='s',caller='PADLOC',status='missing')]
    r=analyze_saved(d)
    assert r['overlap_regions']['coverage_status']=='unassessed'
    assert r['overlap_regions']['union_count']==0
    row=r['inventory_rows']['families'][0]
    assert row['source_counts']['PADLOC']==0
    assert row['source_count_status']['PADLOC']=='unassessed'


def test_hmm_completion_does_not_create_final_members():
    d=base('a');d['calls']=[call('a',kind='hmm_hit')];d['coverage']=[dict(sample_id='s',caller='DefenseFinder',status='missing')]
    r=analyze_saved(d)
    assert r['global_summary']['native_union_CDS']==0


def test_motif_grade_required_and_not_auto_validated():
    d=base('a');d['calls']=[call('a')];d['coverage']=cov('DefenseFinder');d['effector_evidence']=[effect('a',motif_state='supported')]
    with pytest.raises(ValueError,match='motif_grade'):analyze_saved(d)
    d['effector_evidence'][0]['motif_grade']='single_reference_site'
    r=analyze_saved(d)
    assert r['gene_ledger'][0]['motif_grades']==['single_reference_site']


def test_overlapping_family_rows_deduplicate_union_and_keep_unknowns():
    d=base('a');d['calls']=[call('a',family='new A'),call('a',family='new B')];d['coverage']=cov('DefenseFinder')
    r=analyze_saved(d)
    assert r['global_summary']['integrated_CDS']==1
    assert sum(x['N']for x in r['inventory_rows']['families'])==2


def test_catalog_zero_kept_in_scientific_tables_hidden_only_when_complete():
    d=base('a');d['calls']=[call('a')];d['coverage']=cov('DefenseFinder');d['options']={'family_catalog':[{'id':'empty','category':'Unresolved'}]}
    r=analyze_saved(d);empty=next(x for x in r['inventory_rows']['families']if x['id']=='empty')
    assert empty['N']==0 and empty['visible'] is False
    d['coverage'][0]['status']='partial';assert next(x for x in analyze_saved(d)['inventory_rows']['families']if x['id']=='empty')['visible']


def test_identifier_permutation_and_equal_sequences_do_not_change_counts():
    d=base('x','y');d['genes'][1]['protein_sha256']=d['genes'][0]['protein_sha256'];d['calls']=[call('x'),call('y')];d['coverage']=cov('DefenseFinder')
    r=analyze_saved(d)
    renamed=deepcopy(d)
    for arr in [renamed['genes'],renamed['calls']]:
        for row in arr:row['gene_id']={'x':'opaque|α__b','y':'opaque|α__c'}[row['gene_id']]
    assert analyze_saved(renamed)['global_summary']==r['global_summary']
    assert r['global_summary']['integrated_CDS']==2


def test_cross_sample_group_and_duplicate_group_rejected():
    d=base('a');d['calls']=[call('a')];d['coverage']=cov('DefenseFinder');g=group('g',['a']);g['members'][0]['sample_id']='other';d['groups']=[g]
    with pytest.raises(ValueError,match='Cross-sample'):analyze_saved(d)
    g['members'][0].pop('sample_id');d['groups']=[g,g]
    with pytest.raises(ValueError,match='Duplicate group'):analyze_saved(d)


def test_empty_project_truthful_available_lists():
    r=analyze_saved({'samples':[]})
    assert r['gene_ledger']==[] and r['cells']==[]
    assert r['availability']['effector']=='unassessed'
    assert r['global_summary']['integrated_CDS']==0


def test_cas6_processing_does_not_become_effector():
    from defenseviz.analysis.saved import _registry
    rule=next(r for r in _registry()[0]['component_rules']if r['component']=='CRISPR_Cas6')
    d=base('c');d['calls']=[call('c',caller=rule['provider'],model=rule['exact_models'][0])];d['coverage']=cov(rule['provider'])
    assert analyze_saved(d)['global_summary']['effector_CDS']==0


def test_automatic_fgh_requires_same_original_run_and_system():
    d=base('f','g','h');d['coverage']=cov('DefenseFinder')
    d['calls']=[call(x,model='SspBCDE__Ssp'+x.upper(),system='one',run_id='r')for x in ['f','g','h']]
    r=analyze_saved(d)
    assert r['global_summary']['complex_core_CDS']==3
    assert r['global_summary']['individual_effector_CDS']==0
    d['calls'][2]['system_id']='two'
    r=analyze_saved(d)
    assert r['global_summary']['complex_core_CDS']==0
    # Same apparent member list cannot borrow components from another source system.
    d['groups']=[dict(group('all',['f','g','h']),native_system_id='one')]
    with pytest.raises(ValueError,match='exact provider/run/system'):analyze_saved(d)


def test_supplement_cannot_spoof_original_native_group_or_component_transfer():
    d=base('d');d['calls']=[call('d',model='Wadjet__JetD_II',origin='supplement')];d['coverage']=cov('DefenseFinder')
    assert analyze_saved(d)['global_summary']['effector_CDS']==0
    d['groups']=[dict(group('g',['d']),origin='supplement')]
    with pytest.raises(ValueError,match='Supplement cannot'):analyze_saved(d)


def test_saved_decision_availability_is_per_sample():
    d=base('g');d['samples'].append(dict(sample_id='other',label='Other'))
    d['genes'].append(gene('h','other'));d['saved_decisions']=[decision('g')]
    r=analyze_saved(d)
    assert r['availability']['saved_decisions']=='partial'
    assert {x['sample_id']:x['status']for x in r['availability']['saved_decisions_by_sample']}=={'s':'complete','other':'unassessed'}


def test_hmm_endpoint_completion_does_not_complete_final_overlap():
    d=base('g');d['calls']=[call('g',kind='hmm_hit')];d['coverage']=[dict(sample_id='s',caller='DefenseFinder',endpoint='hmm_hits',status='complete')]
    r=analyze_saved(d)
    assert r['overlap_regions']['coverage_status']=='unassessed'


def test_same_group_id_in_different_samples_does_not_merge():
    d=base('a');d['samples'].append(dict(sample_id='t',label='Other'));d['genes'].append(gene('b','t'))
    d['calls']=[call('a'),dict(call('b'),sample_id='t')];d['coverage']=cov('DefenseFinder')+[dict(sample_id='t',caller='DefenseFinder',status='complete')]
    d['groups']=[group('1',['a']),dict(group('1',[]),sample_id='t',members=[{'gene_id':'b','protein_sha256':gene('b','t')['protein_sha256']}])]
    assert analyze_saved(d)['global_summary']['integrated_CDS']==2


def test_general_predictions_only_fallback_not_duplicate_classified_families():
    d=base('a','b');d['calls']=[call('a',family='Wadjet_II',model='Wadjet__JetD_II'),call('a','DefensePredictor',family='General defense'),call('b','DefensePredictor',family='General defense')];d['coverage']=cov('DefenseFinder','DefensePredictor')
    r=analyze_saved(d);fs={v['id']:v for v in r['inventory_rows']['families']}
    assert fs['Wadjet II']['N']==1
    assert fs['General defense only']['member_keys']==[['s','b']]
    assert next(c for c in r['inventory_rows']['categories']if c['id']=='Nuclease')['N']==1


def test_RM_source_counts_not_backfilled_and_prediction_context_labeled():
    d=base('n','d','p');d['calls']=[dict(call('n',family='RM_Type_I'),category='R-M',rm_type='Type I'),call('n','DefensePredictor',family='General defense'),call('p','DefensePredictor',family='General defense')];d['coverage']=cov('DefenseFinder','DefensePredictor');d['saved_decisions']=[decision('d',role='M')]
    r=analyze_saved(d);s=r['sample_summary'][0]
    assert s['rm_source_counts']=={'DefenseFinder':1,'DefensePredictor':1}
    assert s['rm_pool_CDS']==2
    assert next(x for x in r['inventory_rows']['families']if x['id']=='R-M')['N']==2


def test_shared_type_pool_includes_unretained_native_partner():
    d=base('r','m','orphan');d['calls']=[dict(call('r'),category='R-M',rm_type='Type III'),dict(call('m'),category='R-M',rm_type='Type III')];d['coverage']=cov('DefenseFinder')
    d['saved_decisions']=[decision('r',typ='Type III'),decision('m',role='M',typ='Type III',state='not_retained'),decision('orphan',role='M',typ='Type III')]
    t=next(t for t in analyze_saved(d)['type_tables']if t['rm_type']=='Type III')
    assert (t['unique_CDS'],t['pool_CDS'],t['source_counts']['DefenseFinder'])==(2,3,2)
    assert t['role_states']['S']=='not_applicable'


def test_evidence_status_typo_is_error_not_unassessed():
    d=base('a');d['calls']=[call('a')];d['effector_evidence']=[effect('a')];d['effector_evidence'][0]['status']='supportted'
    with pytest.raises(ValueError,match='status'):analyze_saved(d)


def test_completed_empty_retention_vs_missing_file():
    d=base();d['suite']={'capabilities':[{'sample_id':'s','saved_retention':True,'native_candidates':True,'candidate_count':0}]}
    assert analyze_saved(d)['availability']['saved_decisions_by_sample'][0]['status']=='complete'
    d['suite']['capabilities'][0]['saved_retention']=False
    assert analyze_saved(d)['availability']['saved_decisions_by_sample'][0]['status']=='unassessed'


def test_duplicate_core_component_does_not_assign_all_copies_as_required():
    d=base('f','f2','g','h');d['coverage']=cov('DefenseFinder')
    d['calls']=[call(k,model='SspBCDE__Ssp'+c,system='one',run_id='r')for k,c in [('f','F'),('f2','F'),('g','G'),('h','H')]]
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==0


def test_declared_saved_retention_coverage_does_not_become_native_source():
    d=base();d['coverage']=[dict(sample_id='s',caller='DNMB',endpoint='saved_retention',status='complete')]
    r=analyze_saved(d)
    assert r['sources']==[]
    assert r['availability']['saved_decisions_by_sample'][0]['status']=='complete'


def test_complex_explicit_hash_map_cannot_conflict_with_current_protein():
    d=base('f');d['calls']=[call('f')];d['coverage']=cov('DefenseFinder');d['groups']=[group('g',['f'])]
    d['complexes']=[dict(complex_id='c',group_id='g',sample_id='s',status='complete',core_gene_ids=['f'],core_protein_sha256={'f':'0'*64},source_refs=[REF],rule_id='source_review')]
    with pytest.raises(ValueError,match='core protein hash mismatch'):analyze_saved(d)


def test_automatic_complex_requires_no_explicit_unresolved_replicon():
    d=base('f','g','h');d['coverage']=cov('DefenseFinder')
    d['calls']=[call(x,model='SspBCDE__Ssp'+x.upper(),system='one',run_id='r')for x in ['f','g','h']]
    d['groups']=[dict(group('original',['f','g','h']),native_system_id='one',replicon_binding='unassessed')]
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==0
    d['groups'][0]['replicon_binding']='verified'
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==3
    # A historical explicit source group with no newly introduced binding field
    # preserves its original semantics; absence is not a rewritten assessment.
    d['groups'][0].pop('replicon_binding')
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==3
    d['groups'][0]['replicon_binding']='unassessed'
    d['complexes']=[dict(complex_id='independent_review',group_id='original',sample_id='s',status='complete',core_gene_ids=['f','g','h'],source_refs=[REF],rule_id='bound_review',core_protein_sha256={g:gene(g)['protein_sha256']for g in ['f','g','h']})]
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==3


def test_coverage_without_endpoint_preserves_reference_only_observed_kind():
    d=base('a');d['calls']=[call('a','ReferenceDB',kind='reference_hit')];d['coverage']=cov('ReferenceDB')
    r=analyze_saved(d)
    assert r['sources'][0]['kind']=='reference_annotation'
    assert r['global_summary']['native_union_CDS']==0
    d['coverage'][0]['endpoint']='native_calls'
    assert analyze_saved(d)['sources'][0]['kind']=='mixed'


def test_type_pool_recognized_types_and_total_partition():
    d=base('nonrm','iig','unknown');d['coverage']=cov('DefenseFinder')
    d['calls']=[call('nonrm',rm_type='not_applicable'),call('iig',rm_type='IIG',category_override='unused'),call('unknown',rm_type='new_type')]
    d['calls'][1]['category']='R-M';d['calls'][2]['category']='R-M'
    r=analyze_saved(d);rows={x['rm_type']:x for x in r['type_tables']}
    assert rows['Type II']['pool_member_keys']==[['s','iig']]
    assert rows['Unresolved']['pool_member_keys']==[['s','unknown']]
    assert rows['Total']['pool_CDS']==2
    assert set(map(tuple,rows['Total']['pool_member_keys']))==set().union(*(set(map(tuple,x['pool_member_keys']))for x in r['type_tables']if not x['is_total']))
    assert rows['Total']['source_counts']['DefenseFinder']==sum(x['source_counts']['DefenseFinder']for x in r['type_tables']if not x['is_total'])


def test_unassessed_type_cannot_borrow_other_type_decision():
    d=base('a','b');d['coverage']=cov('DefenseFinder')
    d['calls']=[call('a',rm_type='Type I'),call('b',rm_type='Type II')]
    d['saved_decisions']=[decision('a',typ='Type I',role=None,state='unassessed'),decision('b',role='M')]
    row=next(x for x in analyze_saved(d)['type_tables']if x['rm_type']=='Type I')
    assert row['role_states']['R']=='unassessed'
    assert row['context_state']=='unassessed' and row['incomplete_R_context'] is None
    assert row['orphan_M'] is None and row['orphan_state']['M']=='unassessed'


def test_typeiv_orphan_inapplicability_is_explicit():
    d=base('a');d['saved_decisions']=[decision('a',typ='Type IV')]
    row=next(x for x in analyze_saved(d)['type_tables']if x['rm_type']=='Type IV')
    assert row['orphan_M'] is None and row['orphan_S'] is None
    assert row['orphan_state']=={'M':'not_applicable','S':'not_applicable'}


@pytest.mark.parametrize('typ,role',[('Type IV','R'),('Type I','R')])
def test_conflicting_saved_and_reviewed_assignments_rejected(typ,role):
    d=base('a');d['saved_decisions']=[decision('a',typ='Type I',role='M')]
    d['annotations']=[dict(**gene('a'),status='supported',assignment_kind='reviewed_RM_component',family=typ,role=role,rule_id='review',source_refs=[REF])]
    with pytest.raises(ValueError,match='Contradictory saved/reviewed'):analyze_saved(d)


@pytest.mark.parametrize('field,value',[('native_system_id','different'),('run_id','different'),('provider','different')])
def test_native_group_requires_bound_original_membership(field,value):
    d=base('a');d['calls']=[call('a',system='one',run_id='r')];d['coverage']=cov('DefenseFinder')
    g=group('one',['a'],native_system_id='one');g[field]=value;d['groups']=[g]
    with pytest.raises(ValueError,match='Native group member'):analyze_saved(d)


def test_explicit_catalog_keeps_ids_and_receives_display_defaults():
    d=base('f','m');d['calls']=[call('f'),call('m')];d['coverage']=cov('DefenseFinder')
    d['options']={'family_catalog':[{'id':'SspFGH (PT)','category':'Nuclease','label':'SspFGH (PT)'},{'id':'SspABCD (PT)','category':'Modification','label':'SspABCD (PT)'}],
      'family_memberships':[dict(**gene('f'),family_id='SspFGH (PT)',source_refs=[REF]),dict(**gene('m'),family_id='SspABCD (PT)',source_refs=[REF])]}
    rows={r['id']:r for r in analyze_saved(d)['inventory_rows']['families']}
    assert rows['SspFGH (PT)']['label']=='SspE / SspFGH (PT)'
    assert rows['SspFGH (PT)']['display_group']=='Modification'
    assert len(rows['SspFGH (PT)']['band_segments'])==2
    assert rows['SspABCD (PT)']['partner_order']==0
    assert rows['SspFGH (PT)']['partner_order']==1
    d['options']['family_catalog'][0]['display_group']='Explicit custom'
    assert analyze_saved(d)['inventory_rows']['families'][0]['display_group']=='Explicit custom'


@pytest.mark.parametrize('provider,model',[
 ('DefenseFinder','Gabija__GajA'),('PADLOC','GajA_00025'),
 ('DefenseFinder','cas3_I_5'),('DefenseFinder','cas3HD_I_1'),('PADLOC','cd09641'),
 ('DefenseFinder','cas9_II-A_II-B_II-C_3'),('PADLOC','COG3513'),
 ('DefenseFinder','pAgo__pAgo_LongA'),('PADLOC','pAgo_WP_096885432.1')])
def test_named_native_nuclease_annotations_not_motif_or_activity(provider,model):
    d=base('a');d['calls']=[call('a',caller=provider,model=model)];d['coverage']=cov(provider)
    r=analyze_saved(d)
    assert r['global_summary']['individual_effector_CDS']==1
    assert r['global_summary']['motif_CDS']==0 and r['gene_ledger'][0]['activity']=='unassessed'
    d['calls'][0]['evidence_kind']='hmm_hit'
    assert analyze_saved(d)['global_summary']['individual_effector_CDS']==0


@pytest.mark.parametrize('provider,model',[
 ('PADLOC','cd09710'),('DefenseFinder','Gabija__GajB'),('DefenseFinder','pAgo__pAgo_Short'),
 ('PADLOC','PIWI'),('DefenseFinder','cas3_I_5_unreviewed_alias'),('DefenseFinder','BREX__pglZA')])
def test_noncatalytic_or_unspecified_neighbors_not_promoted(provider,model):
    d=base('a');d['calls']=[call('a',caller=provider,model=model)];d['coverage']=cov(provider)
    assert analyze_saved(d)['global_summary']['individual_effector_CDS']==0


def test_new_named_component_respects_explicit_noncatalytic_review():
    d=base('a');d['calls']=[call('a',model='Gabija__GajA')];d['coverage']=cov('DefenseFinder')
    d['annotations']=[dict(**gene('a'),status='supported',assignment_kind='noncatalytic_component',rule_id='specific_nonfunctional_fusion_review',source_refs=[REF])]
    assert analyze_saved(d)['global_summary']['individual_effector_CDS']==0


def test_dnd_complete_core_is_distinct_from_individual_catalysis():
    d=base('f','g','h','b');d['coverage']=cov('DefenseFinder')
    d['calls']=[call(x,model='Dnd_ABCDEFGH__Dpt'+x.upper(),system='one',run_id='r')for x in ['f','g','h']]+[call('b',model='Dnd__DndB',system='one',run_id='r')]
    r=analyze_saved(d)
    assert (r['global_summary']['complex_core_CDS'],r['global_summary']['individual_effector_CDS'],r['global_summary']['associated_CDS'])==(3,0,4)
    assert not next(v for v in r['gene_ledger']if v['gene_id']=='b')['effector']
    d['calls']=[c for c in d['calls']if c['gene_id']!='h']
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==0


@pytest.mark.parametrize('provider,model',[('DefenseFinder','Wadjet__JetD_I'),('PADLOC','JetD1_00014'),('PADLOC','JetD2_00003'),('PADLOC','JetD3_00005')])
def test_official_jetd_component_variants(provider,model):
    d=base('a');d['calls']=[call('a',caller=provider,model=model)];d['coverage']=cov(provider)
    assert analyze_saved(d)['global_summary']['individual_effector_CDS']==1


def test_padloc_exact_dnd_module_components_require_complete_group():
    d=base('f','g','h');d['coverage']=cov('PADLOC')
    d['calls']=[call(g,caller='PADLOC',model=m,system='one',run_id='r')for g,m in [('f','DndF_00003'),('g','DndG_00004'),('h','DndH_00007')]]
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==3
    d['calls'][2]['system_id']='two'
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==0


def test_noncatalytic_veto_cannot_silently_claim_union_membership():
    d=base('a');d['annotations']=[dict(**gene('a'),status='supported',assignment_kind='noncatalytic_component',include_in_union=True,rule_id='motor_review',source_refs=[REF])]
    with pytest.raises(ValueError,match='veto, not an included membership'):analyze_saved(d)


def test_type_and_retained_columns_are_unassessed_without_saved_decisions():
    d=base('a');d['calls']=[call('a',family='RM_Type_I',rm_type='Type I')];d['coverage']=cov('DefenseFinder')
    r=analyze_saved(d);row=next(x for x in r['type_tables']if x['rm_type']=='Type I');sample=r['sample_summary'][0]
    assert row['unique_CDS']==0 and row['counts_status']=='unassessed'
    assert row['pool_CDS']==1 and row['pool_status']=='complete'
    assert sample['retained_CDS']==0 and sample['retained_counts_status']=='unassessed'
    assert sample['retained_ratio_status']=='unassessed'


def test_type_counts_partial_review_cannot_claim_complete_saved_scope():
    d=base('a','b');d['calls']=[call('a',rm_type='Type I'),call('b',rm_type='Type I')];d['coverage']=cov('DefenseFinder')
    d['saved_decisions']=[decision('a',typ='Type I'),decision('b',typ='Type I',state='unassessed')]
    r=analyze_saved(d);row=next(x for x in r['type_tables']if x['rm_type']=='Type I')
    assert row['unique_CDS']==1 and row['counts_status']=='partial'
    assert r['sample_summary'][0]['retained_ratio_status']=='partial'
    d['coverage'].append(dict(sample_id='s',caller='PADLOC',status='missing',endpoint='native_calls'))
    assert next(x for x in analyze_saved(d)['type_tables']if x['rm_type']=='Type I')['pool_status']=='partial'


def test_retained_ratio_unknown_role_is_not_complete():
    d=base('a');d['saved_decisions']=[decision('a',role=None)]
    sample=analyze_saved(d)['sample_summary'][0]
    assert sample['retained_counts_status']=='complete'
    assert sample['retained_RM_counts_status']=='complete'
    assert sample['retained_R_counts_status']=='partial' and sample['retained_ratio_status']=='partial'


def test_complete_saved_zero_can_license_type_zero_only_for_that_sample():
    d=base('a');d['calls']=[call('a',rm_type='Type I')];d['coverage']=cov('DefenseFinder')
    d['saved_decisions']=[decision('a',typ='Type I',state='not_retained')]
    d['samples'].append(dict(sample_id='other',label='Other'))
    d['genes'].append(gene('b','other'))
    d['calls'].append(dict(call('b',rm_type='Type I'),sample_id='other'))
    d['coverage'].append(dict(sample_id='other',caller='DefenseFinder',status='complete'))
    rows={x['sample_id']:x for x in analyze_saved(d)['type_tables']if x['rm_type']=='Type I'}
    assert rows['s']['unique_CDS']==0 and rows['s']['counts_status']=='complete'
    assert rows['other']['unique_CDS']==0 and rows['other']['counts_status']=='unassessed'


def review_group(gid,members):
    return dict(group(gid,members),provider='Review',membership_kind='saved_context',origin='supplement')


def test_source_bound_review_group_projects_current_union_without_retention():
    d=base('native','reviewed','outside');d['calls']=[call('native')];d['coverage']=cov('DefenseFinder')
    d['annotations']=[dict(**gene('reviewed'),status='supported',assignment_kind='reviewed_RM_component',family='Type IV',role='R',include_in_union=True,rule_id='reviewed_component',source_refs=[REF])]
    d['effector_evidence']=[effect('reviewed')];d['groups']=[review_group('r',['native','reviewed','outside'])]
    r=analyze_saved(d)
    assert r['global_summary']['retained_CDS']==0 and r['global_summary']['associated_CDS']==2
    assert r['groups_selected'][0]['membership_scope']=='reviewed_current_integrated_members'
    assert r['groups_selected'][0]['member_keys']==[['s','native'],['s','reviewed']]


def test_review_complete_core_needs_current_union_and_never_automatic_native_rule():
    d=base('f','g','h');d['calls']=[call(x,model='SspBCDE__Ssp'+x.upper(),system='one',run_id='r')for x in ['f','g','h']];d['coverage']=cov('DefenseFinder')
    d['groups']=[review_group('review',['f','g','h'])]
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==0
    d['complexes']=[dict(complex_id='reviewed-core',sample_id='s',group_id='review',status='complete',core_gene_ids=['f','g','h'],source_refs=[REF],rule_id='explicit_complex_review')]
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==3
    d['calls']=d['calls'][:2]
    with pytest.raises(ValueError,match='outside exact group or integrated'):analyze_saved(d)


@pytest.mark.parametrize('missing',['hash','source','retained'])
def test_review_group_authentication_cannot_be_bypassed(missing):
    d=base('a');d['calls']=[call('a')];d['coverage']=cov('DefenseFinder');g=review_group('r',['a']);d['groups']=[g]
    if missing=='hash':g['members'][0].pop('protein_sha256')
    if missing=='source':g['source_refs']=[]
    if missing=='retained':g['members'][0]['retained']=True
    with pytest.raises(ValueError):analyze_saved(d)


def test_unknown_native_run_never_licenses_automatic_complete_core():
    d=base('f','g','h');d['calls']=[call(x,model='SspBCDE__Ssp'+x.upper(),system='one')for x in ['f','g','h']];d['coverage']=cov('DefenseFinder')
    r=analyze_saved(d)
    assert r['global_summary']['complex_core_CDS']==0
    # An explicit group with missing run cannot borrow a known run either.
    for c in d['calls']:c['run_id']='known'
    d['groups']=[dict(group('g',['f','g','h']),run_id=None,native_system_id='one')]
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==0


def test_explicit_incomplete_saved_export_caps_complete_observed_records():
    d=base('a');d['saved_decisions']=[decision('a')]
    d['coverage']=[dict(sample_id='s',caller='DNMB',endpoint='saved_retention',status='partial')]
    r=analyze_saved(d)
    assert r['availability']['saved_decisions_by_sample'][0]['status']=='partial'
    assert r['sample_summary'][0]['retained_counts_status']=='partial'
    assert next(x for x in r['type_tables']if x['rm_type']=='Type II')['counts_status']=='partial'
    d['coverage'].append(dict(sample_id='s',caller='DNMB',endpoint='saved_decisions',status='complete'))
    with pytest.raises(ValueError,match='Conflicting duplicate saved-decision coverage'):analyze_saved(d)


def test_reviewed_only_component_can_join_explicit_review_complex():
    d=base('native','reviewed');d['calls']=[call('native')];d['coverage']=cov('DefenseFinder')
    d['annotations']=[dict(**gene('reviewed'),status='supported',assignment_kind='reviewed_RM_component',family='Type IV',role='R',include_in_union=True,rule_id='component_review',source_refs=[REF])]
    d['groups']=[review_group('context',['native','reviewed'])]
    d['complexes']=[dict(complex_id='reviewed-core',sample_id='s',group_id='context',status='complete',core_gene_ids=['native','reviewed'],source_refs=[REF],rule_id='complex_review')]
    r=analyze_saved(d)
    assert r['global_summary']['complex_core_CDS']==2
    assert r['global_summary']['individual_effector_CDS']==0
    assert r['global_summary']['retained_CDS']==0
    assert r['groups_selected'][0]['membership_scope']=='reviewed_current_integrated_members'


def test_known_group_run_cannot_fill_unknown_member_run():
    d=base('f','g','h');d['coverage']=cov('DefenseFinder')
    d['calls']=[call(x,model='SspBCDE__Ssp'+x.upper(),system='one',run_id='r')for x in ['f','g','h']]
    d['calls'][0].pop('run_id')
    d['groups']=[dict(group('g',['f','g','h']),native_system_id='one')]
    assert analyze_saved(d)['global_summary']['complex_core_CDS']==0
    d['calls'][0]['run_id']='different'
    with pytest.raises(ValueError,match='Native group member conflicts'):analyze_saved(d)

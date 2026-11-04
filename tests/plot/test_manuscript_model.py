from copy import deepcopy
import hashlib

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.plot.publication.manuscript_model import build_manuscript_atlas_model, enrich_manuscript_model


REF={"path":"saved.tsv","sha256":"a"*64,"row":1}


def fixture():
    ids=["r","m","s","other","excluded","iv"]
    d=dict(samples=[dict(sample_id="sample",label="Example strain",genus="Example",species="species",highlight=True,
                        display={"species_color":"#123456","annotation_tags":[{"label":"User group","status":"user_asserted"}]})],
        genes=[dict(sample_id="sample",gene_id=g,protein_sha256=hashlib.sha256(g.encode()).hexdigest())for g in ids],
        calls=[],coverage=[],saved_decisions=[],annotations=[],groups=[],options={"manuscript_atlas":{"group_coverage":{"sample":"complete"}}})
    for source in ["PADLOC","DefenseFinder","DefensePredictor"]:
        d["coverage"].append(dict(sample_id="sample",caller=source,status="complete",endpoint="native_calls"))
        for g in ["r","m","s"]:
            d["calls"].append(dict(sample_id="sample",gene_id=g,caller=source,evidence_kind="final_member",family="RM_Type_I",category="R-M",rm_type="Type I",source_ref=REF))
    d["calls"].append(dict(sample_id="sample",gene_id="other",caller="DefenseFinder",evidence_kind="final_member",family="Other family",category="Nuclease",source_ref=REF))
    d["coverage"].append(dict(sample_id="sample",caller="ReferenceDB",status="complete",endpoint="reference_matches"))
    for g in ["m","excluded"]:
        d["calls"].append(dict(sample_id="sample",gene_id=g,caller="ReferenceDB",evidence_kind="reference_hit",family="Type I",category="R-M",rm_type="Type I",source_ref=REF))
    for g,role in [("r","R"),("m","M"),("s","S"),("excluded","R")]:
        d["saved_decisions"].append(dict(sample_id="sample",gene_id=g,state="not_retained"if g=="excluded"else"retained",family="Type I",role=role,source_ref=REF))
    h={g["gene_id"]:g["protein_sha256"]for g in d["genes"]}
    d["annotations"].append(dict(sample_id="sample",gene_id="iv",protein_sha256=h["iv"],status="supported",assignment_kind="reviewed_RM_component",family="Type IV",role="R",include_in_union=True,rule_id="bound_review",source_refs=[REF]))
    d["groups"].append(dict(sample_id="sample",group_id="context",provider="DNMB",run_id="run",membership_kind="saved_context",family="Type I",context_state="complete",
        members=[dict(gene_id=g,protein_sha256=h[g])for g in ["r","m","s"]],source_refs=[REF]))
    return d


def evidence(d,kind="rm_context_display",gene="r",**kwargs):
    return dict(sample_id="sample",gene_id=gene,protein_sha256=next(g["protein_sha256"]for g in d["genes"]if g["gene_id"]==gene),
        assignment_kind="reviewed_RM_component",status="unassessed",include_in_union=False,rule_id="display_only:"+kind+":"+kwargs.get("group_id",""),evidence_kind=kind,
        source_refs=[REF],**kwargs)


def model(d):
    result=analyze_saved(d)
    return build_manuscript_atlas_model(enrich_manuscript_model(d,result))


def test_original_sources_and_three_source_intersection_are_not_reference_intersection():
    d=fixture();result=analyze_saved(d);original=deepcopy(result);enriched=enrich_manuscript_model(d,result);m=build_manuscript_atlas_model(enriched);row=m["rows"][0]
    assert result==original and {k:v for k,v in enriched.items()if k!="manuscript_atlas"}==original
    assert m["source_order"]==["PADLOC","DefenseFinder","DefensePredictor"]
    assert row["rm_original"]["source_counts"]==dict(PADLOC=3,DefenseFinder=3,DefensePredictor=3)
    assert row["rm_original"]["reference_counts"]=={"ReferenceDB":2}
    assert row["rm_original"]["input_pool_CDS"]==4
    assert row["rm_original"]["intersection_CDS"]==3
    assert result["sample_summary"][0]["rm_intersection_CDS"]==1
    assert row["retained"]==dict(R_CDS=1,total_CDS=3,status="complete")
    assert row["rm_roles"]["IV_R"]["count"]==1


def test_all_saved_compositions_survive_without_any_effector_seed():
    d=fixture();m=model(d);row=m["rows"][0]
    assert row["compositions"]["I_RMS"]["count"]==1
    assert row["compositions"]["I_RMS"]["R_support_all"]is None
    assert row["compositions"]["I_RMS"]["reference_any"]is None
    assert row["rm_roles"]["I_R"]["support_status"]=="unassessed"
    d["saved_decisions"][0]["state"]="not_retained"
    m=model(d)
    assert m["rows"][0]["compositions"]["I_MS"]["count"]==1
    assert m["rows"][0]["compositions"]["I_RMS"]["count"]==0


def test_explicit_supported_and_false_badges_are_separate_from_effector_and_raw_reference():
    d=fixture()
    d["annotations"].append(evidence(d,group_id="context",rm_R_support="supported",rm_reference_link="not_supported"))
    d["annotations"].append(evidence(d,kind="rm_role_display",rm_R_support="supported"))
    r=analyze_saved(d);m=build_manuscript_atlas_model(enrich_manuscript_model(d,r));row=m["rows"][0]
    assert r["global_summary"]["effector_CDS"]==0
    assert row["compositions"]["I_RMS"]["R_support_all"]is True
    assert row["compositions"]["I_RMS"]["reference_any"]is False
    assert row["compositions"]["I_RMS"]["reference_status"]=="complete"
    assert row["rm_roles"]["I_R"]["support_count"]==1
    d["annotations"][-2]["rm_R_support"]="not_supported"
    assert model(d)["rows"][0]["compositions"]["I_RMS"]["R_support_all"]is False


@pytest.mark.parametrize("defect",["hash","source","member","admitting","state"])
def test_display_evidence_is_bound_and_cannot_assign_biology(defect):
    d=fixture();a=evidence(d,group_id="context",rm_R_support="supported",rm_reference_link="not_supported")
    if defect=="hash":a["protein_sha256"]="0"*64
    if defect=="source":a["source_refs"]=[]
    if defect=="member":a["group_id"]="another"
    if defect=="admitting":a["status"]="supported"
    if defect=="state":a["rm_R_support"]="maybe"
    result=analyze_saved(d);d["annotations"].append(a)
    with pytest.raises(ValueError):enrich_manuscript_model(d,result)


def test_methylation_only_context_cannot_inherit_complete_RM_reference_badge():
    d=fixture();d["saved_decisions"][0]["state"]="not_retained"
    d["annotations"].append(evidence(d,gene="m",group_id="context",rm_R_support="not_supported",rm_reference_link="supported"))
    with pytest.raises(ValueError,match="M/MS-only"):model(d)


def test_exact_original_context_dedup_preserves_distinct_boundaries():
    d=fixture();same=deepcopy(d["groups"][0]);same["group_id"]="reciprocal";d["groups"].append(same)
    assert model(d)["rows"][0]["compositions"]["I_RMS"]["count"]==1
    d["groups"][-1]["members"].append(dict(gene_id="excluded",protein_sha256=next(g["protein_sha256"]for g in d["genes"]if g["gene_id"]=="excluded")))
    assert model(d)["rows"][0]["compositions"]["I_RMS"]["count"]==2


def test_masks_are_exact_nonRM_sets_with_explicit_order_and_zero_mask():
    d=fixture();d["options"]["manuscript_atlas"]["source_order"]=["DefenseFinder","PADLOC","DefensePredictor"]
    m=model(d);row=m["rows"][0];masks={x["mask"]:x for x in row["noncanonical"]["masks"]}
    assert row["noncanonical"]["pool_CDS"]==1 and masks["100"]["count"]==1
    assert sum(x["count"]for x in masks.values())==1 and masks["000"]["count"]==0
    assert row["display"]["species_color"]=="#123456" and row["highlight"]is True
    assert all(f["category"]!="R-M"for f in m["families"])


def test_partial_missing_channels_are_not_false_zero():
    d=fixture();d["coverage"][-1]["status"]="partial"
    row=model(d)["rows"][0]
    assert row["rm_original"]["intersection_count_status"]=="complete"
    assert row["rm_original"]["intersection_status"]=="partial"
    d["groups"]=[];d["options"]["manuscript_atlas"].pop("group_coverage")
    assert model(d)["rows"][0]["compositions"]["I_RMS"]["status"]=="unassessed"


def test_identifier_permutation_preserves_every_numeric_endpoint():
    d=fixture();before=model(d);mapping={g["gene_id"]:"renamed/"+str(i)for i,g in enumerate(d["genes"])}
    for table in ["genes","calls","saved_decisions","annotations"]:
        for row in d[table]:row["gene_id"]=mapping[row["gene_id"]]
    for g in d["groups"]:
        for member in g["members"]:member["gene_id"]=mapping[member["gene_id"]]
    after=model(d)
    for field in ["rm_original","retained"]:assert before["rows"][0][field]==after["rows"][0][field]
    assert {k:v["count"]for k,v in before["rows"][0]["rm_roles"].items()}=={k:v["count"]for k,v in after["rows"][0]["rm_roles"].items()}


def test_result_only_fallback_does_not_relabel_all_source_intersection():
    result=analyze_saved(fixture());m=build_manuscript_atlas_model(result)
    assert m["rows"][0]["rm_original"]["intersection_CDS"]is None
    assert m["rows"][0]["compositions"]["I_RMS"]["status"]=="unassessed"


def test_source_membership_mismatch_rejected_before_display():
    d=fixture();result=analyze_saved(d);d["calls"].pop()
    with pytest.raises(ValueError,match="source memberships"):enrich_manuscript_model(d,result)


def test_declared_absence_is_known_zero_not_missing_and_legend_is_declared():
    d=fixture();d['options']['family_catalog']=[dict(id='Other family',category='Nuclease',label='Other family')]
    d['options']['manuscript_atlas'].update(absence_families=['Other family'],absence_label='selected family',
        taxonomy_legend={'genera':[{'key':'Example','label':'E.','color':'#123456'}]})
    assert model(d)['rows'][0]['display']['absence_family_ids']==[]
    d['calls']=[c for c in d['calls']if c['gene_id']!='other']
    m=model(d)
    assert m['rows'][0]['display']['absence_family_ids']==['Other family']
    assert m['taxonomy_legend']['genera'][0]['label']=='E.'
    assert m['taxonomy_legend']['genera'][0]['genomes']==1
    d['coverage'][0]['status']='partial'
    assert model(d)['rows'][0]['display']['absence_family_ids']==[]


def test_family_exclusive_counts_require_the_same_declared_source_scope():
    d=fixture();m=model(d);cell=m['rows'][0]['cells']['Other family']
    assert cell['exclusive_source_counts']==dict(PADLOC=0,DefenseFinder=1,DefensePredictor=0)
    assert cell['exclusive_source_status']=='complete'
    d['calls'].append(dict(d['calls'][-3],gene_id='other',caller='PADLOC',family='Other family',category='Nuclease',rm_type=None))
    m=model(d)
    assert all(v==0 for v in m['rows'][0]['cells']['Other family']['exclusive_source_counts'].values())


def test_known_fused_role_is_component_membership_not_two_unique_genes():
    d=fixture();d['saved_decisions'][0]['role']='RM';d['saved_decisions'][1]['state']='not_retained'
    m=model(d);row=m['rows'][0]
    assert row['rm_roles']['I_R']['count']==row['rm_roles']['I_M']['count']==1
    assert row['compositions']['I_RMS']['count']==1
    assert len(row['compositions']['I_RMS']['member_keys'])==2


def test_conflicting_exact_context_metadata_is_not_row_order_dependent():
    d=fixture();same=deepcopy(d['groups'][0]);same['group_id']='reciprocal';same['context_state']='incomplete';d['groups'].append(same)
    for groups in [d['groups'],list(reversed(d['groups']))]:
        d['groups']=groups
        with pytest.raises(ValueError,match='Contradictory exact saved-context'):model(d)


def test_curated_reference_badge_does_not_create_context_completeness():
    d=fixture();d['groups'][0]['context_state']='unassessed'
    d['annotations'].append(evidence(d,group_id='context',rm_R_support='supported',rm_reference_link='supported'))
    with pytest.raises(ValueError,match='explicitly complete'):model(d)


def test_result_only_preserves_declared_reference_column_and_never_invents_counts():
    result=analyze_saved(fixture());m=build_manuscript_atlas_model(result)
    assert m['reference_sources']==['ReferenceDB']
    assert m['rows'][0]['rm_original']['reference_counts']=={'ReferenceDB':2}
    assert m['rows'][0]['rm_original']['reference_status']=={'ReferenceDB':'complete'}
    result['sample_summary'][0]['rm_source_counts'].pop('ReferenceDB')
    result['sample_summary'][0]['rm_source_count_status'].pop('ReferenceDB')
    m=build_manuscript_atlas_model(result)
    assert m['reference_sources']==['ReferenceDB']
    assert m['rows'][0]['rm_original']['reference_counts']=={'ReferenceDB':None}
    assert m['rows'][0]['rm_original']['reference_status']=={'ReferenceDB':'unassessed'}


def test_result_only_retained_total_preserves_scope_and_prefers_explicit_total():
    result=analyze_saved(fixture());s=result['sample_summary'][0]
    s['retained_RM_CDS']=2
    retained=build_manuscript_atlas_model(result)['rows'][0]['retained']
    assert retained['total_CDS']==3 and retained['total_scope']=='all_saved_retained'
    s.pop('retained_CDS')
    retained=build_manuscript_atlas_model(result)['rows'][0]['retained']
    assert retained['total_CDS']==2 and retained['total_scope']=='saved_retained_RM_only'
    s.pop('retained_RM_CDS');s.pop('retained_ratio_status')
    retained=build_manuscript_atlas_model(result)['rows'][0]['retained']
    assert retained['total_CDS']is None and retained['status']=='unassessed'


def test_context_excluded_retained_member_is_not_reintroduced_or_badged():
    d=fixture();d['groups'][0]['members'][0]['retained']=False
    row=model(d)['rows'][0]
    assert row['compositions']['I_MS']['count']==1
    assert row['compositions']['I_RMS']['count']==0
    d['annotations'].append(evidence(d,group_id='context',rm_R_support='supported'))
    with pytest.raises(ValueError,match='exact retained group member'):model(d)


def test_reciprocal_context_with_conflicting_member_retention_is_rejected():
    d=fixture();other=deepcopy(d['groups'][0]);other['group_id']='reciprocal'
    other['members'][0]['retained']=False;d['groups'].append(other)
    for groups in [d['groups'],list(reversed(d['groups']))]:
        d['groups']=groups
        with pytest.raises(ValueError,match='Contradictory exact saved-context'):model(d)


def test_hmm_only_provider_does_not_enter_native_intersection_or_patterns():
    d=fixture()
    d['calls'].append(dict(sample_id='sample',gene_id='other',caller='ProfileLibrary',evidence_kind='hmm_hit',family='RM_Type_I',source_ref=REF))
    d['coverage'].append(dict(sample_id='sample',caller='ProfileLibrary',endpoint='hmm_hits',status='complete'))
    m=model(d)
    assert m['source_order']==['PADLOC','DefenseFinder','DefensePredictor']
    assert m['rows'][0]['rm_original']['intersection_CDS']==3
    assert m['rows'][0]['rm_original']['intersection_count_status']=='complete'
    assert len(m['rows'][0]['noncanonical']['masks'])==8
    d['options']['manuscript_atlas']['source_order']=['ProfileLibrary']
    with pytest.raises(ValueError,match='original native sources'):model(d)


def test_declared_empty_native_endpoint_is_kept_even_with_only_hmm_observations():
    d=fixture()
    d['calls'].append(dict(sample_id='sample',gene_id='other',caller='NewNative',evidence_kind='hmm_hit',family='RM_Type_I',source_ref=REF))
    d['coverage'].append(dict(sample_id='sample',caller='NewNative',endpoint='native_calls',status='complete'))
    m=model(d);row=m['rows'][0]
    assert 'NewNative' in m['source_order']
    assert row['rm_original']['source_counts']['NewNative']==0
    assert row['rm_original']['source_status']['NewNative']=='complete'
    assert row['rm_original']['intersection_CDS']==0


def test_result_fallback_excludes_profile_only_but_keeps_declared_missing_native():
    d=fixture()
    d['calls'].append(dict(sample_id='sample',gene_id='other',caller='Profiles',evidence_kind='hmm_hit',family='RM_Type_I',source_ref=REF))
    d['coverage'].append(dict(sample_id='sample',caller='AbsentNative',endpoint='native_calls',status='missing'))
    m=build_manuscript_atlas_model(analyze_saved(d))
    assert 'Profiles' not in m['source_order']
    assert 'AbsentNative' in m['source_order']
    assert m['rows'][0]['rm_original']['source_status']['AbsentNative']=='unassessed'


def test_missing_context_type_cannot_be_complete_zero():
    d=fixture();d['groups'][0]['family']=None
    m=model(d)
    assert all(c['status']=='partial' for c in m['rows'][0]['compositions'].values())
    assert m['group_records'][0]['family']is None
    assert any('unresolved type' in w for w in m['warnings'])


def test_group_evidence_unknown_merge_is_order_independent_but_known_conflicts_fail():
    d=fixture()
    a=evidence(d,group_id='context',rm_R_support='supported')
    b=evidence(d,group_id='context',rm_R_support='supported',rm_reference_link='supported')
    b['rule_id']='second_independent_display_record'
    observed=[]
    for pair in [[a,b],[b,a]]:
        dd=deepcopy(d);dd['annotations']+=pair
        observed.append(model(dd)['rows'][0]['compositions']['I_RMS'])
    assert observed[0]==observed[1]
    assert observed[0]['reference_any']is True
    b['rm_reference_link']='not_supported';a['rm_reference_link']='supported'
    for pair in [[a,b],[b,a]]:
        dd=deepcopy(d);dd['annotations']+=pair
        with pytest.raises(ValueError,match='Contradictory atlas group display evidence'):model(dd)


def test_mixed_general_provider_only_projects_its_general_calls_into_rm_pool():
    d=fixture();d['coverage'].append(dict(sample_id='sample',caller='MixedNative',endpoint='native_calls',status='complete'))
    for gene,family in [('m','General defense'),('r','CBASS')]:
        d['calls'].append(dict(sample_id='sample',gene_id=gene,caller='MixedNative',evidence_kind='final_member',family=family,category='Unresolved',source_ref=REF))
    row=model(d)['rows'][0]
    assert row['rm_original']['source_counts']['MixedNative']==1


def test_same_provider_reference_hits_do_not_backfill_its_native_rm_endpoint():
    d=fixture()
    d['calls'].append(dict(sample_id='sample',gene_id='excluded',caller='PADLOC',evidence_kind='reference_hit',family='Type I',category='R-M',rm_type='Type I',source_ref=REF))
    d['coverage'].append(dict(sample_id='sample',caller='PADLOC',endpoint='reference_matches',status='complete'))
    row=model(d)['rows'][0]
    assert row['rm_original']['source_counts']['PADLOC']==3
    assert row['rm_original']['reference_counts']['PADLOC']==1
    d['coverage'].pop()
    row=model(d)['rows'][0]
    assert row['rm_original']['source_status']['PADLOC']=='complete'
    assert row['rm_original']['reference_status']['PADLOC']=='unassessed'


def test_result_only_mixed_provider_preserves_both_channels_without_invented_split():
    d=fixture()
    d['calls'].append(dict(sample_id='sample',gene_id='excluded',caller='PADLOC',evidence_kind='reference_hit',family='Type I',category='R-M',rm_type='Type I',source_ref=REF))
    m=build_manuscript_atlas_model(analyze_saved(d));rm=m['rows'][0]['rm_original']
    assert 'PADLOC' in m['source_order'] and 'PADLOC' in m['reference_sources']
    assert rm['source_counts']['PADLOC']is None and rm['source_status']['PADLOC']=='unassessed'
    assert rm['reference_counts']['PADLOC']is None and rm['reference_status']['PADLOC']=='unassessed'


def test_retained_role_status_is_not_poisoned_by_unassessed_native_pool_members():
    d=fixture()
    d['calls'].append(dict(sample_id='sample',gene_id='other',caller='PADLOC',evidence_kind='final_member',family='RM_Type_I',category='R-M',rm_type='Type I',source_ref=REF))
    r=analyze_saved(d)
    assert next(x for x in r['type_tables']if x['rm_type']=='Type I')['role_states']['R']=='partial'
    row=enrich_manuscript_model(d,r)['manuscript_atlas']['rows'][0]
    assert all(row['rm_roles']['I_'+role]['status']=='complete'for role in ['R','M','S'])
    assert row['rm_roles']['I_R']['status_scope']=='saved_retained_roles'
    d['coverage'].append(dict(sample_id='sample',caller='DNMB',endpoint='saved_retention',status='partial'))
    assert model(d)['rows'][0]['rm_roles']['I_R']['status']=='partial'


def test_unknown_retained_role_or_type_remains_partial_for_retained_display():
    d=fixture();d['saved_decisions'][0]['role']=None
    row=model(d)['rows'][0]
    assert row['rm_roles']['I_R']['status']=='partial'
    assert row['rm_roles']['I_R']['count']==0
    d=fixture();d['saved_decisions'][0]['family']=None
    row=model(d)['rows'][0]
    assert all(row['rm_roles'][name]['status']=='partial'for name in ['I_R','II_R','III_R'])


@pytest.mark.parametrize('gene,field,value',[('s','role','unknown'),('r','role','unknown'),('m','family',None)])
def test_unknown_current_context_member_role_or_type_never_guesses_composition(gene,field,value):
    d=fixture();next(x for x in d['saved_decisions']if x['gene_id']==gene)[field]=value
    m=model(d);row=m['rows'][0]
    assert row['compositions']['I_MS']['count']==0
    assert row['compositions']['I_RMS']['count']==0
    assert all(c['status']=='partial'for c in row['compositions'].values())
    g=m['group_records'][0]
    assert g['composition_class']is None and g['composition_status']=='unassessed'
    assert len(g['member_keys'])==3


def test_excluded_unknown_role_member_does_not_make_current_context_unknown():
    d=fixture();d['groups'][0]['members'][0]['retained']=False;d['saved_decisions'][0]['role']='unknown'
    row=model(d)['rows'][0]
    assert row['compositions']['I_MS']['count']==1
    assert row['compositions']['I_MS']['status']=='complete'

"""Completed native evidence must be bound, scored exactly, and scope-explicit."""
from copy import deepcopy
import csv
from hashlib import sha256
import json
from pathlib import Path

import pytest

from defenseviz.full_evaluation import apply_full_evaluation, _profiles, _divergent, _fresh_hmm_role_ambiguity


def digest(path):return sha256(Path(path).read_bytes()).hexdigest()


def table(path,rows,fields=None):
    fields=fields or list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)


def save(path,obj):Path(path).write_text(json.dumps(obj))


def fixture(tmp_path,roles=('R','M','S')):
    native=tmp_path/'native';hmm=tmp_path/'hmm';native.mkdir();hmm.mkdir()
    library=hmm/'profiles.hmm';library.write_text('synthetic HMM model')
    meta=hmm/'profiles.tsv';table(meta,[dict(profile_name='p',roles='R',type_hints='I',profile_length=100,library_sha256=digest(library))])
    original=tmp_path/'original.gbff';original.write_text('synthetic original genome')
    genes=[];canonical=[];hypotheses=[];ledger=[]
    for i,role in enumerate(roles):
        gid='g'+str(i);seq='M'+'A'*(99+i)
        genes.append(dict(sample_id='s',gene_id=gid,protein_sha256=sha256(seq.encode()).hexdigest(),protein_length=len(seq)))
        canonical.append(dict(locus_tag=gid,translation=seq,contig='c',start=1000*i+1,end=1000*i+300,strand='+'))
        hypothesis=dict(query=gid,family_id='Type I',enzyme_role=role,score_total=75,
            factor_strong=True,factor_moderate=False,factor_profile=False,factor_annotation=False,
            factor_motif=True,factor_context=True,factor_structure=False,factor_alternative_annotation=False,
            factor_conflict=False,factor_partial=False,rease_nuclease_verified=role=='R',hmm_role_ambiguous=False,
            blast_identity='NA',blast_query_coverage='NA',blast_reference_coverage='NA',blast_evalue='NA',
            blast_role_compatible=False,blast_family_compatible='NA',curated_blast_identity='NA',
            curated_blast_query_coverage='NA',curated_blast_reference_coverage='NA',curated_blast_evalue='NA',
            curated_blast_role_compatible=False,curated_blast_family_compatible=False,
            rease_motor_verified=False,mtase_motif_verified=False,rease_motif_hits_raw='NA',
            typei_context_status='complete_rms',typei_context_partners=' | '.join('g'+str(j)+':'+r for j,r in enumerate(roles) if i!=j))
        hypotheses.append(hypothesis);ledger.append(dict(locus_tag=gid,completed=True))
    table(native/'canonical_genes.tsv',canonical);table(native/'native_hypotheses.tsv',hypotheses);table(native/'search_ledger.tsv',ledger)
    query=native/'query.faa';query.write_text(''.join('>'+r['locus_tag']+'\n'+r['translation']+'\n' for r in canonical))
    mapping=hmm/'gene_map.tsv';table(mapping,[dict(sample_id='s',query_id=g['gene_id'],locus_tag=g['gene_id'],aa_length=g['protein_length'],protein_sha256=g['protein_sha256']) for g in genes])
    output=native/'blast.tsv';output.write_text('')
    save(native/'search_commands.json',[dict(argv=['blastp','-query',str(query),'-out',str(output),'-db',str(original)],exit_code=0,query_sha256=digest(query),output_sha256=digest(output))])
    table(hmm/'all_domains.tsv',[],fields=['sample_id','locus_tag'])
    save(hmm/'execution.json',dict(exit_code=0,inputs_unchanged=True,hmm={'exit_code':0},resolution={'exit_code':0},
        signature=dict(sample_id='s',inputs={k:dict(path=str(v),sha256=digest(v)) for k,v in [('hmm',library),('profiles',meta),('fasta',query),('gene_map',mapping)]}),
        artifact_hashes={'all_domains.tsv':digest(hmm/'all_domains.tsv')}))
    save(native/'completion.json',dict(passed=True,full_reference_search=True,selected_policy_applied=False,final_candidates=len(roles),source_hashes={}))
    manifest=dict(input_sha256={str(original):digest(original),str((hmm/'execution.json').resolve()):digest(hmm/'execution.json')},completion_sha256=digest(native/'completion.json'),
        artifacts={p.name:digest(p) for p in native.glob('*') if p.is_file()})
    save(native/'stage_manifest.json',manifest)
    data=dict(samples=[dict(sample_id='s')],genes=genes,calls=[dict(sample_id='s',gene_id='g0',caller='PADLOC')],
        saved_decisions=[dict(sample_id='s',gene_id='g0',state='unassessed',family='Type IV')],
        groups=[dict(group_id='old',sample_id='s',provider='DNMB',membership_kind='saved_context',members=[])],
        coverage=[dict(sample_id='s',endpoint='saved_retention',status='unassessed')],suite={'features':[dict(old=True)]})
    return data,[dict(sample_id='s',native_dir=str(native),hmm_dir=str(hmm))]


def update_native(stages,name,change):
    p=Path(stages[0]['native_dir'])/name
    rows=list(csv.DictReader(p.open(),delimiter='\t'));change(rows);table(p,rows)
    m=p.parent/'stage_manifest.json';v=json.loads(m.read_text());v['artifacts'][name]=digest(p);save(m,v)


def bind_hmm(stages):
    d=Path(stages[0]['native_dir']);h=Path(stages[0]['hmm_dir'])/'execution.json'
    m=json.loads((d/'stage_manifest.json').read_text());m['input_sha256'][str(h.resolve())]=digest(h);save(d/'stage_manifest.json',m)


def test_complete_factors_exact_policy_current_groups_and_input_immutability(tmp_path):
    data,stages=fixture(tmp_path);before=deepcopy(data)
    out=apply_full_evaluation(data,stages)
    assert data==before and out['calls']==before['calls']
    assert len(out['saved_decisions'])==3 and all(r['state']=='retained' for r in out['saved_decisions'])
    assert all(r['score_total']==75 and r['score_min']==r['score_max'] for r in out['saved_decisions'])
    assert out['rm_evaluation']['native_saved_decisions']==before['saved_decisions']
    assert out['suite']['archived_native_context_groups']==before['groups']
    reviews=out['rm_evaluation']['current_context_review']
    assert len(reviews['records'])==1 and reviews['records'][0]['context_state']=='met'
    assert reviews['records'][0]['retained_roles']==['M','R','S']
    assert len(reviews['observed_inventory_coverage'])==2
    assert all(not r['decision']['hypothetical_evaluations'] for r in out['rm_evaluation']['records'])
    json.dumps(out,allow_nan=False)


def test_policy_supported_r_emits_bound_effector_and_current_associated_group(tmp_path):
    from defenseviz.analysis.saved import analyze_saved
    data,stages=fixture(tmp_path);before=deepcopy(data)
    out=apply_full_evaluation(data,stages)
    assert data==before and out['calls']==data['calls']
    evidence=out['effector_evidence'];assert len(evidence)==1
    row=evidence[0]
    assert row['gene_id']=='g0' and row['protein_sha256']==data['genes'][0]['protein_sha256']
    assert row['assignment_kind']=='individual_nuclease_effector'
    assert row['status']=='supported_putative_effector_annotation'
    assert row['policy_sha256']==out['rm_evaluation']['policy_sha256']
    assert row['policy_version']==out['rm_evaluation']['selected_policy']['version']
    assert row['motif_state']=='unassessed' and not row['validated_motif_site_inferred']
    for ref in row['source_refs']:
        assert digest(ref['path'])==ref['sha256']==out['derived_artifact_sha256'][ref['path']]
    out['samples'][0]=dict(out['samples'][0],label='test')
    result=analyze_saved(out)
    assert result['global_summary']['effector_CDS']==1
    assert result['global_summary']['associated_CDS']==3
    assert result['global_summary']['motif_CDS']==0


def test_retained_r_role_without_nuclease_route_is_not_effector(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',))
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(rease_nuclease_verified=False))
    out=apply_full_evaluation(data,stages)
    assert out['saved_decisions'][0]['state']=='retained'
    assert out['saved_decisions'][0]['integrated_R_component_support'] is False
    assert out['effector_evidence']==[]


def test_nuclease_route_without_retention_is_not_effector(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',))
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(factor_context=False,factor_motif=False,score_total=35))
    out=apply_full_evaluation(data,stages)
    assert out['saved_decisions'][0]['state']=='not_retained'
    assert out['saved_decisions'][0]['integrated_R_component_support'] is True
    assert out['effector_evidence']==[]


def add_numeric_r_profile(data,stages,extra=None):
    d=Path(stages[0]['hmm_dir']);receipt=json.loads((d/'execution.json').read_text())
    library=receipt['signature']['inputs']['hmm']['sha256'];sha=data['genes'][0]['protein_sha256']
    metadata=[dict(profile_name='p',roles='R',type_hints='II',profile_length=100,library_sha256=library)]
    domain=dict(sample_id='s',locus_tag='g0',profile_name='p',role_hints='R',type_hints='II',database_sha256=library,
        protein_sha256=sha,resolved_selected=True,off_target=False,mapping_status='unique',sequence_length=100,model_length=100,
        ali_from=1,ali_to=100,hmm_from=1,hmm_to=100,sequence_evalue=1e-40)
    rows=[domain]
    if extra:
        rows.append(dict(domain,profile_name='q',role_hints='M',**extra))
        metadata.append(dict(metadata[0],profile_name='q',roles='M'))
    table(d/'profiles.tsv',metadata);table(d/'all_domains.tsv',rows)
    receipt['signature']['inputs']['profiles']['sha256']=digest(d/'profiles.tsv')
    receipt['artifact_hashes']['all_domains.tsv']=digest(d/'all_domains.tsv');save(d/'execution.json',receipt)
    bind_hmm(stages)
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(family_id='Type II',rease_nuclease_verified=False,
        rease_motif_hits_raw='HNH',factor_motif=False,score_total=55,curated_blast_identity=1,
        curated_blast_query_coverage=1,curated_blast_reference_coverage=1,curated_blast_evalue=0,
        curated_blast_role_compatible=True,curated_blast_family_compatible=True,hmm_role_ambiguous=True,hmm_role_hints='R;M'))


def test_opposing_comparator_nomination_role_does_not_change_dnmb_score(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',));add_numeric_r_profile(data,stages)
    mixed=apply_full_evaluation(data,stages)
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(hmm_role_ambiguous=False,hmm_role_hints='R'))
    native=apply_full_evaluation(data,stages)
    a=mixed['rm_evaluation']['records'][0];b=native['rm_evaluation']['records'][0]
    assert a['decision']==b['decision'] and a['evidence']==b['evidence']
    assert a['decision']['score_total']==75 and a['decision']['integrated_R_component_support']
    assert a['hmm_role_ambiguity']['nomination_union_flag']=='True'
    assert a['hmm_role_ambiguity']['scoring_native_flag'] is False
    assert a['hmm_role_ambiguity']['native_role_hints']==['R']
    assert len(mixed['effector_evidence'])==1 and len(native['effector_evidence'])==1


def test_native_ambiguity_uses_even_weak_resolved_role_hint_without_score_gate(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',))
    add_numeric_r_profile(data,stages,extra=dict(ali_to=10,hmm_to=10,sequence_evalue=.1))
    out=apply_full_evaluation(data,stages);record=out['rm_evaluation']['records'][0]
    assert record['evidence']['hmm_role_ambiguous'] is True
    assert record['decision']['score_total']==55 and not record['decision']['integrated_R_component_support']
    assert out['effector_evidence']==[]


def test_ambiguity_excludes_offtarget_unresolved_and_nonunique_rows():
    base=dict(role_hints='R',resolved_selected=True,off_target=False,mapping_status='unique')
    rows=[(2,base),(3,dict(base,role_hints='M',off_target=True)),
          (4,dict(base,role_hints='M',resolved_selected=False)),(5,dict(base,role_hints='unknown',mapping_status='ambiguous'))]
    assert _fresh_hmm_role_ambiguity(rows)==(False,['R'])
    assert _fresh_hmm_role_ambiguity([])==(False,[])
    assert _fresh_hmm_role_ambiguity([(2,dict(base,role_hints='unknown'))])==(True,['unknown'])


def test_missing_captured_factor_is_rejected_instead_of_false(tmp_path):
    data,stages=fixture(tmp_path)
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(factor_structure='NA'))
    with pytest.raises(ValueError,match='executed boolean'):apply_full_evaluation(data,stages)


@pytest.mark.parametrize('column',['rease_nuclease_verified','rease_motor_verified','mtase_motif_verified',
                                  'rease_motif_hits_raw','curated_blast_identity','blast_role_compatible'])
def test_missing_deployment_column_never_becomes_false_or_null(tmp_path,column):
    data,stages=fixture(tmp_path)
    update_native(stages,'native_hypotheses.tsv',lambda rows:[r.pop(column) for r in rows])
    with pytest.raises(ValueError,match='deployment feature columns missing'):apply_full_evaluation(data,stages)


def test_present_native_na_measurements_keep_frozen_no_hit_semantics(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',))
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(
        rease_nuclease_verified='NA',rease_motor_verified='NA',mtase_motif_verified='NA'))
    out=apply_full_evaluation(data,stages)
    assert out['rm_evaluation']['records'][0]['evidence']['curated_identity'] is None
    assert out['saved_decisions'][0]['integrated_R_component_support'] is False


def test_unbound_hmm_execution_receipt_is_rejected(tmp_path):
    data,stages=fixture(tmp_path);d=Path(stages[0]['native_dir']);h=Path(stages[0]['hmm_dir'])/'execution.json'
    m=json.loads((d/'stage_manifest.json').read_text());m['input_sha256'].pop(str(h.resolve()));save(d/'stage_manifest.json',m)
    with pytest.raises(ValueError,match='does not bind.*HMM execution'):apply_full_evaluation(data,stages)


def test_substituted_hmm_receipt_cannot_rebind_itself(tmp_path):
    data,stages=fixture(tmp_path);h=Path(stages[0]['hmm_dir'])/'execution.json'
    receipt=json.loads(h.read_text());receipt['signature']['other_run']='different same-proteome run';save(h,receipt)
    with pytest.raises(ValueError,match='artifact hash mismatch'):apply_full_evaluation(data,stages)


def test_untranslated_context_member_has_partial_inventory_and_explicit_audit(tmp_path):
    data,stages=fixture(tmp_path);d=Path(stages[0]['native_dir']);h=Path(stages[0]['hmm_dir'])
    data['genes'][2].update(protein_sha256=None,protein_length=0)
    update_native(stages,'canonical_genes.tsv',lambda r:r[2].update(translation='NA'))
    update_native(stages,'native_hypotheses.tsv',lambda r:r[2].update(
        factor_strong=False,factor_motif=False,factor_context=False,score_total=0))
    update_native(stages,'search_ledger.tsv',lambda r:r.pop())
    query=d/'query.faa';lines=query.read_text().splitlines();query.write_text('\n'.join(lines[:4])+'\n')
    commands=json.loads((d/'search_commands.json').read_text());commands[0]['query_sha256']=digest(query);save(d/'search_commands.json',commands)
    m=json.loads((d/'stage_manifest.json').read_text());m['artifacts']['search_commands.json']=digest(d/'search_commands.json');save(d/'stage_manifest.json',m)
    mapping=list(csv.DictReader((h/'gene_map.tsv').open(),delimiter='\t'));table(h/'gene_map.tsv',mapping[:2])
    receipt=json.loads((h/'execution.json').read_text())
    for k,p in [('gene_map',h/'gene_map.tsv'),('fasta',query)]:receipt['signature']['inputs'][k]['sha256']=digest(p)
    save(h/'execution.json',receipt);bind_hmm(stages)
    out=apply_full_evaluation(data,stages);review=out['rm_evaluation']['current_context_review']
    coverage={r['family']:r for r in review['observed_inventory_coverage']}
    assert coverage['Type I']['status']=='partial' and coverage['Type III']['status']=='complete'
    assert len(review['omitted_context_groups'])==1
    omitted=review['omitted_context_groups'][0]
    assert omitted['untranslated_member_keys']==[['s','g2']]
    assert omitted['translated_member_keys']==[['s','g0'],['s','g1']]
    assert len(omitted['source_refs'])==4 and not omitted['orphan_inferred']
    scope=out['rm_evaluation']['candidate_scope'][0]
    assert scope['translated_count']==2 and scope['untranslated_gene_ids']==['g2']
    untranslated=next(r for r in out['rm_evaluation']['records'] if r['gene_id']=='g2')
    assert untranslated['translation_status']=='no_translated_protein' and untranslated['state']=='not_retained'


def test_completed_zero_nominees_retains_true_zero_scope(tmp_path):
    data,stages=fixture(tmp_path);d=Path(stages[0]['native_dir'])
    table(d/'native_hypotheses.tsv',[],fields=['query'])
    table(d/'search_ledger.tsv',[],fields=['locus_tag','completed']);save(d/'search_commands.json',[])
    completion=json.loads((d/'completion.json').read_text());completion['final_candidates']=0;save(d/'completion.json',completion)
    m=json.loads((d/'stage_manifest.json').read_text())
    for name in ('native_hypotheses.tsv','search_ledger.tsv','search_commands.json','completion.json'):m['artifacts'][name]=digest(d/name)
    m['completion_sha256']=digest(d/'completion.json');save(d/'stage_manifest.json',m)
    out=apply_full_evaluation(data,stages)
    assert out['rm_evaluation']['records']==[] and out['rm_evaluation']['candidate_scope'][0]['evaluated_count']==0
    assert all(r['status']=='complete' for r in out['rm_evaluation']['current_context_review']['observed_inventory_coverage'])


def test_captured_factor_score_must_equal_native_total(tmp_path):
    data,stages=fixture(tmp_path)
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(score_total=55))
    with pytest.raises(ValueError,match='factor/score'):apply_full_evaluation(data,stages)


def test_native_artifact_tampering_is_rejected(tmp_path):
    data,stages=fixture(tmp_path)
    p=Path(stages[0]['native_dir'])/'native_hypotheses.tsv';p.write_text(p.read_text().replace('75','74'))
    with pytest.raises(ValueError,match='artifact hash mismatch'):apply_full_evaluation(data,stages)


def test_native_gene_protein_must_match_canonical_hash(tmp_path):
    data,stages=fixture(tmp_path)
    data['genes'][0]['protein_sha256']='a'*64
    with pytest.raises(ValueError,match='protein hash'):apply_full_evaluation(data,stages)


def test_every_translated_nominee_needs_completed_full_search(tmp_path):
    data,stages=fixture(tmp_path)
    update_native(stages,'search_ledger.tsv',lambda r:r[0].update(completed=False))
    with pytest.raises(ValueError,match='lacks completed BLAST'):apply_full_evaluation(data,stages)


def test_prior_only_record_is_out_of_scope_not_biological_negative(tmp_path):
    data,stages=fixture(tmp_path)
    data['saved_decisions'].append(dict(sample_id='s',gene_id='old',state='retained',score_total=80,integrated_R_component_support=True))
    out=apply_full_evaluation(data,stages)
    old=next(r for r in out['saved_decisions'] if r['gene_id']=='old')
    assert old['state']=='out_of_scope' and 'score_total' not in old and 'integrated_R_component_support' not in old
    assert out['rm_evaluation']['candidate_scope'][0]['gene_ids']==['g0','g1','g2']


def test_unknown_role_is_evaluated_without_forced_type_or_negative(tmp_path):
    data,stages=fixture(tmp_path,roles=('NA',))
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(family_id='NA',typei_context_partners='NA'))
    out=apply_full_evaluation(data,stages)
    assert out['saved_decisions'][0]['state']=='retained'
    assert out['saved_decisions'][0]['role'] is None and out['saved_decisions'][0]['rm_type'] is None
    assert out['rm_evaluation']['records'][0]['evidence']['role']=='unknown'


def test_explicit_completed_type_iii_mod_singleton_has_current_m_context(tmp_path):
    data,stages=fixture(tmp_path,roles=('M',))
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(
        family_id='Type III',typeiii_context_status='missing_res_check_neighbors',
        typeiii_context_roles='M',typeiii_context_anchor='g0',typeiii_context_partners='NA'))
    out=apply_full_evaluation(data,stages)
    reviews=out['rm_evaluation']['current_context_review']['records']
    assert len(reviews)==1 and reviews[0]['family']=='Type III'
    assert reviews[0]['retained_roles']==['M'] and reviews[0]['required_roles']==['M','R']
    assert reviews[0]['boundary_source_bound'] and reviews[0]['dispositions_complete']
    assert reviews[0]['context_state']=='missing_retained_partner_in_saved_context'
    assert reviews[0]['orphan_inferred'] is False
    assert reviews[0]['retained_R_absent_in_saved_context'] is None
    groups=[g for g in out['groups'] if g['provider']=='Review']
    assert len(groups)==1 and [m['gene_id'] for m in groups[0]['members']]==['g0']


@pytest.mark.parametrize('change',[
    {'typeiii_context_status':'NA'}, {'typeiii_context_anchor':'other'},
    {'typeiii_context_roles':'NA'}, {'family_id':'Type I'}, {'enzyme_role':'R'}])
def test_missing_or_mismatched_singleton_context_never_invents_m_group(tmp_path,change):
    data,stages=fixture(tmp_path,roles=('M',))
    def change_row(rows):
        rows[0].update(family_id='Type III',typeiii_context_status='missing_res_check_neighbors',
            typeiii_context_roles='M',typeiii_context_anchor='g0',typeiii_context_partners='NA')
        rows[0].update(change)
    update_native(stages,'native_hypotheses.tsv',change_row)
    out=apply_full_evaluation(data,stages)
    assert not out['rm_evaluation']['current_context_review']['records']


def test_structure_points_need_exact_geometry_not_native_flag_only(tmp_path):
    data,stages=fixture(tmp_path)
    update_native(stages,'native_hypotheses.tsv',lambda r:r[0].update(factor_structure=True,score_total=100))
    with pytest.raises(ValueError,match='structure credit'):apply_full_evaluation(data,stages)


def test_no_protein_cannot_retain_despite_native_score(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',))
    data['genes'][0].update(protein_sha256=None,protein_length=0)
    update_native(stages,'canonical_genes.tsv',lambda r:r[0].update(translation='NA'))
    d=Path(stages[0]['hmm_dir']);fasta=d/'empty.faa';fasta.write_text('')
    table(d/'gene_map.tsv',[],fields=['sample_id','query_id','locus_tag','aa_length','protein_sha256'])
    receipt=json.loads((d/'execution.json').read_text())
    receipt['signature']['inputs']['fasta']=dict(path=str(fasta),sha256=digest(fasta))
    receipt['signature']['inputs']['gene_map']['sha256']=digest(d/'gene_map.tsv');save(d/'execution.json',receipt)
    bind_hmm(stages)
    with pytest.raises(ValueError,match='Untranslated'):apply_full_evaluation(data,stages)


def test_completed_hmm_receipt_required(tmp_path):
    data,stages=fixture(tmp_path)
    p=Path(stages[0]['hmm_dir'])/'execution.json';r=json.loads(p.read_text());r['resolution']['exit_code']=1;save(p,r)
    bind_hmm(stages)
    with pytest.raises(ValueError,match='HMM search/resolution'):apply_full_evaluation(data,stages)


def test_terminal_stop_is_removed_only_for_hmm_identity(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',))
    sequence='M'+'A'*99+'*';source_sha=sha256(sequence.encode()).hexdigest()
    data['genes'][0].update(protein_sha256=source_sha,protein_length=len(sequence))
    update_native(stages,'canonical_genes.tsv',lambda r:r[0].update(translation=sequence))
    out=apply_full_evaluation(data,stages)
    assert out['saved_decisions'][0]['protein_sha256']==source_sha
    binding=out['rm_evaluation']['records'][0]['protein_binding']
    assert binding['source_sha256']==binding['canonical_sha256']==source_sha
    assert binding['hmm_sha256']!=source_sha


def test_original_faa_stop_hash_can_differ_from_gbff_only_with_bound_receipt(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',));d=Path(stages[0]['native_dir'])
    source_sha=sha256(('M'+'A'*99+'*').encode()).hexdigest()
    data['genes'][0].update(protein_sha256=source_sha,protein_length=101)
    with pytest.raises(ValueError,match='protein hash'):apply_full_evaluation(data,stages)
    inputs=d.parent/'inputs';inputs.mkdir();p=inputs/'protein_binding.json'
    original=d.parent/'original.gbff'
    save(p,dict(sample_id='s',genbank=str(original),genbank_sha256=digest(original),normalization='uppercase whitespace-free; terminal stop removed for HMM only',source_protein_sha256={'g0':source_sha}))
    manifest=json.loads((d/'stage_manifest.json').read_text());manifest['input_sha256'][str(p.resolve())]=digest(p);save(d/'stage_manifest.json',manifest)
    out=apply_full_evaluation(data,stages)
    assert out['saved_decisions'][0]['protein_sha256']==source_sha


def test_internal_stop_or_amino_acid_change_is_not_terminal_normalization(tmp_path):
    data,stages=fixture(tmp_path,roles=('R',))
    sequence='MA*'+'A'*97
    data['genes'][0].update(protein_sha256=sha256(sequence.encode()).hexdigest(),protein_length=100)
    update_native(stages,'canonical_genes.tsv',lambda r:r[0].update(translation=sequence))
    with pytest.raises(ValueError,match='entire exact canonical'):apply_full_evaluation(data,stages)


def test_candidate_without_hmm_search_cannot_use_complete_scope(tmp_path):
    data,stages=fixture(tmp_path)
    d=Path(stages[0]['hmm_dir']);p=d/'gene_map.tsv'
    rows=list(csv.DictReader(p.open(),delimiter='\t'));table(p,rows[1:])
    receipt=json.loads((d/'execution.json').read_text());receipt['signature']['inputs']['gene_map']['sha256']=digest(p);save(d/'execution.json',receipt)
    bind_hmm(stages)
    with pytest.raises(ValueError,match='entire exact canonical'):apply_full_evaluation(data,stages)


def test_empty_no_hit_blast_still_requires_candidate_in_searched_fasta(tmp_path):
    data,stages=fixture(tmp_path);d=Path(stages[0]['native_dir'])
    query=d/'other_query.faa';query.write_text('>other\nMAAA\n')
    commands=json.loads((d/'search_commands.json').read_text());commands[0]['argv'][2]=str(query);commands[0]['query_sha256']=digest(query);save(d/'search_commands.json',commands)
    manifest=json.loads((d/'stage_manifest.json').read_text());manifest['artifacts']['search_commands.json']=digest(d/'search_commands.json');save(d/'stage_manifest.json',manifest)
    with pytest.raises(ValueError,match='absent from exact completed BLAST'):apply_full_evaluation(data,stages)


def test_exact_profile_union_does_not_double_count_overlaps_or_mix_profiles():
    gene=dict(protein_sha256='a'*64,protein_length=100)
    metadata={name:dict(roles='R',type_hints='II',library_sha256='b'*64,profile_length='100') for name in ('p','q')}
    base=dict(profile_name='p',role_hints='R',type_hints='II',database_sha256='b'*64,protein_sha256='a'*64,
              resolved_selected='true',off_target='false',mapping_status='unique',sequence_length=100,model_length=100,
              ali_from=1,ali_to=40,hmm_from=1,hmm_to=40,sequence_evalue='1e-30')
    rows=[(2,base),(3,dict(base,ali_from=30,ali_to=70,hmm_from=30,hmm_to=70)),(4,dict(base,profile_name='q',ali_from=71,ali_to=100,hmm_from=71,hmm_to=100))]
    profiles,details=_profiles(rows,gene,dict(family_id='Type II',enzyme_role='R'),metadata,'b'*64)
    assert profiles[0].query_coverage==.7 and profiles[1].query_coverage==.3
    assert profiles[0].model_coverage==.7 and details[0]['source_rows']==[2,3]
    with pytest.raises(ValueError,match='library binding'):
        _profiles([(2,dict(base,type_hints='III'))],gene,dict(family_id='Type II',enzyme_role='R'),metadata,'b'*64)


def test_divergent_context_uses_exact_cognate_and_fixed_800bp_bound():
    genes={'r':dict(locus_tag='r',contig='c',start=1,end=300,strand='-'),
           'm':dict(locus_tag='m',contig='c',start=1101,end=1400,strand='+')}
    common=dict(family_id='Type II',curated_blast_family='Type II',curated_blast_role_compatible=True,curated_blast_family_compatible=True,
                partial_status='full_length',evidence_conflict=False,alternative_function_annotation=False,
                curated_blast_identity=.5,curated_blast_query_coverage=.8,curated_blast_reference_coverage=.8,curated_blast_evalue=1e-5)
    rows=[dict(common,query='r',enzyme_role='R',curated_blast_role='R',curated_blast_resolved_subject_id='Example'),
          dict(common,query='m',enzyme_role='M',curated_blast_role='M',curated_blast_resolved_subject_id='M.Example')]
    assert _divergent(rows,genes)[0]=={'r','m'}
    genes['m']['start']=1102
    assert not _divergent(rows,genes)[0]

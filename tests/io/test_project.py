import csv
import hashlib
import json
from pathlib import Path

import pytest
from defenseviz.io import load_project


def table(path, headers, rows=()):
    with path.open('w', newline='', encoding='utf8') as f:
        writer = csv.writer(f, delimiter=',' if path.suffix == '.csv' else '\t')
        writer.writerow(headers)
        writer.writerows(rows)
    return path.name


def project(tmp_path, calls=None, **extra):
    table(tmp_path/'samples.tsv', ['sample_id','label','highlight'], [['S α','菌株 α','true'],['S/β','Isolate β','false']])
    cfg = {'schema_version':1,'samples':'samples.tsv','calls':calls or [],**extra}
    p = tmp_path/'project.json'
    p.write_text(json.dumps(cfg, ensure_ascii=False), encoding="utf-8")
    return p


def normalized(tmp_path, rows=()):
    return table(tmp_path/'calls.tsv',['sample_id','gene_id','caller','family','category','system_id','role','rm_type','evidence_kind'], rows)


def coverage(tmp_path, rows):
    return table(tmp_path/'coverage.tsv',['sample_id','caller','status'],rows)


def test_normalized_unicode_explicit_coverage_and_hashes(tmp_path, monkeypatch):
    calls=normalized(tmp_path,[['S α','g__x','Any caller','New family','Unresolved','','','','final_member']])
    cov=coverage(tmp_path,[['S α','Any caller','complete'],['S/β','Any caller','not_run']])
    p=project(tmp_path,[{'format':'normalized','path':calls}],coverage=cov)
    monkeypatch.chdir(tmp_path.parent)
    result=load_project(p)
    assert result['samples'][0]['label']=='菌株 α'
    assert result['samples'][0]['highlight'] is True
    assert result['calls'][0]['gene_id']=='g__x'
    assert result['coverage'][1]['status']=='not_run'
    assert set(result['inputs']) >= {p.resolve(),tmp_path/'samples.tsv',tmp_path/'calls.tsv',tmp_path/'coverage.tsv'}
    for path in result['inputs']:
        assert result['input_sha256'][str(path)]==hashlib.sha256(path.read_bytes()).hexdigest()


def test_multisample_requires_coverage(tmp_path):
    calls=normalized(tmp_path,[['S α','g1','Tool','X','Unresolved','','','','']])
    p=project(tmp_path,[{'format':'normalized','path':calls}])
    with pytest.raises(ValueError,match='explicit complete coverage'):
        load_project(p)


def test_header_only_complete_and_missing_are_distinct(tmp_path):
    calls=normalized(tmp_path)
    p=project(tmp_path,[{'format':'normalized','path':calls,'sample_id':'S α','caller':'Tool'}])
    result=load_project(p)
    assert result['calls']==[]
    assert result['coverage']==[{'sample_id':'S α','caller':'Tool','status':'complete'}, {'sample_id':'S/β','caller':'Tool','status':'missing'}]


def test_defensefinder_genes_exact_map(tmp_path):
    f=table(tmp_path/'df.tsv',['hit_id','sys_id','type','subtype','gene_name'],[['protein-1','sys1','RM_Type_I','Subtype X','uninterpreted-domain']])
    m=table(tmp_path/'map.tsv',['source_id','gene_id'],[['protein-1','locus:1']])
    p=project(tmp_path,[{'format':'defensefinder','path':f,'sample_id':'S α','id_map':m}])
    r=load_project(p)
    c=r['calls'][0]
    assert c['gene_id']=='locus:1'
    assert c['family']=='Subtype X' and c['category']=='R-M' and c['rm_type']=='Type I'
    assert c['role'] is None and c['evidence_kind']=='final_member'
    assert tmp_path/'map.tsv' in r['inputs']


def test_defensefinder_systems_expand_explicit_members(tmp_path):
    f=table(tmp_path/'systems.tsv',['sys_id','type','subtype','protein_in_syst','genes_count'],[['x','Novel','Novel-A','g1,g2','2']])
    (tmp_path/f).write_text('# macsyfinder metadata\n'+(tmp_path/f).read_text(encoding="utf-8"), encoding="utf-8")
    p=project(tmp_path,[{'format':'defensefinder','path':f,'sample_id':'S α'}])
    r=load_project(p)
    assert [x['gene_id'] for x in r['calls']]==['g1','g2']
    assert {x['system_id'] for x in r['calls']}=={'x'}
    assert {x['category'] for x in r['calls']}=={'Unresolved'}


def test_bad_system_membership_count_rejected(tmp_path):
    f=table(tmp_path/'systems.tsv',['sys_id','type','protein_in_syst','genes_count'],[['x','Novel','g1,g2','3']])
    p=project(tmp_path,[{'format':'defensefinder','path':f,'sample_id':'S α'}])
    with pytest.raises(ValueError,match='genes_count'):
        load_project(p)


def test_padloc_system_numbers_are_contig_scoped(tmp_path):
    f=table(tmp_path/'pad.csv',['system.number','seqid','system','target.name'],[['1','chrA','Mystery','g1'],['1','chrB','Mystery','g2']])
    p=project(tmp_path,[{'format':'padloc','path':f,'sample_id':'S α'}])
    r=load_project(p)
    assert len({x['system_id'] for x in r['calls']})==2
    assert all(x['category']=='Unresolved' for x in r['calls'])


def test_defensepredictor_threshold_is_explicit_and_zero_preserves_coverage(tmp_path):
    f=table(tmp_path/'dp.csv',['locus_tag','mean_log_odds'],[['g1',-8],['g2',1.5]])
    spec={'format':'defensepredictor','path':f,'sample_id':'S α'}
    p=project(tmp_path,[spec])
    with pytest.raises(ValueError,match='explicit.*score_threshold'):
        load_project(p)
    spec['score_threshold']=1.5
    p=project(tmp_path,[spec]); r=load_project(p)
    assert [x['gene_id'] for x in r['calls']]==['g2']
    assert r['calls'][0]['system_id'] is None
    spec['score_threshold']=100
    r=load_project(project(tmp_path,[spec]))
    assert not r['calls']
    assert r['coverage'][0]['status']=='complete'


@pytest.mark.parametrize('score',['nan','inf','bad',''])
def test_defensepredictor_bad_score_not_skipped(tmp_path,score):
    f=table(tmp_path/'dp.csv',['locus_tag','mean_log_odds'],[['g1',score]])
    p=project(tmp_path,[{'format':'defensepredictor','path':f,'sample_id':'S α','score_threshold':0}])
    with pytest.raises(ValueError,match='finite number'):
        load_project(p)


def test_family_mapping_exact_aliases(tmp_path):
    f=table(tmp_path/'df.tsv',['hit_id','sys_id','type','subtype'],[['g1','s1','Something','Version_A']])
    m=table(tmp_path/'families.tsv',['caller','family','category','display_family'],[['DefenseFinder','Version_A','Nuclease','My family']])
    r=load_project(project(tmp_path,[{'format':'defensefinder','path':f,'sample_id':'S α'}],family_map=m))
    assert r['calls'][0]['family']=='My family' and r['calls'][0]['native_family']=='Version_A'
    assert r['calls'][0]['role'] is None


def test_sample_binding_conflict_rejected(tmp_path):
    f=normalized(tmp_path,[['S/β','g1','T','X','Unresolved','','','','']])
    with pytest.raises(ValueError,match='conflicts'):
        load_project(project(tmp_path,[{'format':'normalized','path':f,'sample_id':'S α','caller':'T'}]))


def test_map_collision_and_missing_rejected(tmp_path):
    f=table(tmp_path/'df.tsv',['hit_id','sys_id','type'],[['g1','x','Unknown']])
    m=table(tmp_path/'map.tsv',['source_id','gene_id'],[['g1','L'],['g2','L']])
    spec={'format':'defensefinder','path':f,'sample_id':'S α','id_map':m}
    with pytest.raises(ValueError,match='one-to-one'):
        load_project(project(tmp_path,[spec]))
    table(tmp_path/'map.tsv',['source_id','gene_id'],[['other','L']])
    with pytest.raises(ValueError,match='no entry'):
        load_project(project(tmp_path,[spec]))


def test_duplicate_calls_not_silently_deduplicated(tmp_path):
    f=table(tmp_path/'df.tsv',['hit_id','sys_id','type'],[['g1','x','X'],['g1','x','X']])
    with pytest.raises(ValueError,match='duplicate call'):
        load_project(project(tmp_path,[{'format':'defensefinder','path':f,'sample_id':'S α'}]))


def test_complete_source_conflicting_missing_coverage_rejected(tmp_path):
    f=table(tmp_path/'df.tsv',['hit_id','sys_id','type'])
    cov=coverage(tmp_path,[['S α','DefenseFinder','missing']])
    with pytest.raises(ValueError,match='conflicts'):
        load_project(project(tmp_path,[{'format':'defensefinder','path':f,'sample_id':'S α'}],coverage=cov))


def test_dnmb_evidence_explicit_binding_ignores_previous_decision(tmp_path):
    pth=tmp_path/'e.json'; pth.write_text(json.dumps([{'target_id':'opaque__not_a_sample','evidence':{'family':'Type I'},'decision':{'keep':True}}]), encoding="utf-8")
    m=table(tmp_path/'binding.tsv',['target_id','sample_id','gene_id'],[['opaque__not_a_sample','S α','g8']])
    r=load_project(project(tmp_path,evidence={'format':'dnmb','path':'e.json','id_map':m}))
    assert r['evidence_records'][0]['sample_id']=='S α'
    assert r['evidence_records'][0]['gene_id']=='g8'
    assert 'decision' not in r['evidence_records'][0]


def test_target_only_evidence_requires_binding(tmp_path):
    (tmp_path/'e.json').write_text('[{"target_id":"S α__g1","evidence":{}}]', encoding="utf-8")
    with pytest.raises(ValueError,match='requires id_map'):
        load_project(project(tmp_path,evidence={'format':'dnmb','path':'e.json'}))


def test_duplicate_and_nonfinite_json_rejected(tmp_path):
    p=project(tmp_path)
    p.write_text('{"schema_version":1,"schema_version":1}', encoding="utf-8")
    with pytest.raises(ValueError,match='duplicate JSON key'):
        load_project(p)
    p.write_text('{"schema_version":NaN}', encoding="utf-8")
    with pytest.raises(ValueError,match='non-finite'):
        load_project(p)


@pytest.mark.parametrize('contents',['sample_id\tsample_id\nA\tB\n','sample_id\tlabel\nA\n',''])
def test_malformed_table_rejected(tmp_path,contents):
    p=project(tmp_path)
    (tmp_path/'samples.tsv').write_text(contents, encoding="utf-8")
    with pytest.raises(ValueError):
        load_project(p)


def test_invalid_gene_coordinates_and_hash_rejected(tmp_path):
    g=table(tmp_path/'genes.tsv',['sample_id','gene_id','contig','start','end','strand'],[['S α','g1','c','10','2','+']])
    with pytest.raises(ValueError,match='coordinates'):
        load_project(project(tmp_path,genes=g))


def test_xlsx_explicit_sheet_and_formula_policy(tmp_path):
    openpyxl=pytest.importorskip('openpyxl')
    w=openpyxl.Workbook(); w.active.title='samples'; w.active.append(['sample_id','label']);w.active.append(['A','Alpha']);w.create_sheet('other')
    w.save(tmp_path/'samples.xlsx')
    p=project(tmp_path)
    cfg=json.loads(p.read_text(encoding="utf-8"));cfg['samples']='samples.xlsx';p.write_text(json.dumps(cfg), encoding="utf-8")
    with pytest.raises(ValueError,match='specify sheet'):
        load_project(p)
    cfg['samples']={'path':'samples.xlsx','sheet':'samples'};p.write_text(json.dumps(cfg), encoding="utf-8")
    assert load_project(p)['samples'][0]['sample_id']=='A'
    w['samples']['B2']='=1+1';w.save(tmp_path/'samples.xlsx')
    with pytest.raises(ValueError,match='formulas'):
        load_project(p)


def test_raw_dnmb_features_missing_axes_stay_unassessed(tmp_path):
    f=table(tmp_path/'native.tsv',['locus_tag','rm_type','subunit'],[['g1','Type II','M']])
    r=load_project(project(tmp_path,evidence={'format':'dnmb_features','path':f,'sample_id':'S α'}))
    record=r['evidence_records'][0]
    assert record['evidence'] is None and record['status']=='unassessed'
    assert 'profiles' in record['provenance']['missing_fields']
    assert 'factor_motif' in record['provenance']['missing_fields']


def test_normalized_null_evidence_preserved(tmp_path):
    (tmp_path/'e.json').write_text(json.dumps([{'sample_id':'S α','gene_id':'g1','evidence':None,'provenance':{'reason':'no coordinates'}}]), encoding="utf-8")
    assert load_project(project(tmp_path,evidence='e.json'))['evidence_records'][0]['evidence'] is None


def complete_feature_row():
    from defenseviz.io.project import FEATURE_REQUIRED, FEATURE_FLAGS
    r={k:'' for k in FEATURE_REQUIRED}
    r.update({k:'false' for k in FEATURE_FLAGS})
    r.update(family_id='Type II',enzyme_role='M',sequence_length_aa='100',blast_family_compatible='',
             profiles='[]',divergent_context='false',factor_strong='true',factor_context='true')
    return r


def test_complete_native_factors_map_to_typed_evidence(tmp_path):
    row=complete_feature_row();row['gene_id']='g1'
    f=table(tmp_path/'features.tsv',list(row),[list(row.values())])
    r=load_project(project(tmp_path,evidence={'format':'dnmb_features','path':f,'sample_id':'S α'}))
    e=r['evidence_records'][0]['evidence']
    assert e['old_strong'] is True and e['old_context'] is True
    assert e['profiles']==[] and e['role']=='M'
    from defenseviz.analysis import integrate
    a=integrate(r['samples'],r['calls'],r['evidence_records'],coverage=r['coverage'])
    assert a['decisions'][0]['decision']['score_total']==55


def test_declared_feature_bad_numeric_not_silently_missing(tmp_path):
    row=complete_feature_row();row.update(gene_id='g1',blast_identity='invalid')
    f=table(tmp_path/'features.tsv',list(row),[list(row.values())])
    with pytest.raises(ValueError,match='invalid declared DNMB features'):
        load_project(project(tmp_path,evidence={'format':'dnmb_features','path':f,'sample_id':'S α'}))


def test_declared_feature_missing_boolean_stays_unassessed(tmp_path):
    row=complete_feature_row();row.update(gene_id='g1',factor_motif='NA')
    f=table(tmp_path/'features.tsv',list(row),[list(row.values())])
    r=load_project(project(tmp_path,evidence={'format':'dnmb_features','path':f,'sample_id':'S α'}))
    assert r['evidence_records'][0]['evidence'] is None
    assert 'factor_motif' in r['evidence_records'][0]['provenance']['missing_fields']


def test_extension_json_inputs_are_tracked(tmp_path):
    for name,data in [('registry.json',{'families':[]}),('candidates.json',[]),('rules.json',{'rules':[]})]:
        (tmp_path/name).write_text(json.dumps(data), encoding="utf-8")
    extension={'registry':'registry.json','candidates':'candidates.json','admission_rules':'rules.json','source_root':'.'}
    r=load_project(project(tmp_path,evidence_extension=extension))
    assert r['evidence_extension']['registry_dir']==tmp_path
    assert r['evidence_extension']['source_root']==tmp_path
    assert all(tmp_path/f in r['inputs'] for f in ['registry.json','candidates.json','rules.json'])


def test_sparse_json_native_row_gives_actionable_error(tmp_path):
    (tmp_path/'df.json').write_text(json.dumps([{'hit_id':'g1','sys_id':'s1','type':'Unknown'}, {'hit_id':'g2','type':'Unknown'}]), encoding="utf-8")
    with pytest.raises(ValueError,match='sys_id'):
        load_project(project(tmp_path,[{'format':'defensefinder','path':'df.json','sample_id':'S α'}]))


def test_duplicate_negative_defensepredictor_rows_rejected(tmp_path):
    f=table(tmp_path/'dp.csv',['locus_tag','mean_log_odds'],[['g1',-2],['g1',-3]])
    with pytest.raises(ValueError,match='duplicate DefensePredictor'):
        load_project(project(tmp_path,[{'format':'defensepredictor','path':f,'sample_id':'S α','score_threshold':0}]))


def test_typo_in_source_option_rejected(tmp_path):
    f=normalized(tmp_path)
    with pytest.raises(ValueError,match='unknown source keys'):
        load_project(project(tmp_path,[{'format':'normalized','path':f,'sample_id':'S α','caller':'X','idmap':'wrong.tsv'}]))


def test_shipped_family_catalog_and_user_override(tmp_path):
    f=table(tmp_path/'df.tsv',['hit_id','sys_id','type','subtype'],[['g1','s1','Gabija','Gabija']])
    spec={'format':'defensefinder','path':f,'sample_id':'S α'}
    r=load_project(project(tmp_path,[spec]))
    assert r['calls'][0]['category']=='Nuclease'
    assert r['calls'][0]['native_family']=='Gabija'
    assert r['calls'][0]['category_source']=='builtin_catalog'
    assert any(p.name=='families.json' for p in r['inputs'])
    mapping=table(tmp_path/'override.tsv',['family','category'],[['Gabija','Other mechanisms']])
    r=load_project(project(tmp_path,[spec],family_map=mapping))
    assert r['calls'][0]['category']=='Other mechanisms'
    assert r['calls'][0]['category_source']=='family_map'


def test_hmm_only_source_does_not_imply_native_completeness(tmp_path):
    f=normalized(tmp_path,[['S α','g1','Tool','New family','Unresolved','','','','hmm_hit']])
    r=load_project(project(tmp_path,[{'format':'normalized','path':f,'sample_id':'S α','caller':'Tool'}]))
    assert r['coverage'][0]['status']=='missing'
    assert r['calls'][0]['evidence_kind']=='hmm_hit'
    from defenseviz.analysis import integrate
    result=integrate(r['samples'],r['calls'],r['evidence_records'],coverage=r['coverage'])
    assert result['sample_summary'][0]['caller_counts']['Tool'] is None
    assert len(result['evidence_calls'])==1


def test_reference_only_multisample_table_needs_no_native_completion(tmp_path):
    f=normalized(tmp_path,[['S α','g1','Reference source','New family','Unresolved','','','','reference_hit']])
    r=load_project(project(tmp_path,[{'format':'normalized','path':f}]))
    assert r['coverage']==[]
    from defenseviz.analysis import integrate
    result=integrate(r['samples'],r['calls'],r['evidence_records'],coverage=r['coverage'])
    assert result['callers']==[]
    assert result['evidence_calls'][0]['caller']=='Reference source'


def test_relocated_project_keeps_analysis_provenance_identical(tmp_path):
    import shutil
    one=tmp_path/'one';one.mkdir()
    f=table(one/'df.tsv',['hit_id','sys_id','type'],[['g1','s1','Novel']])
    p=project(one,[{'format':'defensefinder','path':f,'sample_id':'S α'}])
    initial=load_project(p)
    shutil.copytree(one,tmp_path/'two')
    moved=load_project(tmp_path/'two/project.json')
    assert initial['calls']==moved['calls']
    assert initial['calls'][0]['source']=='df.tsv'
    assert initial['inputs']!=moved['inputs']


@pytest.mark.parametrize('category',['Other mechanisms','Unresolved','Nuclease'])
def test_explicit_normalized_category_out_ranks_builtin(tmp_path,category):
    f=normalized(tmp_path,[['S α','g1','DefenseFinder','Gabija',category,'s1','','','final_member']])
    r=load_project(project(tmp_path,[{'format':'normalized','path':f,'sample_id':'S α','caller':'DefenseFinder'}]))
    assert r['calls'][0]['category']==category
    assert r['calls'][0]['category_source']=='explicit'
    assert '_category_explicit' not in r['calls'][0]


def test_boolean_source_id_is_not_missing_default(tmp_path):
    f=table(tmp_path/'df.tsv',['hit_id','sys_id','type'],[['g1','s1','X']])
    with pytest.raises(ValueError,match='caller'):
        load_project(project(tmp_path,[{'format':'defensefinder','path':f,'sample_id':'S α','caller':False}]))


def test_dnmb_existing_record_and_mapping_hash_conflict_rejected(tmp_path):
    (tmp_path/'e.json').write_text(json.dumps([{'target_id':'opaque','sequence_sha256':'a'*64,'evidence':{}}]), encoding="utf-8")
    m=table(tmp_path/'binding.tsv',['target_id','sample_id','gene_id','sequence_sha256'],[['opaque','S α','g1','b'*64]])
    with pytest.raises(ValueError,match='sequence_sha256 conflicts'):
        load_project(project(tmp_path,evidence={'format':'dnmb','path':'e.json','id_map':m}))
    table(tmp_path/'binding.tsv',['target_id','sample_id','gene_id','sequence_sha256'],[['opaque','S α','g1','a'*64]])
    r=load_project(project(tmp_path,evidence={'format':'dnmb','path':'e.json','id_map':m}))
    assert r['evidence_records'][0]['provenance']['sequence_sha256']=='a'*64


def test_native_feature_hash_is_bound_for_extension(tmp_path):
    from defenseviz.analysis.admission import entry_sha256
    from defenseviz.pipeline import evaluate_project
    sequence='ADCEFGHIKLMN';sha=hashlib.sha256(sequence.encode()).hexdigest()
    (tmp_path/'proof.txt').write_text('Synthetic evidence only.\n', encoding="utf-8");source_sha=hashlib.sha256((tmp_path/'proof.txt').read_bytes()).hexdigest()
    model=dict(accession='TEST.1',source='synthetic',model_sha256='1'*64,evidence_kind='nuclease_domain',role='nuclease',family_specific=True,min_model_coverage=.7)
    reference=dict(accession='REF',sequence=sequence,sequence_sha256=sha,numbering_origin='1-based full reference sequence',independent_mapping_methods=['A','B'],catalytic_sites=[dict(position=2,residue='D',role='nuclease',evidence_source='synthetic'),dict(position=4,residue='E',role='nuclease',evidence_source='synthetic')])
    entry=dict(family_id='Synthetic enzyme',component='nuclease',substrate='DNA',domain_models=[model],references=[reference])
    registry=dict(schema_version='1.0',families=[entry])
    hit=dict(accession='TEST.1',source='synthetic',model_sha256='1'*64,query_sequence_sha256=sha,qstart=1,qend=12,hmm_start=1,hmm_end=12,hmm_length=12,target_length=12,independent_evalue=1e-30,qualification=dict(passed=True,rule='synthetic_only',source='proof.txt',source_sha256=source_sha))
    candidate=dict(candidate_id='opaque',sample_id='S α',query_id='g1',sequence=sequence,sequence_sha256=sha,domain_hits=[hit],site_maps=[dict(reference_accession='REF',reference_sequence_sha256=sha,query_sequence_sha256=sha,method=m,positions={'2':2,'4':4},source='proof.txt',source_sha256=source_sha) for m in ['A','B']])
    rule=dict(family_id='Synthetic enzyme',RM_families=['Type II'],roles=['R'],entry_sha256=entry_sha256(entry),model_accession='TEST.1',model_sha256='1'*64,reference_accession='REF',required_anchor_subset=[2,4],source='proof.txt',source_sha256=source_sha)
    for name,data in [('registry.json',registry),('candidates.json',[candidate]),('rules.json',[rule])]:
        (tmp_path/name).write_text(json.dumps(data), encoding="utf-8")
    row=complete_feature_row();row.update(gene_id='g1',sequence_sha256=sha,sequence_length_aa='12',enzyme_role='R',curated_blast_identity='.9',curated_blast_query_coverage='1',curated_blast_reference_coverage='1',curated_blast_evalue='1e-30',curated_blast_role_compatible='true',curated_blast_family_compatible='true')
    f=table(tmp_path/'features.tsv',list(row),[list(row.values())])
    p=project(tmp_path,evidence={'format':'dnmb_features','path':f,'sample_id':'S α'},evidence_extension={'registry':'registry.json','candidates':'candidates.json','admission_rules':'rules.json'})
    data=load_project(p)
    assert data['evidence_records'][0]['provenance']['sequence_sha256']==sha
    result=evaluate_project(data)
    assert result['reference_anchor_admission']['records'][0]['status']=='admitted'
    assert result['decisions'][0]['decision']['score_motif']==20
    row['sequence_sha256']='b'*64;table(tmp_path/'features.tsv',list(row),[list(row.values())])
    with pytest.raises(ValueError,match='sequence mismatch'):
        evaluate_project(load_project(p))


def test_reference_provider_only_becomes_native_caller_with_explicit_coverage(tmp_path):
    f=normalized(tmp_path,[['S α','g1','Custom reference provider','New family','Unresolved','','','','reference_hit']])
    spec={'format':'normalized','path':f,'sample_id':'S α','caller':'Custom reference provider'}
    r=load_project(project(tmp_path,[spec]))
    assert r['coverage']==[]
    cov=coverage(tmp_path,[['S α','Custom reference provider','complete']])
    r=load_project(project(tmp_path,[spec],coverage=cov))
    assert r['coverage'][0]['caller']=='Custom reference provider'
    assert r['coverage'][0]['status']=='complete'


def test_cross_drive_provenance_has_explicit_path_fallback(tmp_path, monkeypatch):
    from defenseviz.io.project import Reader
    reader=Reader(tmp_path/'project.json')
    source=(tmp_path/'data.tsv').resolve()
    def cross_drive(*args, **kwargs):
        raise ValueError('path is on another drive')
    monkeypatch.setattr('defenseviz.io.project.os.path.relpath',cross_drive)
    assert reader.label(source)==source.as_posix()

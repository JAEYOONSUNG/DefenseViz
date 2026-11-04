"""Explicit native artifact lists preserve their memberships and provenance."""
import csv
import hashlib
import json

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.io.suite import SuiteReader, load_suite


def table(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t'); writer.writerow(columns); writer.writerows(rows)
    return path


def fixture(tmp_path, *, second_rows=None, same_replicon=False):
    root=tmp_path/'isolate'; root.mkdir()
    genes=table(root/'genes.tsv',['locus_tag','contig','start','end','direction','translation'],
        [['a','chr',1,12,'+','MAAA'],['b','chr'if same_replicon else'plasmid',1,12,'+','MBBB']])
    mapping=table(root/'map.tsv',['query_id','locus_tag','contig_id'],
        [['q1','a','chr'],['q2','b','chr'if same_replicon else'plasmid']])
    headers=['replicon','hit_id','sys_id','type','subtype','gene_name']
    first=table(root/'one.tsv',headers,[['chr','q1','system1','Wadjet','Wadjet_I','Wadjet__JetD_I']])
    second=table(root/'two.tsv',headers,second_rows if second_rows is not None else
        [['chr'if same_replicon else'plasmid','q2','system1','Wadjet','Wadjet_II','Wadjet__JetD_II']])
    manifest=tmp_path/'manifest.json'
    config=dict(schema_version=1,samples=[dict(sample_id='s',path='isolate',genes='isolate/genes.tsv',
        modules={'defensefinder':{'calls':['isolate/one.tsv','isolate/two.tsv'],'id_map':'isolate/map.tsv'}})])
    manifest.write_text(json.dumps(config))
    return root,manifest,config,first,second


def test_two_native_artifacts_keep_exact_sources_and_original_system_scopes(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path)
    data=load_suite(root,manifest=manifest);result=analyze_saved(data)
    calls=[c for c in data['calls']if c['caller']=='DefenseFinder']
    groups=[g for g in data['groups']if g['provider']=='DefenseFinder']
    assert len(calls)==len(groups)==2
    assert {c['gene_id']for c in calls}=={'a','b'}
    assert {g['replicon_id']for g in groups}=={'chr','plasmid'}
    assert {g['native_system_id']for g in groups}=={'system1'}
    assert len({g['run_id']for g in groups})==2
    assert all(g['replicon_binding']=='verified'for g in groups)
    for path in [first,second]:
        sha=hashlib.sha256(path.read_bytes()).hexdigest()
        assert data['input_sha256'][str(path)]==sha
        call=next(c for c in calls if c['source_ref']['sha256']==sha)
        assert call['source_ref']['row']==2 and call['run_id']==sha
    coverage=[c for c in data['coverage']if c['caller']=='DefenseFinder'and c['endpoint']=='native_calls']
    assert len(coverage)==1 and coverage[0]['status']=='complete'
    detail=next(c for c in data['suite']['coverage']if c['module']=='defensefinder')
    assert detail['artifact_count']==2 and len(detail['selected_sources'])==2
    assert result['global_summary']['native_union_CDS']==2


def test_duplicate_memberships_do_not_inflate_unique_cds(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path,second_rows=
        [['chr','q1','another-system','Wadjet','Wadjet_I','Wadjet__JetD_I']])
    data=load_suite(root,manifest=manifest);result=analyze_saved(data)
    assert len(data['calls'])==2
    assert len(data['groups'])==2
    assert result['global_summary']['native_union_CDS']==1


def test_repeated_identical_artifacts_merge_observations_without_duplicate_groups(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path)
    second.write_bytes(first.read_bytes())
    data=load_suite(root,manifest=manifest)
    assert len(data['calls'])==len(data['groups'])==1
    assert len(data['calls'][0]['source_refs'])==len(data['groups'][0]['source_refs'])==2
    assert str(first)in data['input_sha256'] and str(second)in data['input_sha256']


def test_same_system_id_two_replicons_in_one_artifact_remain_separate(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path)
    with first.open('a')as stream:stream.write('plasmid\tq2\tsystem1\tWadjet\tWadjet_I\tWadjet__JetD_I\n')
    config['samples'][0]['modules']['defensefinder']['calls']=['isolate/one.tsv']
    manifest.write_text(json.dumps(config));data=load_suite(root,manifest=manifest)
    assert len(data['groups'])==2
    assert len({g['group_id']for g in data['groups']})==2
    assert len({g['run_id']for g in data['groups']})==1
    assert {g['replicon_id']for g in data['groups']}=={'chr','plasmid'}


@pytest.mark.parametrize('bad',[[],['isolate/one.tsv','isolate/one.tsv'],['isolate/one.tsv','absent.tsv'],[None]])
def test_bad_artifact_list_fails_without_silently_dropping_a_part(tmp_path,bad):
    root,manifest,config,first,second=fixture(tmp_path)
    config['samples'][0]['modules']['defensefinder']['calls']=bad
    manifest.write_text(json.dumps(config))
    with pytest.raises(ValueError):load_suite(root,manifest=manifest)


def test_shared_id_map_is_exact_for_every_part(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path,second_rows=
        [['plasmid','unknown-id','system1','Wadjet','Wadjet_I','Wadjet__JetD_I']])
    with pytest.raises(ValueError,match='ID map has no entry'):load_suite(root,manifest=manifest)


def test_explicit_replicon_map_conflict_is_rejected(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path,second_rows=
        [['chr','q2','system1','Wadjet','Wadjet_I','Wadjet__JetD_I']])
    with pytest.raises(ValueError,match='replicon contradicts'):load_suite(root,manifest=manifest)


def test_single_path_and_one_element_list_have_same_unique_members(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path)
    config['samples'][0]['modules']['defensefinder']['calls']='isolate/one.tsv'
    manifest.write_text(json.dumps(config));single=load_suite(root,manifest=manifest)
    config['samples'][0]['modules']['defensefinder']['calls']=['isolate/one.tsv']
    manifest.write_text(json.dumps(config));listed=load_suite(root,manifest=manifest)
    assert single['calls']==listed['calls']
    assert single['groups']==listed['groups']
    assert analyze_saved(single)['source_sets']==analyze_saved(listed)['source_sets']


def test_completion_metadata_checks_total_rows_across_parts(tmp_path, monkeypatch):
    root,manifest,config,first,second=fixture(tmp_path)
    monkeypatch.setattr(SuiteReader,'metadata',lambda *a:dict(sentinels={
        'DefenseFinder':dict(sentinel_version=1,alias='DefenseFinder',n_hits=2)}))
    data=load_suite(root,manifest=manifest)
    assert not any(c['code']=='sentinel_live_row_count_mismatch'for c in data['suite']['conflicts'])


def test_redundant_system_file_checks_aggregate_original_memberships(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path)
    table(root/'dnmb_module_defensefinder/defensefinder_systems.tsv',
        ['sys_id','protein_in_syst'],[['system1','q1,q2']])
    data=load_suite(root,manifest=manifest)
    assert not any(c['code']=='genes_systems_membership_mismatch'for c in data['suite']['conflicts'])
    assert len(data['calls'])==2 and len(data['groups'])==2


def test_unresolved_raw_context_is_structured_unknown_not_dropped_zero(tmp_path):
    root,manifest,config,first,second=fixture(tmp_path)
    context=table(root/'dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv',
        ['query','typei_context_status','typei_context_partners'],
        [['a','complete_mrs','missing-member:M:annotation:+1g:20bp']])
    data=load_suite(root,manifest=manifest)
    assert len(data['suite']['contexts'])==1
    unresolved=data['suite']['unresolved_contexts']
    assert len(unresolved)==1
    record=unresolved[0]
    assert record['sample_id']=='s'and record['gene_id']=='a'and record['family']=='Type I'
    assert record['matching_gene_ids']==[]and record['resolved_member_ids']==['a']
    assert record['source_ref']['sha256']==hashlib.sha256(context.read_bytes()).hexdigest()
    assert record['source_ref']['row']==2
    assert not any(g['provider']=='DNMB'for g in data['groups'])

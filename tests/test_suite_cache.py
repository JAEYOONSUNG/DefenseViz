"""Saved DNMB cache evidence is independently bound auxiliary evidence."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from defenseviz.io import suite_cache as module

SEQ='MAKLMNPRSTVWY'


def fixture(tmp_path, *, native_seq=SEQ):
    directory=tmp_path/'sample';m=directory/'dnmb_module_rebasefinder';m.mkdir(parents=True)
    (directory/'.dnmb_module_stage_cache.rds').write_bytes(b'synthetic cache')
    (m/'query_sequences.fasta').write_text('>g\n'+native_seq+'\n')
    (m/'blast_results.txt').write_text(f'g\tref\t100.000\t{len(SEQ)}\t0\t0\t1\t{len(SEQ)}\t1\t{len(SEQ)}\t1e-30\t100\n')
    data=dict(samples=[dict(sample_id='s')],genes=[dict(sample_id='s',gene_id='g',protein_sha256=sha256(SEQ.encode()).hexdigest(),protein_length=len(SEQ))],
              calls=[dict(sample_id='s',gene_id='g',caller='PADLOC')],saved_decisions=[dict(sample_id='s',gene_id='g',state='unassessed')],
              source_metadata=[dict(sample_id='s',source_directory=str(directory))],
              suite=dict(features=[dict(sample_id='s',gene_id='g',native_fields={'hit_label':'ref'})],sources=[]))
    return data,directory,m


def projection(*, seq=SEQ, extra=None):
    values=[dict(kind='metadata',schema_version=1,status='available'),
       dict(kind='table',table='pipeline_result/rm_comprehensive',row_count=2,columns=['locus_tag','translation','rm_type']),
       dict(kind='row',table='pipeline_result/rm_comprehensive',row=1,fields=dict(locus_tag='g',translation=seq,rm_type='Type II',subunit='M',blast_match='ref',blast_qstart=1,blast_qend=len(SEQ))),
       dict(kind='row',table='pipeline_result/rm_comprehensive',row=2,fields=dict(locus_tag='g',translation=seq,rm_type='Type II',subunit='M',blast_match='ref',blast_qstart=1,blast_qend=len(SEQ)))]
    return values+(extra or [])


def mock(monkeypatch, rows):
    monkeypatch.setattr(module.shutil,'which',lambda p:'/bin/echo')
    monkeypatch.setattr(module,'_read_cache',lambda p,r:(rows,{'exit_code':0,'command':['Rscript','reader',str(p)]}))


def test_exact_binding_dedup_provenance_and_native_data_unchanged(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);before=deepcopy(data);mock(monkeypatch,projection())
    out=module.attach_cached_evidence(data,tmp_path)
    assert data==before and out['calls']==data['calls'] and out['saved_decisions']==data['saved_decisions']
    assert out['suite']['features']==data['suite']['features']
    cache=out['suite']['cached_features'];assert len(cache)==1
    assert [x['row'] for x in cache[0]['source_refs']]==[1,2]
    assert cache[0]['protein_sha256']==sha256(SEQ.encode()).hexdigest()
    assert len(out['suite']['cached_blast'])==1
    assert out['suite']['native_cache_evidence']['retention_decisions']==0
    assert all(x['search_is_exhaustive'] is False for x in out['suite']['cached_blast'])
    assert len(out['input_sha256'])==3


def test_conflicting_cached_translation_never_borrows_other_binding(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);mock(monkeypatch,projection(seq=SEQ+'A'))
    out=module.attach_cached_evidence(data,tmp_path)
    assert not out['suite']['cached_features']
    assert out['suite']['cached_blast']  # exact original query FASTA still valid
    assert any(x['reason']=='cached_row_protein_mismatch' for x in out['suite']['native_cache_evidence']['samples'][0]['omitted'])


def test_stale_query_fasta_cannot_borrow_cache_protein_for_raw_blast(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path,native_seq=SEQ+'A');mock(monkeypatch,projection())
    out=module.attach_cached_evidence(data,tmp_path)
    assert out['suite']['cached_features']
    assert not out['suite']['cached_blast']


def test_same_saved_hsp_across_cache_and_file_is_single_measurement(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path)
    hsp=dict(query_id='g',subject_id='ref',pident=100,length=len(SEQ),mismatch=0,gapopen=0,qstart=1,qend=len(SEQ),sstart=1,send=len(SEQ),evalue=1e-30,bitscore=100)
    mock(monkeypatch,projection(extra=[dict(kind='row',table='pipeline_result/blast_raw',row=1,fields=hsp)]))
    out=module.attach_cached_evidence(data,tmp_path)
    hsps=out['suite']['cached_blast'];assert len(hsps)==1
    assert len(hsps[0]['source_refs'])==2
    assert hsps[0]['source_channels']==['blast_results.txt','pipeline_result/blast_raw']


def test_auxiliary_scope_difference_is_visible_not_new_active_candidates(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);mock(monkeypatch,projection(extra=[dict(kind='row',table='hits',row=1,fields={'query':'other','family_id':'Type I'})]))
    out=module.attach_cached_evidence(data,tmp_path)
    summary=out['suite']['native_cache_evidence']['samples'][0]
    assert summary['cached_hit_scope_differs']
    assert len(out['genes'])==1 and len(out['saved_decisions'])==1
    assert any(r['gene_id']=='other' for r in summary['omitted'])


def test_hash_changed_since_import_is_rejected(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);mock(monkeypatch,projection())
    data['input_sha256']={str(d/'.dnmb_module_stage_cache.rds'):'a'*64}
    with pytest.raises(ValueError,match='source changed'):module.attach_cached_evidence(data,tmp_path)


def test_invalid_alignment_does_not_enter_evidence(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);mock(monkeypatch,projection())
    (m/'blast_results.txt').write_text('g\tref\t101\t13\t0\t0\t1\t13\t1\t13\t0\t10\n')
    out=module.attach_cached_evidence(data,tmp_path)
    assert not out['suite']['cached_blast']
    assert any(r['reason']=='invalid_saved_alignment' for r in out['suite']['native_cache_evidence']['samples'][0]['omitted'])


def test_terminal_stop_is_part_of_exact_protein_binding(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path,native_seq=SEQ+'*');mock(monkeypatch,projection(seq=SEQ+'*'))
    out=module.attach_cached_evidence(data,tmp_path)
    assert not out['suite']['cached_features'] and not out['suite']['cached_blast']


def test_no_rscript_does_not_erase_readable_blast_files(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);monkeypatch.setattr(module.shutil,'which',lambda p:None)
    out=module.attach_cached_evidence(data,tmp_path)
    assert out['suite']['native_cache_evidence']['samples'][0]['status']=='reader_unavailable'
    assert out['suite']['cached_blast']
    with pytest.raises(ValueError,match='Configured Rscript'):module.attach_cached_evidence(data,tmp_path,rscript='missing-r')


def test_repeat_attachment_is_explicit_error(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);mock(monkeypatch,projection())
    out=module.attach_cached_evidence(data,tmp_path)
    with pytest.raises(ValueError,match='already attached'):module.attach_cached_evidence(out,tmp_path)


@pytest.mark.skipif(shutil.which('Rscript') is None,reason='Rscript unavailable')
def test_dependency_free_r_reader_real_saved_object(tmp_path):
    path=tmp_path/'cache.rds'
    expression="x<-list(module_results=list(REBASEfinder=list(module_result=list(ok=TRUE,hits=data.frame(query='g',family_id='Type III'),pipeline_result=list(rm_comprehensive=data.frame(locus_tag='g',translation='MAK',blast_qstart=1L,blast_qend=3L))))));saveRDS(x,commandArgs(TRUE)[1])"
    subprocess.run([shutil.which('Rscript'),'--vanilla','-e',expression,str(path)],check=True,capture_output=True)
    rows,command=module._read_cache(path,shutil.which('Rscript'))
    assert rows[0]['schema_version']==1 and command['exit_code']==0
    r=next(r for r in rows if r.get('table')=='pipeline_result/rm_comprehensive' and r['kind']=='row')
    assert r['row']==1 and r['fields']['translation']=='MAK'


@pytest.mark.parametrize('index,value',[(2,'nan'),(10,'inf'),(11,'-1'),(4,'-1'),(1,'')])
def test_bad_numeric_or_subject_hsp_is_quarantined(tmp_path,monkeypatch,index,value):
    data,d,m=fixture(tmp_path);mock(monkeypatch,projection())
    columns=['g','ref','100','13','0','0','1','13','1','13','0','10'];columns[index]=value
    (m/'blast_results.txt').write_text('\t'.join(columns)+'\n')
    out=module.attach_cached_evidence(data,tmp_path)
    assert not out['suite']['cached_blast']


def test_supplemental_hit_uses_its_own_query_fasta(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);mock(monkeypatch,projection())
    sub=m/'supplemental_rebase_blast';sub.mkdir()
    (sub/'supplemental_queries.faa').write_text('>g\n'+SEQ+'A\n')
    (sub/'supplemental_blast_results.txt').write_text((m/'blast_results.txt').read_text().replace('ref','supplemental'))
    out=module.attach_cached_evidence(data,tmp_path)
    assert {r['native_fields']['subject_id'] for r in out['suite']['cached_blast']}=={'ref'}


def test_cache_genome_mismatch_remains_explicit_even_with_exact_protein(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);(d/'test.gbff').write_text('changed genome')
    rows=projection();rows[0]['genbank_files']=[dict(file='test.gbff',md5='0'*32)]
    mock(monkeypatch,rows)
    out=module.attach_cached_evidence(data,tmp_path)
    assert out['suite']['native_cache_evidence']['samples'][0]['genome_binding']=='mismatched'
    assert all(r['genome_binding']=='mismatched' for r in out['suite']['cached_features'])
    assert out['saved_decisions']==data['saved_decisions']


def test_malformed_raw_blast_is_not_silently_dropped(tmp_path,monkeypatch):
    data,d,m=fixture(tmp_path);mock(monkeypatch,projection())
    (m/'blast_results.txt').write_text('g\tref\t100\n')
    with pytest.raises(ValueError,match='12-column'):module.attach_cached_evidence(data,tmp_path)

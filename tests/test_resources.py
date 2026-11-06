import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import pytest

from defenseviz import resources


def profiles(tmp_path):
    hmm=tmp_path/'library.hmm';hmm.write_text('HMMER3/f\nNAME  Restriction\nLENG  100\nGA  20 20;\n//\n')
    meta=tmp_path/'roles.tsv';meta.write_text('profile_name\toff_target\troles\ttype_hints\tprofile_length\tlibrary_sha256\nRestriction\tfalse\tR\tII\t100\t'+hashlib.sha256(hmm.read_bytes()).hexdigest()+'\n')
    return hmm,meta


@pytest.mark.parametrize('fault',['none','wrong_hash','missing_GA','wrong_name','duplicate'])
def test_exact_HMM_role_binding(tmp_path,fault):
    hmm,meta=profiles(tmp_path)
    if fault=='wrong_hash':meta.write_text(meta.read_text().replace(hashlib.sha256(hmm.read_bytes()).hexdigest(),'0'*64))
    if fault=='missing_GA':hmm.write_text(hmm.read_text().replace('GA  20 20;\n',''))
    if fault=='wrong_name':meta.write_text(meta.read_text().replace('Restriction','Another'))
    if fault=='duplicate':hmm.write_text(hmm.read_text()*2)
    if fault=='none':
        result=resources._profile_binding(hmm,meta)
        assert result['profiles']==result['reportable_profiles']==1
        assert result['exact_profile_name_binding']
    else:
        with pytest.raises(ValueError):resources._profile_binding(hmm,meta)


@pytest.mark.parametrize('fault',[
    'missing_hash_column','missing_length_column','empty_hash','empty_length',
    'wrong_length','fractional_length','missing_LENG','duplicate_LENG',
    'zero_LENG','duplicate_column','short_row','long_row',
])
def test_preflight_rejects_metadata_that_full_evaluation_cannot_bind(tmp_path,fault):
    import csv
    hmm,meta=profiles(tmp_path)
    with meta.open() as stream: rows=list(csv.reader(stream,delimiter='\t'))
    header,row=rows
    if fault.startswith('missing_') and fault.endswith('_column'):
        name='library_sha256' if fault=='missing_hash_column' else 'profile_length'
        index=header.index(name);header.pop(index);row.pop(index)
    elif fault=='empty_hash':row[header.index('library_sha256')]=''
    elif fault in {'empty_length','wrong_length','fractional_length'}:
        row[header.index('profile_length')]={'empty_length':'','wrong_length':'99','fractional_length':'100.5'}[fault]
    elif fault in {'missing_LENG','duplicate_LENG','zero_LENG'}:
        replacement={'missing_LENG':'','duplicate_LENG':'LENG  100\nLENG  100\n','zero_LENG':'LENG  0\n'}[fault]
        hmm.write_text(hmm.read_text().replace('LENG  100\n',replacement))
        row[header.index('library_sha256')]=hashlib.sha256(hmm.read_bytes()).hexdigest()
    elif fault=='duplicate_column':header.append('roles');row.append('R')
    elif fault=='short_row':row.pop()
    elif fault=='long_row':row.append('unexpected')
    with meta.open('w',newline='') as stream:csv.writer(stream,delimiter='\t').writerows([header,row])
    with pytest.raises(ValueError):resources._profile_binding(hmm,meta)


def environment(tmp_path,monkeypatch):
    hmm,meta=profiles(tmp_path);cache=tmp_path/'db';cache.mkdir()
    for name in ('rebase_data.rds','rebase_db.fasta','rebase_db.fasta.phr','rebase_db.fasta.pin','rebase_db.fasta.psq'):
        (cache/name).write_text('reference '+name)
    library=tmp_path/'r-library';library.mkdir()
    for name in ('DNMB','DefenseViz'):(library/name).mkdir()
    cfg=dict(schema_version=1,rscript=sys.executable,r_libraries=[str(library)],
        reference_cache=str(cache),hmm_database=str(hmm),profile_metadata=str(meta),
        hmmsearch=sys.executable,blastp=sys.executable,blastdbcmd=sys.executable,cath_resolve_hits=sys.executable)
    config=tmp_path/'engine.json';config.write_text(json.dumps(cfg))
    seen=[]
    def runtime(cfg,work,env):
        seen.append(Path(work))
        value=dict(R_version='R test',R_home='/R',packages=[
            dict(name='DNMB',version='1',path=str(library/'DNMB')),
            dict(name='base',version='1',path='/R/library/base'),
            dict(name='other',version='1',path='/user/R/other')])
        p=work/'r_runtime.json';p.write_text(json.dumps(value))
        return value,{str(p):resources.native.digest(p)}
    def database(prefix,executable,out,env):
        result={'entries':12,'passed':True};out.write_text(json.dumps(result));return result
    monkeypatch.setattr(resources.native,'_capture_r_runtime',runtime)
    monkeypatch.setattr(resources.native,'_validate_reference_index',database)
    monkeypatch.setattr(resources.subprocess,'run',lambda *a,**kw:SimpleNamespace(returncode=0,stdout='tool version1',stderr=''))
    return config,hmm,cache,seen


def test_resource_preflight_reports_observed_versions_origins_and_no_analysis(tmp_path,monkeypatch):
    config,hmm,cache,_=environment(tmp_path,monkeypatch)
    before=resources.native.digest(hmm);out=tmp_path/'readiness'
    result=resources.verify_resources(config,out)
    assert result['ready'] and not result['full_analysis_executed']
    assert result['resource_downloads']==0 and result['reference']['protein_entries']==12
    assert result['hmm']['profiles']==1 and result['R']['version']=='R test'
    assert {p['origin'] for p in result['R']['packages']}=={'configured_library','R_installation','other_installed_library'}
    assert (out/'resource_input_hashes.json').is_file()
    assert result['tools']['blastp']['exit_code']==0
    assert resources.native.digest(hmm)==before
    with pytest.raises(FileExistsError):resources.verify_resources(config,out)


def test_missing_reference_has_actionable_failure_receipt(tmp_path,monkeypatch):
    config,_,cache,_=environment(tmp_path,monkeypatch)
    (cache/'rebase_db.fasta.psq').unlink();out=tmp_path/'failed'
    with pytest.raises(RuntimeError,match='Incomplete REBASE cache.*rebase_db.fasta.psq'):
        resources.verify_resources(config,out)
    assert json.loads((out/'resource_check.json').read_text())['ready'] is False


def test_mutated_resources_never_report_ready(tmp_path,monkeypatch):
    config,hmm,_,_=environment(tmp_path,monkeypatch)
    def altered(*a):hmm.write_text('changed during validation');return {'entries':12}
    monkeypatch.setattr(resources.native,'_validate_reference_index',altered)
    with pytest.raises(RuntimeError,match='changed during validation'):
        resources.verify_resources(config,tmp_path/'failed')


def test_temporary_resource_check_does_not_claim_preserved_logs(tmp_path,monkeypatch):
    config,_,_,seen=environment(tmp_path,monkeypatch)
    result=resources.verify_resources(config)
    assert not result['report_retained'] and 'report_directory' not in result
    assert not seen[0].exists()

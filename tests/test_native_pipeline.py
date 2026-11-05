from pathlib import Path
import hashlib,json,sys
import pytest
from defenseviz.native_pipeline import read_engine_config,engine_config_path,_prepare_sample,_sample_path


def config(tmp_path):
    resources={}
    for key in ('hmm_database','profile_metadata','cath_resolve_hits'):
        p=tmp_path/key;p.write_text('resource');resources[key]=str(p)
    resources['reference_cache']=str(tmp_path)
    resources.update(schema_version=1,rscript=sys.executable,blastp=sys.executable,hmmsearch=sys.executable,cpu=2)
    path=tmp_path/'engine.json';path.write_text(json.dumps(resources));return path,resources


def test_config_paths_resolve_relative_without_shell(tmp_path):
    path,c=config(tmp_path);c['hmm_database']='hmm_database';path.write_text(json.dumps(c))
    result=read_engine_config(path)
    assert result['hmm_database']==str(tmp_path/'hmm_database')
    assert result['_config_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()


def test_relative_executable_is_bound_to_config_directory(tmp_path,monkeypatch):
    path,c=config(tmp_path);tool=tmp_path/'bin'/'blastp';tool.parent.mkdir();tool.write_text('tool');tool.chmod(0o755)
    elsewhere=tmp_path/'elsewhere';(elsewhere/'bin').mkdir(parents=True)
    (elsewhere/'bin'/'blastp').write_text('wrong tool');(elsewhere/'bin'/'blastp').chmod(0o755)
    c['blastp']='bin/blastp';path.write_text(json.dumps(c));monkeypatch.chdir(elsewhere)
    assert read_engine_config(path)['blastp']==str(tool)


@pytest.mark.parametrize('updates',[{'cpu':0},{'cpu':True},{'schema_version':2},{'extra_command':'rm -rf /'},{'hmm_database':'not-there'},{'r_libraries':'bad'}])
def test_invalid_config_fails_closed(tmp_path,updates):
    path,c=config(tmp_path);c.update(updates);path.write_text(json.dumps(c))
    with pytest.raises((ValueError,FileNotFoundError)):read_engine_config(path)


def test_explicit_config_precedes_environment(tmp_path,monkeypatch):
    monkeypatch.setenv('DEFENSEVIZ_ENGINE_CONFIG',str(tmp_path/'env.json'))
    assert engine_config_path(tmp_path/'explicit.json')==tmp_path/'explicit.json'
    assert engine_config_path()==tmp_path/'env.json'


def bundle(sequence='MAC*'):
    return {'genes':[{'sample_id':'a','gene_id':'x','protein_sha256':hashlib.sha256(sequence.encode()).hexdigest()}],
            'suite':{'annotations':[{'sample_id':'a','gene_id':'x','native_fields':{'locus_tag':'x','translation':sequence,'product':'protein','REBASEfinder_old':'stale'}}]}}


def test_preparation_preserves_whole_protein_binding_and_audits_stop(tmp_path):
    directory=tmp_path/'input';directory.mkdir();gb=directory/'genome.gbff';gb.write_text('native genome')
    out=tmp_path/'work';out.mkdir()
    assert _prepare_sample(bundle(),'a',directory,out)==gb
    assert (out/'proteins.faa').read_text()=='>x\nMAC\n'
    assert 'REBASEfinder_old' not in (out/'annotations.tsv').read_text()
    assert json.loads((out/'protein_binding.json').read_text())['source_protein_sha256']['x']==hashlib.sha256(b'MAC*').hexdigest()


def test_preparation_rejects_sequence_mutation(tmp_path):
    data=bundle();data['suite']['annotations'][0]['native_fields']['translation']='MAW'
    with pytest.raises(ValueError,match='binding mismatch'):_prepare_sample(data,'a',tmp_path,tmp_path)


def test_preparation_rejects_omitted_protein(tmp_path):
    data=bundle();data['genes'].append({'sample_id':'a','gene_id':'missing','protein_sha256':'0'*64})
    with pytest.raises(ValueError,match='omits imported proteins'):_prepare_sample(data,'a',tmp_path,tmp_path)


def test_sample_path_uses_source_identity_not_display(tmp_path):
    sample=tmp_path/'misleading label';sample.mkdir();(sample/'genome.gbff').write_text('native')
    module=sample/'dnmb_module_rebasefinder';module.mkdir()
    data={'suite':{'sources':[{'sample_id':'a','resolved_path':str(module/'hits.tsv')}]}}
    assert _sample_path(data,'a',tmp_path)==sample


def test_native_resources_include_fixed_metadata_projection():
    from importlib.resources import files
    source=files('defenseviz').joinpath('data/native/full_rm.R').read_text()
    assert "runtime$r_libraries" in source and '46119' not in source
    assert 'Type III' not in source or "explicit_type" in source
    assert "original_score(strong,moderate,profile" in source
    assert '/Users/' not in source


def test_nomination_union_preserves_qualified_saved_evidence_without_mixed_scores(tmp_path):
    import csv
    from defenseviz.native_pipeline import _nomination_union
    fresh=tmp_path/'fresh.tsv';fresh.write_text('locus_tag\tprofile_name\nnew\tnative_model\n')
    def profile(gid,model,evalue,qcov,mcov):
        return dict(sample_id='a',gene_id=gid,namespace='PADLOC',native_model=model,
                    native_fields={'full.seq.E.value':evalue,'target.coverage':qcov,'hmm.coverage':mcov})
    data={'genes':[dict(sample_id='a',gene_id=g,protein_length=100,protein_sha256='a'*64) for g in ('keep','weak','mixed')],
          'suite':{'profiles':[profile('keep','REase_II_00013',1e-70,.9,.9),
                               profile('weak','REase_II_00013',.1,.9,.9),
                               profile('mixed','REase_II_00013',1e-70,.3,.9),
                               profile('mixed','NEase_II_00007',1e-70,.9,.3),
                               dict(profile('keep','REase_II_00099',1e-90,1,1),sample_id='other')]}}
    output=tmp_path/'union.tsv';audit=_nomination_union(data,'a',fresh,output)
    rows=list(csv.DictReader(output.open(),delimiter='\t'))
    assert {r['locus_tag'] for r in rows}=={'new','keep'}
    assert audit['union_nominees']==2 and audit['saved_qualified_nominees']==1
    assert rows[1]['profile_name']=='REase_II_00013'


def test_full_orchestrator_executes_fresh_stages_and_pins_private_outputs(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from defenseviz import native_pipeline as module
    from defenseviz import full_hmm
    path,cfg=config(tmp_path)
    for name in ('rebase_data.rds','rebase_db.fasta','rebase_db.fasta.phr','rebase_db.fasta.pin','rebase_db.fasta.psq'):
        (tmp_path/name).write_text('fixed resource')
    sample=tmp_path/'sample';sample.mkdir();(sample/'genome.gbff').write_text('native GBFF')
    data=bundle('MAC');data['samples']=[{'sample_id':'a'}]
    data['suite']['samples']=[{'sample_id':'a','directory':str(sample)}]
    data['genes'][0]['protein_length']=3
    seen=[]
    runtime_record=dict(R_version='test R', R_home=str(tmp_path), loaded_DLLs={},
        packages=[dict(name=n,path=str(tmp_path/n),version='1') for n in ('DNMB','DefenseViz','jsonlite','digest')])
    def hmm(args):
        out=Path(args.output);out.mkdir()
        (out/'nominees.tsv').write_text('locus_tag\tprofile_name\n')
        (out/'execution.json').write_text(json.dumps({'exit_code':0,'sample_id':args.sample}))
        seen.append(('hmm',args.sample,Path(args.fasta).read_text()))
    def runtime(cfg,work,env):
        p=work/'r_runtime.json';p.write_text(json.dumps(runtime_record))
        return runtime_record,{str(p):module.digest(p)}
    def validate(prefix,executable,output,env):
        result={'passed':True,'entries':1};output.write_text(json.dumps(result));return result
    def process(argv,**kwargs):
        if 'full_rm.R' in ' '.join(argv):
            out=Path(argv[-2]);seen.append(('native',argv))
            (out/'completion.json').write_text(json.dumps({'passed':True,'full_reference_search':True}))
            (out/'search_commands.json').write_text('[]')
            (out/'runtime_packages.json').write_text(json.dumps(runtime_record))
        return SimpleNamespace(returncode=0,stdout='version',stderr='')
    monkeypatch.setattr(full_hmm,'run',hmm);monkeypatch.setattr(module.subprocess,'run',process)
    monkeypatch.setattr(module,'_capture_r_runtime',runtime)
    monkeypatch.setattr(module,'_resolve_blastdbcmd',lambda cfg:sys.executable)
    monkeypatch.setattr(module,'_validate_reference_index',validate)
    work=tmp_path/'work'
    stages=module.execute_full_pipeline(data,tmp_path,path,work)
    assert seen[0]==('hmm','a','>x\nMAC\n') and seen[1][0]=='native'
    manifest=json.loads((Path(stages[0]['native_dir'])/'stage_manifest.json').read_text())
    assert any(p.endswith('protein_binding.json') for p in manifest['input_sha256'])
    assert str(Path(stages[0]['hmm_dir'])/'execution.json') in manifest['input_sha256']
    assert manifest['native_execution']['sha256']==module.digest(manifest['native_execution']['path'])
    assert 'runtime_validation.json' in manifest['artifacts']
    assert manifest['blast_database']['prefix']==str(work/'cache/db_modules/rebasefinder/cache/rebase_db.fasta')
    assert all(module.digest(p)==v for p,v in manifest['input_sha256'].items())
    assert json.loads((work/'execution.json').read_text())['homology_search_runs']==0
    assert (work/'cache/db_modules/rebasefinder/cache/rebase_data.rds').read_bytes()==(tmp_path/'rebase_data.rds').read_bytes()


def test_failed_full_publication_preserves_logs_without_success_output(tmp_path,monkeypatch):
    from defenseviz import suite_pipeline
    def fail(root,**kwargs):
        work=Path(kwargs['evidence_workspace'])/'full_pipeline';work.mkdir(parents=True)
        (work/'stderr.log').write_text('search failed')
        raise RuntimeError('BLAST exit 2')
    monkeypatch.setattr(suite_pipeline,'prepare_suite',fail)
    with pytest.raises(RuntimeError,match='evidence preserved'):
        suite_pipeline.run_suite(tmp_path,tmp_path/'result',evaluation='full',engine_config='engine.json')
    assert not (tmp_path/'result').exists()
    failed=list(tmp_path.glob('.result.failed-*'));assert len(failed)==1
    assert not json.loads((failed[0]/'failure.json').read_text())['published']
    assert (failed[0]/'evidence/reference/full_pipeline/stderr.log').read_text()=='search failed'


def test_runtime_bootstraps_configured_tools_and_R_libraries(tmp_path,monkeypatch):
    import os,shutil
    from defenseviz.native_pipeline import _runtime_environment
    binary=tmp_path/'tool dir'/('blastp.exe' if os.name=='nt' else 'blastp');binary.parent.mkdir();binary.write_text('executable');binary.chmod(0o755)
    monkeypatch.setenv('PATH','/usr/bin');monkeypatch.setenv('R_LIBS','/existing/R')
    env=_runtime_environment({'blastp':str(binary),'r_libraries':['/private/R']},tmp_path/'cache')
    resolved_binary=shutil.which('blastp',path=env['PATH'])
    assert resolved_binary is not None and Path(resolved_binary).samefile(binary)
    assert env['R_LIBS'].split(os.pathsep)==['/private/R','/existing/R']
    assert env['DNMB_BLASTP']==str(binary)


def test_frozen_source_tamper_fails_before_native_execution(tmp_path):
    from defenseviz.native_pipeline import _verify_bundled_sources
    (tmp_path/'source').mkdir();source=tmp_path/'source'/'policy.R';source.write_text('frozen policy')
    (tmp_path/'provenance.json').write_text(json.dumps({'source_files_sha256':{'policy.R':hashlib.sha256(source.read_bytes()).hexdigest()}}))
    assert str(source) in _verify_bundled_sources(tmp_path)
    source.write_text('changed cutoff')
    with pytest.raises(ValueError,match='frozen provenance'):_verify_bundled_sources(tmp_path)


def test_R_probe_pins_actual_fallback_namespace_assets(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from defenseviz import native_pipeline as module
    packages=[]
    for name in ('DNMB','DefenseViz','jsonlite','digest'):
        p=tmp_path/'default library'/name;(p/'R').mkdir(parents=True)
        (p/'DESCRIPTION').write_text('Version: 1\n');(p/'R'/(name+'.rdb')).write_text('actual namespace')
        packages.append(dict(name=name,path=str(p),version='1'))
    def process(argv,**kwargs):
        Path(argv[3]).write_text(json.dumps(dict(R_version='test R',R_home=str(tmp_path),
            packages=packages,loaded_DLLs={},blastp=sys.executable)))
        return SimpleNamespace(returncode=0,stdout='',stderr='')
    monkeypatch.setattr(module.subprocess,'run',process)
    runtime,hashes=module._capture_r_runtime(dict(rscript=sys.executable,r_libraries=[],blastp=sys.executable),tmp_path,{})
    assert len(runtime['packages'])==4
    assert all(str(Path(p['path'])/'R'/(p['name']+'.rdb')) in hashes for p in packages)


@pytest.mark.parametrize('mutation', ['none','version','library','additional_namespace','DLL','R'])
def test_actual_native_runtime_must_match_preflight(tmp_path,mutation):
    from copy import deepcopy
    from defenseviz.native_pipeline import _verify_native_runtime, digest
    dll=tmp_path/'libdigest.so';dll.write_text('pinned shared library')
    expected=dict(R_version='test R',R_home=str(tmp_path),
        packages=[dict(name=n,path=str(tmp_path/n),version='1') for n in ('DNMB','DefenseViz','jsonlite','digest')],
        loaded_DLLs={'base':'base','digest':str(dll)})
    actual=deepcopy(expected)
    if mutation=='version': actual['packages'][0]['version']='2'
    if mutation=='library': actual['packages'][0]['path']=str(tmp_path/'other-library'/'DNMB')
    if mutation=='additional_namespace': actual['packages'].append(dict(name='surprise',path='/different',version='1'))
    if mutation=='DLL': actual['loaded_DLLs']['digest']=str(tmp_path/'unmeasured.so')
    if mutation=='R': actual['R_version']='different R'
    path=tmp_path/'runtime_packages.json';path.write_text(json.dumps(actual))
    if mutation=='none':
        assert _verify_native_runtime(path,expected,{str(dll):digest(dll)})['passed']
    else:
        with pytest.raises(ValueError,match='different|unpinned'):
            _verify_native_runtime(path,expected,{str(dll):digest(dll)})


@pytest.mark.parametrize('indexed,error',[
    ('ref1\tMAAA\nref2\tMQ*\n',None),
    ('ref1\tMAAT\nref2\tMQ*\n','disagree'),
    ('ref1\tMAAA\n','absent'),
    ('ref1\tMAAA\nref2\tMQ*\nref3\tMQQ\n','absent'),
])
def test_reference_indexes_must_match_every_fasta_entry(tmp_path,indexed,error,monkeypatch):
    import os
    from defenseviz import native_pipeline as module
    from defenseviz.native_pipeline import _validate_reference_index
    fasta=tmp_path/'ref.faa';fasta.write_text('>ref1\nMAAA<>\n>ref2\nMQ*\n')
    executable=tmp_path/'blastdbcmd'
    executable.write_text('synthetic index-stream executable')
    popen=module.subprocess.Popen
    def stream(argv,**kwargs):
        assert argv==[str(executable),'-db',str(fasta),'-entry','all','-outfmt','%t\t%s']
        return popen([sys.executable,'-c','import sys;sys.stdout.write('+repr(indexed)+')'],**kwargs)
    monkeypatch.setattr(module.subprocess,'Popen',stream)
    output=tmp_path/'validation.json'
    if error:
        with pytest.raises(ValueError,match=error):_validate_reference_index(fasta,executable,output,dict(os.environ))
        assert json.loads(output.read_text())['passed'] is False
    else:
        result=_validate_reference_index(fasta,executable,output,dict(os.environ))
        assert result['entries']==2 and result['terminal_markers_removed']==1 and result['passed']

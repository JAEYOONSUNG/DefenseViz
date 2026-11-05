"""Portable local resources preserve identity and reject hidden old-library use."""
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import pytest

spec=importlib.util.spec_from_file_location('runtime_export',Path(__file__).resolve().parents[1]/'scripts/export_runtime_bundle.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def fixture(tmp_path,monkeypatch):
    original=tmp_path/'original';original.mkdir()
    ref=original/'reference';ref.mkdir()
    for name in ('rebase_data.rds','rebase_db.fasta','rebase_db.fasta.phr','rebase_db.fasta.pin','rebase_db.fasta.psq'):
        (ref/name).write_text('reference '+name)
    for name in ('hmm','roles','cath-resolve-hits.bin','cath-wrapper','blastp','blastdbcmd','hmmsearch'):
        (original/name).write_text(name);(original/name).chmod(0o755)
    icu=original/'libicu.68.2.dylib';icu.write_text('exact ICU')
    (original/'libicu.68.dylib').symlink_to(icu.name)
    packages=[]
    for name in ('DNMB','DefenseViz','digest'):
        p=original/'installed'/name;p.mkdir(parents=True)
        (p/'DESCRIPTION').write_text('Package: '+name+'\nVersion: 1\n')
        (p/'R').mkdir();(p/'R'/'package.rdb').write_text('runtime '+name)
        packages.append(dict(name=name,path=str(p),version='1',license='MIT',system=False))
    sources=[]
    for name in ('DNMB','DefenseViz'):
        p=original/'source'/name;p.mkdir(parents=True)
        (p/'DESCRIPTION').write_text('Package: '+name+'\nVersion: 1\n')
        (p/'R').mkdir();(p/'R'/'code.R').write_text('x <- 1')
        sources.append(p)
    cfg=dict(schema_version=1,rscript=sys.executable,r_libraries=[str(original/'installed')],
        reference_cache=str(ref),hmm_database=str(original/'hmm'),profile_metadata=str(original/'roles'),
        cath_resolve_hits=str(original/'cath-wrapper'),blastp=str(original/'blastp'),
        blastdbcmd=str(original/'blastdbcmd'),hmmsearch=str(original/'hmmsearch'),
        runtime_dependencies=[str(original/'cath-resolve-hits.bin'),str(icu),str(original/'libicu.68.dylib')])
    config=original/'engine.json';config.write_text(json.dumps(cfg))
    monkeypatch.setattr(module,'dependency_closure',lambda *a:({'R_version':'R test'},dict(packages=packages),{}))
    return config,sources,original


def test_export_preserves_alias_names_and_relative_resources(tmp_path,monkeypatch):
    cfg,sources,original=fixture(tmp_path,monkeypatch)
    before={str(p):module.sha(p) for p in original.rglob('*') if p.is_file()}
    out=tmp_path/'bundle';result=module.export_bundle(cfg,out,*sources)
    generated=json.loads((out/'engine.json').read_text())
    assert generated['r_libraries']==['r-library']
    assert generated['reference_cache']=='databases/rebase'
    assert (out/'lib/libicu.68.dylib').read_bytes()==(original/'libicu.68.2.dylib').read_bytes()
    assert 'lib/libicu.68.dylib' in generated['runtime_dependencies']
    assert str(original) not in (out/'bin/cath-resolve-hits').read_text()
    assert (out/'sources/DNMB/R/code.R').read_text()=='x <- 1'
    assert (out/'setup/rebuild_r.py').is_file()
    assert result['artifacts']==module.inventory(out)
    assert before=={str(p):module.sha(p) for p in original.rglob('*') if p.is_file()}
    with pytest.raises(FileExistsError):module.export_bundle(cfg,out,*sources)


def test_export_rejects_wrong_source_package_before_writing(tmp_path,monkeypatch):
    cfg,sources,_=fixture(tmp_path,monkeypatch)
    (sources[0]/'DESCRIPTION').write_text('Package: SomethingElse\n')
    with pytest.raises(ValueError,match='Wrong source'):
        module.export_bundle(cfg,tmp_path/'bundle',*sources)
    assert not (tmp_path/'bundle').exists()


def test_tree_export_rejects_directory_symlink_escape(tmp_path):
    source=tmp_path/'source';source.mkdir();external=tmp_path/'external';external.mkdir()
    (source/'data').symlink_to(external,target_is_directory=True)
    with pytest.raises(ValueError,match='Directory symlink'):
        module.copy_tree(source,tmp_path/'copy',{})


@pytest.mark.parametrize('fault',['none','tampered','missing','untracked','old_library'])
def test_relocated_bundle_verifies_identity_and_library_boundary(tmp_path,monkeypatch,fault):
    from defenseviz import native_pipeline
    cfg,sources,original=fixture(tmp_path,monkeypatch)
    source=tmp_path/'before';module.export_bundle(cfg,source,*sources)
    moved=tmp_path/'after';source.rename(moved)
    if fault=='tampered':(moved/'databases/rm-hmm/RM_HMMs.hmm').write_text('changed model')
    if fault=='missing':(moved/'databases/rm-hmm/RM_HMMs.hmm').unlink()
    if fault=='untracked':(moved/'unrecorded').write_text('extra')
    def probe(config,work,env):
        p=(original/'installed/DNMB') if fault=='old_library' else moved/'r-library/DNMB'
        return dict(R_version='R test',R_home='/system/R',packages=[dict(name='DNMB',path=str(p))]),{}
    monkeypatch.setattr(native_pipeline,'_capture_r_runtime',probe)
    monkeypatch.setattr(native_pipeline,'_validate_reference_index',lambda *a:{'entries':3})
    monkeypatch.setattr(module.subprocess,'run',lambda *a,**kw:SimpleNamespace(returncode=0,stdout='CATH test',stderr=''))
    if fault=='none':
        result=module.verify_bundle(moved,tmp_path/'report')
        assert result['passed'] and not result['full_analysis_executed']
    else:
        with pytest.raises(ValueError,match='changed|escaped'):
            module.verify_bundle(moved,tmp_path/'report')


def test_verify_never_writes_into_preserved_bundle(tmp_path):
    with pytest.raises(ValueError,match='outside'):
        module.verify_bundle(tmp_path,tmp_path/'report')

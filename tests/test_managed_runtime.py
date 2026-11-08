import argparse
import json
import os
from pathlib import Path
import subprocess

import pytest

from defenseviz import managed_runtime as m

IMAGE = 'sha256:' + '1' * 64


@pytest.fixture
def managed(tmp_path, monkeypatch):
    if os.name != 'posix':
        pytest.skip('Managed Docker execution uses POSIX paths; Windows users run it inside WSL')
    cache = tmp_path / 'cache with spaces'
    cache.mkdir()
    monkeypatch.setenv('DEFENSEVIZ_CACHE', str(cache))
    monkeypatch.delenv('DEFENSEVIZ_ENGINE_CONFIG', raising=False)
    monkeypatch.delenv('DEFENSEVIZ_MANAGED_CONTAINER', raising=False)
    monkeypatch.setattr('defenseviz.native_pipeline.engine_config_path', lambda explicit=None: explicit)
    monkeypatch.setattr(m, '_docker', lambda: '/tools/docker')
    databases = cache / 'setups/first/databases'
    databases.mkdir(parents=True)
    engine = dict(reference_cache='/resources/rebase', hmm_database='/resources/rm-hmm/RM_HMMs.hmm',
                  profile_metadata='/resources/rm-hmm/profile_roles.tsv', cpu=4)
    m.write_json(databases / 'runtime-engine.json', engine)
    (databases / 'reference.fasta').write_text('>a\nMAA\n')
    profile = dict(schema_version=1, ready=True, kind='docker', image_id=IMAGE, cpu=4,
                   resources=str(databases), resource_sha256=m._inventory(databases),
                   source=dict(content_sha256='a'*64, files={'src/defenseviz/'+p:v for p,v in m._installed_files().items()}))
    m.write_json(cache / 'runtime.json', profile)
    return cache, databases, profile


def test_user_cache_respects_explicit_and_environment(tmp_path, monkeypatch):
    monkeypatch.setenv('DEFENSEVIZ_CACHE', str(tmp_path/'chosen'))
    assert m.cache_directory() == tmp_path/'chosen'
    assert m.cache_directory(tmp_path/'explicit') == tmp_path/'explicit'
    monkeypatch.delenv('DEFENSEVIZ_CACHE')
    monkeypatch.setenv('XDG_CACHE_HOME', str(tmp_path/'xdg'))
    assert m.cache_directory() == tmp_path/'xdg/defenseviz'


def test_corrupt_and_missing_reference_never_reuse(managed):
    cache, databases, _ = managed
    assert m.load_profile()['ready']
    (databases/'reference.fasta').write_text('changed')
    with pytest.raises(ValueError, match='changed or are incomplete'):
        m.load_profile()
    (databases/'reference.fasta').unlink()
    with pytest.raises(ValueError, match='changed or are incomplete'):
        m.load_profile()


def test_profile_cannot_escape_cache(managed, tmp_path):
    cache, _, profile = managed
    profile['resources'] = str(tmp_path)
    with pytest.raises(ValueError, match='outside'):
        m._verify_profile(profile, cache)


@pytest.mark.skipif(os.name != 'posix', reason='Kernel flock is exercised on Linux/macOS/WSL')
def test_lock_releases_after_failure(tmp_path):
    with pytest.raises(RuntimeError, match='example'):
        with m._setup_lock(tmp_path):
            with pytest.raises(RuntimeError, match='Another'):
                with m._setup_lock(tmp_path):
                    pass
            raise RuntimeError('example')
    with m._setup_lock(tmp_path):
        pass


def test_source_build_context_includes_recipe_and_untracked_bytes(tmp_path, monkeypatch):
    actual = Path(m.__file__).resolve().parents[2]
    state = m.source_state(actual)
    assert 'Dockerfile.dnmbsuite.dockerignore' in state['files']
    assert 'src/defenseviz/managed_runtime.py' in state['files']
    assert not any('.egg-info/' in p or '__pycache__/' in p for p in state['files'])
    # A source checkout with missing files cannot silently fall back to a global image.
    with pytest.raises(ValueError, match='complete DefenseViz checkout'):
        m.source_state(tmp_path)


def test_same_setup_is_offline_and_does_not_rebuild(managed, monkeypatch):
    cache, _, profile = managed
    state = dict(profile['source'], path='/checkout', revision='a'*40)
    monkeypatch.setattr(m, 'source_state', lambda p: state)
    calls = []
    def command(argv, record, name, **kw):
        calls.append(name)
        return json.dumps({'OSType':'linux'}) if name == 'docker_info' else json.dumps([{'Id':IMAGE}])
    monkeypatch.setattr(m, '_command', command)
    monkeypatch.setattr('defenseviz.resource_download.fetch_provider_sources', lambda *a, **k: pytest.fail('Reuse attempted a provider download'))
    result = m.setup_managed(source='/checkout', cache=cache)
    assert result['reused'] and result['provider_downloads'] == 0
    assert calls == ['docker_info', 'image_inspect']
    receipt = json.loads(Path(result['setup_receipt']).read_text())
    assert receipt['ready'] and receipt['reused']


@pytest.mark.skipif(os.name != 'posix', reason='Managed setup targets Linux/macOS/WSL')
def test_setup_failure_never_publishes_ready(tmp_path, monkeypatch):
    cache = tmp_path/'new-cache'
    monkeypatch.setattr(m, '_docker', lambda:'/tools/docker')
    monkeypatch.setattr(m, 'source_state', lambda p:dict(content_sha256='a'*64, path='/checkout', revision='a'*40))
    monkeypatch.setattr(m, '_command', lambda *a, **k:json.dumps({'OSType':'linux'}))
    monkeypatch.setattr('defenseviz.resource_download.fetch_provider_sources', lambda *a, **k:(_ for _ in ()).throw(RuntimeError('provider failed')))
    with pytest.raises(RuntimeError, match='provider failed'):
        m.setup_managed(cache=cache)
    assert not (cache/'runtime.json').exists()
    receipts = list((cache/'setups').glob('*/setup.json'))
    assert len(receipts) == 1 and not json.loads(receipts[0].read_text())['ready']
    with m._setup_lock(cache):
        pass


def test_command_failure_retains_exact_exit_and_logs(tmp_path):
    import sys
    argv = [sys.executable, '-c', 'import sys; print("known stdout"); print("known stderr",file=sys.stderr); sys.exit(7)']
    with pytest.raises(RuntimeError, match='exit 7'):
        m._command(argv, tmp_path, 'failure')
    receipt = json.loads((tmp_path/'failure.json').read_text())
    assert receipt['completed'] and receipt['exit_code'] == 7 and receipt['argv'] == argv
    assert (tmp_path/'failure.stdout').read_text() == 'known stdout\n'


def args(root, out):
    return argparse.Namespace(root=root, out=out, manifest=None, supplement=None, rm_evidence=None,
           system_review_rules=None, system_review_evidence=None, layout=None, title=None, strict_inputs=None)


def test_normal_parent_is_one_network_disabled_command(managed, tmp_path, monkeypatch):
    cache, databases, _ = managed
    parent = tmp_path/'균주 folders';parent.mkdir()
    (parent/'sample one').mkdir();(parent/'sample two').mkdir()
    target = tmp_path/'results with spaces'
    calls = []
    def execute(argv, record, name, **kw):
        calls.append(argv)
        result = record/'results';(result/'figures').mkdir(parents=True)
        (result/'receipt.json').write_text('{}')
        (result/'figures/publication.pdf').write_bytes(b'%PDF-1.7\n')
        return '{}'
    monkeypatch.setattr(m, '_command', execute)
    old = m._inventory(databases)
    result = m.managed_analysis(args(parent,target), parent, {'evaluation':'auto'}, [])
    assert len(calls) == 1 and calls[0][calls[0].index('--network')+1] == 'none'
    assert calls[0][calls[0].index('from-suite')+1] == str(parent)
    assert '--manifest' not in calls[0]
    assert '--evaluation' in calls[0] and 'full' in calls[0]
    assert (target/'figures/publication.pdf').exists() and result['provider_downloads'] == 0
    assert m._inventory(databases) == old
    receipt = json.loads(next(tmp_path.glob('.results with spaces.defenseviz-*/managed_run.json')).read_text())
    assert receipt['runtime_profile_sha256'] == m.digest(cache/'runtime.json')


@pytest.mark.parametrize('mode,engine,inside', [('stored',None,False),('full','/custom/engine.json',False),('full',None,True)])
def test_explicit_native_and_stored_precede_managed(tmp_path, monkeypatch, mode, engine, inside):
    monkeypatch.setattr('defenseviz.native_pipeline.engine_config_path',lambda explicit=None:explicit)
    if inside:monkeypatch.setenv('DEFENSEVIZ_MANAGED_CONTAINER','1')
    monkeypatch.setattr(m,'load_profile',lambda:pytest.fail('Should not load managed cache'))
    assert m.managed_analysis(args(tmp_path,tmp_path/'out'),tmp_path,{'evaluation':mode,'engine_config':engine},[]) is None


def test_existing_output_and_input_nesting_are_rejected_without_writes(managed, tmp_path):
    parent = tmp_path/'inputs';parent.mkdir()
    target = parent/'new'/'results'
    with pytest.raises(ValueError,match='outside'):
        m.managed_analysis(args(parent,target),parent,{},[])
    assert not (parent/'new').exists()
    target=tmp_path/'results';target.mkdir();(target/'keep').write_text('original')
    with pytest.raises(FileExistsError):m.managed_analysis(args(parent,target),parent,{},[])
    assert (target/'keep').read_text()=='original'


def test_output_parent_symlink_cannot_write_inside_inputs(managed,tmp_path):
    parent=tmp_path/'input';parent.mkdir()
    alias=tmp_path/'alias';alias.symlink_to(parent,target_is_directory=True)
    with pytest.raises(ValueError,match='outside'):
        m.managed_analysis(args(parent,alias/'result'),parent,{},[])
    assert list(parent.iterdir())==[]


def test_concurrently_created_output_is_preserved(managed,tmp_path,monkeypatch):
    parent=tmp_path/'input';parent.mkdir();target=tmp_path/'output'
    def execute(argv,record,name,**kw):
        result=record/'results';(result/'figures').mkdir(parents=True)
        (result/'receipt.json').write_text('{}');(result/'figures/publication.pdf').write_bytes(b'%PDF')
        target.mkdir()
        return '{}'
    monkeypatch.setattr(m,'_command',execute)
    with pytest.raises(FileExistsError,match='appeared during analysis'):
        m.managed_analysis(args(parent,target),parent,{},[])
    assert target.is_dir() and list(target.iterdir())==[]
    assert len(list(tmp_path.glob('.output.defenseviz-*/results/receipt.json')))==1


def test_missing_full_setup_is_actionable(tmp_path,monkeypatch):
    monkeypatch.setattr('defenseviz.native_pipeline.engine_config_path',lambda explicit=None:None)
    monkeypatch.setattr(m,'load_profile',lambda:None)
    with pytest.raises(RuntimeError,match='defenseviz setup'):
        m.managed_analysis(args(tmp_path,tmp_path/'out'),tmp_path,{'evaluation':'full'},[])


@pytest.mark.parametrize('option', ['rm_evidence', 'reference_cache', 'rscript'])
def test_explicit_saved_or_native_override_does_not_switch_to_managed(tmp_path,monkeypatch,option):
    monkeypatch.setattr(m,'load_profile',lambda:pytest.fail('Explicit override must bypass managed profile'))
    assert m.managed_analysis(args(tmp_path,tmp_path/'out'),tmp_path,{'evaluation':'auto',option:'chosen'},[]) is None


def test_package_upgrade_cannot_silently_run_stale_image(managed,tmp_path,monkeypatch):
    parent=tmp_path/'input';parent.mkdir()
    current=m._installed_files();current['new_feature.py']='0'*64
    monkeypatch.setattr(m,'_installed_files',lambda:current)
    with pytest.raises(ValueError,match='differs from its managed image'):
        m.managed_analysis(args(parent,tmp_path/'out'),parent,{'evaluation':'full'},[])
    assert not (tmp_path/'out').exists()


def test_docker_mount_does_not_interpolate_shell_and_rejects_ambiguous_delimiters(tmp_path):
    path=tmp_path/'$(echo private)'
    assert '$(echo private)' in m._mount(path,'/inputs')[1]
    with pytest.raises(ValueError,match='commas'):
        m._mount(tmp_path/'one,two','/inputs')

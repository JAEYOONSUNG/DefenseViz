import importlib.util
import json
from pathlib import Path
import pytest

spec=importlib.util.spec_from_file_location('native_install',Path(__file__).resolve().parents[1]/'scripts/install_native_runtime.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_default_sources_are_immutable_public_commits():
    assert module.SOURCES['DNMB']['commit']=='976489e9967e9e4dfacc61f7e83e66b1da0f4f2b'
    assert module.SOURCES['DefenseViz']['commit']=='53a1dc4b8c6595945762e21354fbdc86d2b3731a'
    assert all(v['url'].startswith('https://github.com/JAEYOONSUNG/') for v in module.SOURCES.values())
    assert '/inst/' not in module.PATTERNS and '/data/' not in module.PATTERNS


@pytest.mark.parametrize('bad_commit',[False,True])
def test_source_fetch_uses_sparse_checkout_and_verifies_actual_commit(tmp_path,monkeypatch,bad_commit):
    seen=[];dest=tmp_path/'source'
    def run(argv,*a,**kw):
        seen.append(argv)
        if argv[1]=='init':
            (dest/'.git/info').mkdir(parents=True)
            (dest/'DESCRIPTION').write_text('Package: DefenseViz\n')
            (dest/'LICENSE').write_text('MIT License\n')
            (dest/'R').mkdir();(dest/'R/test.R').write_text('x <- 1')
        if 'rev-parse' in argv:return ('0'*40 if bad_commit else module.SOURCES['DefenseViz']['commit'])+'\n'
        return ''
    monkeypatch.setattr(module,'run',run)
    if bad_commit:
        with pytest.raises(ValueError,match='does not match pin'):
            module.fetch_source('DefenseViz',module.SOURCES['DefenseViz'],dest,[],tmp_path/'commands.json')
    else:
        result=module.fetch_source('DefenseViz',module.SOURCES['DefenseViz'],dest,[],tmp_path/'commands.json')
        assert result['working_tree_clean']
        assert any('--filter=blob:none' in a for a in seen)
        assert 'inst' not in (dest/'.git/info/sparse-checkout').read_text()
        assert not any('shell=True' in str(a) for a in seen)


def fixture(tmp_path,monkeypatch,leak=False):
    binary=tmp_path/'bin';binary.mkdir();rscript=binary/'Rscript';rscript.write_text('Rscript');(binary/'R').write_text('R')
    deps=tmp_path/'dependencies';deps.mkdir();calls=[]
    def fetch(name,spec,dest,records,receipt):
        dest.mkdir(parents=True);f=dest/'DESCRIPTION';f.write_text('Package: '+name)
        return dict(**spec,selected_source_sha256={'DESCRIPTION':module.sha(f)})
    def run(argv,*args,**kw):
        argv=list(map(str,argv));calls.append(argv)
        if 'dependency_probe.R' in argv[2]:
            Path(argv[5]).write_text(json.dumps(dict(R_version='R test',R_home='/R',packages=[dict(name='digest',path=str(deps/'digest'))])))
        if 'helper_probe.R' in argv[2]:
            Path(argv[5]).write_text(json.dumps(dict(helpers=[dict(name=n,source_body_matches=True) for n in module.HELPERS],
                packages=[dict(name='DNMB',path=str(deps/'DNMB' if leak else Path(argv[6])/'DNMB'))])))
        return ''
    monkeypatch.setattr(module,'fetch_source',fetch);monkeypatch.setattr(module,'run',run)
    return rscript,deps,calls


def test_install_into_fresh_private_library_and_preserve_dependency_environment(tmp_path,monkeypatch):
    rscript,deps,calls=fixture(tmp_path,monkeypatch)
    out=tmp_path/'new';result=module.install_runtime(out,str(rscript),[str(deps)])
    installs=[a for a in calls if 'INSTALL' in a]
    assert [Path(a[-1]).name for a in installs]==['DefenseViz','DNMB']
    assert all('--library='+str(out/'library') in a for a in installs)
    assert result['automatic_dependency_installations']==0
    assert not result['full_analysis_executed'] and result['helper_source_matches']==19
    assert not list(deps.iterdir())
    with pytest.raises(FileExistsError):module.install_runtime(out,str(rscript),[str(deps)])


def test_old_namespace_fallback_is_not_success(tmp_path,monkeypatch):
    rscript,deps,_=fixture(tmp_path,monkeypatch,leak=True)
    out=tmp_path/'failed'
    with pytest.raises(RuntimeError,match='Old DNMB namespace'):
        module.install_runtime(out,str(rscript),[str(deps)])
    assert not json.loads((out/'failure.json').read_text())['installed']


def test_installer_cannot_write_under_dependency_library(tmp_path,monkeypatch):
    rscript,deps,_=fixture(tmp_path,monkeypatch)
    with pytest.raises(ValueError,match='outside dependency'):
        module.install_runtime(deps/'new',str(rscript),[str(deps)])
    assert not (deps/'new').exists()

import json
from pathlib import Path

from defenseviz import cli


def test_setup_cli_passes_user_cache_source_and_refresh(tmp_path,monkeypatch,capsys):
    seen={}
    def setup(**kwargs):seen.update(kwargs);return {'ready':True}
    monkeypatch.setattr('defenseviz.managed_runtime.setup_managed',setup)
    assert cli.main(['setup','--source',str(tmp_path),'--cache',str(tmp_path/'cache'),'--refresh-databases','--cpu','2'])==0
    assert seen==dict(source=tmp_path,cache=tmp_path/'cache',refresh=True,cpu=2)
    assert json.loads(capsys.readouterr().out)['ready']


def test_parent_command_uses_managed_runtime_without_extra_plot_command(tmp_path,monkeypatch,capsys):
    calls=[]
    def run(args,root,settings,refs):
        calls.append((args,root,settings,refs));return {'completed':True,'figures':'figures/publication.pdf'}
    monkeypatch.setattr('defenseviz.managed_runtime.managed_analysis',run)
    monkeypatch.setattr('defenseviz.suite_pipeline.run_suite',lambda *a,**k: (_ for _ in ()).throw(AssertionError('Should run inside managed image')))
    assert cli.main(['from-suite',str(tmp_path),'--out',str(tmp_path/'out')])==0
    assert len(calls)==1 and calls[0][1]==tmp_path
    assert json.loads(capsys.readouterr().out)['completed']


def test_project_settings_are_resolved_before_managed_dispatch(tmp_path,monkeypatch):
    data=tmp_path/'data';data.mkdir()
    project=tmp_path/'defenseviz_project.json'
    project.write_text(json.dumps(dict(schema_version=1,root='data',evaluation='stored',engine_config='chosen.json')))
    seen=[]
    def run(args,root,settings,refs):seen.append((root,settings,refs));return {'completed':True}
    monkeypatch.setattr('defenseviz.managed_runtime.managed_analysis',run)
    assert cli.main(['from-suite',str(tmp_path),'--out',str(tmp_path/'out')])==0
    assert seen[0][0]==data and seen[0][1]['evaluation']=='stored'
    assert seen[0][1]['engine_config']==tmp_path/'chosen.json'
    assert seen[0][2][0]['path']==str(project)

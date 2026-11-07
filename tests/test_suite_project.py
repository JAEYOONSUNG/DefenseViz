import hashlib
import json
from pathlib import Path

import pytest

from defenseviz import cli
from defenseviz.suite_project import bind_project_inputs, resolve_suite_project


def project(tmp_path, **settings):
    (tmp_path / 'inputs').mkdir(exist_ok=True)
    path = tmp_path / 'defenseviz_project.json'
    path.write_text(json.dumps({'schema_version': 1, 'root': 'inputs', **settings}))
    return path


def test_parent_project_resolves_paths_independently_of_cwd(tmp_path, monkeypatch):
    path = project(tmp_path, manifest='controls/manifest.json', evaluation='stored',
                   external_type_projects=[{'project': 'external/case.json', 'sample_ids': ['case']}])
    monkeypatch.chdir(tmp_path.parent)
    root, settings, refs = resolve_suite_project(tmp_path)
    assert root == tmp_path / 'inputs'
    assert settings['manifest'] == tmp_path / 'controls/manifest.json'
    assert settings['external_type_projects'][0]['project'] == tmp_path / 'external/case.json'
    assert settings['evaluation'] == 'stored'
    assert refs[0]['sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_explicit_options_override_project_but_defaults_do_not(tmp_path):
    project(tmp_path, evaluation='stored', strict_inputs=True, title='Saved comparison')
    _, settings, _ = resolve_suite_project(tmp_path, {'evaluation': None, 'title': None})
    assert settings['evaluation'] == 'stored' and settings['strict_inputs'] is True
    _, settings, _ = resolve_suite_project(tmp_path, {'evaluation': 'full', 'title': 'New comparison'})
    assert settings['evaluation'] == 'full' and settings['title'] == 'New comparison'


def test_no_project_keeps_normal_suite_discovery(tmp_path):
    root, settings, refs = resolve_suite_project(tmp_path)
    assert root == tmp_path and settings == {'evaluation': 'auto', 'strict_inputs': False} and not refs


@pytest.mark.parametrize('settings', [
    {'unknown': 'ignored'}, {'strict_inputs': 'true'}, {'manifest': False},
    {'evaluation': 'guess'}, {'external_type_projects': 'case.json'},
    {'external_type_projects': [{'project': 'case.json', 'sample_ids': ['same', 'same']}]},
    {'external_type_projects': [{'project': 'case.json', 'expected_counts': 7}]},
])
def test_invalid_project_settings_are_rejected(tmp_path, settings):
    project(tmp_path, **settings)
    with pytest.raises(ValueError):
        resolve_suite_project(tmp_path)


def test_project_is_pinned_and_change_is_rejected(tmp_path):
    path = project(tmp_path)
    _, _, refs = resolve_suite_project(path)
    data = {}
    bind_project_inputs(data, refs)
    assert data['input_sha256'][str(path)] == refs[0]['sha256']
    path.write_text(path.read_text() + '\n')
    with pytest.raises(ValueError, match='changed'):
        bind_project_inputs(data, refs)


def test_cli_uses_portable_project_in_one_command(tmp_path, monkeypatch, capsys):
    project(tmp_path, manifest='manifest.json', evaluation='stored', title='Independent comparison',
            layout='layout.json', external_type_projects=['external.json'])
    called = {}
    def run(root, out, **options):
        called.update(root=root, out=out, **options)
        return {'ok': True}
    monkeypatch.setattr('defenseviz.suite_pipeline.run_suite', run)
    assert cli.main(['from-suite', str(tmp_path), '--out', str(tmp_path / 'result')]) == 0
    assert called['root'] == tmp_path / 'inputs'
    assert called['manifest'] == tmp_path / 'manifest.json'
    assert called['evaluation'] == 'stored' and called['title'] == 'Independent comparison'
    assert called['external_type_projects'][0]['project'] == tmp_path / 'external.json'
    assert called['project_refs'][0]['sha256']
    assert json.loads(capsys.readouterr().out) == {'ok': True}

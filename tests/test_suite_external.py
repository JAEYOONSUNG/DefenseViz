import copy
import hashlib
import json
from pathlib import Path

import pytest

from defenseviz import cli
from defenseviz.suite_demo import write_suite_project
from defenseviz.suite_external import prepare_external_type_cases
from defenseviz.suite_pipeline import analyze_suite, prepare_suite, run_suite


def external_project(tmp_path):
    paths = write_suite_project(tmp_path / 'external')
    manifest = json.loads(paths['manifest'].read_text())
    for sample in manifest['samples']:
        sample['sample_id'] = 'external-' + sample['sample_id']
        sample['label'] = 'Éxternal ' + sample['label']
    paths['manifest'].write_text(json.dumps(manifest))
    project = tmp_path / 'external/defenseviz_project.json'
    project.write_text(json.dumps({'schema_version': 1, 'root': 'runs',
                                  'manifest': 'cohort.json', 'evaluation': 'stored'}))
    return project


@pytest.mark.parametrize('external_policy', [False, True])
def test_independent_project_exports_bound_case_without_primary_contamination(tmp_path, monkeypatch, external_policy):
    primary = write_suite_project(tmp_path / 'primary')
    project = external_project(tmp_path)
    if external_policy:
        original_adapter = prepare_external_type_cases
        def evaluated(*args, **kwargs):
            cases, execution, inputs, derived = original_adapter(*args, **kwargs)
            execution[0].update(policy_evaluations=7, policy_decisions=5, scoring_runs=1)
            execution[0]['evaluation_execution'].update(requested_mode='auto', rm_mode='automatic')
            return cases, execution, inputs, derived
        monkeypatch.setattr('defenseviz.suite_external.prepare_external_type_cases', evaluated)
    layout = json.loads(primary['layout'].read_text())
    layout.update(page_composition='letter_composite',
                  type_samples=['isolate-a', 'isolate-b', 'external-isolate-a'])
    primary['layout'].write_text(json.dumps(layout))
    data = prepare_suite(primary['root'], manifest=primary['manifest'],
                         supplement=primary['supplement'], evaluation='stored')
    _, expected = analyze_suite(data)
    root = primary.pop('root')
    result = run_suite(root, tmp_path / 'result', **primary, evaluation='stored',
                       external_type_projects=[{'project': project, 'sample_ids': ['external-isolate-a']}])
    published = json.loads((tmp_path / 'result/analysis.json').read_text())
    receipt = json.loads((tmp_path / 'result/receipt.json').read_text())
    assert result['samples'] == 3
    assert published['global_summary'] == expected['global_summary']
    assert not any(r['sample_id'].startswith('external-') for r in published['gene_ledger'])
    cases = json.loads((tmp_path / 'result/external_type_cases.json').read_text())
    assert len(cases) == 1 and cases[0]['sample_ids'] == ['external-isolate-a']
    from defenseviz.plot.publication.display_scope import model_digest
    assert cases[0]['result_sha256'] == model_digest(cases[0]['result'])
    ref = cases[0]['source_refs'][0]
    assert hashlib.sha256((tmp_path / 'result' / ref['path']).read_bytes()).hexdigest() == ref['sha256']
    assert receipt['hmm_search_runs'] == receipt['homology_search_runs'] == 0
    if external_policy:
        assert receipt['operation'] == 'primary_suite_and_independent_external_evidence_evaluation_and_publication'
        assert receipt['scoring_runs'] == 1
    assert receipt['external_type_case_execution'][0]['analyzed_sample_ids'] == [
        'external-isolate-a', 'external-isolate-b', 'external-isolate-c']
    assert len(list((tmp_path / 'result/figures').glob('publication_page*.svg'))) == 2


def test_external_source_labels_resolve_in_own_project(tmp_path, monkeypatch):
    project = external_project(tmp_path)
    source = tmp_path / 'external/native.tsv'
    source.write_text('external source')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    import defenseviz.suite_pipeline as pipeline
    monkeypatch.setattr(pipeline, 'prepare_suite', lambda *a, **k: {
        'input_sha256': {str(source): digest}, 'suite': {'sources': [
            {'path': 'native.tsv', 'resolved_path': str(source), 'sha256': digest}]}})
    monkeypatch.setattr(pipeline, 'analyze_suite', lambda d: (d, {
        'samples': [{'sample_id': 'external'}], 'source_ref': {'path': 'native.tsv', 'sha256': digest}}))
    cases, _, _, _ = prepare_external_type_cases([{'project': project}], {'primary'}, tmp_path / 'work')
    assert cases[0]['result']['source_ref']['path'] == str(source)


def test_external_selection_rejects_primary_overlap_and_recursive_projects(tmp_path):
    project = external_project(tmp_path)
    with pytest.raises(ValueError, match='duplicate primary'):
        prepare_external_type_cases([{'project': project}], {'external-isolate-a'}, tmp_path / 'work')
    doc = json.loads(project.read_text())
    doc['external_type_projects'] = ['defenseviz_project.json']
    project.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match='recursively'):
        prepare_external_type_cases([{'project': project}], set(), tmp_path / 'work')


def test_inspect_does_not_silently_ignore_invalid_external_project(tmp_path, capsys):
    project = external_project(tmp_path)
    doc = json.loads(project.read_text())
    doc['external_type_projects'] = ['missing.json']
    project.write_text(json.dumps(doc))
    assert cli.main(['inspect-suite', str(project)]) == 2
    captured = capsys.readouterr()
    assert not captured.out and 'existing directory' in captured.err


@pytest.mark.parametrize('option', ['rm_evidence', 'reference_cache', 'engine_config'])
def test_stored_replay_rejects_evaluation_inputs_before_loading(tmp_path, option):
    with pytest.raises(ValueError, match='Stored replay cannot'):
        prepare_suite(tmp_path, evaluation='stored', **{option: tmp_path / 'unused'})

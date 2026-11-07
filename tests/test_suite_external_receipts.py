"""Focused transaction/receipt tests; all external engines are mocked."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from defenseviz.suite_demo import write_suite_project
from defenseviz.suite_pipeline import analyze_suite, prepare_suite, run_suite


def _external_adapter(tmp_path, monkeypatch, *, fail_engine=False):
    """Return independent synthetic evidence and realistic execution metadata."""
    source = tmp_path / 'external' / 'native.tsv'
    source.parent.mkdir()
    source.write_text('Synthetic external source\n', encoding='utf-8')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    source_ref = {'path': str(source), 'sha256': digest}
    independent = {
        'samples': [{'sample_id': 'external-case', 'label': 'Éxternal case'}],
        'sources': [],
        'type_tables': [{'sample_id': 'external-case', 'rm_type': 'Type IV', 'unique_CDS': 1}],
        'gene_ledger': [{'sample_id': 'external-case', 'gene_id': 'external-gene'}],
        'global_summary': {'integrated_CDS': 1},
    }

    def prepare(entries, primary_data, workspace, **kwargs):
        assert entries and 'external-case' not in {s['sample_id'] for s in primary_data['samples']}
        full = Path(workspace) / 'case001' / 'full_pipeline'
        full.mkdir(parents=True)
        log = full / 'engine.log'
        log.write_text('Synthetic external full execution log\n', encoding='utf-8')
        if fail_engine:
            raise RuntimeError('synthetic external engine failure')
        from defenseviz.plot.publication.display_scope import model_digest
        case = {'result': deepcopy(independent), 'sample_ids': ['external-case'],
                'source_refs': [source_ref], 'result_sha256': model_digest(independent)}
        execution = {
            'sample_ids': ['external-case'], 'analyzed_sample_ids': ['external-case'],
            'scope': 'independent type-table examples; excluded from primary aggregates',
            'pre_portable_result_sha256': case['result_sha256'], 'source_refs': [source_ref],
            'evaluation_execution': {'requested_mode': 'full', 'fresh_full_execution': True,
                                     'rm_mode': 'full_native_execution'},
            'elapsed_seconds': 0.001, 'scoring_runs': 1, 'hmm_search_runs': 2,
            'homology_search_runs': 3, 'detector_runs': 1, 'original_defense_tool_runs': 0,
            'policy_evaluations': 7, 'policy_decisions': 5, 'preserved_native_decisions': 2,
            'sequence_feature_evaluations': 11,
            'full_execution': {'hmm_search_runs': 2, 'homology_search_runs': 3,
                               'detector_runs': {'DNMB': 1}},
            'external_processes': [{'argv': ['synthetic-engine', '--sample', 'external-case']}],
        }
        return [case], [execution], {str(source): digest}, {
            str(log): hashlib.sha256(log.read_bytes()).hexdigest()}

    monkeypatch.setattr('defenseviz.suite_external.prepare_external_type_cases', prepare)
    monkeypatch.setattr('defenseviz.suite_pipeline.software_state', lambda: {'source_sha256': 'f' * 64})
    return source


def _renderer(monkeypatch, *, fail=False):
    observed = []

    def render(primary, output, *, external_type_cases=None, **kwargs):
        observed.append(deepcopy(primary))
        assert all(s['sample_id'] != 'external-case' for s in primary['samples'])
        assert all(r['sample_id'] != 'external-case' for r in primary['gene_ledger'])
        assert external_type_cases and external_type_cases[0]['sample_ids'] == ['external-case']
        # Exercise the real envelope/hash binding, without plotting or engines.
        from defenseviz.plot.publication.display_scope import type_case_projection
        displayed, _, bindings = type_case_projection(
            primary, external_type_cases, None, lambda value: None, output.parent)
        assert any(s['sample_id'] == 'external-case' for s in displayed['samples'])
        assert bindings[0]['included_in_primary_aggregates'] is False
        if fail:
            raise ValueError('synthetic downstream rendering failure')
        pdf, svg = output / 'publication.pdf', output / 'publication.svg'
        pdf.write_bytes(b'%PDF-synthetic-receipt-test')
        svg.write_text('<svg/>', encoding='utf-8')
        return [pdf, svg]

    monkeypatch.setattr('defenseviz.plot.publication.render_publication', render)
    return observed


@pytest.mark.parametrize('failure_stage', ['external_engine', 'downstream_render'])
def test_external_full_failure_logs_survive_saved_primary(tmp_path, monkeypatch, failure_stage):
    paths = write_suite_project(tmp_path / 'primary')
    _external_adapter(tmp_path, monkeypatch, fail_engine=failure_stage == 'external_engine')
    _renderer(monkeypatch, fail=failure_stage == 'downstream_render')
    root = paths.pop('root')
    output = tmp_path / 'publication'
    with pytest.raises(RuntimeError, match='full execution evidence preserved'):
        run_suite(root, output, **paths, evaluation='stored',
                  external_type_projects=[{'project': tmp_path / 'external/project.json'}])
    assert not output.exists()
    assert not list(tmp_path.glob('.publication.staging-*'))
    failed = list(tmp_path.glob('.publication.failed-*'))
    assert len(failed) == 1
    log = failed[0] / 'evidence/external_type_cases/case001/full_pipeline/engine.log'
    assert log.read_text(encoding='utf-8') == 'Synthetic external full execution log\n'
    failure = json.loads((failed[0] / 'failure.json').read_text(encoding='utf-8'))
    assert failure['published'] is False and 'preserved' in failure['reason']
    assert not (failed[0] / 'receipt.json').exists()
    assert not (failed[0] / 'evidence/reference/full_pipeline').exists()


def test_saved_primary_external_fresh_receipt_has_scoped_counters(tmp_path, monkeypatch):
    paths = write_suite_project(tmp_path / 'primary')
    source = _external_adapter(tmp_path, monkeypatch)
    observed = _renderer(monkeypatch)
    original = prepare_suite(paths['root'], manifest=paths['manifest'],
                             supplement=paths['supplement'], evaluation='stored')
    _, expected = analyze_suite(original)
    root = paths.pop('root')
    output = tmp_path / 'publication'
    response = run_suite(root, output, **paths, evaluation='stored',
                         external_type_projects=[{'project': tmp_path / 'external/project.json'}])
    receipt = json.loads((output / 'receipt.json').read_text(encoding='utf-8'))
    primary = json.loads((output / 'analysis.json').read_text(encoding='utf-8'))
    native = json.loads((output / 'native_bundle.json').read_text(encoding='utf-8'))
    assert response['samples'] == len(expected['samples'])
    assert primary['global_summary'] == expected['global_summary']
    assert {(r['sample_id'], r['gene_id']) for r in primary['gene_ledger']} == {
        (r['sample_id'], r['gene_id']) for r in expected['gene_ledger']}
    assert native['samples'] == original['samples']
    assert observed and observed[0]['global_summary'] == expected['global_summary']
    assert receipt['operation'] == 'primary_suite_and_independent_external_engine_evaluation_and_publication'
    assert receipt['primary_operation'] == 'saved_suite_result_import_aggregation_and_render'
    assert receipt['evaluation_execution']['requested_mode'] == 'stored'
    assert receipt['evaluation_execution']['fresh_full_execution'] is False
    assert 'primary cohort' in receipt['evaluation_execution_scope']
    assert receipt['full_execution'] == {}
    assert receipt['scoring_runs'] == receipt['detector_runs'] == 1
    assert receipt['hmm_search_runs'] == 2 and receipt['homology_search_runs'] == 3
    assert receipt['sequence_feature_evaluations'] == 11
    assert receipt['policy_evaluations'] == receipt['policy_decisions'] == 0
    assert 'primary cohort' in receipt['policy_counts_scope']
    external = receipt['external_type_case_execution'][0]
    assert external['sample_ids'] == ['external-case']
    assert external['evaluation_execution']['fresh_full_execution'] is True
    assert external['policy_evaluations'] == 7 and external['policy_decisions'] == 5
    assert external['external_processes'][0]['argv'][0] == 'synthetic-engine'
    assert 'excluded from all primary-cohort aggregates' in receipt['scientific_scope']
    assert any(s['name'] == 'independent_external_type_case_analysis'
               and s['sample_ids'] == ['external-case'] for s in receipt['stages'])
    assert not any(s['name'] == 'fresh_DNMB_HMM_BLAST_context_and_feature_execution' for s in receipt['stages'])
    bound = [r for r in receipt['inputs'] if (root / r['path']).resolve() == source]
    assert len(bound) == 1 and bound[0]['sha256'] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert all(hashlib.sha256((output / path).read_bytes()).hexdigest() == sha
               for path, sha in receipt['outputs'].items())


def test_family_scope_summary_qualifies_missing_and_partial_curation():
    from defenseviz.suite_pipeline import _family_scope_capability
    data = {'samples': [{'sample_id': 'a'}, {'sample_id': 'b'}],
            'suite': {'capabilities': [{'sample_id': 'a', 'native_curation_status': 'not_exported'},
                                       {'sample_id': 'b', 'native_curation_status': 'not_exported'}]}}
    status, rows = _family_scope_capability(data, 'auto')
    assert status == 'unavailable_no_native_curation_export'
    assert all(r['bound_curation_records'] == 0 for r in rows)
    data['curation_records'] = [{'sample_id': 'a', 'usable_for_attribution': True}]
    status, rows = _family_scope_capability(data, 'auto')
    assert status == 'partially_available_source_bound_native_curation'
    assert rows[0]['bound_curation_records'] == 1 and rows[1]['bound_curation_records'] == 0
    assert _family_scope_capability(data, 'stored')[0] == 'skipped_stored_replay'


def test_empty_family_attribution_is_exported_with_replay_capability(tmp_path, monkeypatch):
    paths = write_suite_project(tmp_path / 'primary')
    def render(model, out, **kwargs):
        pdf, svg = out / 'page.pdf', out / 'page.svg'
        pdf.write_bytes(b'%PDF-test'); svg.write_text('<svg/>', encoding='utf-8')
        return [pdf, svg]
    monkeypatch.setattr('defenseviz.plot.publication.render_publication', render)
    monkeypatch.setattr('defenseviz.suite_pipeline.software_state', lambda: {'source_sha256': {'test': 'stable'}})
    run_suite(paths.pop('root'), tmp_path / 'result', **paths, evaluation='stored')
    audit = json.loads((tmp_path / 'result/family_attribution.json').read_text(encoding='utf-8'))
    assert audit['records'] == []
    assert audit['summary']['review_status'] == 'skipped_stored_replay'
    assert len(audit['summary']['capabilities']) == 3
    assert len((tmp_path / 'result/tables/family_attribution.tsv').read_text(encoding='utf-8').splitlines()) == 1
    receipt = json.loads((tmp_path / 'result/receipt.json').read_text(encoding='utf-8'))
    assert receipt['primary_operation'] == receipt['operation'] == 'saved_suite_result_import_aggregation_and_render'
    assert 'Replay of saved' in receipt['scientific_scope']
    assert 'automatic stored/sequence evidence evaluation' not in receipt['scientific_scope']


@pytest.mark.parametrize('flags,expected', [([], {}), (['--no-strict-inputs'], {'strict_inputs': False}),
    (['--strict-inputs', '--rscript', 'explicit-R'], {'strict_inputs': True, 'rscript': 'explicit-R'})])
def test_cli_explicit_execution_settings_reach_external_projects(tmp_path, monkeypatch, flags, expected):
    from defenseviz import cli
    paths = write_suite_project(tmp_path / 'primary')
    project = tmp_path / 'primary/defenseviz_project.json'
    project.write_text(json.dumps({'schema_version': 1, 'root': 'runs', 'manifest': 'cohort.json',
        'strict_inputs': True, 'rscript': 'project-R'}), encoding='utf-8')
    observed = []
    def run(*args, **kwargs):
        observed.append(kwargs)
        return {'valid': True}
    monkeypatch.setattr('defenseviz.suite_pipeline.run_suite', run)
    assert cli.main(['from-suite', str(project), '--out', str(tmp_path / 'out'), *flags]) == 0
    assert observed[0]['external_execution_overrides'] == expected
    assert observed[0]['strict_inputs'] is expected.get('strict_inputs', True)
    assert observed[0]['rscript'] == expected.get('rscript', 'project-R')


@pytest.mark.parametrize('overrides,expected', [(None, ('case-R', True)),
    ({'rscript': 'explicit-R', 'strict_inputs': False}, ('explicit-R', False))])
def test_external_projects_keep_own_settings_unless_explicitly_overridden(tmp_path, monkeypatch, overrides, expected):
    from defenseviz.suite_external import prepare_external_type_cases
    write_suite_project(tmp_path / 'external')
    project = tmp_path / 'external/defenseviz_project.json'
    doc = {'schema_version': 1, 'root': 'runs', 'manifest': 'cohort.json', 'evaluation': 'stored'}
    doc.update(rscript='case-R', strict_inputs=True); project.write_text(json.dumps(doc), encoding='utf-8')
    observed = []
    original = prepare_suite
    def prepare(*a, **kw):
        observed.append((kw.get('rscript'), kw.get('strict_inputs')))
        # CLI propagation is under test; avoid invoking a deliberately fake R.
        kw['rscript'] = None; kw['strict_inputs'] = False
        return original(*a, **kw)
    monkeypatch.setattr('defenseviz.suite_pipeline.prepare_suite', prepare)
    prepare_external_type_cases([{'project': project}], {'samples': []}, tmp_path / 'work', execution_overrides=overrides)
    assert observed == [expected]

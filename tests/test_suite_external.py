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
    manifest = json.loads(paths['manifest'].read_text(encoding='utf-8'))
    for sample in manifest['samples']:
        sample['sample_id'] = 'external-' + sample['sample_id']
        sample['label'] = 'Éxternal ' + sample['label']
    paths['manifest'].write_text(json.dumps(manifest), encoding='utf-8')
    # Independent synthetic specimens, not relabelled primary protein sets.
    for genes in (tmp_path / 'external/runs').glob('*/dnmb_module_rebasefinder/rebasefinder_input.tsv'):
        genes.write_text(genes.read_text(encoding='utf-8').replace('\tM', '\tW'), encoding='utf-8')
    project = tmp_path / 'external/defenseviz_project.json'
    project.write_text(json.dumps({'schema_version': 1, 'root': 'runs',
                                  'manifest': 'cohort.json', 'evaluation': 'stored'}), encoding='utf-8')
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
    layout = json.loads(primary['layout'].read_text(encoding='utf-8'))
    layout.update(page_composition='letter_composite',
                  type_samples=['isolate-a', 'isolate-b', 'external-isolate-a'])
    primary['layout'].write_text(json.dumps(layout), encoding='utf-8')
    data = prepare_suite(primary['root'], manifest=primary['manifest'],
                         supplement=primary['supplement'], evaluation='stored')
    _, expected = analyze_suite(data)
    root = primary.pop('root')
    result = run_suite(root, tmp_path / 'result', **primary, evaluation='stored',
                       external_type_projects=[{'project': project, 'sample_ids': ['external-isolate-a']}])
    published = json.loads((tmp_path / 'result/analysis.json').read_text(encoding='utf-8'))
    receipt = json.loads((tmp_path / 'result/receipt.json').read_text(encoding='utf-8'))
    assert result['samples'] == 3
    assert published['global_summary'] == expected['global_summary']
    assert not any(r['sample_id'].startswith('external-') for r in published['gene_ledger'])
    cases = json.loads((tmp_path / 'result/external_type_cases.json').read_text(encoding='utf-8'))
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
    source.write_text('external source', encoding='utf-8')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    import defenseviz.suite_pipeline as pipeline
    monkeypatch.setattr(pipeline, 'prepare_suite', lambda *a, **k: {
        'input_sha256': {str(source): digest}, 'suite': {'sources': [
            {'path': 'native.tsv', 'resolved_path': str(source), 'sha256': digest}]}})
    monkeypatch.setattr(pipeline, 'analyze_suite', lambda d: (d, {
        'samples': [{'sample_id': 'external'}], 'source_ref': {'path': 'native.tsv', 'sha256': digest}}))
    cases, _, _, _ = prepare_external_type_cases([{'project': project}], {'samples': [{'sample_id': 'primary'}]}, tmp_path / 'work')
    assert cases[0]['result']['source_ref']['path'] == str(source)


def test_external_selection_rejects_primary_overlap_and_recursive_projects(tmp_path):
    project = external_project(tmp_path)
    with pytest.raises(ValueError, match='duplicate primary'):
        prepare_external_type_cases([{'project': project}], {'samples': [{'sample_id': 'external-isolate-a'}]}, tmp_path / 'work')
    doc = json.loads(project.read_text(encoding='utf-8'))
    doc['external_type_projects'] = ['defenseviz_project.json']
    project.write_text(json.dumps(doc), encoding='utf-8')
    with pytest.raises(ValueError, match='recursively'):
        prepare_external_type_cases([{'project': project}], {'samples': []}, tmp_path / 'work')


def test_inspect_does_not_silently_ignore_invalid_external_project(tmp_path, capsys):
    project = external_project(tmp_path)
    doc = json.loads(project.read_text(encoding='utf-8'))
    doc['external_type_projects'] = ['missing.json']
    project.write_text(json.dumps(doc), encoding='utf-8')
    assert cli.main(['inspect-suite', str(project)]) == 2
    captured = capsys.readouterr()
    assert not captured.out and 'existing directory' in captured.err


@pytest.mark.parametrize('option', ['rm_evidence', 'reference_cache', 'engine_config'])
def test_stored_replay_rejects_evaluation_inputs_before_loading(tmp_path, option):
    with pytest.raises(ValueError, match='Stored replay cannot'):
        prepare_suite(tmp_path, evaluation='stored', **{option: tmp_path / 'unused'})


def _identity(sample_id, *, protein='a', source=None, module='suite', role='protein_sequences'):
    data = {'samples': [{'sample_id': sample_id, 'label': 'Same display label'}],
            'genes': [{'sample_id': sample_id, 'gene_id': 'renamable', 'protein_sha256': protein * 64}],
            'input_sha256': {}, 'input_roles': {}}
    if source is not None:
        path = str(source.resolve())
        data['input_sha256'][path] = hashlib.sha256(source.read_bytes()).hexdigest()
        data['input_roles'][path] = [{'sample_id': sample_id, 'module': module, 'kind': role}]
    return data


@pytest.mark.parametrize('axis', ['accession', 'genbank', 'fasta', 'proteins', 'specimen_path'])
def test_relabelled_specimens_rejected_on_bound_identity(tmp_path, axis):
    from defenseviz.suite_external import cohort_identities, reject_cohort_overlap
    first, second = tmp_path / 'first.dat', tmp_path / 'second.dat'
    first.write_bytes(b'bound specimen sequence'); second.write_bytes(first.read_bytes())
    primary = _identity('primary', protein='a')
    external = _identity('external', protein='b')
    if axis == 'accession':
        primary['samples'][0]['sample_id'] = 'GCF_000000001.1'
        external['samples'][0]['sample_id'] = 'external_GCF_000000001.2_case'
    elif axis in {'genbank', 'fasta'}:
        module = 'genbank' if axis == 'genbank' else 'suite'
        primary = _identity('primary', protein='a', source=first, module=module)
        external = _identity('external', protein='b', source=second, module=module)
    elif axis == 'proteins':
        external['genes'][0]['protein_sha256'] = 'a' * 64
        external['genes'][0]['gene_id'] = 'different-locus'
    else:
        primary = _identity('primary', protein='a', source=first, role='native_calls')
        external = _identity('external', protein='b', source=first, role='native_calls')
    with pytest.raises(ValueError, match='duplicate primary'):
        reject_cohort_overlap(cohort_identities(external), cohort_identities(primary))


def test_shared_resource_hashes_and_one_shared_protein_do_not_imply_specimen_overlap(tmp_path):
    from defenseviz.suite_external import cohort_identities, reject_cohort_overlap
    shared = tmp_path / 'reference.fasta'; shared.write_bytes(b'common database')
    a, b = _identity('a'), _identity('b')
    a['genes'].append({'sample_id': 'a', 'gene_id': 'a2', 'protein_sha256': 'b' * 64})
    b['genes'].append({'sample_id': 'b', 'gene_id': 'b2', 'protein_sha256': 'c' * 64})
    for data in (a, b):
        data['input_sha256'][str(shared)] = hashlib.sha256(shared.read_bytes()).hexdigest()
        # Resource/config inputs have no sample-specific specimen role.
    reject_cohort_overlap(cohort_identities(b), cohort_identities(a))
    a['genes'] = b['genes'] = []
    reject_cohort_overlap(cohort_identities(b), cohort_identities(a))


def test_external_overlap_guard_precedes_full_engine_execution(tmp_path, monkeypatch):
    import defenseviz.native_pipeline as native
    monkeypatch.setattr(native, 'engine_config_path', lambda value: tmp_path / 'engine.json')
    monkeypatch.setattr(native, 'read_engine_config', lambda value: {'rscript': 'Rscript'})
    monkeypatch.setattr('defenseviz.io.suite.load_suite', lambda *a, **k: _identity('duplicate'))
    monkeypatch.setattr(native, 'execute_full_pipeline', lambda *a, **k: pytest.fail('engine must not start'))
    def guard(data):
        raise ValueError('duplicate specimen preflight')
    with pytest.raises(ValueError, match='duplicate specimen preflight'):
        prepare_suite(tmp_path, evaluation='full', evidence_workspace=tmp_path / 'work', _input_guard=guard)


def test_undisplayed_external_sample_must_also_be_independent(tmp_path):
    primary = write_suite_project(tmp_path / 'primary')
    project = external_project(tmp_path)
    genes = tmp_path / 'external/runs/isolate-a/dnmb_module_rebasefinder/rebasefinder_input.tsv'
    genes.write_text(genes.read_text(encoding='utf-8').replace('\tW', '\tM'), encoding='utf-8')
    data = prepare_suite(primary['root'], manifest=primary['manifest'], evaluation='stored')
    with pytest.raises(ValueError, match='imported_protein_multiset'):
        prepare_external_type_cases([{'project': project, 'sample_ids': ['external-isolate-b']}], data, tmp_path / 'work')


def _small_renderer(monkeypatch, mutate=None):
    models = []
    def render(model, out, **kwargs):
        models.append(copy.deepcopy(model))
        if mutate:
            mutate()
        pdf, svg = out / 'page.pdf', out / 'page.svg'
        pdf.write_bytes(b'%PDF-test'); svg.write_text('<svg/>', encoding='utf-8')
        return [pdf, svg]
    monkeypatch.setattr('defenseviz.plot.publication.render_publication', render)
    monkeypatch.setattr('defenseviz.suite_pipeline.software_state', lambda: {'source_sha256': {'test': 'stable'}})
    return models


def test_external_and_wrapper_hashes_cannot_change_primary_analysis(tmp_path, monkeypatch):
    from defenseviz.suite_project import resolve_suite_project, ANALYSIS_OPTIONS
    paths = write_suite_project(tmp_path / 'primary')
    external = external_project(tmp_path)
    project = tmp_path / 'primary/defenseviz_project.json'
    doc = {'schema_version': 1, 'root': 'runs', 'manifest': 'cohort.json',
           'supplement': 'review.json', 'layout': 'layout.json', 'evaluation': 'stored'}
    models = _small_renderer(monkeypatch)
    receipts = []
    for index, attached in enumerate((False, True)):
        if attached:
            doc['external_type_projects'] = [str(external)]
            paths['layout'].write_text(json.dumps({'type_samples': ['external-isolate-a']}), encoding='utf-8')
        project.write_text(json.dumps(doc), encoding='utf-8')
        root, settings, refs = resolve_suite_project(project)
        run_suite(root, tmp_path / f'result{index}',
                  **{k: v for k, v in settings.items() if k in ANALYSIS_OPTIONS},
                  layout=settings['layout'], project_refs=refs,
                  external_type_projects=settings.get('external_type_projects'))
        receipts.append(json.loads((tmp_path / f'result{index}/receipt.json').read_text(encoding='utf-8')))
    assert (tmp_path / 'result0/analysis.json').read_bytes() == (tmp_path / 'result1/analysis.json').read_bytes()
    assert models[0] == models[1]
    assert receipts[0]['input_scopes']['primary_scientific'] == receipts[1]['input_scopes']['primary_scientific']
    assert receipts[0]['input_scopes']['workflow_and_presentation'] != receipts[1]['input_scopes']['workflow_and_presentation']
    assert not receipts[0]['input_scopes']['external'] and receipts[1]['input_scopes']['external']
    external_paths = {r['path'] for r in receipts[1]['input_scopes']['external']}
    assert external_paths <= {r['path'] for r in receipts[1]['inputs']}
    execution = receipts[1]['external_type_case_execution'][0]
    assert 'pre_portable_result_sha256' in execution and 'result_sha256' not in execution
    assert execution['rendered_result_sha256']


def test_separate_external_inputs_still_rehashed_after_render(tmp_path, monkeypatch):
    paths = write_suite_project(tmp_path / 'primary'); external = external_project(tmp_path)
    _small_renderer(monkeypatch, lambda: external.write_text('{}', encoding='utf-8'))
    with pytest.raises(ValueError, match='Input changed during suite publication'):
        run_suite(paths.pop('root'), tmp_path / 'result', **paths, evaluation='stored',
                  external_type_projects=[{'project': external}])
    assert not (tmp_path / 'result').exists()


def test_empty_sequence_exports_do_not_prove_same_specimen(tmp_path):
    from defenseviz.suite_external import cohort_identities, reject_cohort_overlap
    a, b = tmp_path / 'a.faa', tmp_path / 'b.faa'
    a.write_bytes(b''); b.write_bytes(b'')
    primary = _identity('a', protein='a', source=a)
    external = _identity('b', protein='b', source=b)
    primary['genes'] = external['genes'] = []
    left, right = cohort_identities(primary), cohort_identities(external)
    assert left[0]['identity_status'] == 'unassessed_no_assembly_or_sequence_binding'
    reject_cohort_overlap(right, left)

"""Source binding and atomic publication safety, independent of plot geometry."""
import copy
import hashlib
import json
from pathlib import Path
import sys
from types import ModuleType

import pytest

from defenseviz.supplement import apply_supplement
from defenseviz.suite_pipeline import run_suite


def binding_fixture(tmp_path):
    source = tmp_path / 'native.tsv'
    source.write_text('gene\tvalue\na\tretained\n', encoding='utf-8')
    review = tmp_path / 'review.txt'
    review.write_text('Synthetic explicit component review', encoding='utf-8')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    protein = hashlib.sha256(b'MADEUPPROTEIN').hexdigest()
    data = {'samples': [{'sample_id': 's', 'label': 'Synthetic'}],
            'genes': [{'sample_id': 's', 'gene_id': 'a', 'protein_sha256': protein}],
            'inputs': [source], 'input_sha256': {str(source): sha(source)},
            'input_roles': {str(source): [{'sample_id': 's', 'native_result': True, 'kind': 'native_calls'}]},
            'calls': [], 'coverage': [], 'annotations': []}
    doc = {'schema_version': 1, 'input_sha256': {'native.tsv': sha(source)},
           'annotations': [{'sample_id': 's', 'gene_id': 'a', 'protein_sha256': protein,
                            'status': 'supported', 'assignment_kind': 'reviewed_RM_component',
                            'family': 'Type IV', 'role': 'R', 'include_in_union': True,
                            'rule_id': 'synthetic-rule',
                            'source_refs': [{'path': 'review.txt', 'sha256': sha(review)}]}]}
    supplement = tmp_path / 'supplement.json'
    supplement.write_text(json.dumps(doc), encoding='utf-8')
    return data, doc, supplement


def test_supplement_binds_native_artifact_and_current_protein(tmp_path):
    data, _, path = binding_fixture(tmp_path)
    before = copy.deepcopy(data)
    result = apply_supplement(data, path)
    assert data == before
    assert len(result['annotations']) == 1
    assert len(result['inputs']) == 3
    for raw, expected in result['input_sha256'].items():
        assert hashlib.sha256(Path(raw).read_bytes()).hexdigest() == expected
    assert result['annotations'][0]['source_refs'][0]['path'] == str(tmp_path / 'review.txt')


@pytest.mark.parametrize('mutation,match', [
    ('native', 'binding'), ('protein', 'protein hash'), ('source', 'source hash'),
    ('unbound', 'bindings'), ('unknown', 'Unknown'), ('unknown_gene', 'native gene'),
])
def test_stale_or_unbound_supplements_are_rejected(tmp_path, mutation, match):
    data, doc, path = binding_fixture(tmp_path)
    row = doc['annotations'][0]
    if mutation == 'native': doc['input_sha256']['native.tsv'] = '0' * 64
    if mutation == 'protein': row['protein_sha256'] = '0' * 64
    if mutation == 'source': row['source_refs'][0]['sha256'] = '0' * 64
    if mutation == 'unbound': doc['input_sha256'] = {}
    if mutation == 'unknown': doc['hardcoded_count'] = 500
    if mutation == 'unknown_gene': row['gene_id'] = 'not-a'
    path.write_text(json.dumps(doc), encoding='utf-8')
    with pytest.raises(ValueError, match=match): apply_supplement(data, path)
    assert not data['annotations']


def test_complex_review_binds_every_core_gene(tmp_path):
    data, doc, path = binding_fixture(tmp_path)
    protein = data['genes'][0]['protein_sha256']
    refs = doc['annotations'][0]['source_refs']
    doc['groups'] = [{'sample_id': 's', 'group_id': 'context', 'membership_kind': 'saved_context',
                      'members': [{'gene_id': 'a', 'protein_sha256': protein}], 'source_refs': refs}]
    doc['complexes'] = [{'sample_id': 's', 'group_id': 'context', 'complex_id': 'core',
                         'core_gene_ids': ['a'], 'status': 'complete', 'source_refs': refs}]
    path.write_text(json.dumps(doc), encoding='utf-8')
    with pytest.raises(ValueError, match='core_protein_sha256'): apply_supplement(data, path)
    doc['complexes'][0]['core_protein_sha256'] = {'a': protein}
    path.write_text(json.dumps(doc), encoding='utf-8')
    assert len(apply_supplement(data, path)['complexes']) == 1


def fake_stages(monkeypatch, data, *, fail=False, mutate=False):
    """Keep transaction tests independent of scientific/renderer implementations."""
    import defenseviz.suite_pipeline as pipeline
    monkeypatch.setattr(pipeline, 'prepare_suite', lambda *a, **kw: copy.deepcopy(data))
    analysis = ModuleType('defenseviz.analysis.saved')
    analysis.analyze_saved = lambda d: {'samples': d['samples'], 'global_summary': {'integrated_CDS': 1},
                                       'inventory_rows': {'categories': [], 'families': []},
                                       'warnings': [], 'availability': {'effector': 'unassessed'},
                                       'system_review': {'units': [], 'gene_reviews': [], 'rule_evaluations': [],
                                                         'pending_queue': [], 'detection_accounting': [],
                                                         'sample_coverage': [], 'quarantine': [],
                                                         'accounting': {}, 'provenance': {}},
                                       'gene_ledger': [], 'priority_ledger': [], 'cells': [],
                                       'rm_components': [], 'type_tables': [], 'groups_selected': []}
    plot = ModuleType('defenseviz.plot.publication')
    def render(result, out, **kwargs):
        # A stand-in tests transaction/manifest handling, not PDF validity.
        path = out / 'publication.pdf'
        path.write_bytes(b'%PDF-transaction-fixture')
        svg = out / 'publication.svg'
        svg.write_text('<svg/>', encoding='utf-8')
        if fail: raise ValueError('synthetic renderer failure')
        if mutate: Path(data['inputs'][0]).write_text('changed', encoding='utf-8')
        return [path, svg]
    plot.render_publication = render
    model = ModuleType('defenseviz.plot.publication.manuscript_model')
    model.enrich_manuscript_model = lambda bundle, result: result
    monkeypatch.setitem(sys.modules, analysis.__name__, analysis)
    monkeypatch.setitem(sys.modules, plot.__name__, plot)
    monkeypatch.setitem(sys.modules, model.__name__, model)


def test_suite_transaction_records_exact_inputs_and_outputs(tmp_path, monkeypatch):
    data, _, _ = binding_fixture(tmp_path)
    fake_stages(monkeypatch, data)
    summary = run_suite(tmp_path, tmp_path / 'result', argv=['test', 'from-suite'])
    assert summary['samples'] == 1
    receipt = json.loads((tmp_path / 'result/receipt.json').read_text(encoding='utf-8'))
    assert receipt['detector_runs'] == 0 and receipt['threshold_fitting'] is False
    assert receipt['operation'] == 'saved_suite_result_import_aggregation_and_render'
    assert receipt['argv'] == ['test', 'from-suite']
    for name, expected in receipt['outputs'].items():
        assert hashlib.sha256((tmp_path / 'result' / name).read_bytes()).hexdigest() == expected
    assert (tmp_path / 'result/tables/priority_ledger.tsv').read_text().startswith('sample_id\tgene_id\tselection_route')
    assert 'unassessed' in (tmp_path / 'result/tables/availability.tsv').read_text()
    assert str(tmp_path) not in (tmp_path / 'result/native_bundle.json').read_text()
    with pytest.raises(FileExistsError): run_suite(tmp_path, tmp_path / 'result')


def test_automatic_preparation_keeps_original_import_separate(tmp_path):
    from defenseviz.io.suite import load_suite
    from defenseviz.suite_demo import write_suite_project
    from defenseviz.suite_pipeline import prepare_suite
    paths = write_suite_project(tmp_path / 'inputs')
    original = load_suite(paths['root'], manifest=paths['manifest'])
    prepared = prepare_suite(paths['root'], manifest=paths['manifest'])
    assert prepared['_native_import_bundle'] == original
    assert 'rm_evaluation' not in prepared['_native_import_bundle']
    assert 'rm_evaluation' in prepared
    # Later scientific review must not leak into the original source snapshot.
    prepared['saved_decisions'].append({'sample_id': 'synthetic-review'})
    prepared['suite'].setdefault('reference_alignments', []).append({'derived': True})
    assert prepared['_native_import_bundle'] == original


def test_reference_information_is_exported_but_not_a_figure_input(tmp_path, monkeypatch):
    data, _, _ = binding_fixture(tmp_path)
    fake_stages(monkeypatch, data)
    plot = sys.modules['defenseviz.plot.publication']
    render = plot.render_publication
    consumed = []
    def inspect(result, *args, **kwargs):
        consumed.append(copy.deepcopy(result))
        return render(result, *args, **kwargs)
    monkeypatch.setattr(plot, 'render_publication', inspect)
    run_suite(tmp_path, tmp_path / 'result')
    analyzed = json.loads((tmp_path / 'result/analysis.json').read_text())
    assert analyzed['system_information'] == json.loads((tmp_path / 'result/system_information.json').read_text())
    assert consumed and 'system_information' not in consumed[0]
    assert consumed[0] == {key: value for key, value in analyzed.items() if key != 'system_information'}


def test_native_blank_fqn_does_not_abort_publication(tmp_path):
    import csv
    from defenseviz.suite_demo import write_suite_project
    paths = write_suite_project(tmp_path / 'inputs')
    for path in paths['root'].rglob('defensefinder_best_solution_genes.tsv'):
        with path.open(encoding='utf-8', newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            fields = [*reader.fieldnames, 'model_fqn']
            rows = list(reader)
        with path.open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    # No supplement: its native-artifact binding intentionally refers to the
    # unmodified synthetic input. This checks the direct native import route.
    out = tmp_path / 'publication'
    result = run_suite(paths['root'], out, manifest=paths['manifest'])
    assert result['samples'] == 3
    info = json.loads((out / 'system_information.json').read_text(encoding='utf-8'))
    assert any(row['provider'] == 'DefenseFinder' and row['native_model_fqn'] == ''
               and 'native_fqn_field_present_but_missing' in row['reason_codes']
               for row in info['detections'])
    assert info['accounting']['unaccounted_detections'] == 0
    assert (out / 'figures/publication.pdf').stat().st_size > 1000


@pytest.mark.parametrize('reason', ['failure', 'source_changed'])
def test_failed_suite_publication_leaves_no_partial_output(tmp_path, monkeypatch, reason):
    data, _, _ = binding_fixture(tmp_path)
    fake_stages(monkeypatch, data, fail=reason == 'failure', mutate=reason == 'source_changed')
    with pytest.raises(ValueError): run_suite(tmp_path, tmp_path / 'result')
    assert not (tmp_path / 'result').exists()
    assert not list(tmp_path.glob('.result.staging-*'))


def test_suite_cli_options_are_explicit():
    from defenseviz.cli import parser
    args = parser().parse_args(['from-suite', 'native', '--out', 'result', '--strict-inputs', '--supplement', 'review.json', '--rm-evidence', 'evidence.json'])
    assert args.root == Path('native') and args.strict_inputs
    assert args.supplement == Path('review.json')
    assert args.rm_evidence == Path('evidence.json')


def context_supplement_fixture(tmp_path):
    data, doc, path = binding_fixture(tmp_path)
    members = [{'gene_id': 'a', 'protein_sha256': data['genes'][0]['protein_sha256']}]
    data['groups'] = [{'sample_id': 's', 'group_id': 'old', 'membership_kind': 'saved_context',
                       'provider': 'DNMB', 'members': members, 'source_refs': doc['annotations'][0]['source_refs']}]
    doc['groups'] = [{'sample_id': 's', 'group_id': 'current', 'membership_kind': 'saved_context',
                      'members': members, 'supersedes_group_ids': ['old'],
                      'source_refs': doc['annotations'][0]['source_refs']}]
    path.write_text(json.dumps(doc))
    return data, doc, path


def test_explicit_current_context_replaces_only_bound_saved_context(tmp_path):
    data, _, path = context_supplement_fixture(tmp_path)
    before = copy.deepcopy(data)
    result = apply_supplement(data, path)
    assert data == before
    assert [g['group_id'] for g in result['groups']] == ['current']
    assert result['suite']['superseded_context_groups'][0]['group'] == data['groups'][0]
    assert result['groups'][0]['provider'] == 'Review'


@pytest.mark.parametrize('mutation', ['unknown', 'native_system', 'outside_boundary', 'duplicate', 'complex'])
def test_current_context_cannot_hide_unrelated_or_still_referenced_groups(tmp_path, mutation):
    data, doc, path = context_supplement_fixture(tmp_path)
    if mutation == 'unknown': doc['groups'][0]['supersedes_group_ids'] = ['absent']
    if mutation == 'native_system': data['groups'][0]['membership_kind'] = 'native_system'
    if mutation == 'outside_boundary': data['groups'][0]['members'] = []
    if mutation == 'duplicate':
        second = copy.deepcopy(doc['groups'][0]); second['group_id'] = 'other'
        doc['groups'].append(second)
    if mutation == 'complex': data['complexes'] = [{'sample_id': 's', 'group_id': 'old'}]
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError): apply_supplement(data, path)


def test_suite_reevaluation_is_reported_separately_from_engine_run(tmp_path, monkeypatch):
    data, _, _ = binding_fixture(tmp_path)
    data['rm_evaluation'] = {'operation': 'fresh_evaluation_of_stored_RM_evidence', 'scoring_runs': 1, 'records': []}
    fake_stages(monkeypatch, data)
    run_suite(tmp_path, tmp_path / 'result', rm_evidence=tmp_path / 'evidence.json')
    receipt = json.loads((tmp_path / 'result/receipt.json').read_text())
    assert receipt['operation'] == 'saved_suite_evidence_reevaluation_aggregation_and_render'
    assert receipt['scoring_runs'] == 1 and receipt['detector_runs'] == 0
    assert (tmp_path / 'result/rm_evaluation.json').is_file()
    analyzed = json.loads((tmp_path / 'result/analysis.json').read_text())
    assert analyzed['provenance']['scoring_runs'] == 1


def test_import_warnings_reach_cli_and_receipt(tmp_path, monkeypatch, capsys):
    from defenseviz.cli import main
    data, _, _ = binding_fixture(tmp_path)
    data['suite'] = {'warnings': ['RDS metadata unavailable', 'unbound structure identifier'],
                     'commands': [{'argv': ['/recorded/Rscript', '--vanilla', '/recorded/helper.R'],
                                   'exit_code': 0, 'elapsed_s': 0.1}]}
    fake_stages(monkeypatch, data)
    assert main(['inspect-suite', str(tmp_path)]) == 0
    inspected = json.loads(capsys.readouterr().out)
    assert set(data['suite']['warnings']) <= set(inspected['warnings'])
    summary = run_suite(tmp_path, tmp_path / 'result')
    receipt = json.loads((tmp_path / 'result/receipt.json').read_text())
    assert set(data['suite']['warnings']) <= set(summary['warnings'])
    assert receipt['warnings'] == summary['warnings']
    assert receipt['external_processes'] == data['suite']['commands']


def test_implementation_change_during_publication_is_rejected(tmp_path, monkeypatch):
    import defenseviz.suite_pipeline as pipeline
    data, _, _ = binding_fixture(tmp_path)
    fake_stages(monkeypatch, data)
    states = iter([{'source_sha256': {'analysis.py': 'before'}},
                   {'source_sha256': {'analysis.py': 'after'}}])
    monkeypatch.setattr(pipeline, 'software_state', lambda: next(states))
    with pytest.raises(ValueError, match='implementation changed'):
        run_suite(tmp_path, tmp_path / 'result')
    assert not (tmp_path / 'result').exists()
    assert not list(tmp_path.glob('.result.staging-*'))


def test_supplement_requires_native_result_binding_for_each_sample(tmp_path):
    data, doc, path = binding_fixture(tmp_path)
    data['input_roles'][str(tmp_path / 'native.tsv')][0]['native_result'] = False
    with pytest.raises(ValueError, match='native result artifacts'): apply_supplement(data, path)
    data['input_roles'][str(tmp_path / 'native.tsv')][0].update(native_result=True, sample_id='other')
    with pytest.raises(ValueError, match='no bound native result artifact'): apply_supplement(data, path)


@pytest.mark.parametrize('spoof', ['saved_decisions', 'native_system', 'native_provider', 'retained_member', 'origin'])
def test_reviews_cannot_spoof_native_decisions_or_systems(tmp_path, spoof):
    data, doc, path = binding_fixture(tmp_path)
    row = doc['annotations'][0]
    group = {'sample_id': 's', 'group_id': 'review-group', 'membership_kind': 'saved_context',
             'members': [{'gene_id': 'a', 'protein_sha256': row['protein_sha256']}],
             'source_refs': row['source_refs']}
    if spoof == 'saved_decisions': doc['saved_decisions'] = [{'state': 'retained'}]
    if spoof == 'native_system': group['membership_kind'] = 'native_system'
    if spoof == 'native_provider': group['provider'] = 'DefenseFinder'
    if spoof == 'retained_member': group['members'][0]['retained'] = True
    if spoof == 'origin': row['origin'] = 'native'
    doc['groups'] = [group]
    path.write_text(json.dumps(doc), encoding='utf-8')
    with pytest.raises(ValueError): apply_supplement(data, path)


@pytest.mark.parametrize('bad', ['missing', 'wrong_type', 'nonfinite'])
def test_invalid_analysis_cannot_be_published(tmp_path, monkeypatch, bad):
    data, _, _ = binding_fixture(tmp_path)
    fake_stages(monkeypatch, data)
    module = sys.modules['defenseviz.analysis.saved']
    original = module.analyze_saved
    def analyze(d):
        result = original(d)
        if bad == 'missing': del result['priority_ledger']
        if bad == 'wrong_type': result['type_tables'] = {}
        if bad == 'nonfinite': result['global_summary']['bad'] = float('nan')
        return result
    module.analyze_saved = analyze
    with pytest.raises(ValueError): run_suite(tmp_path, tmp_path / 'result')
    assert not (tmp_path / 'result').exists()
    assert not list(tmp_path.glob('.result.staging-*'))


def test_direct_native_cli_preserves_distinct_endpoints_and_editable_outputs(tmp_path, capsys):
    from defenseviz.suite_demo import write_suite_project
    from defenseviz.cli import main
    paths = write_suite_project(tmp_path / 'inputs')
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (tmp_path / 'inputs').rglob('*') if p.is_file()}
    output = tmp_path / 'published'
    code = main(['from-suite', str(paths['root']), '--manifest', str(paths['manifest']),
                 '--supplement', str(paths['supplement']), '--out', str(output)])
    assert code == 0, capsys.readouterr().err
    result = json.loads((output / 'analysis.json').read_text(encoding='utf-8'))
    summary = result['global_summary']
    assert result['manuscript_atlas']['schema_version'] == 'manuscript-atlas-v1'
    assert len(result['manuscript_atlas']['rows']) == len(result['samples'])
    assert result['manuscript_atlas']['provenance']['scientific_assignments_changed'] is False
    assert summary['native_union_CDS'] == 11
    assert summary['source_union_CDS'] == 13  # reference-only includes an unretained observation
    assert summary['integrated_CDS'] == 12  # saved retained orphan enters DNMB, raw-only hit does not
    assert summary['individual_effector_CDS'] == 2
    assert summary['complex_core_CDS'] == 3
    assert summary['effector_CDS'] == 5 and summary['priority_CDS'] == 6
    assert summary['motif_CDS'] == 1
    assert len(result['priority_ledger']) == len({(r['sample_id'], r['gene_id']) for r in result['priority_ledger']})
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == sha for p, sha in before.items())
    assert (output / 'figures/publication.pdf').stat().st_size > 1000
    svg = next((output / 'figures').glob('*.svg')).read_text(encoding='utf-8')
    assert '<text' in svg and '<image' not in svg
    receipt = json.loads((output / 'receipt.json').read_text(encoding='utf-8'))
    assert receipt['exit_code'] == 0 and receipt['detector_runs'] == 0
    assert receipt['options']['strict_inputs'] is False
    assert 'conflicts' in receipt and 'external_processes' in receipt

"""Display semantics for type subsets, tool-call pools, fraction scope and Letter exports."""
from copy import deepcopy
import hashlib
import json

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.plot.publication import render_publication, validate_publication_model, _config
from defenseviz.plot.publication.atlas import _type_cell, type_pages
from defenseviz.plot.publication.composite import composite_pages, LETTER
from defenseviz.plot.publication.display_scope import type_case_projection, model_digest
from defenseviz.plot.publication.drawing import Canvas, Fonts, STROKE_FLOOR_PT, stroke_artists
from tests.analysis.test_saved import base, call, cov, decision
from tests.plot.test_composite_publication import combined, external_result, envelope, five_type_rows


def ledger(path):
    return json.loads((path/'publication_layout.json').read_text(encoding='utf-8'))


def row(**values):
    data = dict(sample_id='s', rm_type='Type I', unique_CDS=3, pool_CDS=3, R=2, M=1, S=1, A=0,
                within_pool_role_states={k: 'complete' for k in 'RMSA'},
                context_applicability='saved_I_III_context', context_state='complete', incomplete_R_context=0,
                orphan_M=0, orphan_S=0, orphan_state={'M': 'complete', 'S': 'complete'})
    data.update(values)
    return data


def test_three_subset_states_and_structural_inapplicability_are_distinct():
    assessed_zero = row()
    assessed_n = row(incomplete_R_context=1, orphan_M=1, orphan_S=1)
    not_assessed = row(context_state='unassessed', incomplete_R_context=None, orphan_M=None, orphan_S=None,
                       orphan_state={'M': 'unassessed', 'S': 'unassessed'})
    for field, zero, some, unknown in (('R', '2', '2 (1)', '2 (—)'), ('M', '1', '1 (1)', '1 (—)'), ('S', '1', '1 (1)', '1 (—)')):
        assert _type_cell(assessed_zero, field)['display'] == zero
        assert _type_cell(assessed_n, field)['display'] == some
        assert _type_cell(not_assessed, field)['display'] == unknown
        assert {_type_cell(r, field)['subset_assessment'] for r in (assessed_zero, assessed_n)} == {'assessed'}
        assert _type_cell(not_assessed, field)['subset_assessment'] == 'not_assessed'
    # Type IV R context is structurally outside the I-III rule: always bare.
    type_iv = row(rm_type='Type IV', R=2, context_applicability='not_assessed_by_I_III_rule',
                  context_state='unassessed', incomplete_R_context=None)
    assert _type_cell(type_iv, 'R')['display'] == '2'
    assert _type_cell(type_iv, 'R')['subset_assessment'] == 'not_applicable'
    # A zero count never receives a subset mark; unavailable stays a bare dash.
    assert _type_cell(row(M=0, orphan_M=None, orphan_state={'M': 'unassessed'}), 'M')['display'] == '0'
    unavailable = row(within_pool_role_states={'R': 'unassessed', 'M': 'complete', 'S': 'complete', 'A': 'complete'})
    assert _type_cell(unavailable, 'R')['display'] == '—'
    inapplicable = row(rm_type='Type II', S=None, orphan_state={'M': 'complete', 'S': 'not_applicable'},
                       within_pool_role_states={'R': 'complete', 'M': 'complete', 'S': 'not_applicable', 'A': 'not_applicable'})
    assert _type_cell(inapplicable, 'S')['display'] == '—'


def _mark_states(data):
    rows = [r for r in data['type_tables'] if r['sample_id'] == data['samples'][0]['sample_id']]
    rows[0].update(M=1, orphan_M=None, orphan_state={'M': 'unassessed', 'S': 'complete'})
    rows[1].update(M=2, orphan_M=1, orphan_state={'M': 'complete', 'S': 'not_applicable'})
    return rows[0]['rm_type'], rows[1]['rm_type']


@pytest.mark.parametrize('kind', ['Type I', 'Type II', 'Type III', 'Total', 'Type IV'])
@pytest.mark.parametrize('missing', [True, False])
def test_legacy_missing_context_applicability_does_not_imply_assessed_zero(kind, missing):
    data = five_type_rows(combined(1))
    target = next(r for r in data['type_tables'] if r['rm_type'] == kind)
    target.update(context_state='unassessed', incomplete_R_context=None)
    if missing:
        target.pop('context_applicability', None)
    else:
        target['context_applicability'] = None
    # This is an accepted direct-render/legacy model, not an invalid row.
    validate_publication_model(data)
    cell = _type_cell(target, 'R')
    assert cell['display'] == str(target['R']) + ('' if kind == 'Type IV' else ' (—)')
    assert cell['subset_assessment'] == ('not_applicable' if kind == 'Type IV' else 'not_assessed')


@pytest.mark.parametrize('composition', ['separate', 'letter_composite'])
def test_both_renderers_mark_not_assessed_and_define_it_without_plus_or_question_marks(tmp_path, composition):
    data = five_type_rows(combined(3))
    unknown_type, assessed_type = _mark_states(data)
    sid = data['samples'][0]['sample_id']
    render_publication(data, tmp_path, page_composition=composition, sections=['atlas', 'type_benchmarks'], type_samples=[sid])
    report = ledger(tmp_path)
    cells = {(c['rm_type'], c['field']): c for p in report['pages'] for c in p['cells'] if c.get('kind') == 'type_cell'}
    assert cells[(unknown_type, 'M')]['display'] == '1 (—)'
    assert cells[(unknown_type, 'M')]['subset_assessment'] == 'not_assessed'
    assert cells[(assessed_type, 'M')]['display'] == '2 (1)'
    texts = [t['text'] for p in report['pages'] for t in p['texts'] if p['section'] in ('type_benchmark', 'atlas_composite')]
    assert '(—): not assessed' in texts
    numeric = [t['text'] for p in report['pages'] for t in p['texts'] if t['role'] == 'type_numeric']
    assert not any(mark in t for t in numeric for mark in ('+', '?', 'unassessed'))
    caption = (tmp_path/'publication_caption.md').read_text(encoding='utf-8')
    assert '(—) after a positive count means context or orphan status was not assessed' in caption
    assert 'Type IV R context lies outside the Type I–III context rule' in caption


def test_tool_columns_count_any_family_calls_in_the_assigned_pool(tmp_path):
    data = base('iv')
    data['calls'] = [call('iv', caller='PADLOC', family='PDC-S (pooled)')]
    data['coverage'] = cov('PADLOC')
    data['saved_decisions'] = [decision('iv', role='R', typ='Type IV')]
    result = analyze_saved(data)
    typed = next(r for r in result['type_tables'] if r['rm_type'] == 'Type IV')
    # The accepted endpoint is CDS inclusion: a non-R-M PADLOC call in the
    # DNMB-assigned pool is counted and must be labelled as such.
    assert typed['source_counts']['PADLOC'] == 1
    render_publication(result, tmp_path, sections=['type_benchmarks'], type_samples=['s'])
    page = ledger(tmp_path)['pages'][0]
    pad = next(c for c in page['cells'] if c.get('kind') == 'type_cell' and c['rm_type'] == 'Type IV' and c['field'] == 'PADLOC')
    assert (pad['value'], pad['display']) == (1, '1')
    header = [t for t in page['texts'] if t['text'].startswith('Tool calls in pool')]
    assert [t['text'] for t in header] == ['Tool calls in pool (any family)']
    assert 'Shared-pool CDSs' not in {t['text'] for t in page['texts']}
    # Even with one tool column the full label stays inside its header band.
    band = [c['bounds_mm'] for c in page['cells'] if c.get('kind') == 'type_cell'
            and c['rm_type'] == 'Type IV' and c['field'] in ('Type', 'Pool', 'PADLOC', 'DNMB')]
    left, right = min(b[0] for b in band), max(b[0]+b[2] for b in band)
    x, _, w, _ = header[0]['bounds_mm']
    assert left <= x and x+w <= right
    caption = (tmp_path/'publication_caption.md').read_text(encoding='utf-8')
    assert 'in any defense family' in caption and 'PDC-S or DMS_other' in caption
    assert 'rather than agreement on R-M type' in caption


@pytest.mark.parametrize('scope,label', [('retained_I_III', 'I–III CDSs'), ('assigned_I_IV', 'I–IV CDSs')])
def test_fraction_header_names_its_scope_without_changing_values(tmp_path, scope, label):
    data = combined(3)
    data['gene_ledger'] = [dict(sample_id=s['sample_id'], gene_id=f'g{i}', retained=True,
                                rm_assignments=[dict(rm_type='Type I', role='R' if i == 0 else 'M')])
                           for s in data['samples'] for i in range(2)]
    render_publication(data, tmp_path, rm_fraction_scope=scope, sections=['atlas'])
    page = ledger(tmp_path)['pages'][0]
    texts = [t['text'] for t in page['texts']]
    assert {'Assigned', label, '(R/total)'} <= set(texts)
    # The former singular heading line is gone; 7.3 pt 'CDS' sublabels remain.
    assert not [t for t in page['texts'] if t['text'] == 'CDS' and t['font_pt'] == pytest.approx(9.24)]
    if scope == 'retained_I_III':
        assert [(r['R_CDS'], r['total_CDS']) for r in page['rm_fraction_records']] == [(1, 2)]*3
    assert not page['validation']['errors']


def test_external_case_is_named_in_caption_and_cued_on_its_page(tmp_path):
    data = combined(3)
    external = external_result()
    render_publication(data, tmp_path/'figures', external_type_cases=[envelope(tmp_path, external)],
                       page_composition='letter_composite', sections=['atlas', 'type_benchmarks'],
                       type_samples=['isolate-000', 'independent/external'])
    caption = (tmp_path/'figures/publication_caption.md').read_text(encoding='utf-8')
    assert 'External type cases (tan headers): Independentus species DSM 999.' in caption
    assert 'not included in the atlas, overlap, family or workbook totals' in caption
    assert 'Primary cohort' not in caption and 'target strain' not in caption
    report = json.loads((tmp_path/'figures/publication_layout.json').read_text(encoding='utf-8'))
    texts = [t['text'] for p in report['pages'] for t in p['texts']]
    assert 'separately analyzed external case; not in cohort totals' in texts
    cells = [c for p in report['pages'] for c in p['cells'] if c.get('kind') == 'type_cell']
    assert {c['external_case'] for c in cells if c['sample_id'] == 'independent/external'} == {True}
    assert {c['external_case'] for c in cells if c['sample_id'] == 'isolate-000'} == {False}


def _render_pages(data, projection, selection, font='Arial', sections=('atlas', 'type_benchmarks')):
    fonts = Fonts(font)
    config = _config('Synthetic', dict(page_composition='letter_composite', sections=list(sections)))
    type_config = {**config, 'type_samples': selection}
    return composite_pages(data, projection, fonts, config, type_config,
                           lambda w, h, s: Canvas(w, h, fonts, s, 'synthetic'))


def test_letter_strokes_have_a_physical_floor_and_metrics_are_recorded(tmp_path):
    data = five_type_rows(combined(52))
    render_publication(data, tmp_path, page_composition='letter_composite')
    report = ledger(tmp_path)
    manifest = json.loads((tmp_path/'publication_manifest.json').read_text(encoding='utf-8'))
    assert manifest['letter_export_stroke_floor_pt'] == STROKE_FLOOR_PT
    assert [m['page'] for m in manifest['physical_metrics']] == list(range(1, report['page_count']+1))
    for page, metric in zip(report['pages'], manifest['physical_metrics']):
        assert metric == dict(page=page['page_id'], section=page['section'], **page['physical_metrics'])
        assert metric['min_text_pt'] == pytest.approx(min(t['font_pt'] for t in page['texts']))
        if 'export_transform' in page:
            assert page['size_mm'] == list(LETTER)
            assert metric['min_stroke_pt'] >= STROKE_FLOOR_PT-1e-9
            assert metric['export_scale'] == page['export_transform']['scale']
    caption = (tmp_path/'publication_caption.md').read_text(encoding='utf-8')
    assert f'kept at {STROKE_FLOOR_PT} pt or thicker' in caption
    assert 'small for print at Letter size' in caption
    first = manifest['physical_metrics'][0]
    assert f"smallest text {first['min_text_pt']:.2f} pt" in caption


def test_live_strokes_are_floored_on_every_letter_page():
    data = five_type_rows(combined(6))
    for sample in data['samples']:
        sample['highlight'] = True
    pages = _render_pages(data, data, None)
    lettered = [p for p in pages if 'export_transform' in p.ledger]
    assert lettered
    for page in lettered:
        widths = [w for _, w in stroke_artists(page.ax)]
        assert widths and min(widths) >= STROKE_FLOOR_PT-1e-9
        page.validate()
        assert not page.ledger['validation']['errors']
        page.fig.clear()


def test_compact_type_panel_preserves_lower_legend_outside_its_horizontal_span():
    from defenseviz.plot.publication.composite import _type_plans, _merged_layout, _resize_source
    fonts = Fonts('Arial')
    config = _config('Synthetic', dict(page_composition='letter_composite', type_samples=['isolate-000']))
    data = five_type_rows(combined(1))
    plan = _type_plans(data, fonts, config)[0]
    canvas = Canvas(447, 340, fonts, 'atlas', 'synthetic')
    canvas.text(10, 150, 'Existing left panel', 8)
    canvas.text(400, 320, 'Taxonomy legend', 8)
    existing = canvas.validate()['texts']
    bottom = max(t['bounds_mm'][1]+t['bounds_mm'][3] for t in existing)
    placed = _merged_layout(canvas, plan, fonts)
    assert placed['height'] > bottom
    _resize_source(canvas, placed['width'], placed['height'])
    assert not canvas.validate()['validation']['errors']
    canvas.fig.clear()


def test_continuation_type_pages_use_their_own_letter_area_and_say_continued():
    data = five_type_rows(combined(7))
    for sample in data['samples']:
        sample['highlight'] = True
    pages = _render_pages(data, data, None)
    first = pages[0]
    continuation = [p for p in pages if p.ledger['section'] == 'type_benchmark']
    assert first.ledger['section'] == 'atlas_composite' and first.ledger['type_continued'] is False
    assert continuation and all(p.ledger['type_continued'] for p in continuation)
    shown = [sid for p in [first]+continuation for sid in p.ledger['type_sample_ids']]
    assert shown == [s['sample_id'] for s in data['samples']]
    for page in continuation:
        assert page.ledger['export_transform']['scale'] >= 0.85
        assert page.ledger['size_mm'] == list(LETTER)
        titles = [t['text'] for t in page.measure() if t['text'].startswith('CDS inclusion')]
        assert titles == ['CDS inclusion by assigned R-M type (continued)']
        page.validate()
        assert not page.ledger['validation']['errors']
    for page in pages:
        page.fig.clear()


def test_long_labels_keep_a_generic_continuation_fallback():
    data = five_type_rows(combined(5, long_sample=True))
    for sample in data['samples']:
        sample['highlight'] = True
        sample['label'] = sample['label']+' with an unusually long strain designation for wrapping tests'
        sample.pop('display', None)
    pages = _render_pages(data, data, None, 'DejaVu Sans', sections=['type_benchmarks'])
    shown = [sid for p in pages for sid in p.ledger.get('type_sample_ids', [])]
    assert shown == [s['sample_id'] for s in data['samples']]
    assert pages[0].ledger['type_continued'] is False
    assert all(p.ledger['type_continued'] for p in pages[1:])
    for page in pages:
        page.validate()
        assert not page.ledger['validation']['errors']
        page.fig.clear()


def test_merged_panel_c_uses_measured_ink_by_horizontal_span():
    data = five_type_rows(combined(52))
    selection = [s['sample_id'] for s in data['samples'] if s['highlight']]
    pages = _render_pages(data, data, selection)
    first = pages[0]
    placement = first.ledger['type_placement']
    assert placement['method'] == 'measured_ink_by_horizontal_span'
    native = first.ledger['source_layout']
    texts = native['texts']
    title = next(t for t in texts if t['text'] == 'CDS inclusion by assigned R-M type')
    heading_texts = [t for t in texts if t['anchor_mm'][1] >= placement['title_anchor_mm']-1e-6]
    above = [t for t in texts if t['anchor_mm'][1] < placement['title_anchor_mm']-1e-6]
    right = max(t['bounds_mm'][0]+t['bounds_mm'][2] for t in heading_texts if t['bounds_mm'][1] < placement['tables_top_mm'])
    overlap = [t for t in above if t['bounds_mm'][0] < right and t['bounds_mm'][0]+t['bounds_mm'][2] > 3]
    assert title['bounds_mm'][1] >= max(t['bounds_mm'][1]+t['bounds_mm'][3] for t in overlap)+1.5
    assert not first.ledger['validation']['errors']
    for page in pages:
        page.fig.clear()


def _renamed(data, sample_id, label):
    for key in ('samples', 'cells', 'type_tables', 'sample_summary'):
        for item in data[key]:
            item['sample_id'] = sample_id
    data['samples'][0].update(label=label)
    data['samples'][0].pop('display', None)
    return data


def test_direct_render_rejects_a_known_assembly_alias_but_not_a_shared_label(tmp_path):
    from tests.plot.test_publication import fixture
    primary = _renamed(fixture(samples=1, families=1), 'GCF_000832905.1', 'Heyndrickxia coagulans DSM 1')
    alias = _renamed(fixture(samples=1, families=1), 'external_GCF_000832905.1_DSM1', 'B. coagulans DSM 1')
    with pytest.raises(ValueError, match='must not duplicate primary or other external specimens'):
        render_publication(primary, tmp_path/'figures', external_type_cases=[envelope(tmp_path, alias)],
                           sections=['type_benchmarks'])
    assert not (tmp_path/'figures').exists()
    # Display labels are not identity evidence: an identical label with a
    # different explicit assembly accession remains a separate case.
    other = _renamed(fixture(samples=1, families=1), 'external_GCF_000999999.1_X', 'Heyndrickxia coagulans DSM 1')
    _, _, bindings = type_case_projection(primary, [envelope(tmp_path, other)], None, validate_publication_model, tmp_path)
    assert bindings[0]['identity_check']['status'] == 'no_shared_assembly_accession'
    assert bindings[0]['identity_check']['external_assembly_accessions'] == ['GCF_000999999']


def test_missing_assembly_metadata_is_recorded_as_unchecked(tmp_path):
    from tests.plot.test_publication import fixture
    primary = fixture(samples=1, families=1)
    external = external_result()
    _, _, bindings = type_case_projection(primary, [envelope(tmp_path, external)], None, validate_publication_model, tmp_path)
    check = bindings[0]['identity_check']
    assert check['status'] == 'partly_unchecked_missing_assembly_metadata'
    assert set(check['samples_without_assembly_metadata']) == {primary['samples'][0]['sample_id'], 'independent/external'}
    assert check['external_assembly_accessions'] == []

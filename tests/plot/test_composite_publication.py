"""Synthetic contracts for optional page composition, scope and external cases."""
from copy import deepcopy
import hashlib
import json

import pytest

from defenseviz.plot.publication import render_publication, validate_publication_model, _config
from defenseviz.plot.publication.display_scope import fraction_projection, type_case_projection, model_digest
from defenseviz.plot.publication.composite import composite_pages
from defenseviz.plot.publication.drawing import Canvas, Fonts
from tests.plot.test_publication import fixture
from tests.plot.test_manuscript_atlas import synthetic_result


def combined(samples=3, **kwargs):
    data = fixture(samples=samples, families=3)
    template = synthetic_result(samples=samples, **kwargs)
    old = [s['sample_id'] for s in data['samples']]
    data.update(template)
    mapping = dict(zip(old, [s['sample_id'] for s in data['samples']]))
    for key in ('cells', 'type_tables', 'sample_summary'):
        for row in data[key]:
            row['sample_id'] = mapping[row['sample_id']]
    for row in data['type_tables']:
        if row['rm_type'] != 'Total': row['rm_type'] = 'Type '+row['rm_type']
        row['pool_CDS'] = row['unique_CDS']
        row['pool_status'] = 'complete'
    return data


def envelope(tmp_path, result):
    path = tmp_path/'external.json'
    content = json.dumps(result, ensure_ascii=False).encode('utf-8')
    path.write_bytes(content)
    return dict(result=result, source_refs=[dict(path=path.name, sha256=hashlib.sha256(content).hexdigest())], result_sha256=model_digest(result))


def external_result():
    data = fixture(samples=1, families=1)
    for key in ('samples', 'cells', 'type_tables', 'sample_summary'):
        for row in data[key]:
            row['sample_id'] = 'independent/external'
    data['samples'][0].update(label='Independentus species DSM 999', display={
        'label_runs': [{'text':'Independentus species', 'italic':True}, {'text':'DSM 999', 'italic':False}]})
    return data


def test_external_projection_requires_actual_matching_file_and_never_merges_primary(tmp_path):
    data = fixture()
    before = deepcopy(data)
    external = external_result()
    bound = envelope(tmp_path, external)
    projection, selected, refs = type_case_projection(data, [bound], None, validate_publication_model, tmp_path)
    assert data == before
    assert selected == ['sample/1', 'independent/external']
    assert set(projection) == {'samples', 'sources', 'type_tables'}
    assert len(projection['type_tables']) == len(data['type_tables'])+len(external['type_tables'])
    assert refs[0]['included_in_primary_aggregates'] is False
    assert refs[0]['result_sha256'] == model_digest(external)
    bad = deepcopy(bound)
    bad['result']['type_tables'][0]['R'] = 2
    with pytest.raises(ValueError, match='result SHA256'):
        type_case_projection(data, [bad], None, validate_publication_model, tmp_path)
    bad.pop('result_sha256')
    with pytest.raises(ValueError, match='bind the separately analyzed'):
        type_case_projection(data, [bad], None, validate_publication_model, tmp_path)
    (tmp_path/'external.json').write_bytes(b'{}')
    with pytest.raises(ValueError, match='source SHA256 mismatch'):
        type_case_projection(data, [bound], None, validate_publication_model, tmp_path)


def test_external_sample_identity_collision_rejected(tmp_path):
    data = fixture(samples=1)
    with pytest.raises(ValueError, match='overlaps'):
        type_case_projection(data, [envelope(tmp_path, data)], None, validate_publication_model, tmp_path)


def test_strict_fraction_uses_cds_identities_and_retained_flags_not_aggregate_numbers():
    data = fixture(samples=1, families=1)
    sid = data['samples'][0]['sample_id']
    data['sample_summary'][0].update(current_retained_counts_status='complete', assigned_R_CDS=999, assigned_RM_CDS=1000)
    rows = [('a', True, 'Type I', 'R'), ('b', True, 'Type II', 'M'), ('c', False, 'Type IV', 'R'),
            ('d', True, 'Type IV', 'R'), ('e', False, 'Type I', 'R')]
    data['gene_ledger'] = [dict(sample_id=sid, gene_id=gene, retained=kept, rm_assignments=[dict(rm_type=kind, role=role)]) for gene,kept,kind,role in rows]
    before = deepcopy(data)
    scoped = fraction_projection(data, 'retained_I_III')[sid]
    assert (scoped['R_CDS'], scoped['total_CDS']) == (1, 2)
    assert scoped['member_keys'] == [[sid, 'a'], [sid, 'b']]
    assert scoped['R_member_keys'] == [[sid, 'a']]
    assert data == before
    assert fraction_projection(data, 'assigned_I_IV') is None
    data['gene_ledger'].append(deepcopy(data['gene_ledger'][0]))
    with pytest.raises(ValueError, match='Duplicate CDS'):
        fraction_projection(data, 'retained_I_III')


@pytest.mark.parametrize('samples,kwargs', [(1, {'categories':[]}), (6, {}), (52, {}), (80, {'long_sample':True})])
def test_letter_composition_preserves_rows_type_values_and_physical_bounds(samples, kwargs):
    data = combined(samples, **kwargs)
    before = deepcopy(data)
    config = _config('Synthetic composition', dict(page_composition='letter_composite', sections=['atlas','type_benchmarks']))
    fonts = Fonts()
    def factory(w,h,section):
        return Canvas(w,h,fonts,section,'synthetic')
    pages = composite_pages(data, data, fonts, config, config, factory)
    assert data == before
    strain_ids = [cell['sample_id'] for page in pages for cell in page.ledger['cells'] if cell['kind']=='strain']
    assert set(strain_ids) == {s['sample_id'] for s in data['samples']}
    from collections import Counter
    assert len(set(Counter(strain_ids).values())) == 1
    assert all(page.ledger['export_transform']['scale'] >= .6 for page in pages)
    assert all(page.ledger['size_mm'] == [279.4, 215.9] for page in pages)
    cells = [cell for page in pages for cell in page.ledger['cells'] if cell['kind']=='type_cell' and cell['field']=='DNMB']
    selected = [s['sample_id'] for s in data['samples'] if s['highlight']]
    expected = {(row['sample_id'],row['rm_type']):row['unique_CDS'] for row in data['type_tables'] if row['sample_id'] in selected}
    assert {(c['sample_id'],c['rm_type']):c['value'] for c in cells} == expected
    for page in pages:
        page.validate()
        scale = page.ledger['export_transform']['scale']
        assert len(page.ledger['bars']) == len(page.ledger['source_layout']['bars'])
        for bar, native in zip(page.ledger['bars'], page.ledger['source_layout']['bars']):
            for field in ('track_width_mm', 'filled_width_mm'):
                if field in native:
                    assert bar[field] == pytest.approx(native[field] * scale)
            if native.get('track_width_mm') and 'filled_width_mm' in native:
                assert bar['filled_width_mm'] / bar['track_width_mm'] == pytest.approx(
                    native['filled_width_mm'] / native['track_width_mm'])
        for glyph in page.ledger['glyphs']:
            if glyph['outer_radius_mm']:
                assert glyph['inner_radius_mm']**2/glyph['outer_radius_mm']**2 == pytest.approx(glyph['inner_area_fraction'])
        assert not page.ledger['validation']['errors']
        page.fig.clear()
    if samples <= 52 and fonts.base == 'Arial':
        assert pages[0].ledger['section'] == 'atlas_composite'


def test_letter_svg_is_live_text_external_rows_scoped_and_long_unicode_styles_preserved(tmp_path):
    data = combined(3)
    external = external_result()
    external['samples'][0].update(label='Independentus äß α species isolate designation with literal $x_1$ and ἀλφα', display={'label_runs':[
        {'text':'Independentus äß α species', 'italic':True},
        {'text':'isolate designation with literal $x_1$ and ἀλφα', 'italic':False}]})
    bound = envelope(tmp_path, external)
    render_publication(data, tmp_path/'figures', external_type_cases=[bound], page_composition='letter_composite',
                       sections=['atlas','type_benchmarks'], type_samples=['isolate-000','isolate-002','independent/external'])
    report = json.loads((tmp_path/'figures/publication_layout.json').read_text())
    assert report['input_model_sha256'] == model_digest(data)
    assert report['external_type_bindings'][0]['included_in_primary_aggregates'] is False
    texts = [t for page in report['pages'] for t in page['texts'] if t['role']=='type_sample_label' and t['source_ref']['sample_id']=='independent/external']
    assert any(t['italic'] for t in texts) and any(not t['italic'] for t in texts)
    assert ' '.join(t['text'] for t in texts) == external['samples'][0]['label']
    svg = (tmp_path/'figures/publication_page001.svg').read_text()
    assert '<text' in svg and 'data:image' not in svg
    assert all(not p['validation']['errors'] for p in report['pages'])


def five_type_rows(data):
    """Synthetic I–IV+Total rows, with explicit structural S applicability."""
    rows=[]
    for sample in data['samples']:
        seed=next(r for r in data['type_tables'] if r['sample_id']==sample['sample_id'])
        for kind,n,r,m,s,a in [('Type I',3,1,1,1,0),('Type II',2,1,1,None,0),
                             ('Type III',2,1,1,None,0),('Type IV',2,1,0,0,1),('Total',9,4,3,1,1)]:
            row=deepcopy(seed)
            row.update(rm_type=kind,unique_CDS=n,pool_CDS=n,pool_status='complete',R=r,M=m,S=s,A=a,
                role_states={role:'not_applicable' if role=='S' and s is None else 'complete' for role in 'RMSA'},
                source_counts={source:min(n,2) for source in row['source_counts']},
                incomplete_R_context=None if kind=='Type IV' else 1 if kind in ('Type I','Type III') else 2 if kind=='Total' else 0,
                context_applicability='not_applicable' if kind=='Type IV' else 'applicable',
                orphan_M=1 if kind in ('Type III','Total') else 0,orphan_S=0)
            rows.append(row)
    data['type_tables']=rows
    return data


@pytest.mark.parametrize('font_family',['Arial','DejaVu Sans'])
def test_three_five_row_type_cases_fit_52row_letter_frame_with_styled_long_case_label(tmp_path,font_family):
    data=five_type_rows(combined(52))
    independent=five_type_rows(external_result())
    taxon='Independentus thermosyntheticus subspecies exemplaris'
    strain='DSM 999 collection accession SYNTHETIC0001'
    independent['samples'][0].update(label=taxon+' '+strain,
        display={'label_runs':[dict(text=taxon,italic=True),dict(text=strain,italic=False)]})
    projection,selection,_=type_case_projection(data,[envelope(tmp_path,independent)],None,validate_publication_model,tmp_path)
    fonts=Fonts(font_family)
    config=_config('Synthetic five-row type cases',dict(page_composition='letter_composite',sections=['atlas','type_benchmarks']))
    type_config={**config,'type_samples':selection}
    pages=composite_pages(data,projection,fonts,config,type_config,
        lambda w,h,s:Canvas(w,h,fonts,s,'synthetic-five-row'))
    if fonts.base=='Arial':assert len(pages)==1 and pages[0].ledger['section']=='atlas_composite'
    for page in pages:page.validate()
    cells=[cell for page in pages for cell in page.ledger['cells'] if cell['kind']=='type_cell']
    assert len({cell['sample_id'] for cell in cells})==3
    assert {cell['rm_type'] for cell in cells}=={'Type I','Type II','Type III','Type IV','Total'}
    assert len([cell for cell in cells if cell['field']=='S' and cell['state']=='not_applicable'])==6
    assert all(cell['display']=='—' for cell in cells if cell['field']=='S' and cell['state']=='not_applicable')
    assert all(page.ledger['export_transform']['scale']>=.6 for page in pages)
    assert len({cell['numeric_font_pt'] for cell in cells})==1
    assert any(cell['display']=='1 (1)' for cell in cells if cell['field']=='M')
    for page in pages:page.fig.clear()

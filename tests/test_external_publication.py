"""External cohorts keep typography and assessment scope explicit."""
import copy
import json

import pytest

from defenseviz.plot.publication import render_publication
from defenseviz.plot.publication.atlas import _type_case_label
from defenseviz.plot.publication.drawing import Canvas, Fonts
from matplotlib.text import Text
from tests.plot.test_publication import fixture


def _pages(path):
    return json.loads((path / 'publication_layout.json').read_text())['pages']


def test_long_taxon_and_strain_have_visible_gap_at_original_size(tmp_path):
    data=fixture(samples=1,families=1)
    taxon='Exemplothermobacterium hypotheticalis subsp. synthetica'
    strain='EX-001'
    data['samples'][0].update(label=taxon+' '+strain,
        display={'taxon_label':taxon,'strain_label':strain})
    render_publication(data,tmp_path,sections=['type_benchmarks'])
    page=_pages(tmp_path)[0]
    parts=[t for t in page['texts'] if t['role']=='type_sample_label']
    assert [t['text'] for t in parts]==[taxon,strain]
    a,b=(t['bounds_mm'] for t in parts)
    assert b[0]-(a[0]+a[2])>=1.1
    assert all(t['font_pt']==9 for t in parts)
    assert parts[0]['italic'] and not parts[1]['italic']
    assert not page['validation']['errors']


def test_comparison_headers_reserve_visible_space(tmp_path):
    render_publication(fixture(),tmp_path,sections=['comparison'])
    for page in _pages(tmp_path):
        for genome in (t for t in page['texts'] if t['text']=='Number of\ngenomes'):
            dnmb=next(t for t in page['texts'] if t['text']=='DNMB' and t['anchor_mm'][1]==genome['anchor_mm'][1])
            a,b=genome['bounds_mm'],dnmb['bounds_mm']
            assert b[0]-(a[0]+a[2])>=1.5
        assert not page['validation']['errors']


@pytest.mark.parametrize('effector,rm,observed',[
    ('complete','complete',False),
    ('partial','complete',True),
    ('complete','unassessed',True),
    ('partial','partial',True),
])
def test_priority_heading_qualifies_evidence_separately_from_source_coverage(tmp_path,effector,rm,observed):
    data=fixture(samples=1,families=1)
    data['global_summary']['metric_status']={'effector':effector}
    data['availability']={'saved_decisions':'partial',
        'current_saved_decisions_by_sample':[{'sample_id':'sample/0','status':rm}]}
    before=copy.deepcopy(data)
    render_publication(data,tmp_path,sections=['comparison'])
    assert data==before
    page=_pages(tmp_path)[0]
    selected=next(r for r in page['overlaps'] if r['endpoint']=='Prioritized CDSs')
    assert selected['union_count']==data['priority_overlap_regions']['union_count']
    assert selected['set_sizes']=={
        source:sum(row['count'] for row in data['priority_overlap_regions']['regions'] if row['mask'][i]=='1')
        for i,source in enumerate(data['priority_overlap_regions']['source_order'])}
    assert selected['display_heading']==('Observed priorities' if observed else 'Prioritized CDSs')
    assert page['selection_assessment']['complete'] is (not observed)
    assert page['selection_assessment']['rm_origin']=='availability.current_saved_decisions_by_sample'
    assert data['priority_overlap_regions']['coverage_status']=='complete'
    assert not page['validation']['errors']


def test_unknown_selection_scope_is_qualified_in_intersection_fallback(tmp_path):
    data=fixture(samples=1,families=1,sources=5)
    render_publication(data,tmp_path,sections=['comparison'],overlap_mode='intersection')
    page=_pages(tmp_path)[0]
    assert any(t['text']=='Observed priorities' for t in page['texts'])
    assert page['selection_assessment']['complete'] is False
    assert sum(r['priority'] for p in _pages(tmp_path) for r in p.get('intersections',[]))==data['priority_overlap_regions']['union_count']
    assert all(not p['validation']['errors'] for p in _pages(tmp_path))


def test_type_label_spacing_uses_actual_renderer_advance():
    # The live Agg renderer can advance small hinted text farther than the
    # unhinted outline used by the publication ink ledger.
    c=Canvas(150,40,Fonts(),'type_benchmark',1)
    taxon='Exemplothermobacterium hypotheticalis subsp. synthetica'
    strain='EX-001'
    sample={'sample_id':'synthetic','label':taxon+' '+strain,
            'display':{'taxon_label':taxon,'strain_label':strain}}
    _type_case_label(c,4,10,142,10,sample,sample['label'])
    c.fig.canvas.draw()
    renderer=c.fig.canvas.get_renderer()
    assert len(c.text_artists)==2
    a,b=[Text.get_window_extent(artist,renderer) for artist,_ in c.text_artists]
    assert (b.x0-a.x1)*25.4/c.fig.dpi>=1.1
    c.fig.clear()

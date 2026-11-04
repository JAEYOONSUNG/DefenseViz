"""Visible-ink row sizing without smaller fonts or omitted observations."""
from copy import deepcopy
import json
import xml.etree.ElementTree as ET

import pytest

from defenseviz.plot.publication import render_publication


def dataset(families=49, categories=8):
    sources=['PADLOC','DefenseFinder','DefensePredictor','REBASE']
    def row(i,prefix):
        return dict(id=f'{prefix}/{i}',label=f'{prefix} (gyp) {i}',category='Nuclease',
                    N=1192,G=52,associated_CDS=703,effector_CDS=293,
                    metric_status=dict(inventory='partial',associated='partial',effector='partial'),
                    source_counts={s:879 for s in sources},
                    source_count_status={s:'partial' for s in sources})
    overlap=dict(source_order=sources,regions=[dict(mask=f'{i:04b}',count=1) for i in range(1,16)],
                 union_count=15,coverage_status='complete')
    return dict(samples=[dict(sample_id=f'isolate/{i}',label=f'Isolate {i}') for i in range(52)],
                sources=[dict(id=s,kind='reference_annotation' if s=='REBASE' else 'native_final_member') for s in sources],
                inventory_rows=dict(categories=[row(i,'Category') for i in range(categories)],
                                    families=[row(i,'Family') for i in range(families)]),
                overlap_regions=overlap,priority_overlap_regions=deepcopy(overlap),
                cells=[],type_tables=[],sample_summary=[])


def rendered(data,path,**config):
    before=deepcopy(data)
    render_publication(data,path,sections=['comparison'],**config)
    assert data==before
    # The ledger is UTF-8 even when the platform's default text codec is not.
    result=json.loads((path/'publication_layout.json').read_text(encoding='utf-8'))
    assert all(not p['validation']['errors'] for p in result['pages'])
    return result


def test_fifty_seven_rows_fit_letter_at_unchanged_font_with_descenders_and_partial_commas(tmp_path):
    data=dataset();data['inventory_rows']['families'][3]['category']='Unresolved'
    result=rendered(data,tmp_path)
    assert result['page_count']==1
    page=result['pages'][0]
    assert page['size_mm']==[215.9,279.4]
    rows=page['cells']
    assert len(rows)==57
    assert max(r['bounds_mm'][1]+r['bounds_mm'][3] for r in rows)<264
    assert all(2.6<r['bounds_mm'][3]<2.75 for r in rows)
    numbers=[t for t in page['texts'] if t['role']=='numeric']
    assert all(t['font_pt']==7.8 and t['alignment']=='right' and t['bounds_kind']=='glyph_ink' for t in numbers)
    assert any(t['text']=='1,192+' for t in numbers)
    assert all(b['bounds_mm'][2:]==[12.,1.4] for b in page['bars'])
    for r in rows:
        matching=[t for t in page['texts'] if (t.get('source_ref') or {}).get('row_id')==r['row_id']]
        assert matching
        assert {round(t['anchor_mm'][1],9) for t in matching}=={round(r['text_anchor_y_mm'],9)}
        top=r['bounds_mm'][1];bottom=top+r['bounds_mm'][3]
        assert all(t['bounds_mm'][1]>=top+.059 for t in matching)
        assert all(t['bounds_mm'][1]+t['bounds_mm'][3]<=bottom-.059 for t in matching)
    headers=[t for t in page['texts'] if t['role'] in ('total_header','effector_header')]
    for total,effector in zip(headers[::2],headers[1::2]):
        assert effector['bounds_mm'][0]-(total['bounds_mm'][0]+total['bounds_mm'][2])>=2.5
    assert any('REB: REBASE' in t['text'] and 'DF: DefenseFinder' in t['text'] for t in page['texts'] if t['role']=='source_abbreviations')
    assert '#eef1f3' in (tmp_path/'publication_page001.svg').read_text()
    # All numbers keep independent, common right-aligned column anchors.
    assert all(r['numeric_anchors']==rows[0]['numeric_anchors'] for r in rows)


def test_multiline_fallback_unicode_and_literal_svg_remain_editable_and_fit(tmp_path):
    data=dataset(4,1)
    label='Óñgj αβ $x_{y}$ & literal\nMultiple lines with descenders gypq and further qualification'
    data['inventory_rows']['families'][0]['label']=label
    result=rendered(data,tmp_path,font_family='Definitely not installed')
    assert result['font_resolved']!='Definitely not installed'
    assert result['font_files'] and all(f['sha256'] for f in result['font_files'])
    rows=[r for p in result['pages'] for r in p['cells']]
    long=next(r for r in rows if r['row_id']=='Family/0')
    assert len(long['label_lines'])>2 and long['bounds_mm'][3]>6
    assert all(t['font_pt']==7.8 for p in result['pages'] for t in p['texts'] if t['role']=='numeric')
    xml=ET.parse(tmp_path/'publication_page001.svg').getroot()
    assert not xml.findall('.//{http://www.w3.org/2000/svg}image')
    strings=[''.join(t.itertext()) for t in xml.findall('.//{http://www.w3.org/2000/svg}text')]
    assert all(line in strings for line in long['label_lines'])
    assert '$x_{y}$' in ' '.join(strings)
    assert 'Óñgj' in ' '.join(strings)


def test_arbitrary_long_families_paginate_without_losing_rows_or_smaller_fonts(tmp_path):
    data=dataset(90,3)
    for r in data['inventory_rows']['families'][::5]:
        r['label']='An unknown family with a long explanatory label containing gypq descenders'
    result=rendered(data,tmp_path)
    assert result['page_count']>1
    rows=[r for p in result['pages'] for r in p['cells']]
    expected={(table,r['id']) for table,items in data['inventory_rows'].items() for r in items}
    assert {(r['table'],r['row_id']) for r in rows}==expected
    assert {(r['table'],r['row_id'],source) for r in rows for source in r['source_order']}=={(table,rid,source) for table,rid in expected for source in ['DefenseFinder','PADLOC','DefensePredictor']}
    assert all(t['font_pt']==7.8 for p in result['pages'] for t in p['texts'] if t['role']=='numeric')
    assert all(r['source_order']==[s for s in ['DefenseFinder','PADLOC','DefensePredictor'] if s in r['source_order']] for r in rows)


def test_unassessed_rows_and_explicit_larger_fonts_are_not_compressed_to_fit(tmp_path):
    data=dataset(49,8)
    for r in data['inventory_rows']['families']:
        r['metric_status']={key:'unassessed' for key in r['metric_status']}
        r['N']=r['associated_CDS']=r['effector_CDS']=0
        r['source_counts']={s:0 for s in r['source_counts']}
        r['source_count_status']={s:'unassessed' for s in r['source_counts']}
    result=rendered(data,tmp_path,body_font_size=9)
    assert result['page_count']>1
    rows=[r for p in result['pages'] for r in p['cells']]
    assert len({(r['table'],r['row_id']) for r in rows})==57
    numbers=[t for p in result['pages'] for t in p['texts'] if t['role']=='numeric']
    assert all(t['font_pt']==9 for t in numbers)
    assert any(t['text']=='NA' for t in numbers)


def test_single_unrenderable_label_errors_instead_of_clipping_or_omission(tmp_path):
    data=dataset(1,0)
    data['inventory_rows']['families'][0]['label']='An extremely long unclassified family description '*300
    with pytest.raises(ValueError,match='cannot fit a single page'):
        render_publication(data,tmp_path,sections=['comparison'])
    assert not list(tmp_path.iterdir())

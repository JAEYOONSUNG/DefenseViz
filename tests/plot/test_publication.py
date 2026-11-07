import copy
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from defenseviz.plot.publication import render_publication, validate_publication_model
from defenseviz.plot.publication.drawing import Fonts
from defenseviz.plot.publication.overlap import venn_fits


def fixture(samples=2, families=5, sources=4):
    names=["PADLOC","DefenseFinder","DefensePredictor","REBASE"][:sources]
    names += [f"Novel source {i}" for i in range(max(0,sources-4))]
    ss=[{"id":s,"label":s,"kind":"reference_annotation" if s=="REBASE" else "native_final_member","coverage_status":"complete"} for s in names]
    sample_rows=[{"sample_id":f"sample/{i}","label":f"Synthetic species isolate {i}","highlight":i==1} for i in range(samples)]
    metrics={"inventory":"observed","associated":"observed","effector":"observed","motif":"unassessed"}
    fs=[];cells=[]
    for i in range(families):
        row={"id":f"family/{i}","label":f"Synthetic family {i}","category":"Modification" if i%2==0 else "Nuclease",
             "category_tags":["Modification","Nuclease"] if i==0 else [],"N":10,"G":samples,
             "source_counts":{s:8-j%3 for j,s in enumerate(names)},"source_count_status":{s:"complete" for s in names},
             "associated_CDS":6,"effector_CDS":3,"motif_CDS":0,"metric_status":metrics.copy(),"counts_status":"observed"}
        fs.append(row)
        for sample in sample_rows:
            cells.append(dict(copy.deepcopy(row),sample_id=sample["sample_id"],family_id=row["id"]))
    masks=[format(i,f"0{len(names)}b") for i in range(1,2**len(names))] if names else []
    masks=masks[:70]
    before={"source_order":names,"regions":[{"mask":m,"count":1} for m in masks],"union_count":len(masks),"coverage_status":"complete"}
    after={"source_order":names,"regions":[{"mask":m,"count":1} for m in masks if m[-1:]=="1"],"union_count":sum(m[-1:]=="1" for m in masks),"coverage_status":"complete"}
    type_rows=[]
    for sample in sample_rows:
        for kind,n,r,m,s,a in (("I",3,1,1,1,0),("IV",2,1,0,0,1),("Total",5,2,1,1,1)):
            type_rows.append({"sample_id":sample["sample_id"],"rm_type":kind,"unique_CDS":n,
                "counts_status":"complete",
                "source_counts":{source:min(n,2) for source in names},"source_count_status":{source:"complete" for source in names},
                "R":r,"M":m,"S":s,"A":a,"role_states":{role:"observed" for role in "RMSA"},
                "context_applicability":"not_applicable" if kind=="IV" else "applicable",
                "context_state":"unassessed" if kind=="IV" else "complete",
                "incomplete_R_context":None if kind=="IV" else 1,"orphan_M":0,"orphan_S":0,"orphan_state":{"M":"observed","S":"observed"}})
    return {"schema_version":1,"samples":sample_rows,"sources":ss,"inventory_rows":{"families":fs,"categories":copy.deepcopy(fs[:2])},
            "cells":cells,"type_tables":type_rows,"overlap_regions":before,"priority_overlap_regions":after,
            "sample_summary":[{"sample_id":s["sample_id"],"integrated_CDS":15,
                "source_counts":{name:8 for name in names},"source_count_status":{name:"complete" for name in names},"metric_status":{"inventory":"observed"}} for s in sample_rows],
            "availability":{},"global_summary":{},"endpoint_definitions":{}}


def ledger(path):
    return json.loads((path/"publication_layout.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize('style',['manuscript','compact'])
def test_public_assignment_labels_preserve_saved_result_values(tmp_path,style):
    data=fixture(samples=1,families=1)
    before=copy.deepcopy(data)
    render_publication(data,tmp_path,sections=['atlas'],atlas_style=style)
    assert data==before
    labels=[t['text'] for p in ledger(tmp_path)['pages'] for t in p['texts']]
    assert any('assigned' in text.lower() for text in labels)
    assert all('retained' not in text.lower() for text in labels)


@pytest.mark.parametrize('numerator,denominator,status,expected',[
    (4,9,'complete','4/9'),(0,0,'complete','0/0'),
    (4,9,'partial','4/9'),(None,9,'unassessed','—')])
def test_compact_assignment_fraction_prefers_explicit_assigned_scope(tmp_path,numerator,denominator,status,expected):
    data=fixture(samples=1,families=1)
    data['sample_summary'][0].update(retained_R_CDS=2,retained_RM_CDS=6,retained_ratio_status='complete',
        assigned_R_CDS=numerator,assigned_RM_CDS=denominator,assigned_counts_status=status,
        rm_source_counts=data['sample_summary'][0]['source_counts'])
    before=copy.deepcopy(data)
    render_publication(data,tmp_path,sections=['atlas'],atlas_style='compact')
    assert data==before
    values=[t['text'] for p in ledger(tmp_path)['pages'] for t in p['texts']
            if (t.get('source_ref') or {}).get('field')=='retained_R_fraction']
    assert values==[expected]


def test_comparison_right_alignment_split_columns_and_source_order(tmp_path):
    data=fixture();before=copy.deepcopy(data)
    render_publication(data,tmp_path,sections=["comparison"])
    assert data==before
    g=ledger(tmp_path)
    rows=[r for p in g["pages"] for r in p["cells"] if r["kind"]=="comparison_row"]
    assert len(rows)==7
    assert all(r["source_order"]==["DefenseFinder","PADLOC","DefensePredictor"] for r in rows)
    for page in g["pages"]:
        nums=[t for t in page["texts"] if t["role"]=="numeric"]
        assert all(t["alignment"]=="right" and t["font_pt"]==7.8 for t in nums)
        for t in nums:
            assert "(" not in t["text"]
        assert not page["validation"]["errors"]
        for bar in page["bars"]:
            assert bar["bounds_mm"][2:]==[12.0,1.4]
            assert bar["filled_width_mm"]==pytest.approx(12*bar["value"]/bar["denominator"])
    assert g["pages"][0]["overlaps"][0]["source_order"]==data["overlap_regions"]["source_order"]


def test_nested_glyph_area_dual_bands_and_pink_highlight(tmp_path):
    render_publication(fixture(),tmp_path,sections=["atlas"])
    g=ledger(tmp_path);page=g["pages"][0]
    assert page["bands"][0]["secondary"]==["Nuclease"]
    for glyph in page["glyphs"]:
        assert glyph["inner_radius_mm"]**2/glyph["outer_radius_mm"]**2==pytest.approx(.3)
    assert any(c.get("highlight") for c in page["cells"])
    raw=(tmp_path/"publication_page001.svg").read_text(encoding="utf-8")
    assert "#fbe0eb" in raw and "#7953a2" in raw


def test_type_parentheses_equal_font_and_no_activity_claim(tmp_path):
    render_publication(fixture(),tmp_path,sections=["type_benchmarks"])
    page=ledger(tmp_path)["pages"][0]
    assert page['sample_ids']==['sample/1']
    assert page['sample_selection']=='highlighted_default'
    numeric=[t for t in page["texts"] if t["role"]=="type_numeric"]
    assert {t["font_pt"] for t in numeric}=={9.4}
    assert all(t['alignment']=='right' for t in numeric)
    assert any(t["text"]=="1 (1)" for t in numeric)
    assert all("Active" not in t["text"] for t in page["texts"])
    assert all(c["included_subset"] is None for c in page["cells"] if c["field"]=="R" and c["rm_type"]=="IV")


def test_unknown_evidence_is_dash_not_zero(tmp_path):
    data=fixture(samples=1,families=1)
    for row in data["inventory_rows"]["families"]+data["inventory_rows"]["categories"]+data["cells"]:
        row["metric_status"]["associated"]="unassessed"
        row["metric_status"]["effector"]="unassessed"
        row["associated_CDS"]=row["effector_CDS"]=0
        row["source_count_status"]["DefenseFinder"]="missing"
        row["source_counts"]["DefenseFinder"]=0
    render_publication(data,tmp_path,sections=["comparison"])
    nums=[t for p in ledger(tmp_path)["pages"] for t in p["texts"] if t["role"]=="numeric"]
    for t in nums:
        if t["source_ref"]["field"] in ("associated","effector","DefenseFinder"):
            assert t["text"]=="—"


def test_partial_coverage_is_audited_and_missing_zero_family_retained(tmp_path):
    data=fixture(samples=1,families=1)
    row=data["inventory_rows"]["families"][0]
    row.update(N=0,associated_CDS=0,effector_CDS=0,source_counts={s:0 for s in row["source_counts"]})
    row["metric_status"]={k:"unassessed" for k in row["metric_status"]}
    data["inventory_rows"]["categories"][0]["source_count_status"]["DefenseFinder"]="partial"
    render_publication(data,tmp_path,sections=["comparison"])
    pages=ledger(tmp_path)["pages"]
    assert any(c.get("table")=="families" for p in pages for c in p["cells"])
    partial=[t for p in pages for t in p["texts"] if t["role"]=="numeric"
             and t["source_ref"]["field"]=="DefenseFinder" and t["source_ref"]["status"]=="partial"]
    assert partial and all(t["text"]==str(t["source_ref"]["value"]) for t in partial)
    missing=[t for p in pages for t in p["texts"] if t["role"]=="numeric"
             and t["source_ref"]["table"]=="families" and t["source_ref"]["field"]=="DNMB"]
    assert missing and all(t["text"]=="—" and t["source_ref"]["status"]=="unassessed" for t in missing)


def test_more_than_four_sources_exact_intersections_and_all_zero_mask(tmp_path):
    data=fixture(samples=1,families=2,sources=5)
    data["priority_overlap_regions"]["regions"].append({"mask":"00000","count":2})
    data["priority_overlap_regions"]["union_count"]+=2
    render_publication(data,tmp_path,sections=["comparison"])
    pages=ledger(tmp_path)["pages"]
    regions=[r for p in pages for r in p.get("intersections",[])]
    assert sum(r["original"] for r in regions)==data["overlap_regions"]["union_count"]
    assert sum(r["priority"] for r in regions)==data["priority_overlap_regions"]["union_count"]
    assert next(r for r in regions if r["mask"]=="00000")["priority"]==2
    shown={s for p in pages for c in p["cells"] for s in c.get("source_order",[])}
    assert shown=={s["id"] for s in data["sources"] if s["kind"]!="reference_annotation"}


def test_many_families_and_samples_paginate_without_dropping_identities(tmp_path):
    data=fixture(samples=7,families=13)
    render_publication(data,tmp_path,sections=["atlas"],max_rows_per_page=3,max_families_per_page=5)
    pages=ledger(tmp_path)["pages"]
    identities={(g["sample_id"],g["family_id"]) for p in pages for g in p["glyphs"]}
    assert identities=={(s["sample_id"],f["id"]) for s in data["samples"] for f in data["inventory_rows"]["families"]}
    assert len(pages)==9


def test_svg_live_literal_text_and_font_fallback_recorded(tmp_path):
    data=fixture(samples=1,families=1)
    data["samples"][0]["label"]=r"$literal_{$ & alpha α"
    render_publication(data,tmp_path,sections=["atlas"],font_family="Not An Installed Font")
    g=ledger(tmp_path)
    assert g["font_resolved"]!="Not An Installed Font"
    assert g["font_files"] and all(r["sha256"] for r in g["font_files"])
    xml=ET.parse(tmp_path/"publication_page001.svg").getroot()
    assert not xml.findall(".//{http://www.w3.org/2000/svg}image")
    texts=["".join(t.itertext()) for t in xml.findall(".//{http://www.w3.org/2000/svg}text")]
    assert data["samples"][0]["label"] in texts


def test_invalid_counts_and_masks_refused_before_any_output(tmp_path):
    data=fixture();data["cells"][0]["effector_CDS"]=12
    with pytest.raises(ValueError,match="effector"):render_publication(data,tmp_path)
    assert not list(tmp_path.iterdir())
    data=fixture();data["overlap_regions"]["union_count"]+=1
    with pytest.raises(ValueError,match="union count"):render_publication(data,tmp_path)


def test_empty_dataset_readable_and_no_overwrite(tmp_path):
    data=fixture(samples=0,families=0,sources=0)
    render_publication(data,tmp_path,sections=["comparison"])
    assert ledger(tmp_path)["page_count"]>=1
    with pytest.raises(FileExistsError):render_publication(data,tmp_path,sections=["comparison"])


def test_venn_fit_depends_on_physical_number_boxes():
    data=fixture();fonts=Fonts()
    assert venn_fits((data["overlap_regions"],data["priority_overlap_regions"]),fonts)
    large=copy.deepcopy(data["overlap_regions"])
    for r in large["regions"]:r["count"]=10**12
    assert not venn_fits((large,large),fonts)


def test_large_counts_keep_three_columns_and_common_anchors(tmp_path):
    data=fixture(samples=1,families=4)
    data['inventory_rows']['categories'][0]['N']=1192
    data['inventory_rows']['families'][0]['label']='A long synthetic family label (pooled) with additional qualification'
    render_publication(data,tmp_path,sections=['comparison'])
    pages=ledger(tmp_path)['pages']
    rows=[r for p in pages for r in p['cells'] if r['kind']=='comparison_row']
    assert all(r['source_order']==['DefenseFinder','PADLOC','DefensePredictor'] for r in rows)
    assert all(r['numeric_anchors']==rows[0]['numeric_anchors'] for r in rows)
    assert any(r['bounds_mm'][3]>5 for r in rows)
    assert any(r['bounds_mm'][3]<4 for r in rows)
    assert all('(pooled)' not in t['text'] for p in pages for t in p['texts'])


def test_saved_context_and_unassessed_effector_glyph_states(tmp_path):
    data=fixture(samples=1,families=1)
    summary=data['sample_summary'][0]
    summary.update(rm_source_counts=summary['source_counts'],rm_source_count_status=summary['source_count_status'],
                   rm_pool_CDS=15,rm_intersection_CDS=2,retained_R_CDS=1,retained_RM_CDS=3,
                   complete_saved_contexts=2,complete_context_state='observed_saved_assessment')
    data['cells'][0]['metric_status']['effector']='unassessed'
    render_publication(data,tmp_path,sections=['atlas'],atlas_style='compact')
    texts=ledger(tmp_path)['pages'][0]['texts']
    assert next(t['text'] for t in texts if (t.get('source_ref') or {}).get('field')=='complete_saved_contexts')=='2'
    glyph=ledger(tmp_path)['pages'][0]['glyphs'][0]
    assert glyph['total']==10 and glyph['outer_radius_mm']>0
    assert glyph['effector_status']=='unassessed' and glyph['inner_radius_mm']==0
    assert not any(t['text']=='?' for t in texts)


def test_scope_mismatch_rejected_and_partial_scope_disclosed(tmp_path):
    data=fixture(samples=2,families=1)
    for key in ['overlap_regions','priority_overlap_regions']:
        data[key].update(sample_ids=['sample/0'],excluded_sample_ids=['sample/1'],coverage_status='partial')
    data['priority_overlap_regions']['excluded_selected_count']=3
    render_publication(data,tmp_path/'ok',sections=['comparison'])
    assert any('Shared coverage: 1/2 genomes' in t['text'] for p in ledger(tmp_path/'ok')['pages'] for t in p['texts'])
    data['priority_overlap_regions']['sample_ids']=['sample/1']
    with pytest.raises(ValueError,match='same sample scope'):render_publication(data,tmp_path/'bad')


def test_null_inventory_and_explicit_type_context_applicability(tmp_path):
    data=fixture(samples=1,families=1)
    for row in data['inventory_rows']['families']:
        row['N']=None;row['metric_status']['inventory']='unassessed'
    iv=next(r for r in data['type_tables'] if r['rm_type']=='IV')
    iv['context_applicability']='not_assessed_by_I_III_rule';iv['context_state']='unassessed'
    render_publication(data,tmp_path,sections=['type_benchmarks','comparison'])
    texts=[t for p in ledger(tmp_path)['pages'] for t in p['texts']]
    r=next(t for t in texts if t['role']=='type_numeric' and (t.get('source_ref') or {}).get('type')=='IV' and t['source_ref']['field']=='R')
    assert r['text']=='1'


def test_invalid_typed_pool_and_alias_collision_refused(tmp_path):
    data=fixture(samples=1,families=1)
    data['type_tables'][0]['pool_CDS']=1
    with pytest.raises(ValueError,match='shared typed pool'):render_publication(data,tmp_path)
    with pytest.raises(ValueError,match='aliases collide'):
        render_publication(fixture(),tmp_path,aliases={'PADLOC':'S','DefenseFinder':'S'})
    with pytest.raises(ValueError,match='at least270'):
        render_publication(fixture(),tmp_path,atlas_page_size=[120,120])


def test_compact_atlas_respects_explicit_roles_and_display_partners(tmp_path):
    data=fixture(samples=1,families=3)
    first,middle,last=data['inventory_rows']['families']
    first.update(partner_group='synthetic pair',partner_order=0,display_group='Modification')
    last.update(partner_group='synthetic pair',partner_order=1,display_group='Modification',category='Nuclease',
                band_segments=[{'category':'Modification','fraction':.5},{'category':'Nuclease','fraction':.5}])
    data['rm_components']=[]
    for r in data['type_tables']:
        r['is_total']=r['rm_type']=='Total'
        if r['rm_type']=='I':r['role_states']['A']='not_applicable'
        if r['rm_type']=='IV':r['role_states'].update(M='not_applicable',S='not_applicable')
        for role in 'RMSA':
            data['rm_components'].append({'sample_id':r['sample_id'],'rm_type':r['rm_type'],
                                         'role':'Aux' if role=='A' else role,'count':r[role]})
    render_publication(data,tmp_path,sections=['atlas','type_benchmarks'],atlas_style='compact')
    page=ledger(tmp_path)['pages'][0]
    assert page['family_ids']==[first['id'],last['id'],middle['id']]
    components={(r['rm_type'],r['role']) for r in page['cells'] if r.get('field')=='component'}
    assert components=={('I','R'),('I','M'),('I','S'),('IV','R'),('IV','Aux')}
    assert next(b for b in page['bands'] if b['family_id']==last['id'])['segments']==last['band_segments']
    sample=next(t for t in page['texts'] if (t.get('source_ref') or {}).get('sample_id')=='sample/0' and t['role']=='label')
    assert sample['bounds_mm'][0]>max(g['bounds_mm'][0] for g in page['glyphs'])
    assert any(c.get('rm_type')=='Total' for c in ledger(tmp_path)['pages'][1]['cells'])


def test_actual_native_demo_headings_fit_and_category_order_is_stable(tmp_path):
    from defenseviz.suite_demo import write_suite_project
    from defenseviz.suite_pipeline import prepare_suite
    from defenseviz.analysis.saved import analyze_saved
    from defenseviz.plot.publication import DEFAULTS
    from defenseviz.plot.publication.comparison import ordered_families,visible_rows

    paths=write_suite_project(tmp_path/'native')
    data=analyze_saved(prepare_suite(paths['root'],manifest=paths['manifest'],supplement=paths['supplement']))
    # Discovery order is not a display contract; reversing inputs must not reorder known categories.
    data['inventory_rows']['families'].reverse()
    render_publication(data,tmp_path/'figures',sections=['atlas','comparison'])
    atlas=next(p for p in ledger(tmp_path/'figures')['pages'] if p['section']=='atlas')
    assert atlas['section_headings']['mode']=='group_anchored_staggered'
    assert not atlas['validation']['errors']
    headings=atlas['section_headings']['headings']
    assert headings
    assert all(h['group_span_mm'][0] <= h['badge_x'] <= h['group_span_mm'][1] for h in headings)
    assert {int(t['text']) for t in atlas['texts'] if t['role']=='stage_badge'}=={h['stage_id'] for h in headings}
    assert all(h['text_bounds_mm'] and h['badge_bounds_mm'] for h in headings)
    category_sequence=list(dict.fromkeys(b['primary'] for b in atlas['bands']))
    known=[x for x in DEFAULTS['category_order'] if x in category_sequence]
    assert category_sequence[:len(known)]==known
    rows=visible_rows(data['inventory_rows']['families'])
    preserved=ordered_families(rows,{'category_order':[]})
    assert preserved[0].get('display_group',preserved[0]['category'])==rows[0].get('display_group',rows[0]['category'])
    assert sum(r['count'] for r in data['overlap_regions']['regions'])==13


@pytest.mark.parametrize('coverage,display',[('complete','1'),('partial','1')])
def test_missing_retention_is_dash_and_typed_pool_records_coverage(tmp_path,coverage,display):
    from defenseviz.analysis.saved import analyze_saved
    data=analyze_saved({'samples':[{'sample_id':'s','label':'s','highlight':True}],
        'calls':[{'sample_id':'s','gene_id':'a','caller':'DefenseFinder','family':'RM_Type_I',
                  'category':'R-M','rm_type':'Type I','system_id':'x'}],
        'coverage':[{'sample_id':'s','caller':'DefenseFinder','status':coverage,'endpoint':'native_calls'}]})
    render_publication(data,tmp_path,sections=['atlas','type_benchmarks'],atlas_style='compact')
    pages=ledger(tmp_path)['pages']
    type_texts=[t for p in pages for t in p['texts'] if t['role']=='type_numeric']
    assert {t['text'] for t in type_texts if t['source_ref']['field']=='DNMB'}=={'—'}
    pool=next(t['text'] for t in type_texts if t['source_ref']['field']=='Pool' and t['source_ref']['type']=='Type I')
    assert pool==display
    retained=[t for p in pages for t in p['texts'] if (t.get('source_ref') or {}).get('field')=='retained_R_fraction']
    assert retained and {t['text'] for t in retained}=={'—'}


def test_actual_native_missing_rebase_is_unassessed_overlap_not_zero(tmp_path):
    import shutil
    from defenseviz.suite_demo import write_suite_project
    from defenseviz.suite_pipeline import prepare_suite
    from defenseviz.analysis.saved import analyze_saved
    paths=write_suite_project(tmp_path/'native')
    for module in (tmp_path/'native/runs').glob('*/dnmb_module_rebasefinder'):
        shutil.rmtree(module)
    data=analyze_saved(prepare_suite(paths['root'],manifest=paths['manifest']))
    assert data['global_summary']['native_union_CDS']>0
    assert data['overlap_regions']['coverage_status']=='unassessed'
    render_publication(data,tmp_path/'figures',sections=['atlas','comparison'],atlas_style='compact')
    pages=ledger(tmp_path/'figures')['pages']
    overlaps=[o for p in pages for o in p.get('overlaps',[])]
    assert overlaps and all(o['mode']=='unassessed' and o['reported_counts'] is None for o in overlaps)
    texts=[t['text'] for p in pages for t in p['texts']]
    assert 'Source overlap: unavailable' in texts
    assert 'No observed source memberships' not in texts
    assert not any(t['role']=='venn_count' for p in pages for t in p['texts'])
    atlas=next(p for p in pages if p['section']=='atlas')
    shared=[t for t in atlas['texts'] if (t.get('source_ref') or {}).get('field')=='rm_intersection_CDS']
    pools=[t for t in atlas['texts'] if (t.get('source_ref') or {}).get('field')=='rm_pool_CDS']
    assert shared and {t['text'] for t in shared}=={'—'}
    assert pools and all(t['text'].isdecimal() for t in pools)
    pool_cells=[c for c in atlas['cells'] if c.get('field')=='rm_pool_CDS']
    assert pool_cells and all(c['status']=='partial' for c in pool_cells)


def test_only_documented_positive_context_subsets_get_parentheses(tmp_path):
    data=fixture(samples=1,families=1)
    row=data['type_tables'][0]
    row.update(context_state='partial',incomplete_R_context=0,orphan_M=1,orphan_state={'M':'unassessed','S':'complete'})
    render_publication(data,tmp_path,sections=['type_benchmarks'])
    texts=[t for p in ledger(tmp_path)['pages'] for t in p['texts'] if t['role']=='type_numeric' and t['source_ref']['type']=='I']
    assert next(t['text'] for t in texts if t['source_ref']['field']=='R')=='1'
    assert next(t['text'] for t in texts if t['source_ref']['field']=='M')=='1'


def test_actual_native_partial_rebase_scope_stays_with_its_sample(tmp_path):
    import csv
    from defenseviz.suite_demo import write_suite_project
    from defenseviz.suite_pipeline import prepare_suite
    from defenseviz.analysis.saved import analyze_saved
    paths=write_suite_project(tmp_path/'native')
    source=tmp_path/'native/runs/isolate-a/dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv'
    with source.open() as stream:fields=next(csv.reader(stream,delimiter='\t'))
    with source.open('a') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t')
        writer.writerow({'query':'unmapped-synthetic-reference','curation_keep':'FALSE',
                         'final_family':'Type I','final_role':'R','blast_match':'Synthetic unmatched candidate'})
    data=analyze_saved(prepare_suite(paths['root'],manifest=paths['manifest']))
    render_publication(data,tmp_path/'figures',sections=['atlas'],atlas_style='compact')
    values={(t['source_ref']['sample_id'],t['source_ref']['field']):t['text']
            for p in ledger(tmp_path/'figures')['pages'] for t in p['texts']
            if (t.get('source_ref') or {}).get('field') in ('rm_pool_CDS','rm_intersection_CDS')}
    cells={(c['sample_id'],c['field']):c for p in ledger(tmp_path/'figures')['pages']
           for c in p['cells'] if c.get('field') in ('rm_pool_CDS','rm_intersection_CDS')}
    for field in ('rm_pool_CDS','rm_intersection_CDS'):
        assert values['isolate-a',field]==str(cells['isolate-a',field]['value'])
        assert cells['isolate-a',field]['status']=='partial'
        assert values['isolate-b',field]==str(cells['isolate-b',field]['value'])
        assert cells['isolate-b',field]['status']=='complete'
        assert values['isolate-c',field]=='0'

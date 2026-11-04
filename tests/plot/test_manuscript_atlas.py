"""Independent synthetic regressions for the approved manuscript visual grammar.

All identifiers, samples and counts are generated here. Physical coordinates
describe the drawing template; no manuscript observations are test fixtures.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import math

from matplotlib.colors import to_rgba
from matplotlib.patches import Circle, Rectangle
import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.plot.publication import DEFAULTS, _config
from defenseviz.plot.publication.drawing import Canvas, Fonts
from defenseviz.plot.publication.manuscript_atlas import manuscript_atlas_pages
from defenseviz.plot.publication.manuscript_model import enrich_manuscript_model


SOURCES = ["PADLOC", "DefenseFinder", "DefensePredictor"]
MASKS = ["100", "010", "001", "110", "101", "011", "111"]
CATEGORIES = [("Modification", 7), ("Nuclease", 17), ("CRISPR-Cas", 1),
              ("Abortive infection", 11), ("Other mechanisms", 4), ("Unresolved", 8)]
EXPECTED_COLUMNS = [("I_R", "role"), ("I_M", "role"), ("I_S", "role"),
                    ("I_MS", "context"), ("I_RMS", "context"),
                    ("II_R", "role"), ("II_M", "role"), ("III_R", "role"),
                    ("III_M", "role"), ("III_M", "context"), ("III_MR", "context"),
                    ("IV_R", "role"), ("IV_A", "role")]


def synthetic_result(samples=3, categories=CATEGORIES, long_family=False, long_sample=False):
    """A declared display model with distinguishable source/count channels."""
    families = []
    for category, count in categories:
        for index in range(count):
            fid = f"{category}-{index:02d}"
            label = f"Family {len(families)+1:02d}"
            if long_family:
                label += " with an extended descriptive synthetic annotation label"
            families.append(dict(id=fid, label=label, category=category, N=8*samples,
                                 counts_status="complete", metric_status={"inventory":"complete"}))
    rows = []
    for index in range(samples):
        sid = f"isolate-{index:03d}"
        label = f"Syntheticus example isolate {index:03d}"
        if long_sample:
            label = f"Syntheticus example subspecies with an unusually long strain designation specimen {index:03d}"
        sample = dict(sample_id=sid, label=label, genus="Syntheticus", species="example",
                      highlight=index in (0, 2), display={"transformable": index == 0})
        rm_roles = {key: dict(count=1, status="complete", support_count=0, support_status="complete")
                    for key, kind in EXPECTED_COLUMNS if kind == "role"}
        compositions = {key: dict(count=1, status="complete", R_support_all=key in {"I_RMS", "III_MR"},
                                 R_support_status="complete", reference_any=key == "I_RMS",
                                 reference_status="complete")
                        for key, kind in EXPECTED_COLUMNS if kind == "context"}
        masks = [dict(mask=mask, count=n, status="complete") for n, mask in enumerate(MASKS, 1)]
        counts = {source: sum(m["count"] for m in masks if m["mask"][j] == "1")
                  for j, source in enumerate(SOURCES)}
        cells = {f["id"]: dict(sample_id=sid, family_id=f["id"], N=8, associated_CDS=4,
                              effector_CDS=2, counts_status="complete",
                              metric_status={"inventory":"complete", "associated":"complete", "effector":"complete"})
                 for f in families}
        rows.append(dict(sample, rm_original=dict(source_counts=dict(zip(SOURCES, (3, 4, 5))),
            source_status={s:"complete" for s in SOURCES}, reference_counts={"REBASE":10},
            reference_status={"REBASE":"complete"}, input_pool_CDS=12, pool_status="complete",
            intersection_CDS=2, intersection_status="complete", intersection_sources=SOURCES),
            rm_roles=rm_roles, compositions=compositions, retained=dict(R_CDS=2, total_CDS=6, status="complete"),
            noncanonical=dict(source_counts=counts, source_status={s:"complete" for s in SOURCES},
                pool_CDS=28, pool_status="complete", intersection_CDS=7, intersection_status="complete",
                masks=[dict(mask="000", count=0, status="complete"), *masks]), cells=cells))
    model = dict(schema_version="manuscript-atlas-v1", source_order=SOURCES.copy(), reference_sources=["REBASE"],
                 samples=[{k:deepcopy(v) for k,v in r.items() if k in {"sample_id","label","genus","species","highlight","display"}}
                          for r in rows], rows=rows, families=families, group_records=[],
                 provenance={"operation":"synthetic_display_fixture", "scientific_assignments_changed":False})
    return {"samples":model["samples"], "manuscript_atlas":model}


def render_pages(result, **options):
    config = deepcopy(DEFAULTS)
    config.update(options)
    fonts = Fonts()
    def factory(width, height, section):
        return Canvas(width, height, fonts, section, "synthetic-atlas")
    return manuscript_atlas_pages(result, fonts, config, factory)


def texts(page):
    return [meta["text"] for _, meta in page.text_artists]


def rectangles_at(page, x, y, w, h):
    return [p for p in page.ax.patches if isinstance(p, Rectangle)
            and all(math.isclose(a, b, abs_tol=1e-7) for a,b in
                    zip((p.get_x(),p.get_y(),p.get_width(),p.get_height()),(x,y,w,h)))]


def test_48_synthetic_families_reproduce_all_fourteen_template_panel_bounds():
    result = synthetic_result()
    before = deepcopy(result)
    pages = render_pages(result)
    assert result == before and len(pages) == 1
    page = pages[0]
    # Independently pinned physical template positions, not copied observations.
    expected = {
        "RM_callers":(6.04875,28.72875), "REBASE":(29.62125,37.18125),
        "RM_pool":(38.07375,55.92375), "RM_integrated":(56.81625,107.21625),
        "Retained_summary":(108.10875,119.86875), "NonRM_callers":(124.06875,146.74875),
        "cases":(150.94875,169.84875), "NonRM_summary":(170.74125,185.23125),
        "Modification":(189.43125,211.48125), "Nuclease":(215.68125,269.23125),
        "Cas":(273.43125,276.58125), "Abortive infection":(277.47375,312.12375),
        "Other mechanisms":(313.01625,325.61625), "Unresolved":(326.50875,351.70875)}
    assert list(page.ledger["panel_bounds"]) == list(expected)
    for name, (left, right) in expected.items():
        observed = page.ledger["panel_bounds"][name]
        assert (observed["left"], observed["right"]) == pytest.approx((left, right))
        assert (observed["top"], observed["bottom"]) == pytest.approx((90, 99.15))
    assert page.width == 447
    page.validate()


def test_thirteen_rm_columns_preserve_four_group_compositions_and_rotation():
    page = render_pages(synthetic_result(samples=1))[0]
    symbols = page.ledger["rm_symbols"]
    assert [(x["field"], x["kind"]) for x in symbols] == EXPECTED_COLUMNS
    assert sum(x["count_unit"] == "groups" for x in symbols) == 4
    assert sum(x["count_unit"] == "CDSs" for x in symbols) == 9
    assert next(x for x in symbols if x["field"] == "I_RMS")["record"]["reference_any"] is True
    rotated = [(m["text"],m["rotation"]) for _,m in page.text_artists if m["rotation"]]
    assert ("MS",55) in rotated and ("RMS",55) in rotated and ("MR",55) in rotated
    assert ("Cas",90) in rotated


def test_regular_major_gaps_and_emphasis_frames_match_original_physical_units():
    page = render_pages(synthetic_result(samples=1))[0]
    p = page.ledger["panel_bounds"]
    for left, right in [("Retained_summary","NonRM_callers"),("NonRM_callers","cases"),
                        ("NonRM_summary","Modification"),("Modification","Nuclease"),("Nuclease","Cas")]:
        assert p[right]["left"]-p[left]["right"] == pytest.approx(4.2)
    for left, right in [("RM_callers","REBASE"),("REBASE","RM_pool"),("RM_pool","RM_integrated"),
                        ("RM_integrated","Retained_summary"),("cases","NonRM_summary"),
                        ("Cas","Abortive infection"),("Abortive infection","Other mechanisms"),("Other mechanisms","Unresolved")]:
        assert p[right]["left"]-p[left]["right"] == pytest.approx(.8925)
    for name in ("RM_pool","Retained_summary"):
        frame=p[name]
        patches=rectangles_at(page,frame["left"],frame["top"],frame["right"]-frame["left"],3.05)
        assert any(math.isclose(x.get_linewidth(),1.1625) and x.get_edgecolor()==to_rgba("#287C78") for x in patches)


def test_source_minibars_and_three_distinct_fraction_denominators():
    page = render_pages(synthetic_result(samples=1))[0]
    bars = page.ledger["bars"]
    source = {(b["source_ref"]["panel"],b["source_ref"]["source"]):b for b in bars if b["kind"]=="atlas_source"}
    assert {(key,b["value"]) for key,b in source.items()} >= {
        (("RM_callers","PADLOC"),3), (("RM_callers","DefenseFinder"),4),
        (("RM_callers","DefensePredictor"),5), (("REBASE","REBASE"),10),
        (("NonRM_callers","PADLOC"),17), (("NonRM_callers","DefenseFinder"),19),
        (("NonRM_callers","DefensePredictor"),21)}
    measured=page.measure()
    for panel in {key[0] for key in source}:
        ratios=[b["bounds_mm"][2]/b["value"] for key,b in source.items() if key[0]==panel]
        assert ratios==pytest.approx([ratios[0]]*len(ratios))
    for b in source.values():
        assert b["bounds_mm"][2]>0
        assert rectangles_at(page,*b["bounds_mm"])
        labels=[t for t in measured if t.get("source_ref")==b["source_ref"]]
        assert len(labels)==1 and labels[0]["text"]==str(b["value"])
        assert labels[0]["bounds_mm"][0]-(b["bounds_mm"][0]+b["bounds_mm"][2])>.1
    fractions = {b["source_ref"]["endpoint"]:b for b in bars if b["kind"]=="atlas_fraction"}
    assert {key:(v["value"],v["denominator"]) for key,v in fractions.items()} == {
        "intersection":(2,12), "retained":(2,6), "noncanonical":(7,28)}
    assert fractions["intersection"]["bounds_mm"][2] == 5
    assert fractions["retained"]["bounds_mm"][2] == 3.25
    assert fractions["noncanonical"]["bounds_mm"][2] == 3.7
    for endpoint,width,color in [("intersection",5*2/12,"#68A19C"),
                                 ("retained",3.25*2/6,"#124F4C"),
                                 ("noncanonical",3.7*7/28,"#7953A2")]:
        x,y,_,h=fractions[endpoint]["bounds_mm"]
        assert any(p.get_facecolor()==to_rgba(color) for p in rectangles_at(page,x,y,width,h))


def test_source_track_and_numeric_label_stay_separate_at_the_cohort_maximum():
    result=synthetic_result(samples=2,categories=[("Nuclease",1)])
    for row,count in zip(result["manuscript_atlas"]["rows"],(1,73)):
        row["rm_original"]["source_counts"]["PADLOC"]=count
        row["rm_original"].update(input_pool_CDS=100,intersection_CDS=1)
    page=render_pages(result)[0]
    page.validate()
    bars=[b for b in page.ledger["bars"] if b["kind"]=="atlas_source"
          and b["source_ref"]["source"]=="PADLOC" and b["source_ref"]["panel"]=="RM_callers"]
    assert [b["value"] for b in bars]==[1,73]
    assert bars[1]["bounds_mm"][2]/bars[0]["bounds_mm"][2]==pytest.approx(73)
    for b in bars:
        label=next(t for t in page.ledger["texts"] if t.get("source_ref")==b["source_ref"])
        assert label["bounds_mm"][0]-(b["bounds_mm"][0]+b["bounds_mm"][2])>.1


def test_six_exclusive_body_columns_and_seven_mask_key_preserve_all_source_patterns():
    page = render_pages(synthetic_result(samples=1))[0]
    cells = [x for x in page.ledger["cells"] if x["kind"]=="exclusive"]
    assert [(x["mask"],x["count"]) for x in cells] == list(zip(MASKS[:-1],range(1,7)))
    assert not any(x["mask"]=="111" for x in cells)
    for cell,radius in zip(cells,(.4032,.5712,.5712,.7644,.7644,.7644)):
        left,top,width,height=cell["bounds_mm"]
        dots=[q for q in page.ax.patches if isinstance(q,Circle)
              and math.isclose(q.center[0],left+width/2,abs_tol=1e-7)
              and math.isclose(q.center[1],top+height/2,abs_tol=1e-7)]
        assert len(dots)==1 and dots[0].radius==pytest.approx(radius)
        assert dots[0].get_facecolor()==to_rgba("#E78E51")
    p = page.ledger["panel_bounds"]; bottom=p["cases"]["bottom"]
    xs=[p["cases"]["left"]+(i+.5)*3.15 for i in range(6)]+[p["NonRM_summary"]["left"]+1.8]
    ys=[bottom+32.4+4*i for i in range(3)]
    for mask,x in zip(MASKS,xs):
        for bit,y in zip(mask,ys):
            dots=[q for q in page.ax.patches if isinstance(q,Circle)
                  and math.isclose(q.center[0],x,abs_tol=1e-7) and math.isclose(q.center[1],y,abs_tol=1e-7)]
            assert len(dots)==1
            assert dots[0].radius == pytest.approx(.65 if bit=="1" else .3)
    assert "n/pool" in texts(page)


@pytest.mark.parametrize("count,expected",[(1,.4032),(2,.5712),(3,.5712),(4,.7644),(7,.7644),(8,.966),(15,.966),(16,1.176)])
def test_family_outer_radius_and_effector_area_follow_original_count_bins(count,expected):
    result=synthetic_result(samples=1,categories=[("Nuclease",1)])
    cell=next(iter(result["manuscript_atlas"]["rows"][0]["cells"].values()))
    cell.update(N=count,associated_CDS=count,effector_CDS=1)
    result["manuscript_atlas"]["families"][0]["N"]=count
    page=render_pages(result)[0]
    glyph=page.ledger["glyphs"][0]
    assert glyph["outer_radius_mm"]==pytest.approx(expected)
    assert glyph["inner_radius_mm"]==pytest.approx(expected/math.sqrt(count))
    assert glyph["inner_area_fraction"]==pytest.approx(1/count)
    assert not glyph["activity_inferred"]


def test_default_atlas_preserves_five_cards_and_manual_two_of_three_target_rows():
    assert DEFAULTS["atlas_style"]=="manuscript"
    page=render_pages(synthetic_result())[0]
    assert [x["number"] for x in page.ledger["manuscript_workflow"]["cards"]]==[1,2,3,4,5]
    assert "MethREfinder" in texts(page) and "Type IV candidates" in texts(page)
    rows=[x for x in page.ledger["cells"] if x["kind"]=="strain"]
    assert [x["highlight"] for x in rows]==[True,False,True]
    for row in rows:
        x,y,w,h=row["bounds_mm"]
        target=rectangles_at(page,x-.6825,y,58,h)
        assert any(p.get_facecolor()==to_rgba("#F3BCD2") for p in target)==row["highlight"]
    assert "Target strain" not in texts(page) and "Transformable" not in texts(page)


def test_absence_shading_requires_explicit_family_and_complete_zero():
    result=synthetic_result(samples=4)
    rows=result["manuscript_atlas"]["rows"];fid=result["manuscript_atlas"]["families"][0]["id"]
    for row,state,value in zip(rows,("complete","partial","unassessed","complete"),(0,0,None,8)):
        row["display"]["absence_family_ids"]=[fid]
        cell=row["cells"][fid];cell.update(N=value,counts_status=state,associated_CDS=0,effector_CDS=0)
        cell["metric_status"]["inventory"]=state
    result["manuscript_atlas"]["families"][0].update(N=8,counts_status="partial",metric_status={"inventory":"partial"})
    page=render_pages(result)[0]
    matching=[g for g in page.ledger["glyphs"] if g["family_id"]==fid]
    shades=[]
    for glyph in matching:
        patches=rectangles_at(page,*glyph["bounds_mm"])
        shades.append(any(p.get_facecolor()==to_rgba("#D0CFE5") for p in patches))
    assert shades==[True,False,False,False]


def test_transformable_marker_requires_explicit_true_and_does_not_follow_highlight():
    result=synthetic_result()
    rows=result["manuscript_atlas"]["rows"]
    for row,value in zip(rows,(False,"yes",None)):row["display"]["transformable"]=value
    page=render_pages(result)[0]
    assert "Transformable" not in texts(page)
    assert not any(isinstance(p,Rectangle) and math.isclose(p.get_width(),1.3125)
                   and p.get_facecolor()==to_rgba("#526A75") for p in page.ax.patches)


def test_source_bound_three_native_intersection_excludes_reference_from_numerator():
    ref={"path":"synthetic_saved.tsv","sha256":"a"*64}
    ids=["a","b","c","reference-only"]
    bundle={"samples":[{"sample_id":"s","label":"Syntheticus example"}],
            "genes":[dict(sample_id="s",gene_id=g,protein_sha256=hashlib.sha256(g.encode()).hexdigest())for g in ids],
            "calls":[],"coverage":[],"saved_decisions":[],"options":{}}
    for provider,genes in [("PADLOC",ids[:3]),("DefenseFinder",ids[:2]),("DefensePredictor",ids[:3]),("REBASE",["b","reference-only"])]:
        reference=provider=="REBASE"
        bundle["coverage"].append(dict(sample_id="s",caller=provider,status="complete",endpoint="reference_matches"if reference else"native_calls"))
        for gene in genes:bundle["calls"].append(dict(sample_id="s",gene_id=gene,caller=provider,family="Type I",rm_type="Type I",category="R-M",
            evidence_kind="reference_hit"if reference else"final_member",source_ref=ref))
    for gene,role in zip(ids[:3],("R","M","S")):
        bundle["saved_decisions"].append(dict(sample_id="s",gene_id=gene,state="retained",family="Type I",role=role,source_ref=ref))
    analyzed=analyze_saved(bundle);result=enrich_manuscript_model(bundle,analyzed)
    row=result["manuscript_atlas"]["rows"][0]
    assert row["rm_original"]["intersection_sources"]==SOURCES
    assert row["rm_original"]["intersection_CDS"]==2
    assert row["rm_original"]["reference_counts"]=={"REBASE":2}
    assert analyzed["sample_summary"][0]["rm_intersection_CDS"]==1
    page=render_pages(result)[0]
    bar=next(x for x in page.ledger["bars"] if x["kind"]=="atlas_fraction" and x["source_ref"]["endpoint"]=="intersection")
    assert (bar["value"],bar["denominator"])==(2,4)


def test_complete_zero_families_hidden_but_unknown_zero_families_retained():
    result=synthetic_result(samples=1,categories=[("Nuclease",3)])
    families=result["manuscript_atlas"]["families"]
    families[0].update(N=0,counts_status="complete",metric_status={"inventory":"complete"})
    families[1].update(N=0,counts_status="unassessed",metric_status={"inventory":"unassessed"})
    for family,state in zip(families[:2],("complete","unassessed")):
        cell=result["manuscript_atlas"]["rows"][0]["cells"][family["id"]]
        cell.update(N=0,associated_CDS=0,effector_CDS=0,counts_status=state)
        cell["metric_status"]["inventory"]=state
    pages=render_pages(result)
    shown={g["family_id"] for p in pages for g in p.ledger["glyphs"]}
    assert families[0]["id"] not in shown and families[1]["id"] in shown and families[2]["id"] in shown


@pytest.mark.parametrize("samples,max_rows,max_families",[(7,3,20),(65,20,64)])
def test_genome_and_family_pagination_preserves_every_sample_family_pair(samples,max_rows,max_families):
    result=synthetic_result(samples=samples)
    pages=render_pages(result,max_rows_per_page=max_rows,max_families_per_page=max_families)
    observed=Counter((g["sample_id"],g["family_id"])for page in pages for g in page.ledger["glyphs"])
    expected={(row["sample_id"],family["id"])for row in result["manuscript_atlas"]["rows"]for family in result["manuscript_atlas"]["families"]}
    assert set(observed)==expected and set(observed.values())=={1}
    assert len(pages)>1
    for page in pages:
        assert len([r for r in page.ledger["cells"]if r["kind"]=="strain"])<=max_rows
        assert len({g["family_id"]for g in page.ledger["glyphs"]})<=max_families
        assert page.width<=DEFAULTS["atlas_page_size"][0] and page.height<=DEFAULTS["atlas_page_size"][1]
        page.validate()


@pytest.mark.parametrize("long_family,long_sample",[(False,False),(True,False),(False,True)])
def test_short_and_long_synthetic_labels_fit_without_text_collisions(long_family,long_sample):
    result=synthetic_result(samples=3,long_family=long_family,long_sample=long_sample)
    pages=render_pages(result)
    for page in pages:
        assert page.width<=DEFAULTS["atlas_page_size"][0] and page.height<=DEFAULTS["atlas_page_size"][1]
        page.validate()


def test_explicit_taxon_and_strain_labels_determine_the_required_page_width():
    result=synthetic_result(samples=1)
    display={"taxon_label":"Syntheticus example",
             "strain_label":"subspecies with an unusually long strain designation specimen 000"}
    result["manuscript_atlas"]["rows"][0]["display"].update(display)
    result["manuscript_atlas"]["samples"][0]["display"].update(display)
    pages=render_pages(result)
    assert display["strain_label"] in texts(pages[0])
    assert pages[0].width<=DEFAULTS["atlas_page_size"][0]
    pages[0].validate()


@pytest.mark.parametrize("size",[(446.99,365),(506,199.99),(300,145)])
def test_manuscript_page_limits_reject_sizes_below_the_fixed_template(size):
    with pytest.raises(ValueError,match="manuscript atlas_page_size.*447.*200"):
        _config("Synthetic atlas",{"atlas_style":"manuscript","atlas_page_size":size})


def test_compact_and_manuscript_minimum_sizes_are_distinct():
    compact=_config("Synthetic atlas",{"atlas_style":"compact","atlas_page_size":[270,145]})
    manuscript=_config("Synthetic atlas",{"atlas_style":"manuscript","atlas_page_size":[447,200]})
    assert compact["atlas_page_size"]==[270,145]
    assert manuscript["atlas_page_size"]==[447,200]
    page=render_pages(synthetic_result(samples=1),atlas_page_size=[447,200])[0]
    assert page.width<=447 and page.height<=200
    page.validate()


def test_family_labels_beyond_the_page_limit_fail_with_actionable_guidance():
    result=synthetic_result(samples=1,categories=[("Nuclease",1)])
    result["manuscript_atlas"]["families"][0]["label"]="Extended synthetic family annotation "*30
    with pytest.raises(ValueError,match="exceed atlas_page_size.*increase.*shorter display labels"):
        render_pages(result)


def test_empty_inventory_keeps_the_schematic_without_fabricated_data_rows():
    pages=render_pages(synthetic_result(samples=0,categories=[]))
    assert len(pages)==1 and len(pages[0].ledger["manuscript_workflow"]["cards"])==5
    assert pages[0].ledger["glyphs"]==[]
    assert not [c for c in pages[0].ledger["cells"]if c["kind"]=="strain"]
    pages[0].validate()


def _texts_in_cell(page, bounds):
    """Locate the actual text artists in one rendered body cell."""
    x,y,w,h=bounds
    return [meta["text"] for _,meta in page.text_artists
            if x-1e-7 <= meta["anchor_mm"][0] <= x+w+1e-7
            and y-1e-7 <= meta["anchor_mm"][1] <= y+h+1e-7]


@pytest.mark.parametrize("count",[0,2])
def test_partial_family_role_context_and_exclusive_cells_remain_visibly_distinct(count):
    result=synthetic_result(samples=3,categories=[("Nuclease",1)])
    model=result["manuscript_atlas"];family=model["families"][0];fid=family["id"]
    family.update(N=count*3,counts_status="partial",metric_status={"inventory":"partial"})
    for row,state in zip(model["rows"],("complete","partial","unassessed")):
        cell=row["cells"][fid]
        cell.update(N=count,associated_CDS=min(count,1),effector_CDS=min(count,1),counts_status=state)
        cell["metric_status"].update(inventory=state,associated="complete",effector="complete")
        row["rm_roles"]["I_R"].update(count=count,status=state)
        row["compositions"]["I_MS"].update(count=count,status=state)
        next(m for m in row["noncanonical"]["masks"] if m["mask"]=="100").update(count=count,status=state)
    before=deepcopy(result);page=render_pages(result)[0]
    assert result==before
    panel=page.ledger["panel_bounds"]["RM_integrated"]
    family_cells={g["sample_id"]:g["bounds_mm"] for g in page.ledger["glyphs"]}
    mask_cells={c["sample_id"]:c["bounds_mm"] for c in page.ledger["cells"]
                if c["kind"]=="exclusive" and c["mask"]=="100"}
    for row_index,(row,state) in enumerate(zip(model["rows"],("complete","partial","unassessed"))):
        sid=row["sample_id"];y=panel["top"]+row_index*3.05
        # These are the established template's first role and first group cells.
        bounds={"role":[panel["left"],y,3.15,3.05],
                "context":[panel["left"]+9.45,y,4.725,3.05],
                "family":family_cells[sid],"exclusive":mask_cells[sid]}
        for kind,box in bounds.items():
            labels=_texts_in_cell(page,box)
            if state=="unassessed":assert "NA" in labels,(kind,labels)
            elif state=="partial":
                assert ("0+" in labels if count==0 else any("+" in v for v in labels)),(kind,labels)
            else:
                assert not any("+" in v or v=="NA" for v in labels),(kind,labels)
                if count==0:assert not labels,(kind,labels)
            if count==0:
                x,yy,w,h=box
                assert not any(isinstance(p,Circle) and math.isclose(p.center[0],x+w/2,abs_tol=1e-7)
                    and math.isclose(p.center[1],yy+h/2,abs_tol=1e-7) for p in page.ax.patches)
    page.validate()


@pytest.mark.parametrize("states,expected",[
    (("complete","complete"),"PAD only (5 CDSs)"),
    (("complete","partial"),"PAD only (5+ CDSs)"),
    (("partial","partial"),"PAD only (5+ CDSs)"),
    (("complete","unassessed"),"PAD only (NA CDSs)"),
    (("unassessed","unassessed"),"PAD only (NA CDSs)"),
])
def test_rotated_exclusive_total_labels_respect_all_sample_coverage(states,expected):
    result=synthetic_result(samples=2,categories=[("Nuclease",1)])
    for row,state,count in zip(result["manuscript_atlas"]["rows"],states,(2,3)):
        next(m for m in row["noncanonical"]["masks"] if m["mask"]=="100").update(count=count,status=state)
    page=render_pages(result)[0]
    labels=[meta for _,meta in page.text_artists if meta["text"].startswith("PAD only (")]
    assert len(labels)==1 and labels[0]["text"]==expected
    assert labels[0]["rotation"]==90
    page.validate()


@pytest.mark.parametrize("states,expected_value",[
    (("complete","complete"),"5"),
    (("complete","partial"),"5+"),
    (("partial","partial"),"5+"),
    (("complete","unassessed"),"NA"),
    (("unassessed","unassessed"),"NA"),
])
def test_members_outside_shown_sources_are_disclosed_beneath_the_mask_key(states,expected_value):
    result=synthetic_result(samples=2,categories=[("Nuclease",1)])
    for row,state,count in zip(result["manuscript_atlas"]["rows"],states,(2,3)):
        next(m for m in row["noncanonical"]["masks"] if m["mask"]=="000").update(count=count,status=state)
        row["noncanonical"]["pool_CDS"]+=count
    before=deepcopy(result);page=render_pages(result)[0]
    assert result==before
    disclosures=[meta for _,meta in page.text_artists if meta["text"].startswith("Outside shown sources:")]
    assert len(disclosures)==1
    assert disclosures[0]["text"]==f"Outside shown sources: {expected_value} CDSs"
    key_label=next(meta for _,meta in page.text_artists if meta["text"]=="n/pool")
    assert disclosures[0]["anchor_mm"][1]>key_label["anchor_mm"][1]
    fractions=[b for b in page.ledger["bars"] if b["kind"]=="atlas_fraction"
               and b["source_ref"]["endpoint"]=="noncanonical"]
    assert [b["denominator"] for b in fractions]==[30,31]
    assert [b["value"] for b in fractions]==[7,7]
    page.validate()


def test_complete_zero_outside_shown_sources_adds_no_unnecessary_note():
    page=render_pages(synthetic_result(samples=1,categories=[("Nuclease",1)]))[0]
    assert not any(label.startswith("Outside shown sources:") for label in texts(page))


def _declared_taxonomy_result():
    result=synthetic_result(samples=2,categories=[("Nuclease",1)])
    result["manuscript_atlas"]["taxonomy_legend"]={
        "genera":[{"key":"Syntheticus","label":"Syntheticus","color":"#123456"}],
        "species":[{"key":"example","label":"S. example","color":"#C45123"}]}
    return result


@pytest.mark.parametrize("explicit_matching",[False,True])
def test_declared_taxonomy_colors_match_actual_row_strips_and_legend_swatches(explicit_matching):
    result=_declared_taxonomy_result();model=result["manuscript_atlas"]
    if explicit_matching:
        for sample in model["samples"]+model["rows"]:
            sample["display"].update(genus_color="#123456",species_color="#C45123")
    before=deepcopy(result);page=render_pages(result)[0]
    assert result==before
    top=page.ledger["panel_bounds"]["RM_callers"]["top"]
    for x,color in [(2.35,"#123456"),(4.03,"#C45123")]:
        for row in range(2):
            strips=rectangles_at(page,x,top+row*3.05,1.365,3.05)
            assert len(strips)==1 and strips[0].get_facecolor()==to_rgba(color)
        swatches=[p for p in page.ax.patches if isinstance(p,Rectangle)
                  and math.isclose(p.get_width(),2.2) and math.isclose(p.get_height(),2.2)
                  and p.get_facecolor()==to_rgba(color)]
        assert len(swatches)==1
    page.validate()


@pytest.mark.parametrize("field",["genus_color","species_color"])
def test_contradictory_sample_and_declared_taxonomy_colors_are_rejected(field):
    result=_declared_taxonomy_result()
    # A later sample must not escape validation because the first row agrees.
    model=result["manuscript_atlas"]
    for sample in (model["samples"][-1],model["rows"][-1]):
        sample["display"][field]="#FF0000"
    before=deepcopy(result)
    with pytest.raises(ValueError,match="(?i)(taxonomy.*colo|contradict.*colo|colo.*conflict)"):
        render_pages(result)
    assert result==before


def test_later_explicit_taxonomy_color_fills_earlier_unset_rows_and_null_legend_color():
    result=synthetic_result(samples=3,categories=[("Nuclease",1)])
    model=result["manuscript_atlas"]
    # This is the default projection when the first observed sample had no
    # explicit color: the key exists, but its legend color is still unset.
    model["taxonomy_legend"]={
        "genera":[{"key":"Syntheticus","label":"Syntheticus","color":None}],
        "species":[{"key":"example","label":"S. example","color":None}]}
    for sample in (model["samples"][1],model["rows"][1]):
        sample["display"].update(genus_color="#123456",species_color="#C45123")
    before=deepcopy(result);page=render_pages(result)[0]
    assert result==before
    top=page.ledger["panel_bounds"]["RM_callers"]["top"]
    for x,color in [(2.35,"#123456"),(4.03,"#C45123")]:
        for row in range(3):
            strips=rectangles_at(page,x,top+row*3.05,1.365,3.05)
            assert len(strips)==1 and strips[0].get_facecolor()==to_rgba(color)
        swatches=[p for p in page.ax.patches if isinstance(p,Rectangle)
                  and math.isclose(p.get_width(),2.2) and math.isclose(p.get_height(),2.2)
                  and p.get_facecolor()==to_rgba(color)]
        assert len(swatches)==1
    page.validate()


def test_conflicting_later_taxonomy_colors_rejected_when_first_sample_color_is_unset():
    for field in ("genus_color","species_color"):
        result=synthetic_result(samples=3,categories=[("Nuclease",1)])
        model=result["manuscript_atlas"]
        model["taxonomy_legend"]={
            "genera":[{"key":"Syntheticus","label":"Syntheticus","color":None}],
            "species":[{"key":"example","label":"S. example","color":None}]}
        for index,color in [(1,"#123456"),(2,"#654321")]:
            for sample in (model["samples"][index],model["rows"][index]):
                sample["display"][field]=color
        before=deepcopy(result)
        with pytest.raises(ValueError,match="(?i)(taxonomy.*colo|contradict.*colo|colo.*conflict)"):
            render_pages(result)
        assert result==before

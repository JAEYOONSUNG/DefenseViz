"""Actual editable titles, circular badges and group anchors on small atlases."""
from copy import deepcopy
from itertools import combinations
import xml.etree.ElementTree as ET

import matplotlib
from matplotlib.patches import Circle
import pytest

from tests.plot.test_manuscript_atlas import CATEGORIES, render_pages, synthetic_result
from defenseviz.plot.publication.manuscript_atlas import _boxes_overlap


TITLES = {
    1: ["R-M evidence and curation"],
    2: ["Methylation-associated", "evidence"],
    3: ["Non-canonical", "DNA modification", "system"],
    4: ["Nuclease-associated", "evidence"],
    5: ["Complementary defense priorities"],
}
COLORS = {1:"#287C78", 2:"#829764", 3:"#B29BCB", 4:"#D9878C", 5:"#52626A"}
COMPACT = [("Modification",4),("Nuclease",6),("CRISPR-Cas",1),
           ("Abortive infection",5),("Unresolved",6)]


def assert_actual_headers(page, expected):
    ledger = page.ledger["section_headings"]
    headings = {h["stage_id"]:h for h in ledger["headings"]}
    assert set(headings) == set(expected)
    texts = page.measure()
    for sid in expected:
        h=headings[sid]
        assert h["title"] == " ".join(TITLES[sid])
        left,right=h["group_span_mm"]
        assert left <= h["badge_x"] <= right
        assert h["group_keys"]
        assert left == page.ledger["panel_bounds"][h["group_keys"][0]]["left"]
        assert right == page.ledger["panel_bounds"][h["group_keys"][-1]]["right"]
        badges=[p for p in page.ax.patches if p.get_gid()==f"atlas-stage-badge-{sid}"]
        assert len(badges)==1 and isinstance(badges[0],Circle)
        assert badges[0].radius == 2.1
        assert badges[0].center == (h["badge_x"],h["y"])
        title=[t for t in texts if t["role"]=="stage_title" and t["source_ref"]["stage_id"]==sid]
        number=[t for t in texts if t["role"]=="stage_badge" and t["source_ref"]["stage_id"]==sid]
        assert [t["text"] for t in title] == TITLES[sid]
        assert [t["text"] for t in number] == [str(sid)]
        assert all(t["font_pt"]==12 and t["bold"] for t in title)
        assert all(left<=t["bounds_mm"][0]<right for t in title)
        assert number[0]["font_pt"]==10.5 and number[0]["color"]==COLORS[sid]
        assert [t["bounds_mm"] for t in title] == h["text_bounds_mm"]
    # Circle strokes and full title ink must not touch another numbered group.
    for a,b in combinations(headings.values(),2):
        for aa in [a["badge_bounds_mm"],*a["text_bounds_mm"]]:
            for bb in [b["badge_bounds_mm"],*b["text_bounds_mm"]]:
                assert not _boxes_overlap(aa,bb,.15)
    # Inspect real line artists including their physical stroke, not only the
    # bracket ledger. Text-only validation cannot detect a strike-through.
    bracket_lines=[line for line in page.ax.lines
                   if (line.get_gid() or "").startswith("atlas-stage-bracket-5-")]
    assert len(bracket_lines)==(3 if 5 in expected else 0)
    for line in bracket_lines:
        xs=list(line.get_xdata());ys=list(line.get_ydata())
        pad=line.get_linewidth()*25.4/72/2+.05
        stroke=[min(xs)-pad,min(ys)-pad,max(xs)-min(xs)+2*pad,max(ys)-min(ys)+2*pad]
        assert not any(_boxes_overlap(stroke,t["bounds_mm"],0) for t in texts)
    return headings


@pytest.mark.parametrize("samples",[2,3])
def test_compact_inventory_keeps_all_literal_titles_and_circular_svg_badges(samples,tmp_path):
    result=synthetic_result(samples,categories=COMPACT);before=deepcopy(result)
    page=render_pages(result)[0]
    assert result==before
    headings=assert_actual_headers(page,range(1,6))
    assert page.ledger["section_headings"]["mode"]=="group_anchored_staggered"
    assert any(h["y"]!=64 for h in headings.values())
    assert headings[4]["anchor_mode"]=="left"
    assert headings[4]["badge_x"]==pytest.approx(headings[4]["group_span_mm"][0]+2.1)
    bracket=page.ledger["section_headings"]["bracket"]
    assert bracket["segments_mm"][0][1] < min(t["bounds_mm"][1] for t in page.measure() if t["text"]=="Abortive infection" and t["rotation"]==90)
    page.validate()
    path=tmp_path/"atlas.svg"
    with matplotlib.rc_context({"svg.fonttype":"none"}):page.fig.savefig(path,format="svg")
    svg=ET.parse(path).getroot();ns={"s":"http://www.w3.org/2000/svg"}
    for sid in range(1,6):
        badge=svg.find(f".//*[@id='atlas-stage-badge-{sid}']")
        assert badge is not None and badge.find("s:path",ns) is not None
        numeral=svg.find(f".//*[@id='atlas-stage-number-{sid}']")
        assert [t.text for t in numeral.findall(".//s:text",ns)] == [str(sid)]
        for i,line in enumerate(TITLES[sid]):
            title=svg.find(f".//*[@id='atlas-stage-title-{sid}-{i}']")
            assert [t.text for t in title.findall(".//s:text",ns)] == [line]


@pytest.mark.parametrize("samples",[3,51])
def test_full_width_original_anchors_are_independent_of_sample_count(samples):
    page=render_pages(synthetic_result(samples,categories=CATEGORIES))[0]
    headings=assert_actual_headers(page,range(1,6))
    assert page.ledger["section_headings"]["mode"]=="original_group_anchors"
    assert page.ledger["section_headings"]["body_offset_mm"]==0
    assert [headings[i]["badge_x"] for i in range(1,6)] == pytest.approx(
        [8.14875,126.16875,191.53125,237.18125,290.43125])
    assert all(h["y"]==64 for h in headings.values())
    assert all(p["top"]==90 for p in page.ledger["panel_bounds"].values())
    bracket=page.ledger["section_headings"]["bracket"]
    assert bracket["segments_mm"][0][1]==73
    assert bracket["segments_mm"][1][3]==80


def test_single_column_groups_measure_ink_and_absent_groups_do_not_get_fake_columns():
    page=render_pages(synthetic_result(2,categories=[("Modification",1),("Nuclease",1),("Unresolved",1)]))[0]
    heads=assert_actual_headers(page,range(1,6));page.validate()
    assert all(heads[i]["anchor_mode"]=="stacked_left" for i in (3,4,5))
    absent=render_pages(synthetic_result(2,categories=[("Nuclease",1)]))[0]
    assert_actual_headers(absent,[1,2,4]);absent.validate()
    assert {h["stage_id"] for h in absent.ledger["section_headings"]["absent"]}=={3,5}
    assert all(h["reason"]=="no_displayed_column_group" for h in absent.ledger["section_headings"]["absent"])
    assert "Modification" not in absent.ledger["panel_bounds"]
    assert "Unresolved" not in absent.ledger["panel_bounds"]
    assert len(absent.ledger["glyphs"])==2


def test_fraction_numerator_and_denominator_use_their_own_declared_status():
    result=synthetic_result(3,categories=COMPACT)
    for i,row in enumerate(result["manuscript_atlas"]["rows"]):
        row["retained"].update(status="partial",R_status="partial" if i==1 else "complete",
                               total_status="partial" if i==2 else "complete")
        row["rm_original"].update(intersection_status="partial",intersection_count_status="complete",
                                  pool_status="partial" if i==2 else "complete")
    page=render_pages(result)[0];page.validate()
    ts=page.ledger["texts"]
    def labels(endpoint,part):
        return [t["text"] for t in ts if t["source_ref"]=={"endpoint":endpoint,"part":part}]
    assert labels("retained","numerator")==["2","2+","2"]
    assert labels("retained","denominator")==["6","6","6+"]
    assert labels("intersection","numerator")==["2"]*3
    assert labels("intersection","denominator")==["12","12","12+"]
    bars=[b for b in page.ledger["bars"] if b["source_ref"]=={"endpoint":"retained"}]
    assert [(b["status"],b["denominator_status"]) for b in bars]==[("complete","complete"),("partial","complete"),("complete","partial")]
    cells=[x for x in page.ledger["cells"] if x.get("field")=="rm_intersection_CDS"]
    assert [x["status"] for x in cells]==["complete"]*3


def test_unknown_fraction_numerator_preserves_known_denominator_without_zero_fill():
    from defenseviz.plot.publication.drawing import Canvas, Fonts
    from defenseviz.plot.publication.manuscript_atlas import _fraction, _fraction_layout
    records=[(None,11,"unassessed","complete"),(4,None,"complete","unassessed"),
             (None,11,"unassessed","partial"),(None,None,"unassessed","unassessed")]
    fonts=Fonts();layout=_fraction_layout(fonts,records,"retained")
    page=Canvas(50,30,fonts,"fraction","independent-scope")
    for i,(n,d,state,ds) in enumerate(records):
        _fraction(page,8,5+i*5,layout["width_mm"],n,d,state,"retained",11,ds,layout)
    page.validate()
    assert [t["text"] for t in page.ledger["texts"]]==["NA","/","11","4","/","NA","NA","/","11+","NA"]
    bars=page.ledger["bars"]
    assert len(bars)==3
    assert bars[0]["value"] is None and bars[0]["status"]=="unassessed"
    assert bars[0]["denominator"]==11 and bars[0]["denominator_status"]=="complete"
    assert bars[0]["bounds_mm"][2]==pytest.approx(3.25) and bars[0]["filled_width_mm"]==0
    assert bars[1]["value"]==4 and bars[1]["denominator"] is None
    assert bars[1]["bounds_mm"][2]==0 and bars[1]["filled_width_mm"]==0


@pytest.mark.parametrize("current_gap",[False,True])
def test_analysis_to_actual_atlas_preserves_current_scope_and_real_partial_marks(current_gap):
    from tests.test_atlas_rm_scope import bundle, project
    data=bundle()
    if current_gap:
        data["saved_decisions"][1]["state"]="unassessed"
    result,row=project(data)
    assert result["availability"]["saved_decisions_by_sample"][0]["status"]=="partial"
    page=render_pages(result,aliases={"NativeA":"A","NativeB":"B","NativeC":"C","Reference":"REF"})[0];page.validate()
    ts=page.ledger["texts"]
    labels=[t["text"] for t in ts if (t.get("source_ref") or {}).get("endpoint")=="retained"]
    assert labels==(["1+","/","2+"] if current_gap else ["1","/","3"])
    intersections=[t["text"] for t in ts if (t.get("source_ref") or {}).get("endpoint")=="intersection"]
    assert intersections==["3","/","3"]


def test_result_only_fallback_does_not_suppress_genuine_partial_retention():
    from tests.test_atlas_rm_scope import bundle
    from defenseviz.analysis.saved import analyze_saved
    data=bundle();del data["rm_evaluation"]
    result=analyze_saved(data)
    assert "manuscript_atlas" not in result
    page=render_pages(result,aliases={"NativeA":"A","NativeB":"B","NativeC":"C","Reference":"REF"})[0];page.validate()
    labels=[t["text"] for t in page.ledger["texts"] if (t.get("source_ref") or {}).get("endpoint")=="retained"]
    assert labels==["1+","/","3+"]

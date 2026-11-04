"""Physical fraction slots reserve complete labels, including partial markers."""
from copy import deepcopy

import pytest

from defenseviz.plot.publication import DEFAULTS
from defenseviz.plot.publication.drawing import Canvas, Fonts, number
from defenseviz.plot.publication.manuscript_atlas import (
    NUM, _fraction, _fraction_layout, _layout, manuscript_atlas_pages,
)


def fraction_canvas(records, kind="retained"):
    fonts=Fonts()
    layout=_fraction_layout(fonts,records,kind)
    canvas=Canvas(layout["width_mm"]+20,10+len(records)*4,fonts,"fraction","test-fraction")
    scale=max([1]+[d for _,d,_,_ in records if d is not None])
    for i,(n,d,state,denstate) in enumerate(records):
        _fraction(canvas,10,6+i*4,layout["width_mm"],n,d,state,kind,scale,
                  denstate=denstate,layout=layout)
    canvas.validate()
    return canvas,layout


def test_partial_markers_expand_labels_without_sacrificing_original_retained_track():
    records=[(4,11,"partial","partial"),(3,9,"partial","partial"),(2,8,"complete","complete")]
    canvas,layout=fraction_canvas(records)
    assert layout["width_mm"] > 11.76
    assert layout["track_width_mm"] == pytest.approx(3.25)
    bars=canvas.ledger["bars"]
    assert len({b["track_width_mm"] for b in bars}) == 1
    assert [b["bounds_mm"][2]/b["denominator"] for b in bars] == pytest.approx([layout["track_width_mm"]/11]*3)
    assert [b["filled_width_mm"]/b["value"] for b in bars] == pytest.approx([layout["track_width_mm"]/11]*3)
    texts=canvas.ledger["texts"]
    assert [t["text"] for t in texts] == ["4+","/","11+","3+","/","9+","2","/","8"]
    for part in ("numerator","denominator"):
        slots=[t for t in texts if t["source_ref"]["part"]==part]
        assert all(t["alignment"] == "right" for t in slots)
        assert len({t["anchor_mm"][0] for t in slots}) == 1
        assert all(t["font_pt"] == NUM for t in slots)
        assert all(t["bold"] == (part=="numerator") for t in slots)
    for bar,num in zip(bars,[t for t in texts if t["source_ref"]["part"]=="numerator"]):
        assert num["bounds_mm"][0] - (bar["bounds_mm"][0]+bar["bounds_mm"][2]) > .12


@pytest.mark.parametrize("kind",["retained","intersection","noncanonical"])
@pytest.mark.parametrize("state,denstate",[("complete","complete"),("partial","partial"),("complete","partial"),("partial","complete")])
def test_two_to_four_digit_counts_with_separate_denominator_states(kind,state,denstate):
    records=[(12,34,state,denstate),(123,456,state,denstate),(1234,5678,state,denstate),(0,5678,state,denstate)]
    canvas,layout=fraction_canvas(records,kind)
    assert layout["track_width_mm"] >= .5-1e-8
    assert canvas.ledger["validation"]["errors"] == []
    assert all(t["font_pt"]==NUM for t in canvas.ledger["texts"])
    assert [t["text"] for t in canvas.ledger["texts"] if t["source_ref"]["part"]=="denominator"] == [number(d,denstate) for _,d,_,_ in records]
    for t in canvas.ledger["texts"]:
        x,y,w,h=t["bounds_mm"]
        assert x >= 10-.08 and x+w <= 10+layout["width_mm"]+.08


def test_complete_small_counts_keep_approved_panel_and_track_dimensions():
    plans={}
    for kind,nominal,track in [("retained",11.76,3.25),("intersection",17.85,5.),("noncanonical",14.49,3.7)]:
        _,plans[kind]=fraction_canvas([(2,6,"complete","complete")],kind)
        assert plans[kind]["width_mm"]==pytest.approx(nominal)
        assert plans[kind]["track_width_mm"]==pytest.approx(track)
    before=_layout([],['PADLOC','DefenseFinder','DefensePredictor'],['REBASE'],{})[0]
    after=_layout([],['PADLOC','DefenseFinder','DefensePredictor'],['REBASE'],{},plans)[0]
    assert after==before


def test_four_digit_partial_counts_expand_column_without_hiding_bar():
    records=[(1234,5678,"partial","partial")]
    _,layout=fraction_canvas(records)
    assert layout["width_mm"]>11.76
    assert layout["track_width_mm"]==pytest.approx(3.25)
    panels=_layout([],['PADLOC'],[],{},dict(retained=layout))[0]
    assert panels["Retained_summary"]["right"]-panels["Retained_summary"]["left"]==pytest.approx(layout["width_mm"])
    assert panels["NonRM_callers"]["left"]-panels["Retained_summary"]["right"]==pytest.approx(4.2)
    c=Canvas(40,20,Fonts(),"fraction","too-narrow")
    with pytest.raises(ValueError,match="expand the atlas fraction column"):
        _fraction(c,5,10,11.76,*records[0][:3],"retained",denstate="partial")


def test_missing_endpoint_and_zero_counts_remain_explicit():
    c,layout=fraction_canvas([(None,None,"unassessed",None),(0,0,"complete","complete"),(0,20,"partial","unassessed")])
    assert [t["text"] for t in c.ledger["texts"]]==["NA","0","/","0","0+","/","NA"]
    assert c.ledger["bars"][0]["bounds_mm"][2]==0


def atlas_result():
    rows=[]
    sources=["PADLOC","DefenseFinder","DefensePredictor"]
    for i,(n,d) in enumerate([(4,11),(3,9),(2,8)]):
        rows.append(dict(sample_id=f"synthetic-{i}",label=f"Synthetic isolate {i}",
            retained=dict(R_CDS=n,total_CDS=d,status="partial"),
            rm_original=dict(source_counts=dict(zip(sources,(20,24,30))),source_status={s:"complete" for s in sources},
                             reference_counts={"REBASE":14},reference_status={"REBASE":"complete"},
                             intersection_CDS=12,input_pool_CDS=34,intersection_status="partial",pool_status="partial"),
            noncanonical=dict(source_counts={s:34 for s in sources},source_status={s:"complete" for s in sources},intersection_CDS=34,pool_CDS=34,
                              intersection_status="partial",pool_status="partial",
                              masks=[dict(mask=format(j,"03b"),count=34 if j==7 else 0,status="complete")for j in range(8)]),
            rm_roles={},compositions={},cells={}))
    return dict(samples=[dict(sample_id=r["sample_id"],label=r["label"])for r in rows],
                manuscript_atlas=dict(schema_version="manuscript-atlas-v1",rows=rows,
                    source_order=sources,reference_sources=["REBASE"],families=[],group_records=[],provenance={}))


def test_full_atlas_uses_column_reservations_across_page_chunks_without_input_mutation():
    result=atlas_result();before=deepcopy(result)
    fonts=Fonts();config=deepcopy(DEFAULTS);config["max_rows_per_page"]=2
    pages=manuscript_atlas_pages(result,fonts,config,lambda w,h,s:Canvas(w,h,fonts,s,"synthetic"))
    assert result==before and len(pages)==2
    for page in pages:
        page.validate()
    assert pages[0].ledger["fraction_layout"]==pages[1].ledger["fraction_layout"]
    tracks=[bar["track_width_mm"] for page in pages for bar in page.ledger["bars"] if bar["source_ref"].get("endpoint")=="retained"]
    assert len(set(tracks))==1


def test_large_counts_respect_existing_page_size_limit_instead_of_shrinking_fonts():
    result=atlas_result()
    for row in result["manuscript_atlas"]["rows"]:
        row["retained"].update(R_CDS=1234,total_CDS=5678)
    fonts=Fonts();config=deepcopy(DEFAULTS);config["atlas_page_size"]=[100,500]
    with pytest.raises(ValueError,match="atlas_page_size"):
        manuscript_atlas_pages(result,fonts,config,lambda w,h,s:Canvas(w,h,fonts,s,"synthetic"))

"""Schematic four-set Venns with an exact, counted intersection fallback."""
from __future__ import annotations

import math
from matplotlib.patches import Ellipse
from matplotlib.transforms import Affine2D

from .drawing import INK, TEAL, alias, source_color

# Data-independent four-ellipse geometry. No dataset memberships or counts live here.
ELLIPSES = ((.35,.4,.72,.43,140),(.45,.5,.72,.43,140),
            (.544,.5,.72,.43,40),(.644,.4,.72,.43,40))
ANCHORS = {"0001":(.8952066,.5273268),"0010":(.6321199,.7153146),
 "0011":(.7693998,.5915631),"0100":(.3623601,.7151546),
 "0101":(.6835940,.2840049),"0110":(.497,.6460481),
 "0111":(.6780047,.4375889),"1000":(.0987934,.5273268),
 "1001":(.497,.1711642),"1010":(.3084060,.2857649),
 "1011":(.3606801,.2545119),"1100":(.2249821,.5957034),
 "1101":(.6333199,.2545119),"1110":(.3149553,.4361391),
 "1111":(.4970181,.3532812)}
BOUNDS=(.0415336997,.9524663003,.1159691186,.7840308814)


def regions(overlap):
    return {r["mask"]:r["count"] for r in overlap.get("regions", [])}


def totals(overlap):
    order=overlap.get("source_order",[]);r=regions(overlap)
    return [sum(n for mask,n in r.items() if mask[i]=="1") for i in range(len(order))]


def scope_note(before,after):
    excluded=before.get("excluded_sample_ids",[])
    selected=after.get("excluded_selected_count")
    if excluded or before.get("coverage_status") not in (None,"complete"):
        included=len(before.get("sample_ids",[]));total=included+len(excluded)
        return f"Shared coverage: {included}/{total} genomes"+(f" · {selected:,} priorities outside scope" if selected else "")
    return ""


def _quadratic_extrema(box, ellipse):
    x1,y1,x2,y2=box;cx,cy,w,h,angle=ellipse
    a=math.radians(angle);c,s=math.cos(a),math.sin(a)
    aa=(c/(w/2))**2+(s/(h/2))**2
    bb=c*s*((w/2)**-2-(h/2)**-2)
    dd=(s/(w/2))**2+(c/(h/2))**2
    def q(x,y):
        dx,dy=x-cx,y-cy
        return aa*dx*dx+2*bb*dx*dy+dd*dy*dy
    candidates=[(x,y) for x in (x1,x2) for y in (y1,y2)]
    for x in (x1,x2): candidates.append((x,min(y2,max(y1,cy-bb*(x-cx)/dd))))
    for y in (y1,y2): candidates.append((min(x2,max(x1,cx-bb*(y-cy)/aa)),y))
    if x1<=cx<=x2 and y1<=cy<=y2: candidates.append((cx,cy))
    values=[q(x,y) for x,y in candidates]
    return min(values),max(values)


def venn_fits(overlaps, fonts, width=87., height=37., pt=7.5):
    if any(len(o.get("source_order",[]))!=4 for o in overlaps): return False
    if overlaps[0]["source_order"]!=overlaps[1]["source_order"]: return False
    sx=width/(BOUNDS[1]-BOUNDS[0]);sy=height/(BOUNDS[3]-BOUNDS[2])
    for mask,(x,y) in ANCHORS.items():
        tw=max(fonts.width(str(regions(o).get(mask,0)),pt,True) for o in overlaps)+.2
        th=fonts.height("0123456789",pt,True)+.2
        box=(x-tw/2/sx,y-th/2/sy,x+tw/2/sx,y+th/2/sy)
        for bit,e in zip(mask,ELLIPSES):
            lo,hi=_quadratic_extrema(box,e)
            if bit=="1" and hi>=.997: return False
            if bit=="0" and lo<=1.003: return False
    return True


def _priority_heading(overlap, complete_label="Prioritized CDSs"):
    assessment=overlap.get("selection_assessment")
    return "Observed priorities" if assessment is not None and not assessment.get("complete",False) else complete_label


def draw_venn(canvas, overlap, x, y, config, *, heading):
    width,height=87.,37.;sx=width/(BOUNDS[1]-BOUNDS[0]);sy=height/(BOUNDS[3]-BOUNDS[2])
    ox=x-BOUNDS[0]*sx;oy=y+BOUNDS[3]*sy
    order=overlap["source_order"];counts=totals(overlap)
    union=sum(regions(overlap).values())
    display_heading=_priority_heading(overlap,heading) if heading=="Prioritized CDSs" else heading
    canvas.text(x+width/2,y-8,display_heading+f"  {union:,}",8.6,bold=True,ha="center")
    for i,(name,n) in enumerate(zip(order,counts)):
        canvas.text(x+(i+.5)*width/4,y-3,alias(name,config)+f" {n:,}",7.5,bold=True,
                    ha="center",color=source_color(name,config),role="source_label")
    for name,(cx,cy,w,h,angle) in zip(order,ELLIPSES):
        shade=source_color(name,config)
        e=Ellipse((cx,cy),w,h,angle=angle,facecolor=shade,edgecolor=shade,alpha=.14,linewidth=.6)
        e.set_transform(Affine2D().scale(sx,-sy).translate(ox,oy)+canvas.ax.transData)
        canvas.ax.add_patch(e)
        outline=Ellipse((cx,cy),w,h,angle=angle,fill=False,edgecolor=shade,linewidth=.45)
        outline.set_transform(e.get_data_transform());canvas.ax.add_patch(outline)
    vals=regions(overlap)
    records=[]
    for mask,(nx,ny) in ANCHORS.items():
        xx,yy=ox+nx*sx,oy-ny*sy;n=vals.get(mask,0)
        canvas.text(xx,yy,str(n),7.5,bold=True,ha="center",role="venn_count",
                    ref={"endpoint":heading,"mask":mask,"value":n})
        records.append({"mask":mask,"count":n,"anchor_mm":[xx,yy]})
    canvas.ledger.setdefault("overlaps",[]).append({"mode":"venn4","source_order":order,
        "regions":records,"set_sizes":dict(zip(order,counts)),"union_count":union,
        "bounds_mm":[x,y,width,height],"same_geometry":True,"areas_schematic":True,"endpoint":heading,
        "display_heading":display_heading,"selection_assessment":overlap.get("selection_assessment"),
        "normalized_ellipses":ELLIPSES,"transform":{"sx":sx,"sy":sy,"ox":ox,"oy":oy}})
    outside=vals.get("0000",0)
    if outside:
        canvas.text(x+width/2,y+height+3,f"Outside original sources: {outside:,}",7,ha="center")


def draw_flow(canvas, before, after, config):
    """Return false if an exact intersection-table page is needed."""
    canvas.ledger["selection_assessment"]=after.get("selection_assessment")
    if before.get("coverage_status")=="unassessed" or after.get("coverage_status")=="unassessed":
        canvas.rect(11,28,194,43,"#F7FAF9","#C8D3D2",.4)
        canvas.text(108,39,"Source overlap: NA",10,bold=True,ha="center",role="missing_overlap")
        canvas.text(108,49,"No common complete source coverage",8,ha="center")
        excluded=len(before.get("excluded_sample_ids",[]))
        canvas.text(108,59,f"0 of {excluded} genomes have all configured source outputs",7.5,ha="center")
        canvas.ledger.setdefault("overlaps",[]).append({"mode":"unassessed","source_order":before.get("source_order",[]),
            "reported_counts":None,"reason":"no_common_complete_source_scope","sample_ids":before.get("sample_ids",[]),
            "excluded_sample_ids":before.get("excluded_sample_ids",[]),"bounds_mm":[11,28,194,43]})
        return True
    note=scope_note(before,after)
    if note:canvas.text(canvas.width-8.45,21,note,7,ha="right",role="coverage_note")
    labels_fit=all(canvas.fonts.width(alias(s,config)+f" {n:,}",7.5,True)<=19.75
                   for o in (before,after) for s,n in zip(o.get("source_order",[]),totals(o)))
    if config.get("overlap_mode","auto")=="intersection" or not labels_fit or not venn_fits((before,after),canvas.fonts):
        order=before.get("source_order",[]);b,a=regions(before),regions(after)
        masks=sorted((m for m in set(b)|set(a) if b.get(m,0) or a.get(m,0)),key=lambda m:(-m.count("1"),m))
        if config.get("overlap_mode")!="intersection" and len(order)<=6 and len(masks)<=15:
            # Exact compact matrix keeps a small source set on the same comparison page.
            canvas.text(11,27,"Membership",8,bold=True)
            canvas.text(171,27,"Original",7.8,ha="right",bold=True)
            canvas.text(204,27,_priority_heading(after,"Prioritized"),7.8,ha="right",bold=True)
            xs=[62+j*16 for j in range(len(order))]
            if xs and xs[-1]>142:xs=[42+j*16 for j in range(len(order))]
            for name,x in zip(order,xs):canvas.text(x,27,canvas.fonts.wrap(alias(name,config),14,7.5,True),7.5,bold=True,ha="center",color=source_color(name,config))
            for i,m in enumerate(masks):
                yy=32+i*2.9
                for name,x,bit in zip(order,xs,m):canvas.circle(x,yy,.65,source_color(name,config) if bit=="1" else "#DCE3E3")
                canvas.text(171,yy,f"{b.get(m,0):,}",7.5,ha="right",role="numeric")
                canvas.text(204,yy,f"{a.get(m,0):,}",7.5,ha="right",role="numeric")
            canvas.ledger["source_order"]=order
            canvas.ledger["intersections"]=[{"mask":m,"original":b.get(m,0),"priority":a.get(m,0)} for m in masks]
            if not masks:canvas.text(11,40,"No observed source memberships",7.8)
            return True
        return False
    draw_venn(canvas,before,11,34,config,heading="Original annotations")
    draw_venn(canvas,after,118,34,config,heading="Prioritized CDSs")
    canvas.text(107.5,46,"DNMB",7.8,ha="center",bold=True,color=TEAL)
    canvas.text(107.5,50,"prioritization",7,ha="center",color=TEAL)
    canvas.ax.annotate("",xy=(116,56),xytext=(100,56),arrowprops={"arrowstyle":"-|>","color":TEAL,"lw":.7})
    canvas.text(11,75,"Areas schematic",7,color="#6F8185")
    if not regions(after).get("0000",0):
        canvas.text(161.5,75,config.get("selection_label","R subunits + nuclease effectors"),7,ha="center",color="#6F8185")
    return True


def intersection_pages(before, after, fonts, config, page_factory):
    """Complete exact counted-mask fallback, including selected zero masks."""
    order=before.get("source_order",[])
    if after.get("source_order",[])!=order:
        raise ValueError("Before/after overlap source orders must agree")
    b,a=regions(before),regions(after)
    masks=sorted(set(b)|set(a),key=lambda m:(-sum(int(c) for c in m),m))
    masks=[m for m in masks if b.get(m,0) or a.get(m,0)]
    # A source column is never silently dropped; wide source sets use a wider page.
    width=max(215.9,40+len(order)*16+50)
    pages=[]
    for start in range(0,max(1,len(masks)),45):
        c=page_factory(width,279.4,"overlap")
        c.text(8.45,10,"DNMB prioritization: source intersections",11,bold=True)
        c.text(8.45,17,scope_note(before,after) or "Exact shared memberships; original source definitions preserved",7,role="coverage_note")
        for j,name in enumerate(order):
            label=fonts.wrap(alias(name,config),15,7.5,True)
            c.text(23+j*16,29,label,7.5,ha="center",bold=True,color=source_color(name,config))
        xcount=width-39
        c.text(xcount,29,"Original",7.5,ha="right",bold=True)
        c.text(width-9,29,_priority_heading(after,"Prioritized"),7.5,ha="right",bold=True)
        c.ledger["selection_assessment"]=after.get("selection_assessment")
        page_masks=masks[start:start+45]
        for i,m in enumerate(page_masks):
            yy=38+i*4.8
            if i%2==0:c.rect(8.45,yy-2.4,width-16.9,4.8,"#F6F8F8")
            positions=[23+j*16 for j,v in enumerate(m) if v=="1"]
            if positions:c.line(min(positions),yy,max(positions),yy,"#9CAAAA",.5)
            for j,v in enumerate(m):
                c.circle(23+j*16,yy,.85,source_color(order[j],config) if v=="1" else "#DFE4E4")
            c.text(xcount,yy,f"{b.get(m,0):,}",7.8,ha="right",role="numeric",ref={"mask":m,"endpoint":"original"})
            c.text(width-9,yy,f"{a.get(m,0):,}",7.8,ha="right",role="numeric",ref={"mask":m,"endpoint":"priority"})
            c.ledger.setdefault("intersections",[]).append({"mask":m,"original":b.get(m,0),"priority":a.get(m,0)})
        c.text(8.45,269,"Filled dots: source membership. An all-open row is outside all original source sets.",7)
        c.ledger["source_order"]=order
        pages.append(c)
    return pages

"""The approved five-stage/landscape composition, driven by saved identities.

The constants below are physical drawing dimensions, never cohort quantities.
Source memberships, group composition, glyph areas and labels come from the
manuscript model. No annotation or biological activity is inferred here.
"""
from __future__ import annotations

from collections import OrderedDict
from itertools import combinations
import math

from .drawing import alias, number, status_kind, stable_color, display_label
from .comparison import ordered_families, visible_rows, metric_state

TEAL = "#287C78"
PURPLE = "#7953A2"
REFERENCE = "#7B3294"
INNER = "#65408C"
ORANGE = "#DA6D25"
GRID = "#E3E6E7"
PINK = "#FBE0EB"
TARGET = "#F3BCD2"
PALETTE = {"Modification":"#B29BCB", "Nuclease":"#D9878C", "Cas":"#A8BF72",
           "Abortive infection":"#7FAED0", "Other mechanisms":"#C6AD75", "Unresolved":"#D2D4D6",
           "Former R-M residual":"#AAB2B8"}
CELL = 3.15
ROW = 3.05
GAP = .8925
MAJOR = 4.2
BODY_TOP = 90.
NUM = 8.82
ROLE_COLUMNS = [
    ("I_R","role","I","R",1), ("I_M","role","I","M",1), ("I_S","role","I","S",1),
    ("I_MS","context","I","MS",1.5), ("I_RMS","context","I","RMS",2),
    ("II_R","role","II","R",1.5), ("II_M","role","II","M",1),
    ("III_R","role","III","R",1), ("III_M","role","III","M",1),
    ("III_M","context","III","M",1.5), ("III_MR","context","III","MR",1.5),
    ("IV_R","role","IV","R",1), ("IV_A","role","IV","A",1)]


def radius(n):
    return .84 * (.48 if n<2 else .68 if n<4 else .91 if n<8 else 1.15 if n<16 else 1.4)


def _text(c,x,y,label,pt=NUM,**kw):
    from .manuscript_text import draw_manuscript_text
    label=str(label)
    if "\n" in label:
        lines=label.splitlines();step=pt*25.4/72*1.12
        artists=[draw_manuscript_text(c,x,y+(i-(len(lines)-1)/2)*step,line,pt,**kw) for i,line in enumerate(lines)]
        return artists[-1]
    return draw_manuscript_text(c,x,y,label,pt,**kw)


def _category(f):
    cat=f.get("display_group",f.get("category","Unresolved"))
    return {"CRISPR-Cas":"Cas","Other":"Other mechanisms","Former R-M residual":"Unresolved"}.get(cat,cat)


def _color(cat,config):
    cat={"CRISPR-Cas":"Cas","Other":"Other mechanisms"}.get(cat,cat)
    return config.get("category_colors",{}).get(cat,PALETTE.get(cat,stable_color(cat)))


def _cell(c,x,y,w,h,fill=None):
    c.rect(x,y,w,h,fill or "none",GRID,.165)


def _value(c,x,y,w,value,state,ref,blank_zero=False):
    label=number(value,state)
    if blank_zero and value==0 and status_kind(state)=="complete":label=""
    if label:_text(c,x+w-.65,y,label,ha="right",role="numeric",ref=ref)


def _mini_bar(c,x,y,w,value,state,scale,color,ref,reserve_width=None):
    if value and status_kind(state) in ("complete","partial"):
        # Keep one linear scale while reserving measured digits and a clear gap.
        capacity=min(w-3.4,w-1.35-(reserve_width if reserve_width is not None else c.fonts.width(str(max(scale,value)),NUM)))
        width=max(0.,capacity)*value/max(scale,value)
        c.rect(x+.5,y-.675,width,1.35,color)
        c.ledger["bars"].append({"bounds_mm":[x+.5,y-.675,width,1.35],"value":value,
          "denominator":scale,"status":status_kind(state),"source_ref":ref,"kind":"atlas_source"})
    _value(c,x,y,w,value,state,ref,True)


def _fraction_layout(fonts,records,kind):
    """One measured label reservation/linear bar track for the whole column."""
    nominal,preferred,left = {"intersection":(17.85,5.,.5),
                              "retained":(11.76,3.25,.4),
                              "noncanonical":(14.49,3.7,.4)}[kind]
    labels=[(number(n,state),number(den,denstate or state))
            for n,den,state,denstate in records
            if (n is not None and status_kind(state) in ("complete","partial"))
            or (den is not None and status_kind(denstate or state) in ("complete","partial"))]
    numerator=max((fonts.width(n,NUM,kind=="retained") for n,_ in labels),default=0.)
    denominator=max((fonts.width(d,NUM) for _,d in labels),default=0.)
    slash=fonts.width("/",NUM)
    # Keep the approved font and physical column width for ordinary counts.
    # A four-digit partial fraction can require a wider column; never hide its
    # bar or let the text spill into adjacent evidence columns to fit it.
    # Retained bars keep their original visual scale: make room for wider
    # labels by widening the column, never by reducing that 3.25 mm track.
    right=.7;gap=.18;track_gap=.24;minimum_track=preferred if kind=="retained" else .5
    reserved=numerator+denominator+slash+2*gap
    width=max(nominal,left+minimum_track+track_gap+reserved+right)
    track=preferred if kind=="retained" else min(preferred,width-left-right-reserved-track_gap)
    return dict(width_mm=width,track_width_mm=track,numerator_width_mm=numerator,
                denominator_width_mm=denominator,slash_width_mm=slash,
                left_padding_mm=left,right_padding_mm=right,text_gap_mm=gap,
                track_gap_mm=track_gap,font_pt=NUM)


def _fraction(c,x,y,w,n,den,state,kind,scale=None,denstate=None,layout=None):
    denstate=state if denstate is None else denstate
    num_known=n is not None and status_kind(state) in ("complete","partial")
    den_known=den is not None and status_kind(denstate) in ("complete","partial")
    if not num_known and not den_known:
        _text(c,x+w-.7,y,"—",ha="right",role="numeric")
        return
    layout=layout or _fraction_layout(c.fonts,[(n,den,state,denstate)],kind)
    if w+1e-8<layout["width_mm"]:
        raise ValueError("Fraction counts exceed available column width; expand the atlas fraction column")
    track=layout["track_width_mm"];left=x+layout["left_padding_mm"]
    if kind=="intersection":
        base="#C8DDDA";fill="#68A19C"
        bw=track if den_known and den else 0;fw=track*n/den if num_known and den_known and den else 0
    elif kind=="retained":
        base=TEAL;fill="#124F4C"
        bw=track*den/max(scale or den,den,1) if den_known else 0
        fw=track*n/max(scale or den,den,1) if num_known and den_known else 0
    else:
        base="#F2D6C3";fill=PURPLE
        bw=track if den_known and den else 0;fw=track*n/den if num_known and den_known and den else 0
    if bw:c.rect(left,y-.725,bw,1.45,base)
    if fw:c.rect(left,y-.725,fw,1.45,fill)
    c.ledger["bars"].append({"bounds_mm":[left,y-.725,bw,1.45],"value":n,"denominator":den,
        "status":status_kind(state),"denominator_status":status_kind(denstate or state),
        "source_ref":{"endpoint":kind},"kind":"atlas_fraction",
        "track_width_mm":track,"filled_width_mm":fw,"scale":scale if kind=="retained" else den})
    den_right=x+w-layout["right_padding_mm"]
    slash_center=den_right-layout["denominator_width_mm"]-layout["text_gap_mm"]-layout["slash_width_mm"]/2
    num_right=slash_center-layout["slash_width_mm"]/2-layout["text_gap_mm"]
    first_artist=len(c.text_artists)
    for i,(label,anchor,ha) in enumerate(zip((number(n,state),"/",number(den,denstate or state)),
                                           (num_right,slash_center,den_right),("right","center","right"))):
        _text(c,anchor,y,label,ha=ha,bold=kind=="retained" and i==0,role="numeric",
              owner=[x,y-ROW/2,w,ROW],ref={"endpoint":kind,"part":("numerator","separator","denominator")[i]})
    # Commas in large counts extend below the usual capital-height anchor.
    # Fit the real ink bounds inside the same row, moving the whole fraction
    # baseline together without reducing the approved type size.
    artists=c.text_artists[first_artist:]
    renderer=c.fig.canvas.get_renderer();bounds=[]
    for artist,_ in artists:
        box=artist.get_window_extent(renderer)
        bounds.append((c.ax.transData.inverted().transform((box.x0,box.y1))[1],
                       c.ax.transData.inverted().transform((box.x1,box.y0))[1]))
    low=y-ROW/2+.04-min(a for a,_ in bounds)
    high=y+ROW/2-.04-max(b for _,b in bounds)
    if low>high:
        raise ValueError("Fraction text exceeds the atlas row height; increase the row size")
    shift=min(max(0.,low),high)
    if shift:
        for artist,meta in artists:
            artist.set_y(artist.get_position()[1]+shift)
            meta["anchor_mm"][1]+=shift;meta["baseline_mm"][1]+=shift


def _family_glyph(c,x,y,row,family,config):
    n=row.get("N");a=row.get("associated_CDS");e=row.get("effector_CDS")
    state=metric_state(row,"inventory");ast=metric_state(row,"associated");est=metric_state(row,"effector")
    ref={"sample_id":row.get("sample_id"),"family_id":family["id"]}
    outer=inner=0.
    if n is None or status_kind(state)=="unassessed":
        _text(c,x,y,"—",6.7,ha="center",role="missing",ref=ref)
    elif n:
        outer=radius(n);category=family.get("category","Unresolved")
        c.circle(x,y,outer,_color(category,config))
        if a and status_kind(ast) in ("complete","partial"):
            c.circle(x,y,outer+.1,None,"white",.7875)
            c.circle(x,y,outer+.1,None,PURPLE,.4125)
        if e is not None and status_kind(est) in ("complete","partial"):
            inner=outer*math.sqrt(e/n)
            if inner:c.circle(x,y,inner,INNER)
    c.ledger["glyphs"].append(dict(ref,total=n,associated=a,effector=e,
      inventory_status=state,associated_status=ast,effector_status=est,
      outer_radius_mm=outer,inner_radius_mm=inner,inner_area_fraction=e/n if n and inner else 0,
      bounds_mm=[x-CELL/2,y-ROW/2,CELL,ROW],category=family.get("category"),activity_inferred=False))


def _role_glyph(c,x,y,w,record,key,kind):
    n=record.get("count");state=record.get("status")
    if n is None or status_kind(state)=="unassessed":
        # Status marks need visible separation even in adjacent one-CDS cells.
        mark_pt=min(6.7,6.7*(CELL-.45)/c.fonts.width("—",6.7))
        _text(c,x+w/2,y,"—",mark_pt,ha="center");return
    if not n:
        if kind=="context":_text(c,x+w/2,y,"0",6.7,ha="center",role="numeric")
        return
    r=radius(n);cx=x+w/2;numeric=key in ("II_R","IV_R") or kind=="context"
    if kind=="context":
        # M/MS-only composition is not evidence of biological orphanhood.
        methylation_only=key in ("I_MS","III_M");cx-=.9
        c.circle(cx,y,r,"#D8CAE5" if methylation_only else "#BCD9D4","#81669C" if methylation_only else "#292929",.6)
        if record.get("R_support_all") is True and not methylation_only:c.circle(cx,y,r+.254,None,TEAL,.7125)
        if record.get("reference_any") is True:c.circle(cx,y,.34,REFERENCE)
        _text(c,x+w/2+1.05,y,number(n,state),ha="center",role="numeric")
    else:
        color=TEAL if key.endswith("_R") else "#8EADA8" if key.endswith("_M") else "#B4C8C3"
        if key.startswith("IV"):color="#54666D" if key.endswith("R") else "#89969A"
        if key=="II_R":cx-=.9
        from matplotlib.colors import to_rgba
        c.circle(cx,y,r,to_rgba(color,.25 if numeric else .78))
        if record.get("support_count",0) and status_kind(record.get("support_status")) in ("complete","partial"):
            c.circle(cx,y,r+.104,None,TEAL,.7)
        if numeric:_text(c,x+w/2+(1.05 if key=="II_R" else 0),y,number(n,state),ha="center",role="numeric")
    c.ledger.setdefault("rm_symbols",[]).append({"field":key,"kind":kind,"count":n,
        "record":record,"center_mm":[cx,y],"radius_mm":r,"count_unit":"groups" if kind=="context" else "CDSs"})


def _mask_order(n):
    return ["".join("1" if j in selected else "0" for j in range(n))
            for k in range(1,n+1) for selected in combinations(range(n),k)]


def _remap_masks(row,old,new):
    return {"".join(item["mask"][old.index(s)] for s in new):item
            for item in row.get("noncanonical",{}).get("masks",[]) if len(item["mask"])==len(old)}


def _layout(families,sources,refs,config,fractions=None):
    fractions=fractions or {}
    panels=OrderedDict();x=6.04875
    def add(name,width,gap=GAP):
        nonlocal x
        panels[name]={"left":x,"right":x+width};x+=width+gap
    add("RM_callers",len(sources)*7.56)
    if refs:add("REBASE",len(refs)*7.56)
    add("RM_pool",fractions.get("intersection",{}).get("width_mm",17.85))
    add("RM_integrated",sum(z[4]*CELL for z in ROLE_COLUMNS))
    add("Retained_summary",fractions.get("retained",{}).get("width_mm",11.76),MAJOR)
    add("NonRM_callers",len(sources)*7.56,MAJOR)
    masks=_mask_order(len(sources));exclusive=masks[:-1]
    add("cases",len(exclusive)*CELL)
    add("NonRM_summary",fractions.get("noncanonical",{}).get("width_mm",14.49),MAJOR)
    groups=OrderedDict()
    for f in families:groups.setdefault(_category(f),[]).append(f)
    fpos={}
    for cat,fs in groups.items():
        width=len(fs)*CELL
        add(cat,width,MAJOR if cat in ("Modification","Nuclease") else GAP)
        for i,f in enumerate(fs):fpos[f["id"]]=panels[cat]["left"]+i*CELL
    data_right=max(p["right"] for p in panels.values());sampleleft=data_right+2.6825
    return panels,groups,fpos,masks,sampleleft,max(447.,data_right+91.)


def _badge(c,x,y,n,color):
    from matplotlib.colors import to_rgba
    c.circle(x,y,2.1,to_rgba(color,.15),color,.6)
    c.ax.patches[-1].set_gid(f"atlas-stage-badge-{n}")
    artist=_text(c,x,y,n,10.5,bold=True,ha="center",color=color,
                 role="stage_badge",ref={"stage_id":n})
    artist.set_gid(f"atlas-stage-number-{n}")


def _stage_specs(panels,groups):
    """Approved titles, attached to observed column groups rather than slots."""
    definitions=[
        (1,["RM_callers","REBASE","RM_pool","RM_integrated","Retained_summary"],
         ["R-M evidence and curation"],TEAL,2.1),
        (2,["NonRM_callers","cases","NonRM_summary"],
         ["Methylation-associated","evidence"],"#829764",2.1),
        (3,["Modification"],["Non-canonical","DNA modification","system"],PALETTE["Modification"],2.1),
        (4,["Nuclease"],["Nuclease-associated","evidence"],PALETTE["Nuclease"],21.5),
        (5,[k for k in groups if k not in ("Modification","Nuclease")],
         ["Complementary defense priorities"],"#52626A",17.)]
    stages=[];absent=[]
    for stage_id,keys,lines,color,offset in definitions:
        keys=[k for k in keys if k in panels and panels[k]["right"]>panels[k]["left"]]
        if not keys:
            absent.append({"stage_id":stage_id,"title":" ".join(lines),
                           "reason":"no_displayed_column_group"});continue
        left=panels[keys[0]]["left"];right=panels[keys[-1]]["right"]
        # A clipped right-edge badge puts the entire title over the next group.
        # Preserve the original anchor only while its leading text remains in
        # this group; otherwise use a left anchor. A one-column group stacks
        # the title below its badge so both still start over that actual group.
        badge_x=left+offset;anchor_mode="original"
        if badge_x+3.7>=right-.5:
            badge_x=left+min(2.1,(right-left)/2);anchor_mode="left"
        title_x=badge_x+3.7;title_dy=0.
        if title_x>=right-.5:
            title_x=left+.2;title_dy=5.5;anchor_mode="stacked_left"
        stages.append(dict(stage_id=stage_id,title=" ".join(lines),lines=lines,color=color,
                           group_keys=keys,group_span_mm=[left,right],badge_x=badge_x,y=64.,
                           title_x=title_x,title_offset_y_mm=title_dy,anchor_mode=anchor_mode))
    return stages,absent


def _draw_stages(c,stages):
    for stage in stages:
        sid=stage["stage_id"];x=stage["badge_x"];y=stage["y"]
        _badge(c,x,y,sid,stage["color"])
        for i,line in enumerate(stage["lines"]):
            artist=_text(c,stage["title_x"],y+stage["title_offset_y_mm"]+i*4.7,line,12,bold=True,
                         color="black" if sid==5 else stage["color"],role="stage_title",
                         ref={"stage_id":sid,"line":i})
            artist.set_gid(f"atlas-stage-title-{sid}-{i}")


def _boxes_overlap(a,b,gap=.6):
    x,y,w,h=a;u,v,s,t=b
    return min(x+w,u+s)+gap>max(x,u) and min(y+h,v+t)+gap>max(y,v)


def _header_plan(fonts,width,panels,groups,fpos,sources,refs,config):
    """Keep original anchors unless measured ink requires a second tier.

    Family/sample counts do not select a layout. Circle bounds participate in
    the same clearance test as title ink, including narrow observed groups.
    """
    from .drawing import Canvas
    stages,absent=_stage_specs(panels,groups)
    probe=Canvas(width,200,fonts,"header_measure","header_measure")
    _headers(probe,panels,groups,fpos,sources,refs,config,
             plan={"stages":stages,"absent":absent,"body_offset_mm":0.},record=False)
    measured=probe.measure()
    boxes={s["stage_id"]:[[s["badge_x"]-2.1,61.9,4.2,4.2]] for s in stages}
    for t in measured:
        if t["role"]=="stage_title":boxes[t["source_ref"]["stage_id"]].append(t["bounds_mm"])
    settled=[]
    for stage in stages:
        original=boxes[stage["stage_id"]];shift=0.
        while True:
            moved=[[x,y+shift,w,h] for x,y,w,h in original]
            collisions=[(a,b) for a in moved for b in settled if _boxes_overlap(a,b)]
            if not collisions:break
            shift+=max(b[1]+b[3]+1.-a[1] for a,b in collisions)
        stage["y"]+=shift;stage["tier_offset_mm"]=shift
        stage["layout_boxes_mm"]=moved;settled.extend(moved)
    # Only lower leaf headers/grid. The b title and stage anchors remain tied
    # to their own group; this avoids a detached full-page numbered key.
    leaves=[t["bounds_mm"] for t in measured
            if t["role"] not in ("stage_title","stage_badge") and t["anchor_mm"][1]>60]
    remaining=[panels[k] for k in groups if k not in ("Modification","Nuclease")]
    # A tall rotated category title must also stay below the b-title/workflow
    # floor even if its narrow x interval misses every stage-title box.
    dy=max(0.,61.5-min((b[1] for b in leaves),default=61.5))
    while True:
        moved=[[x,y+dy,w,h] for x,y,w,h in leaves]
        collisions=[(a,b) for a in moved for b in settled if _boxes_overlap(a,b)]
        if not collisions:break
        dy+=max(b[1]+b[3]+.8-a[1] for a,b in collisions)
    bracket=None
    if remaining:
        ll=remaining[0]["left"];rr=remaining[-1]["right"]
        in_span=lambda b: b[0]+b[2]>ll-.1 and b[0]<rr+.1
        # The bracket is an independent drawing object, not a text box. Its
        # horizontal stroke must sit below overlapping stage ink and above
        # every leaf label, including tall rotated category names. Moving it
        # with the leaf labels would preserve an existing strike-through.
        bracket_y=max(73.,max((b[1]+b[3]+1. for b in settled if in_span(b)),default=73.))
        dy=max(dy,bracket_y+.8-min((b[1] for b in leaves if in_span(b)),default=bracket_y+.8))
        bracket={"stage_id":5,"bounds_mm":[ll,bracket_y,rr-ll,80.+dy-bracket_y],
                 "segments_mm":[[ll,bracket_y,rr,bracket_y],
                                [ll,bracket_y,ll,80.+dy],[rr,bracket_y,rr,80.+dy]],
                 "clearance_mm":.8}
    probe.fig.clear()
    return {"stages":stages,"absent":absent,"body_offset_mm":dy,
            "bracket":bracket,
            "mode":"group_anchored_staggered" if dy or any(s["tier_offset_mm"] for s in stages) else "original_group_anchors"}


def _headers(c,panels,groups,fpos,sources,refs,config,dy=0,plan=None,record=True):
    # Leaf labels may move down; the panel title retains its original anchor.
    def text(x,y,*a,**kw):return _text(c,x,y+dy,*a,**kw)
    _text(c,3,58.2,"b",14,bold=True)
    _text(c,10,58.2,config.get("atlas_title") or "Comparative defense landscape across genomes",13,bold=True)
    if plan is None:
        stages,absent=_stage_specs(panels,groups)
        plan={"stages":stages,"absent":absent,"body_offset_mm":dy,"mode":"original_group_anchors"}
    _draw_stages(c,plan["stages"])
    remaining=[panels[k] for k in groups if k not in ("Modification","Nuclease")]
    if remaining:
        ll=remaining[0]["left"];rr=remaining[-1]["right"]
        bracket=plan.get("bracket") or {"segments_mm":[[ll,73+dy,rr,73+dy],
            [ll,73+dy,ll,80+dy],[rr,73+dy,rr,80+dy]]}
        for i,segment in enumerate(bracket["segments_mm"]):
            c.line(*segment,"#AAB7BC",.5625)
            c.ax.lines[-1].set_gid(f"atlas-stage-bracket-5-{i}")
    p=panels
    text((p["RM_callers"]["left"]+p["RM_callers"]["right"])/2,80,"R-M CDSs in",9.24,bold=True,ha="center")
    if refs:
        text((p["REBASE"]["left"]+p["REBASE"]["right"])/2,80,"Hits",9.24,bold=True,ha="center")
    for name,keys in (("RM_callers",sources),("REBASE",refs),("NonRM_callers",sources)):
        if not keys:continue
        for i,source in enumerate(keys):text(p[name]["left"]+(i+.5)*7.56,85,alias(source,config),8.82,ha="center")
    pp=p["RM_pool"];text((pp["left"]+pp["right"])/2,80,"Intersection",9.24,bold=True,ha="center")
    text((pp["left"]+pp["right"])/2,85,"(n/R-M pool)",9.6,ha="center")
    pp=p["RM_integrated"];text((pp["left"]+pp["right"])/2,72.4,"R-M subunits and contexts",10,bold=True,ha="center")
    x=pp["left"];spans=OrderedDict()
    for key,kind,typ,label,units in ROLE_COLUMNS:
        w=units*CELL;spans.setdefault(typ,[x,x]);spans[typ][1]=x+w
        if kind=="context":c.rect(x,80.9+dy,w,8.8,"#F6F8F7")
        angle=55 if kind=="context" and len(label)>1 else 0
        text(x+w/2,85.5 if angle else 87,label,8.82,ha="center",rotation=angle)
        x+=w
    for typ,(left,right) in spans.items():
        text((left+right)/2,75.5,typ,8.8,bold=True,ha="center")
        c.line(left+.3,77.5+dy,right-.3,77.5+dy,TEAL,.525)
    x=pp["left"]
    for typ in ("I","II","III","IV"):
        cols=[z for z in ROLE_COLUMNS if z[2]==typ]
        for kind in ("role","context"):
            selected=[z for z in cols if z[1]==kind]
            if selected:
                width=sum(z[4]*CELL for z in selected)
                if typ!="IV":text(x+width/2,79.5,"CDS" if kind=="role" else "Groups",7.3,ha="center")
                x+=width
    pp=p["Retained_summary"]
    for yy,label,bold,size in ((76.4,"Assigned",True,9.24),(80.6,"CDS",True,9.24),(84.9,"(R/total)",False,9.6)):
        text((pp["left"]+pp["right"])/2,yy,label,size,bold=bold,ha="center")
    pp=p["NonRM_callers"]
    text((pp["left"]+pp["right"])/2,76,"Non-canonical",9.24,bold=True,ha="center")
    text((pp["left"]+pp["right"])/2,80,"Annotated CDSs",9.24,bold=True,ha="center")
    cx=(p["cases"]["left"]+p["NonRM_summary"]["right"])/2
    text(cx,76,"Exclusive CDSs",9.24,bold=True,ha="center")
    text(cx,80,"R-M-excluded union",9.6,ha="center")
    pp=p["NonRM_summary"];text((pp["left"]+pp["right"])/2,85,"(n/pool)",9.6,ha="center")
    for cat,fs in groups.items():
        pp=p[cat];cx=(pp["left"]+pp["right"])/2
        label="Other" if cat=="Other mechanisms" else cat
        if cat=="Cas" or c.fonts.width(label,9.6,True)>pp["right"]-pp["left"]:
            text(cx,85.8,label,9.6,bold=True,rotation=90,ha="left")
        else:text(cx,80.8,label,9.6,bold=True,ha="center")
    for name in ("RM_callers","REBASE","RM_pool","Retained_summary","NonRM_callers","cases","NonRM_summary"):
        if name not in p:continue
        pp=p[name];col="#C8DDDA" if name=="RM_pool" else PURPLE if name=="NonRM_summary" else ORANGE if name in ("cases","NonRM_callers") else TEAL
        c.rect(pp["left"],87.1+dy,pp["right"]-pp["left"],1.7,col)
    for cat,fs in groups.items():
        for f in fs:
            xx=fpos[f["id"]];segments=f.get("band_segments")
            primary=_color(cat,config)
            c.rect(xx,87.1+dy,CELL,1.7,primary)
            tags=f.get("category_tags",[])
            if isinstance(tags,str):tags=[tags]
            secondary=[t for t in tags if t!=cat]
            if f.get("category")!=cat and f.get("category") in PALETTE:secondary.insert(0,f["category"])
            if segments and len(segments)>1:
                yy=87.1+dy;total=sum(z["fraction"] for z in segments)
                for z in segments:
                    h=1.7*z["fraction"]/total;c.rect(xx,yy,CELL,h,_color(z["category"],config));yy+=h
            elif secondary:c.rect(xx,87.95+dy,CELL,.85,_color(secondary[0],config))
            overlay_secondary=f.get("secondary_category")
            if overlay_secondary:
                color=_color(overlay_secondary,config);yy=88.02+dy
                if f.get("secondary_style") in ("striped","proposed_context"):
                    c.rect(xx,yy,CELL,.78,"white")
                    for offset in (.15,1.25,2.35):c.rect(xx+offset*CELL/3.2,yy,.65*CELL/3.2,.78,color)
                else:c.rect(xx,yy,CELL,.78,color)
            c.ledger["bands"].append({"family_id":f["id"],"primary":cat,"secondary":secondary,
                "bounds_mm":[xx,87.1+dy,CELL,1.7]})
        for left,right in zip(fs,fs[1:]):
            if left.get("partner_group") and left.get("partner_group")==right.get("partner_group"):
                c.line(fpos[left["id"]]+.15,86.7+dy,fpos[right["id"]]+CELL-.15,86.7+dy,"#9AA9AA",.375)
    if record:
        measured=c.measure();records=[]
        for stage in plan["stages"]:
            sid=stage["stage_id"]
            records.append({**stage,"badge_bounds_mm":[stage["badge_x"]-2.1,stage["y"]-2.1,4.2,4.2],
                "text_bounds_mm":[t["bounds_mm"] for t in measured if t["role"]=="stage_title" and t["source_ref"]["stage_id"]==sid],
                "font_pt":12,"badge_radius_mm":2.1})
        c.ledger["section_headings"]={"mode":plan["mode"],"heading_count":len(records),
            "body_offset_mm":dy,"headings":records,"absent":plan["absent"],
            "bracket":plan.get("bracket")}


def _strain_parts(row):
    label=str(row.get("label",row["sample_id"]));display=row.get("display",{})
    taxon=display.get("taxon_label");strain=display.get("strain_label")
    if taxon is None:
        if not row.get("genus") and not row.get("species"):
            return "", label
        bits=label.split();cut=1 if len(bits)<2 or bits[1] in ("sp.","genomosp.") else 2
        taxon=" ".join(bits[:cut]);strain=" ".join(bits[cut:])
    return str(taxon), str(strain or "")


def _strain(c,x,y,row):
    taxon,strain=_strain_parts(row)
    if not taxon:
        _text(c,x,y,strain,NUM,bold=bool(row.get("highlight")))
        return
    bold=bool(row.get("highlight"));t=_text(c,x,y,taxon,NUM,bold=bold,italic=True)
    # Exact italic advance is measured with the selected font, not a character count.
    from matplotlib.textpath import TextToPath
    width=TextToPath().get_text_width_height_descent(taxon,t.get_fontproperties(),False)[0]*25.4/72
    if strain:_text(c,x+width+.8,y,strain,NUM,bold=bold)


def _absence_cell(row,fid):
    cell=row.get("cells",{}).get(fid,{})
    return fid in row.get("display",{}).get("absence_family_ids",[]) and cell.get("N")==0 and status_kind(metric_state(cell,"inventory"))=="complete"


def _marker_selected(row,marker):
    return marker is not None and row.get("display",{}).get(marker["field"]) is True


def _taxonomy_legend_plan(fonts,model,rows,groups,fpos,config,taxonomy_colors,available):
    """Keep the original two-column key when it fits measured italic text.

    Long entries widen their own column. If two columns cannot fit, stack
    entries at the same readable size and reserve their actual vertical space.
    """
    from matplotlib import font_manager
    from matplotlib.textpath import TextToPath
    genus=OrderedDict();species=OrderedDict()
    for row in rows:
        disp=row.get("display",{});g=row.get("genus");s=row.get("species")
        if g:genus.setdefault(g,{"label":disp.get("genus_label",g),"color":disp.get("genus_color",stable_color(g))})
        if s:species.setdefault((g,s),{"label":disp.get("species_label",s),"color":disp.get("species_color",stable_color(s))})
    declared=model.get("taxonomy_legend",{})
    if declared.get("genera"):genus=OrderedDict((z["key"],z) for z in declared["genera"] if z["key"] in genus)
    if declared.get("species"):
        present_species={r.get("species") for r in rows}
        species=OrderedDict((z["key"],z) for z in declared["species"] if z["key"] in present_species)
    sections=[]
    for title,entries,field in (("Genus",genus,"genus"),("Species",species,"species")):
        items=[]
        for key,entry in entries.items():
            color_key=key[-1] if isinstance(key,tuple) else key
            runs=entry.get("label_runs")
            if runs is None:
                runs=([{"text":entry["abbreviation"],"italic":True},{"text":"=","italic":False},
                       {"text":entry["label"],"italic":True}] if title=="Genus" and entry.get("abbreviation")
                      else [{"text":entry["label"],"italic":entry.get("italic",True)}])
            items.append(dict(key=str(key),runs=runs,color=taxonomy_colors[field].get(color_key) or entry.get("color") or stable_color(str(color_key))))
        sections.append((title,items,3.,5. if title=="Genus" else 6.))
    extras=[]
    if groups.get("Other mechanisms"):extras.append(("Other",_color("Other mechanisms",config)))
    if config.get("highlight_label") and any(r.get("highlight") for r in rows):extras.append((config["highlight_label"],TARGET))
    if groups.get("Unresolved"):extras.append(("Mechanism not assigned yet",_color("Unresolved",config)))
    if any(_absence_cell(r,fid) for r in rows for fid in fpos):extras.append(("No "+str(model.get("absence_label","saved family")),"#D0CFE5"))
    marker=config.get("sample_marker")
    if any(_marker_selected(r,marker) for r in rows):extras.append((marker["label"],marker["color"]))
    sections.append(("Annotations",[dict(key=label,runs=[dict(text=label,italic=False)],color=color) for label,color in extras],4.,4.))
    plan=dict(sections=[],bottom_mm=3.4);y=3.4
    for title,items,step,after in sections:
        if not items:continue
        for item in items:
            width=0.
            for run in item["runs"]:
                prop=font_manager.FontProperties(family=fonts.pick(run["text"]),size=7.8,
                    style="italic" if run.get("italic") else "normal")
                width+=TextToPath().get_text_width_height_descent(run["text"],prop,False)[0]*25.4/72+.65
            # Four mm separates each swatch/text; two mm protects the next
            # swatch from italic overhang and gives a visible inter-column gap.
            item["width_mm"]=4.+width+2.
        columns=2 if len(items)>1 else 1
        widths=[max(40.,max(item["width_mm"] for item in items[j::columns])) for j in range(columns)]
        if sum(widths)>available:
            columns=1;widths=[max(item["width_mm"] for item in items)]
        if sum(widths)>available:
            raise ValueError("Taxonomy legend labels exceed atlas_page_size; increase the page limit or provide shorter display labels")
        for i,item in enumerate(items):
            item.update(x_mm=sum(widths[:i%columns]),y_mm=y+4+(i//columns)*step)
        plan["sections"].append(dict(title=title,y_mm=y,columns=columns,column_widths_mm=widths,items=items))
        y+=after+math.ceil(len(items)/columns)*step
        plan["bottom_mm"]=max(plan["bottom_mm"],items[-1]["y_mm"]+2.)
    return plan


def _legend(c,model,panels,groups,fpos,masks,sources,rows,bottom,sampleleft,config,taxonomy_colors,taxonomy_plan):
    y=bottom+3.4
    _text(c,6,y,"Coding",8.8,bold=True);_text(c,6,y+3,"sequences (CDSs)",8.8,bold=True)
    for i,n in enumerate((1,2,4,8,16)):
        cy=y+7+3*i;c.circle(7,cy,radius(n)/1.05,ORANGE)
        _text(c,10,cy,("1","2-3","4-7","8-15","16+")[i],7.8)
    _text(c,36,y,"Subunit groups",8.8,bold=True)
    reference_label="Any REBASE-linked group" if model.get("reference_sources")==["REBASE"] else "Any reference-linked group"
    for i,label in enumerate(("One symbol + group count","All groups R-supported",reference_label,"M(S)-only groups")):
        cy=y+4+3*i
        if i==2:c.circle(37,cy,.34,REFERENCE)
        else:
            c.circle(37,cy,.7,"#D8CAE5" if i==3 else "#BCD9D4","#81669C" if i==3 else "#292929",.6)
            if i==1:c.circle(37,cy,.95,None,TEAL,.7125)
        if i==0:_text(c,40.5,cy,"(n)",7.7,ha="center")
        _text(c,45,cy,label,7.7)
    _text(c,36,y+19,"Family evidence",8.8,bold=True)
    for i,label in enumerate(("Putative nuclease-associated CDSs","Effector CDSs")):
        cy=y+23+4*i;c.circle(37,cy,1.1,None,PURPLE,.45)
        if i:c.circle(37,cy,.635,INNER)
        _text(c,40,cy,label,7.7)
    # Seven source masks are a legend, never mistaken for additional genomes.
    ys=[bottom+32.4+4*i for i in range(len(sources))]
    if sources:
        for source,cy in zip(sources,ys):_text(c,panels["cases"]["left"]+.55,cy,alias(source,config),6.7,ha="right")
        for i,mask in enumerate(masks):
            xx=panels["NonRM_summary"]["left"]+1.8 if i==len(masks)-1 else panels["cases"]["left"]+(i+.5)*CELL
            active=[ys[k] for k,bit in enumerate(mask) if bit=="1"];color=PURPLE if i==len(masks)-1 else ORANGE
            if len(active)>1:c.line(xx,min(active),xx,max(active),color,.6)
            for bit,cy in zip(mask,ys):c.circle(xx,cy,.65 if bit=="1" else .3,color if bit=="1" else "#E3E7E8")
        _text(c,panels["NonRM_summary"]["left"]+5.5,ys[-1],"n/pool",8.4,bold=True)
        outside=[_remap_masks(row,model["source_order"],sources).get("0"*len(sources),{}) for row in model["rows"]]
        if any(item.get("count",0) for item in outside):
            from .manuscript_model import _intersection_state
            state=_intersection_state(item.get("status","unassessed") for item in outside)
            label=number(sum(item.get("count") or 0 for item in outside),state)
            _text(c,panels["cases"]["left"],ys[-1]+5,"Outside shown sources: "+label+" CDSs",7.7)
    from matplotlib.textpath import TextToPath
    c.ledger["taxonomy_legend_layout"]={**taxonomy_plan,"origin_mm":[sampleleft,bottom]}
    for section in taxonomy_plan["sections"]:
        _text(c,sampleleft,bottom+section["y_mm"],section["title"],8.8 if section["title"]=="Genus" else 8.1,bold=True)
        for item in section["items"]:
            xx=sampleleft+item["x_mm"];yy=bottom+item["y_mm"]
            c.rect(xx,yy-1.2,2.2,2.2,item["color"]);xx+=4
            for run in item["runs"]:
                artist=_text(c,xx,yy,run["text"],7.8,italic=bool(run.get("italic")),role="taxonomy_legend",
                    ref={"section":section["title"],"key":item["key"]})
                xx+=TextToPath().get_text_width_height_descent(run["text"],artist.get_fontproperties(),False)[0]*25.4/72+.65


def manuscript_atlas_pages(result,fonts,config,page_factory):
    from .manuscript_model import build_manuscript_atlas_model
    from .manuscript_workflow import draw_manuscript_workflow
    model=build_manuscript_atlas_model(result)
    from .display_scope import fraction_projection
    scoped=fraction_projection(result,config.get("rm_fraction_scope","assigned_I_IV"))
    if scoped is not None:
        for row in model["rows"]:row["retained"]=scoped[row["sample_id"]]
    sources=list(model.get("source_order",[]))
    desired=config.get("atlas_source_order",["PADLOC","DefenseFinder","DefensePredictor"])
    sources=[s for s in desired if s in sources]+[s for s in sources if s not in desired]
    refs=list(model.get("reference_sources",[]))
    from matplotlib.colors import to_rgba
    taxonomy_colors={field:{item["key"]:item.get("color") for item in model.get("taxonomy_legend",{}).get(plural,[])}
                     for field,plural in (("genus","genera"),("species","species"))}
    for row in model.get("rows",[]):
        for field in taxonomy_colors:
            explicit=row.get("display",{}).get(field+"_color")
            declared=taxonomy_colors[field].get(row.get(field))
            if explicit and declared and to_rgba(explicit)!=to_rgba(declared):
                raise ValueError("Conflicting taxonomy colors between sample display and legend")
            if explicit and not declared and row.get(field):taxonomy_colors[field][row[field]]=explicit
    families=ordered_families(visible_rows(model.get("families",[])),config)
    families=[f for f in families if f.get("category")!="R-M"]
    allrows=model.get("rows",[])
    fractions={}
    for kind,record_key,nkey,dkey,skey,dskey in (
        ("intersection","rm_original","intersection_CDS","input_pool_CDS","intersection_count_status","pool_status"),
        ("retained","retained","R_CDS","total_CDS","R_status","total_status"),
        ("noncanonical","noncanonical","intersection_CDS","pool_CDS","intersection_status","pool_status")):
        records=[row.get(record_key,{}) for row in allrows]
        fallback="intersection_status" if kind=="intersection" else "status"
        fractions[kind]=_fraction_layout(fonts,[(r.get(nkey),r.get(dkey),
            r.get(skey,r.get(fallback)),r.get(dskey,r.get(fallback))) for r in records],kind)
    # Reserve adjacent inter-panel whitespace for headings as well as numbers.
    # A wider fallback font must not overlap the separate reference heading.
    fractions["intersection"]["width_mm"] = max(
        fractions["intersection"]["width_mm"],
        fonts.width("Intersection",9.24,True)-2*GAP+.4,
        fonts.width("(n/R-M pool)",9.6)-2*GAP+.4)
    # Keep the Assigned heading clear of IV with wider fallback fonts.
    fractions["retained"]["width_mm"] = max(
        fractions["retained"]["width_mm"],
        fonts.width("Assigned",9.24,True)+fonts.width("IV",8.8,True)-2*CELL-2*GAP+.8)
    family_pages=[];block=[]
    maxwidth=config["atlas_page_size"][0]
    labelwidth=max((sum(fonts.width(part,NUM,bool(r.get("highlight"))) for part in _strain_parts(r))+2 for r in allrows),default=58.)
    for f in families:
        trial=block+[f]
        trial_layout=_layout(trial,sources,refs,config,fractions)
        required_width=max(trial_layout[-1],trial_layout[-2]+labelwidth+8)
        if block and (required_width>maxwidth or len(trial)>config["max_families_per_page"]):
            family_pages.append(block);block=[]
        block.append(f)
    family_pages.append(block)
    pages=[]
    for fs in family_pages:
        panels,groups,fpos,masks,sampleleft,width=_layout(fs,sources,refs,config,fractions)
        heading_plan=_header_plan(fonts,width,panels,groups,fpos,sources,refs,config)
        dy=heading_plan["body_offset_mm"]
        width=max(width,sampleleft+labelwidth+8)
        maxrows=max(1,min(config["max_rows_per_page"],int((config["atlas_page_size"][1]-BODY_TOP-dy-55)/ROW)))
        chunks=[];cursor=0
        # Budget the legend for each page's taxa, not the entire cohort. A
        # tall one-column key can reduce rows without shrinking any labels.
        while cursor<len(allrows):
            count=min(maxrows,len(allrows)-cursor)
            while count>1:
                trial=_taxonomy_legend_plan(fonts,model,allrows[cursor:cursor+count],groups,fpos,config,taxonomy_colors,width-sampleleft-6)
                if BODY_TOP+dy+count*ROW+max(50.,trial["bottom_mm"]+5.)<=config["atlas_page_size"][1]:break
                count-=1
            chunks.append(allrows[cursor:cursor+count]);cursor+=count
        if not chunks:chunks=[[]]
        for rows in chunks:
            top=BODY_TOP+dy;bottom=top+len(rows)*ROW
            longlabel=max((fonts.width(display_label(f),8.3) for f in fs),default=0.)
            taxonomy_plan=_taxonomy_legend_plan(fonts,model,rows,groups,fpos,config,taxonomy_colors,width-sampleleft-6)
            taxonomy_height=taxonomy_plan["bottom_mm"]+5.
            height=max(160.,bottom+max(50.,longlabel+7,taxonomy_height))
            if width>maxwidth or height>config["atlas_page_size"][1]:
                raise ValueError("Manuscript labels/evidence columns exceed atlas_page_size; increase the page limit or provide shorter display labels")
            c=page_factory(width,height,"atlas")
            c.ledger["atlas_style"]="manuscript"
            c.ledger["rm_fraction_scope"]=config.get("rm_fraction_scope","assigned_I_IV")
            c.ledger["rm_fraction_records"]=[dict(sample_id=row["sample_id"],**row["retained"]) for row in rows]
            c.ledger["panel_bounds"]={k:dict(v,top=top,bottom=bottom) for k,v in panels.items()}
            c.ledger["fraction_layout"]=fractions
            c.ledger["model_provenance"]=model.get("provenance",{})
            draw_manuscript_workflow(c,result,config,[0,0,447.,55.53])
            _headers(c,panels,groups,fpos,sources,refs,config,dy,heading_plan)
            _text(c,sampleleft,85+dy,"Strains",9.45,bold=True)
            scale_rm=max([25]+[n for r in allrows for n in r.get("rm_original",{}).get("source_counts",{}).values() if isinstance(n,int)])
            scale_ref=max([40]+[n for r in allrows for n in r.get("rm_original",{}).get("reference_counts",{}).values() if isinstance(n,int)])
            scale_nc=max([80]+[n for r in allrows for n in r.get("noncanonical",{}).get("source_counts",{}).values() if isinstance(n,int)])
            scale_ret=max([1]+[r.get("retained",{}).get("total_CDS") or 0 for r in allrows])
            reservations={}
            for panel,rootkey,countkey,statuskey,scale in (
                ("RM_callers","rm_original","source_counts","source_status",scale_rm),
                ("REBASE","rm_original","reference_counts","reference_status",scale_ref),
                ("NonRM_callers","noncanonical","source_counts","source_status",scale_nc)):
                widths=[fonts.width(str(scale),NUM)]
                for row in allrows:
                    record=row.get(rootkey,{})
                    for source,count in record.get(countkey,{}).items():
                        if count and status_kind(record.get(statuskey,{}).get(source)) in ("complete","partial"):
                            widths.append(fonts.width(number(count,record[statuskey][source]),NUM))
                reservations[panel]=max(widths)
            for ri,row in enumerate(rows):
                sid=row["sample_id"];y=top+ri*ROW;cy=y+ROW/2;display=row.get("display",{})
                c.ledger["cells"].append({"kind":"strain","sample_id":sid,"highlight":bool(row.get("highlight")),"bounds_mm":[sampleleft,y,58,ROW]})
                c.rect(panels["RM_callers"]["left"],y,sampleleft+58-panels["RM_callers"]["left"],ROW,PINK if row.get("highlight") else "white")
                if row.get("highlight"):c.rect(sampleleft-.6825,y,58,ROW,TARGET)
                for k,field in enumerate(("genus","species")):
                    if row.get(field):c.rect(2.35+k*1.68,y,1.365,ROW,display.get(field+"_color") or taxonomy_colors[field].get(row[field]) or stable_color(row[field]))
                marker=config.get("sample_marker")
                if _marker_selected(row,marker):
                    bounds=[sampleleft-2.6825,y,1.3125,ROW]
                    c.rect(*bounds,marker["color"])
                    c.ledger.setdefault("sample_annotations",[]).append({"sample_id":sid,**marker,"bounds_mm":bounds})
                _strain(c,sampleleft,cy,row)
                raw=row.get("rm_original",{});nc=row.get("noncanonical",{})
                for name,keys,data,status,scale,color in (
                    ("RM_callers",sources,raw.get("source_counts",{}),raw.get("source_status",{}),scale_rm,"#8EADA8"),
                    ("REBASE",refs,raw.get("reference_counts",{}),raw.get("reference_status",{}),scale_ref,"#B4C8C3"),
                    ("NonRM_callers",sources,nc.get("source_counts",{}),nc.get("source_status",{}),scale_nc,"#AAB7BC")):
                    if not keys:continue
                    for j,source in enumerate(keys):
                        x=panels[name]["left"]+j*7.56;_cell(c,x,y,7.56,ROW)
                        _mini_bar(c,x,cy,7.56,data.get(source),status.get(source),scale,color,{"sample_id":sid,"source":source,"panel":name},reservations[name])
                pp=panels["RM_pool"];_cell(c,pp["left"],y,pp["right"]-pp["left"],ROW)
                _fraction(c,pp["left"],cy,pp["right"]-pp["left"],raw.get("intersection_CDS"),raw.get("input_pool_CDS"),raw.get("intersection_count_status",raw.get("intersection_status")),"intersection",denstate=raw.get("pool_status"),layout=fractions["intersection"])
                for field,value,state in (("rm_pool_CDS",raw.get("input_pool_CDS"),raw.get("pool_status")),("rm_intersection_CDS",raw.get("intersection_CDS"),raw.get("intersection_count_status",raw.get("intersection_status")))):
                    c.ledger["cells"].append({"kind":"atlas_evidence","sample_id":sid,"field":field,"value":value,"status":state,"bounds_mm":[pp["left"],y,pp["right"]-pp["left"],ROW]})
                x=panels["RM_integrated"]["left"]
                for key,kind,typ,label,units in ROLE_COLUMNS:
                    w=units*CELL;fill=PINK if row.get("highlight") else "#FAF8FC" if key in ("I_MS","III_M") and kind=="context" else "#F6F8F7" if kind=="context" else None
                    _cell(c,x,y,w,ROW,fill);rec=row.get("compositions" if kind=="context" else "rm_roles",{}).get(key,{})
                    _role_glyph(c,x,cy,w,rec,key,kind);x+=w
                    if key in ("I_RMS","II_M","III_MR","IV_A"):c.line(x,y,x,y+ROW,"#8C999E",.4875)
                pp=panels["Retained_summary"];ret=row.get("retained",{});_cell(c,pp["left"],y,pp["right"]-pp["left"],ROW)
                _fraction(c,pp["left"],cy,pp["right"]-pp["left"],ret.get("R_CDS"),ret.get("total_CDS"),ret.get("R_status",ret.get("status")),"retained",scale_ret,denstate=ret.get("total_status",ret.get("status")),layout=fractions["retained"])
                rowmasks=_remap_masks(row,model.get("source_order",[]),sources)
                for j,mask in enumerate(masks[:-1]):
                    x=panels["cases"]["left"]+j*CELL;_cell(c,x,y,CELL,ROW)
                    item=rowmasks.get(mask,{});n=item.get("count");state=item.get("status")
                    if n is None or status_kind(state)=="unassessed":_text(c,x+CELL/2,cy,"—",6.7,ha="center")
                    elif n:
                        c.circle(x+CELL/2,cy,radius(n),"#E78E51")
                    c.ledger["cells"].append({"kind":"exclusive","sample_id":sid,"mask":mask,"count":n,"status":state,"bounds_mm":[x,y,CELL,ROW]})
                pp=panels["NonRM_summary"];_cell(c,pp["left"],y,pp["right"]-pp["left"],ROW)
                _fraction(c,pp["left"],cy,pp["right"]-pp["left"],nc.get("intersection_CDS"),nc.get("pool_CDS"),nc.get("intersection_status"),"noncanonical",denstate=nc.get("pool_status"),layout=fractions["noncanonical"])
                for f in fs:
                    x=fpos[f["id"]];cell=row.get("cells",{}).get(f["id"],{})
                    absent=_absence_cell(row,f["id"])
                    _cell(c,x,y,CELL,ROW,"#D0CFE5" if absent else None)
                    _family_glyph(c,x+CELL/2,cy,cell,f,config)
                    selected=f.get("highlight_source")
                    if selected and status_kind(cell.get("exclusive_source_status")) in ("complete","partial") and cell.get("exclusive_source_counts",{}).get(selected,0)>0:
                        color=config.get("atlas_link_colors",{}).get(selected,"#D56B25" if selected=="PADLOC" else "#347DA7")
                        c.rect(x+.02,cy-.85,.16,1.7,color)
                if ri==len(rows)-1 or (row.get("species"),row.get("genus"))!=(rows[ri+1].get("species"),rows[ri+1].get("genus")):
                    c.line(2.35,y+ROW,sampleleft+58,y+ROW,"#9EA8AC",.315)
            for name,p in panels.items():
                c.rect(p["left"],top,p["right"]-p["left"],len(rows)*ROW,"none",TEAL if name in ("RM_pool","Retained_summary") else "#8C999E",1.1625 if name in ("RM_pool","Retained_summary") else .3375)
            for f in fs:
                source=f.get("highlight_source");color=config.get("atlas_link_colors",{}).get(source,"#D56B25" if source=="PADLOC" else "#347DA7") if source else "black"
                _text(c,fpos[f["id"]]+CELL/2,bottom+1.7,display_label(f),8.3,ha="right",rotation=90,color=color)
            for j,mask in enumerate(masks[:-1]):
                label=" + ".join(alias(s,config) for s,b in zip(sources,mask) if b=="1")
                if mask.count("1")==1:label+=" only"
                source=next((s for s,b in zip(sources,mask) if b=="1"),None)
                color="black"
                if mask.count("1")==1 and source in ("PADLOC","DefenseFinder"):
                    from .manuscript_model import _intersection_state
                    items=[_remap_masks(r,model["source_order"],sources).get(mask,{}) for r in allrows]
                    count=sum(item.get("count") or 0 for item in items)
                    status=_intersection_state(item.get("status","unassessed") for item in items)
                    label+=f" ({number(count,status)} CDSs)";color="#D56B25" if source=="PADLOC" else "#347DA7"
                _text(c,panels["cases"]["left"]+(j+.5)*CELL,bottom+1.7,label,7.7,ha="right",rotation=90,color=color)
            x=panels["RM_integrated"]["left"]
            for key,kind,typ,label,units in ROLE_COLUMNS:
                if typ=="IV":_text(c,x+units*CELL/2,bottom+1.7,"R: nuclease" if key=="IV_R" else "A: auxiliary",7.7,ha="right",rotation=90)
                x+=units*CELL
            _legend(c,model,panels,groups,fpos,masks,sources,rows,bottom,sampleleft,config,taxonomy_colors,taxonomy_plan)
            c.ledger["workflow"]={"operation":"saved_result_import_and_aggregation","new_engine_runs":False,"activity_assessed":False}
            pages.append(c)
    return pages

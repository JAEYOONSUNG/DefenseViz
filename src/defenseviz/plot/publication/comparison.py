"""Uniform row-union bars and independently aligned comparison columns."""
from __future__ import annotations

from collections import Counter

from matplotlib import font_manager

from .drawing import (TEAL, ASSOCIATED, EFFECTOR, alias, category_color,
                      number, source_color, status_kind, display_label)
from .overlap import draw_flow, intersection_pages
from .manuscript_text import (MM_PER_POINT, _cap_height, _literal_path,
                              draw_manuscript_text)


def _line_ink(fonts, value, pt):
    """Visible top/bottom relative to the live text's cap-height anchor."""
    family=fonts.pick(value)
    path=font_manager.findfont(font_manager.FontProperties(family=family),fallback_to_default=False)
    outline=_literal_path(path,value)
    if outline is None:return (0.,0.)
    bounds=outline.get_extents()
    factor=pt*MM_PER_POINT/100
    cap=_cap_height(path)*pt*MM_PER_POINT/2
    return cap-bounds.y1*factor,cap-bounds.y0*factor


def _row_metrics(row, geom, fonts, config):
    """One consistent baseline for ordinary rows; grow for actual extra ink.

    The reference line includes capitals, descenders, comma and partial marks.
    Every actual label/number is measured too, including its selected fallback
    font. Multiline labels retain independent baselines and expand only their
    own row. Font size and the validation tolerances are unchanged.
    """
    pt=config["body_font_size"]
    label_start=max(1.7,len(_tag_colors(row,config))*.8+.4)
    label=fonts.wrap(display_label(row),geom["label_width"]-label_start-1.3,pt)
    lines=label.split("\n")
    spans=[_line_ink(fonts,line,pt) for line in lines]
    step=max(pt*MM_PER_POINT*1.12,max((b-a for a,b in spans),default=0)+.12)
    offsets=[(i-(len(lines)-1)/2)*step for i in range(len(lines))]
    extents=[_line_ink(fonts,"Agjpq 0123456789,+()-",pt)]
    extents.extend((a+d,b+d) for (a,b),d in zip(spans,offsets))
    values=[(row.get("G",row.get("genome_count")),metric_state(row,"inventory")),
            (row.get("N"),metric_state(row,"inventory")),
            (row.get("associated_CDS"),metric_state(row,"associated")),
            (row.get("effector_CDS"),metric_state(row,"effector"))]
    values.extend((row.get("source_counts",{}).get(col["source"]),
                   row.get("source_count_status",{}).get(col["source"])) for col in geom["source_cols"])
    extents.extend(_line_ink(fonts,number(value,state),pt) for value,state in values)
    top=min(a for a,b in extents);bottom=max(b for a,b in extents)
    return {"height":bottom-top+.12,"anchor_offset":.06-top,
            "label_lines":lines,"line_offsets":offsets,"label_start":label_start}


def source_order(result, config):
    names=[s["id"] for s in result.get("sources",[])]
    # Reference homology belongs in the overlap, not native annotation-count columns.
    native=[s["id"] for s in result.get("sources",[]) if s.get("kind")!="reference_annotation"]
    requested=config.get("source_order",["DefenseFinder","PADLOC","DefensePredictor"])
    return [s for s in requested if s in native]+[s for s in native if s not in requested]


def metric_state(row, key):
    return row.get("metric_status",{}).get(key,row.get("counts_status") if key=="inventory" else None)


def visible_rows(rows):
    # Unassessed all-zero rows are not hidden as biological absences.
    return [r for r in rows if (r.get("N") or 0)>0 or status_kind(metric_state(r,"inventory"))!="complete"]


def ordered_families(rows, config=None):
    """Respect explicit display grouping and companion order without reclassifying."""
    groups={}
    for row in rows:
        group=row.get("display_group",row.get("category","Unresolved"))
        groups.setdefault(group,{}).setdefault(row.get("partner_group") or row["id"],[]).append(row)
    order=(config or {}).get("category_order",[])
    keys=[key for key in order if key in groups]+[key for key in groups if key not in order]
    ordered=[r for key in keys for pair in groups[key].values()
             for r in sorted(pair,key=lambda r:r.get("partner_order",0))]
    # Distinct declared pools must not acquire indistinguishable display names.
    # Resolve only presentation collisions; keep every identity/member/count.
    for _ in range(len(ordered)+1):
        counts=Counter(display_label(r) for r in ordered)
        if max(counts.values(),default=0)<2:return ordered
        resolved=[dict(r,_display_label=str(r["id"])) if counts[display_label(r)]>1 else r for r in ordered]
        if [display_label(r) for r in resolved]==[display_label(r) for r in ordered]:
            raise ValueError("Distinct publication rows need unique identities")
        ordered=resolved
    raise ValueError("Publication row labels could not be disambiguated")


def _tag_colors(row, config):
    tags=row.get("category_tags") or row.get("mechanism_tags") or [row.get("category","Unresolved")]
    if isinstance(tags,str):tags=[tags]
    primary=row.get("display_group",row.get("category","Unresolved"))
    return [category_color(t,config) for t in dict.fromkeys([primary,row.get("category",primary),*tags])]


def layout(width, fonts, rows, sources, config):
    """Reserve count banks by actual maximum text width, not cell-edge spacing."""
    left=8.45;available=width-2*left
    body=config["body_font_size"]
    number_width=max(7.,max((fonts.width(number(r.get("N"),metric_state(r,"inventory")),body) for r in rows),default=0)+1.2)
    for r in rows:
        for source in sources:
            number_width=max(number_width,fonts.width(number(r.get("source_counts",{}).get(source),r.get("source_count_status",{}).get(source)),body)+1.2)
    source_width=12+1.2+number_width+8.8
    label_width=32.;genome_width=15.;dnmb_width=max(12.,number_width+3.)
    nuc_width=max(45.,12+1.2+number_width+5.+number_width+4.2,
                  18.9+number_width+fonts.width("Effector CDSs",7.5)+3.,
                  fonts.width("Putative nuclease-associated CDSs",7.8,True)+.8)
    required=label_width+genome_width+dnmb_width+nuc_width+len(sources)*source_width
    if required>available and sources:
        # Use spare space after each number bank before adding source pages.
        # Neither the 12-mm tracks nor text/count widths are reduced.
        reduction=min(2.,(required-available)/len(sources))
        source_width-=reduction;required-=reduction*len(sources)
    if required>available+.01:
        return None
    label_width+=available-required
    cursor=left+label_width
    genome_right=cursor+genome_width-2
    cursor+=genome_width
    dnmb_right=cursor+dnmb_width-2
    cursor+=dnmb_width
    nuc_x=cursor+1.5;nuc_total=nuc_x+12+1.2+number_width
    nuc_effector=cursor+nuc_width-4.2
    cursor+=nuc_width
    source_cols=[]
    for s in sources:
        track=cursor+1.5
        source_cols.append({"source":s,"track":track,"number":track+12+1.2+number_width,"left":cursor,"width":source_width})
        cursor+=source_width
    return {"left":left,"width":available,"label_width":label_width,"genome_right":genome_right,
            "dnmb_right":dnmb_right,"nuc_x":nuc_x,"nuc_total":nuc_total,"nuc_effector":nuc_effector,
            "source_cols":source_cols,"number_width":number_width,"body_pt":body,"nuc_width":nuc_width}


def draw_table(c, rows, sources, config, y, name, *, row_height=2.7, header=True, layout_rows=None):
    geom=layout(c.width,c.fonts,rows if layout_rows is None else layout_rows,sources,config)
    if geom is None:raise ValueError("Source columns do not fit publication page; paginate them first")
    x=geom["left"];w=geom["width"];pt=geom["body_pt"]
    header_h=9.5 if header else 0
    c.line(x,y,x+w,y,"#ACBDBE",.4)
    if header:
        c.text(x+1.8,y+4.5,"Category" if name=="categories" else "Family / pool",8,bold=True)
        c.text(geom["genome_right"]-4,y+4.5,"Number of\ngenomes",7.5,ha="center",bold=True)
        c.text(geom["dnmb_right"]-4.5,y+4.5,"DNMB",8,ha="center",bold=True,color=TEAL)
        c.text(geom["nuc_x"]+geom["nuc_width"]/2-1.5,y+2.2,"Putative nuclease-associated CDSs",7.8,ha="center",bold=True,color=TEAL)
        c.text(geom["nuc_total"],y+6.5,"Total",7.5,ha="right",color="#6F8185",role="total_header")
        c.text(geom["nuc_effector"],y+6.5,"Effector CDSs",7.5,ha="right",color="#6F8185",role="effector_header")
        for col in geom["source_cols"]:
            label=c.fonts.wrap(alias(col["source"],config),col["width"]-1,8,True)
            c.text((col["track"]+col["number"])/2,y+4.5,label,8,ha="center",bold=True,color=source_color(col["source"],config))
    y+=header_h
    c.line(x,y,x+w,y,"#ACBDBE",.4)
    cursor=y
    for i,r in enumerate(rows):
        rh=row_height[i] if isinstance(row_height,list) else row_height
        metrics=_row_metrics(r,geom,c.fonts,config)
        if rh+1e-8<metrics["height"]:
            raise ValueError("Comparison family label needs a larger row height")
        top=cursor;yy=top+rh/2;cursor+=rh
        text_y=top+metrics["anchor_offset"]+(rh-metrics["height"])/2
        bounds=[x,top,w,rh]
        shade="#EEF1F3" if r.get("category")=="Unresolved" else ("#FAFBFB" if i%2 else "white")
        c.rect(x,top,w,rh,shade)
        c.rect(geom["nuc_x"]-1.5,top,geom["nuc_width"],rh,"#F4F9F6")
        colors=_tag_colors(r,config)
        for j,color in enumerate(colors):c.rect(x+j*.8,top+.1,.65,rh-.2,color)
        for line,offset in zip(metrics["label_lines"],metrics["line_offsets"]):
            draw_manuscript_text(c,x+metrics["label_start"],text_y+offset,line,pt,
                color="#24393B",owner=[x,top,geom["label_width"],rh],
                ref={"row_id":r["id"],"field":"label"})
        ref={"table":name,"row_id":r["id"]}
        for val,state,anchor,field in ((r.get("G",r.get("genome_count")),metric_state(r,"inventory"),geom["genome_right"],"genomes"),
                  (r.get("N"),metric_state(r,"inventory"),geom["dnmb_right"],"DNMB")):
            draw_manuscript_text(c,anchor,text_y,
                number(val,state),pt,ha="right",color="#24393B",owner=bounds,
                role="numeric",ref=dict(ref,field=field,value=val,status=state))
        a,e=r.get("associated_CDS"),r.get("effector_CDS")
        ast,est=metric_state(r,"associated"),metric_state(r,"effector")
        segments=None
        if status_kind(est) in ("complete","partial") and e is not None and a is not None:
            segments=[(e,EFFECTOR),(a-e,ASSOCIATED)]
        c.bar(geom["nuc_x"],yy,a,r.get("N"),ast,ASSOCIATED,segments=segments,ref=dict(ref,field="associated"))
        for val,state,anchor,field in ((a,ast,geom["nuc_total"],"associated"),(e,est,geom["nuc_effector"],"effector")):
            draw_manuscript_text(c,anchor,text_y,
                number(val,state),pt,ha="right",color="#24393B",role="numeric",
                owner=bounds,ref=dict(ref,field=field,value=val,status=state))
        for col in geom["source_cols"]:
            source=col["source"];val=r.get("source_counts",{}).get(source);state=r.get("source_count_status",{}).get(source)
            c.bar(col["track"],yy,val,r.get("N"),state,source_color(source,config),ref=dict(ref,field=source))
            draw_manuscript_text(c,col["number"],text_y,
                number(val,state),pt,ha="right",color="#24393B",role="numeric",
                owner=bounds,ref=dict(ref,field=source,value=val,status=state))
        c.line(x,top+rh,x+w,top+rh)
        c.ledger["cells"].append({"kind":"comparison_row","table":name,"row_id":r["id"],"bounds_mm":bounds,
            "count":r.get("N"),"genomes":r.get("G",r.get("genome_count")),"associated":a,"effector":e,
            "row_height_basis":"visible_glyph_ink","label_lines":metrics["label_lines"],
            "text_anchor_y_mm":text_y,
            "source_counts":{s:r.get("source_counts",{}).get(s) for s in sources},"source_order":sources,
            "numeric_anchors":{ "DNMB":geom["dnmb_right"],"Total":geom["nuc_total"],"Effector CDSs":geom["nuc_effector"],
              **{col["source"]:col["number"] for col in geom["source_cols"]}},"source_ref":ref})
    return cursor


def _row_height(rows, geom, fonts, config):
    return max((_row_metrics(r,geom,fonts,config)["height"] for r in rows),default=0.)


def comparison_pages(result, fonts, config, page_factory):
    before=result.get("overlap_regions",{});after=result.get("priority_overlap_regions",{})
    categories=ordered_families(visible_rows(result.get("inventory_rows",{}).get("categories",[])),config)
    families=ordered_families(visible_rows(result.get("inventory_rows",{}).get("families",[])),config)
    sources=source_order(result,config)
    legend_sources=list(dict.fromkeys([*before.get("source_order",[]),*after.get("source_order",[]),*sources]))
    # Whole source groups repeat across independent source-column pages.
    chunks=[];current=[]
    for source in sources:
        if current and layout(215.9,fonts,categories+families,current+[source],config) is None:
            chunks.append(current);current=[]
        current.append(source)
    chunks.append(current)
    pages=[]
    first=page_factory(215.9,279.4,"comparison")
    first.text(8.45,8,config["title"],12,bold=True)
    first.text(8.45,14,f"{len(result.get('samples',[])):,} genomes · Coding sequences (CDSs)",7.5,color="#6F8185")
    first.text(8.45,21,"a  DNMB prioritization",9,bold=True)
    flow=draw_flow(first,before,after,config)
    if not flow:
        # Drop the unfilled overview page and create explicit intersection pages.
        first.fig.clear()
        pages.extend(intersection_pages(before,after,fonts,config,page_factory))
        first=None
    for block_index,source_chunk in enumerate(chunks):
        c=first if block_index==0 and first is not None else page_factory(215.9,279.4,"comparison")
        top=80 if c is first else 19
        if c is not first:c.text(8.45,9,config["title"]+" - inventory",11,bold=True)
        pending=[("categories",categories),("families",families)]
        for name,allrows in pending:
            rows=list(allrows);part=0
            while rows or (part==0 and not allrows):
                geom=layout(215.9,fonts,categories+families,source_chunk,config)
                if geom is None:raise ValueError("A source count cannot fit the configured comparison page")
                heights=[_row_height([r],geom,fonts,config) for r in rows]
                room=264-top-3.5-9.5;capacity=0;used=0.
                for rh in heights:
                    if used+rh>room:break
                    used+=rh;capacity+=1
                if not allrows:capacity=1
                if capacity<1:
                    if top==19:raise ValueError("A comparison label cannot fit a single page at the requested font size")
                    _footer(c,config,legend_sources);pages.append(c);c=page_factory(215.9,279.4,"comparison")
                    c.text(8.45,9,config["title"]+" - continued",11,bold=True);top=19;continue
                take=rows[:capacity]
                title="Defense categories" if name=="categories" else "Defense families"
                c.text(8.45,top,("b  " if name=="categories" else "c  ")+title+(" (continued)" if part else ""),9,bold=True)
                end=draw_table(c,take,source_chunk,config,top+3.5,name,row_height=heights[:capacity],layout_rows=categories+families)
                if not allrows:
                    c.text(10,top+17,"No observed inventory supplied",7.8)
                    end=top+20
                rows=rows[len(take):];part+=1;top=end+3.5
                if not rows:break
                _footer(c,config,legend_sources);pages.append(c);c=page_factory(215.9,279.4,"comparison")
                c.text(8.45,9,config["title"]+" - continued",11,bold=True);top=19
        _footer(c,config,legend_sources);pages.append(c)
    return pages


def _footer(c,config,sources=()):
    c.rect(8.45,266.1,2,1.4,ASSOCIATED);c.text(11.4,266.8,"Putative nuclease-associated CDSs",7.0)
    c.rect(65,266.1,2,1.4,EFFECTOR);c.text(68,266.8,"Effector CDSs",7.0)
    c.text(c.width-8.45,266.8,"Rows overlap",7.0,ha="right")
    c.text(8.45,271,"Track: row union · NA: unassessed · +: partial · Activity untested",7.0)
    labels=[f"{alias(s,config)}: {s}" for s in sources if alias(s,config)!=s]
    if labels:
        text="; ".join(labels)
        if c.fonts.width(text,7.0)>c.width-16.9:
            text="Source abbreviations: see publication caption"
        c.text(8.45,275.3,text,7.0,role="source_abbreviations")

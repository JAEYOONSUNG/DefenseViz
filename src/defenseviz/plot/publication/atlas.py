"""Cohort-independent workflow, nested-CDS atlas and type benchmark pages."""
from __future__ import annotations

import math
from itertools import product

from .drawing import TEAL, PURPLE, GRID, alias, category_color, source_color, number, status_kind, display_label
from .comparison import source_order, metric_state, visible_rows, ordered_families


def _chunks(values, n):
    return [values[i:i+n] for i in range(0,len(values),n)] or [[]]


def workflow(c, result, config, y=19):
    c.text(8.45,y,"a  Saved evidence and DNMB prioritization",10,bold=True)
    availability=result.get("availability",{})
    stages=[("Source annotations","Original final calls and reference matches", "Saved outputs"),
            ("Saved DNMB decisions","Retained annotations and explicit reviewed components",availability.get("saved_decisions","unassessed")),
            ("Effector evidence","Individual effectors and reviewed complex cores",availability.get("effector","unassessed")),
            ("Review priorities","R annotations + effector CDSs; activity untested", "Annotation-based")]
    gap=5.;w=(c.width-16.9-3*gap)/4
    for i,(heading,body,state) in enumerate(stages):
        x=8.45+i*(w+gap)
        c.rect(x,y+4,w,20,"#F7FAF9",TEAL,.45)
        c.text(x+2,y+8,heading,8,bold=True,color=TEAL)
        label=c.fonts.wrap(body,w-4,7.2)
        c.text(x+2,y+15,label,7.2,owner=[x,y+4,w,20])
        c.text(x+2,y+22,str(state).replace("_"," "),7,color="#6F8185")
        if i<3:c.ax.annotate("",xy=(x+w+gap-.5,y+14),xytext=(x+w+.5,y+14),arrowprops={"arrowstyle":"-|>","color":TEAL,"lw":.6})
    c.ledger["workflow"]={"operation":"saved_result_import_and_aggregation","new_engine_runs":False,
                           "capabilities":result.get("availability",{}),"activity_assessed":False}


def _outer_radius(count):
    # Stable presentation bins, not biological thresholds.
    return .45 if count<=1 else .7 if count<=3 else 1.0 if count<=7 else 1.3 if count<=15 else 1.55


def section_headings(fonts,config,groups,prefix,width):
    """Use anchored headings only when their real group spans can contain them."""
    titles=config.get("atlas_section_titles",{})
    records=[{"section":"R-M","text":titles.get("R-M","1 R-M evidence and curation"),
              "left":8.45,"width":max(1.,prefix-2),"color":TEAL}]
    for category,left,right in groups:
        if category in ("Modification","Nuclease"):
            records.append({"section":category,"text":titles.get(category,"3 Non-canonical DNA modification system" if category=="Modification" else "4 Nuclease-associated evidence"),
                            "left":left,"width":right-left,"color":category_color(category,config)})
    others=[(left,right) for category,left,right in groups if category not in ("R-M","Modification","Nuclease")]
    if others:
        records.append({"section":"complementary","text":titles.get("complementary","5 Complementary defense priorities"),
                        "left":others[0][0],"width":others[-1][1]-others[0][0],"color":"#24393B"})
    for r in records:
        r["group_span_mm"]=[r["left"],r["width"]]
        r["label"]=fonts.wrap(r["text"],max(1.,r["width"]-1),7.8,True)
    anchored=all(len(r["label"].splitlines())<=2 and fonts.width(r["label"],7.8,True)<=r["width"]-1 for r in records)
    if not anchored:
        gap=4.;slot=(width-16.9-gap*(len(records)-1))/len(records)
        for i,r in enumerate(records):
            r["left"]=8.45+i*(slot+gap);r["width"]=slot
            r["label"]=fonts.wrap(r["text"],slot-1,7.8,True)
    height=max(7.,max(fonts.height(r["label"],7.8,True) for r in records)+1)
    return records,height,"anchored" if anchored else "section_key"


def draw_glyph(c,x,y,row,family,config,bounds):
    n=row.get("N");state=metric_state(row,"inventory")
    a,e=row.get("associated_CDS"),row.get("effector_CDS")
    ast,est=metric_state(row,"associated"),metric_state(row,"effector")
    ref={"sample_id":row.get("sample_id"),"family_id":row.get("family_id",row.get("id"))}
    if status_kind(state)=="unassessed" or n is None:
        c.text(x,y,"NA",7,ha="center",owner=bounds,role="missing",ref=ref)
        radius=inner=0
    elif n==0:
        c.circle(x,y,.16,"#D7DFDF");radius=inner=0
    else:
        radius=_outer_radius(n)
        c.circle(x,y,radius,category_color(family.get("category","Unresolved"),config))
        if status_kind(ast) in ("complete","partial") and a:
            c.circle(x,y,radius+.10,None,PURPLE,.55)
        inner=radius*math.sqrt(e/n) if e is not None and status_kind(est) in ("complete","partial") else 0
        if inner:c.circle(x,y,inner,PURPLE)
        if "unassessed" in (status_kind(ast),status_kind(est)):
            c.text(x,y,"?",7,ha="center",role="missing_evidence",ref=ref)
        if "partial" in (status_kind(state),status_kind(ast),status_kind(est)):
            c.text(x+1.65,y-1.25,"+",6.5,ha="center",role="partial",ref=ref)
    c.ledger["glyphs"].append(dict(ref,bounds_mm=bounds,total=n,associated=a,effector=e,
        inventory_status=state,associated_status=ast,effector_status=est,outer_radius_mm=radius,
        inner_radius_mm=inner,inner_area_fraction=(e/n if n and inner else 0),
        category=family.get("category"),outer_size="count_bin",activity_inferred=False))


def atlas_pages(result,fonts,config,page_factory):
    families=ordered_families(visible_rows(result.get("inventory_rows",{}).get("families",[])),config)
    samples=result.get("samples",[])
    samplemap={s["sample_id"]:s for s in result.get("sample_summary",[])}
    cellmap={(r["sample_id"],r.get("family_id",r.get("id"))):r for r in result.get("cells",[])}
    types={(r["sample_id"],r["rm_type"]):r for r in result.get("type_tables",[])}
    componentmap={(r["sample_id"],r["rm_type"],r["role"]):r for r in result.get("rm_components",[])}
    components=list(dict.fromkeys((r["rm_type"],r["role"]) for r in result.get("rm_components",[])))
    compact=[]
    for typ,role in components:
        rr=[r for r in result.get("type_tables",[]) if r["rm_type"]==typ]
        if any(r.get("is_total",str(typ).lower() in ("total","all")) for r in rr):continue
        if rr and all(r.get("unique_CDS")==0 and all(status_kind(s) in ("complete","not_applicable") for s in r.get("role_states",{}).values()) for r in rr):continue
        if rr and all(status_kind(r.get("role_states",{}).get("A" if role=="Aux" else role))=="not_applicable" for r in rr):continue
        compact.append((typ,role))
    components=compact
    native=source_order(result,config)
    allsources=native+[s["id"] for s in result.get("sources",[]) if s["id"] not in native]
    maxwidth,maxheight=config["atlas_page_size"]
    has_rm=any("rm_source_counts" in s for s in samplemap.values())
    # Source and role columns are explicit observed endpoints, never family-row sums.
    descriptors=[]
    if has_rm:
        descriptors += [{"key":"rm_source_counts","source":source,"label":alias(source,config),"group":"Original R-M","width":8.} for source in allsources]
        descriptors += [{"key":"rm_pool_CDS","label":"Pool","group":"R-M evidence","width":10.},
                        {"key":"rm_intersection_CDS","label":"Shared","group":"R-M evidence","width":11.},
                        {"key":"retained_R_fraction","label":"R / total","group":"DNMB retained","width":14.}]
        if any("complete_saved_contexts" in r for r in samplemap.values()):
            descriptors.append({"key":"complete_saved_contexts","label":"Complete","group":"Saved contexts","width":14.})
    else:
        descriptors += [{"key":"source_counts","source":source,"label":alias(source,config),"group":"Original sources","width":10.} for source in native]
        descriptors += [{"key":"integrated_CDS","label":"DNMB","group":"Integrated","width":13.}]
    descriptors += [{"key":"component","type":typ,"role":role,"label":"A" if role=="Aux" else role,"group":str(typ),"width":6.} for typ,role in components]
    if has_rm:
        descriptors += [{"key":"noncanonical_source_counts","source":source,"label":alias(source,config),"group":"Non-canonical CDSs","width":8.} for source in native]
    # Keep every descriptor; exceptionally many source/role columns use separate pages.
    descriptor_pages=[];current=[];used=0.
    for desc in descriptors:
        if current and used+desc["width"]>maxwidth-140:
            descriptor_pages.append(current);current=[];used=0
        current.append(desc);used+=desc["width"]
    descriptor_pages.append(current)
    pages=[]
    for ds in descriptor_pages:
        prefix=sum(d["width"] for d in ds)+4
        familywidth=4.3
        family_pages=[];current=[];used=0.;last=None
        for f in families:
            cat=f.get("display_group",f.get("category","Unresolved"));extra=familywidth+(1.5 if current and cat!=last else 0)
            if current and (used+extra>maxwidth-22-prefix-64 or len(current)>=config["max_families_per_page"]):
                family_pages.append(current);current=[];used=0.;extra=familywidth
            current.append(f);used+=extra;last=cat
        family_pages.append(current)
        for fs in family_pages:
            fpositions=[];cursor=8.45+prefix;last=None;groups=[]
            for f in fs:
                category=f.get("display_group",f.get("category","Unresolved"))
                if last is not None and category!=last:cursor+=1.5
                if category!=last:groups.append([category,cursor,cursor])
                fpositions.append(cursor);cursor+=familywidth;groups[-1][2]=cursor;last=category
            width=max(270.,min(maxwidth,cursor+64+8.45))
            headings,headingheight,headingmode=section_headings(fonts,config,groups,prefix,width)
            header_offset=max(0.,headingheight-7.)
            descriptor_spans=[]
            for descriptor in ds:
                if descriptor_spans and descriptor_spans[-1][0]==descriptor["group"]:
                    descriptor_spans[-1][1]+=descriptor["width"]
                else:descriptor_spans.append([descriptor["group"],descriptor["width"]])
            group_height=max((fonts.height(fonts.wrap(group,span-.7,7,True),7,True) for group,span in descriptor_spans),default=0.)
            group_extra=max(0.,group_height-6.)
            header_offset+=group_extra
            labeldepth=max(23.,min(80.,max((fonts.width(display_label(f),7) for f in fs),default=10)+4))
            bodytop=82.+header_offset;room=maxheight-bodytop-labeldepth-23
            slicelist=[];rows=[];used=0.
            for sample in samples:
                lab=fonts.wrap(sample.get("label",sample["sample_id"]),60,7.8,sample.get("highlight",False))
                rh=max(3.6,fonts.height(lab,7.8,sample.get("highlight",False))+.4)
                if rows and (used+rh>room or len(rows)>=config["max_rows_per_page"]):
                    slicelist.append(rows);rows=[];used=0
                if rh>room:raise ValueError("Sample label cannot fit one atlas page")
                rows.append((sample,lab,rh));used+=rh
            if rows or not samples:slicelist.append(rows)
            for sample_rows in slicelist:
                height=min(maxheight,max(145.,bodytop+sum(r[2] for r in sample_rows)+labeldepth+23))
                c=page_factory(width,height,"atlas")
                c.text(8.45,8,config["title"],12,bold=True)
                workflow(c,result,config)
                c.text(8.45,54,"b  Comparative defense landscape",10,bold=True)
                for h in headings:
                    c.text(h["left"],57+headingheight/2,h["label"],7.8,bold=True,color=h["color"],
                           owner=[h["left"],57,h["width"],headingheight],role="section_heading",ref={"section":h["section"]})
                c.ledger["section_headings"]={"mode":headingmode,"height_mm":headingheight,"records":headings}
                sampleleft=width-8.45-62
                c.text(sampleleft+1.5,75+header_offset,"Samples",8,bold=True)
                x=8.45;previous=None;group_start=x
                for di,d in enumerate(ds):
                    group=d["group"]
                    if group!=previous:
                        if previous is not None:
                            label=fonts.wrap(previous,x-group_start-.7,7,True)
                            c.text((group_start+x)/2,68+header_offset-group_extra/2,label,7,bold=True,ha="center")
                            c.line(group_start,72+header_offset,x-.6,72+header_offset,TEAL,.4)
                        group_start=x;previous=group
                    d=dict(d,x=x)
                    c.text(x+d["width"]/2,76+header_offset,d["label"],7.5,bold=True,ha="center",color=source_color(d["source"],config) if d.get("source") else TEAL)
                    x+=d["width"]
                if ds:
                    label=fonts.wrap(previous,x-group_start-.7,7,True)
                    c.text((group_start+x)/2,68+header_offset-group_extra/2,label,7,bold=True,ha="center");c.line(group_start,72+header_offset,x-.6,72+header_offset,TEAL,.4)
                for category,left,right in groups:
                    label=fonts.wrap(category,right-left,7,True)
                    # Very narrow groups use the complete category label in the legend.
                    if len(label.splitlines())<=2:c.text((left+right)/2,68+header_offset,label,7,bold=True,ha="center",color=category_color(category,config))
                for f,xx in zip(fs,fpositions):
                    primary=f.get("display_group",f.get("category","Unresolved"))
                    tags=f.get("category_tags") or f.get("mechanism_tags") or []
                    tags=[tags] if isinstance(tags,str) else tags
                    secondary=list(dict.fromkeys(t for t in [f.get("category",primary),*tags] if t!=primary))
                    segments=f.get("band_segments") or [{"category":tag,"fraction":1/(1+len(secondary))} for tag in [primary,*secondary]]
                    total=sum(segment["fraction"] for segment in segments)
                    if total<=0 or any(segment["fraction"]<=0 for segment in segments):raise ValueError("Category band fractions must be positive")
                    bandheight=2.5-.3*(len(segments)-1);yy=77.+header_offset
                    if bandheight<=0:raise ValueError("Too many category bands for an atlas column")
                    for segment in segments:
                        h=bandheight*segment["fraction"]/total
                        c.rect(xx,yy,familywidth,h,category_color(segment["category"],config));yy+=h+.3
                    c.ledger["bands"].append({"family_id":f["id"],"primary":primary,"secondary":secondary,"segments":segments,"bounds_mm":[xx,77+header_offset,familywidth,2.5]})
                y=bodytop
                for sample,label,rh in sample_rows:
                    sid=sample["sample_id"];yc=y+rh/2
                    c.rect(8.45,y,width-16.9,rh,"#FBE0EB" if sample.get("highlight") else "white")
                    if sample.get("highlight"):c.rect(sampleleft,y,62,rh,"#F3BCD2")
                    c.text(sampleleft+1.5,yc,label,7.8,bold=sample.get("highlight",False),owner=[sampleleft,y,62,rh],ref={"sample_id":sid})
                    summary=samplemap.get(sid,{});x=8.45
                    for d in ds:
                        field=d["key"];state=None;value=None
                        if field=="component":
                            rec=componentmap.get((sid,d["type"],d["role"]),{})
                            type_rec=types.get((sid,d["type"]),{})
                            role="A" if d["role"]=="Aux" else d["role"]
                            state=type_rec.get("role_states",{}).get(role)
                            value=rec.get("count")
                            if value is None and role in type_rec:value=type_rec[role]
                        elif d.get("source"):
                            value=summary.get(field,{}).get(d["source"])
                            statusfield=field.replace("counts","count_status")
                            state=summary.get(statusfield,summary.get("source_count_status",{})).get(d["source"])
                        elif field=="retained_R_fraction":
                            numerator=summary.get("retained_R_CDS");denominator=summary.get("retained_RM_CDS")
                            state=summary.get("retained_ratio_status")
                            value=(number(numerator,summary.get("retained_R_counts_status",state))+"/"+
                                   number(denominator,summary.get("retained_RM_counts_status",state))) if numerator is not None and denominator is not None and status_kind(state)!="unassessed" else "NA"
                        elif field in ("rm_pool_CDS","rm_intersection_CDS"):
                            value=summary.get(field)
                            state=summary.get("rm_pool_status" if field=="rm_pool_CDS" else "rm_intersection_status")
                        else:
                            value=summary.get(field)
                            state=summary.get("complete_context_state",summary.get("complete_saved_contexts_status")) if field=="complete_saved_contexts" else summary.get("metric_status",{}).get("inventory")
                        display=value if field=="retained_R_fraction" else number(value,state)
                        if field=="component":c.rect(x,y,d["width"],rh,"#CCD3D8" if status_kind(state)=="not_applicable" else "#F0F2F3")
                        c.text(x+d["width"]-.8,yc,display,7.8,ha="right",owner=[x,y,d["width"],rh],role="numeric",ref={"sample_id":sid,"field":field,"source":d.get("source"),"rm_type":d.get("type"),"role":d.get("role"),"value":value,"status":state})
                        c.line(x,y,x,y+rh)
                        c.ledger["cells"].append({"kind":"atlas_evidence","sample_id":sid,"field":field,"source":d.get("source"),"rm_type":d.get("type"),"role":d.get("role"),"value":value,"status":state,"bounds_mm":[x,y,d["width"],rh]})
                        x+=d["width"]
                    for f,xx in zip(fs,fpositions):
                        row=cellmap.get((sid,f["id"]),{"sample_id":sid,"family_id":f["id"],"N":None})
                        draw_glyph(c,xx+familywidth/2,yc,row,f,config,[xx,y,familywidth,rh]);c.line(xx,y,xx,y+rh)
                    c.line(8.45,y+rh,width-8.45,y+rh)
                    c.ledger["cells"].append({"kind":"sample_row","sample_id":sid,"highlight":sample.get("highlight",False),"bounds_mm":[8.45,y,width-16.9,rh]})
                    y+=rh
                for f,xx in zip(fs,fpositions):
                    label=fonts.wrap(display_label(f),labeldepth-3,7)
                    c.text(xx+familywidth/2,y+3,label,7,rotation=90,ha="right",role="family_label",ref={"family_id":f["id"]})
                foot=height-14
                c.circle(10,foot,1,None,PURPLE);c.text(13,foot,"Putative nuclease-associated CDSs",7.2)
                c.circle(76,foot,1,None,PURPLE);c.circle(76,foot,.6,PURPLE);c.text(79,foot,"Effector CDSs (inner area / total)",7.2)
                c.text(8.45,height-8,"Count bins: 1, 2-3, 4-7, 8-15, 16+ · NA / ?: unassessed · +: observed lower bound · Activity untested",7)
                c.ledger["family_ids"]=[f["id"] for f in fs];c.ledger["source_order"]=allsources
                pages.append(c)
    return pages

def _type_display(r,field):
    subset=None
    if field=="Type":return str(r["rm_type"]),r["rm_type"],"complete",None
    if field=="Pool":value=r.get("pool_CDS");state=r.get("pool_status")
    elif field=="DNMB":value=r.get("unique_CDS");state=r.get("counts_status")
    elif field in ("R","M","S","A"):value=r.get(field);state=r.get("role_states",{}).get(field)
    else:value=r.get("source_counts",{}).get(field);state=r.get("source_count_status",{}).get(field)
    display=number(value,state)
    if field=="R":
        subset=r.get("incomplete_R_context")
        assessed=r.get("context_state")
        applicable=r.get("context_applicability") not in (None,"not_applicable","TypeIV_not_applicable","not_assessed_by_I_III_rule")
    elif field in ("M","S"):
        subset=r.get("orphan_"+field);assessed=r.get("orphan_state",{}).get(field);applicable=True
    else:assessed=None;applicable=False
    if value and status_kind(state) in ("complete","partial"):
        if applicable and status_kind(assessed)=="unassessed":display+=" (NA)"
        elif applicable and (subset or status_kind(assessed)=="partial"):
            display+=" ("+number(subset,assessed)+")"
    return display,value,state,subset


def type_pages(result,fonts,config,page_factory):
    rows=result.get("type_tables",[])
    samples={s["sample_id"]:s for s in result.get("samples",[])}
    selected=config.get("type_samples")
    caseids=[s for s in samples if any(r["sample_id"]==s for r in rows) and (selected is None or s in selected)]
    if selected is None:
        highlighted=[sid for sid in caseids if samples[sid].get("highlight")]
        if highlighted:caseids=highlighted
    if not caseids:return []
    sources=source_order(result,config);pages=[]
    for source_slice in _chunks(sources,3):
        fields=["Type"]+(["Pool"] if any("pool_CDS" in r for r in rows) else [])+source_slice+["DNMB","R","M","S","A"]
        widths={field:max(13.,fonts.width(alias(field,config),8,True)+3,
                          max((fonts.width(_type_display(r,field)[0],9.4,True)+3 for r in rows),default=0)) for field in fields}
        casew=sum(widths.values());gap=6
        perpage=max(1,min(3,int((config["atlas_page_size"][0]-16.9+gap)/(casew+gap))))
        for cases in _chunks(caseids,perpage):
            width=max(215.9,len(cases)*casew+(len(cases)-1)*gap+16.9)
            maxrows=max(sum(r["sample_id"]==sid for r in rows) for sid in cases)
            labels={sid:fonts.wrap(samples[sid].get("label",sid),casew-4,9,True) for sid in cases}
            strainh=max(7,max(fonts.height(label,9,True)+2 for label in labels.values()))
            bodytop=23+strainh+4+6
            height=max(85,bodytop+maxrows*6+20)
            c=page_factory(width,height,"type_benchmark")
            c.ledger["sample_ids"]=cases
            c.ledger["sample_selection"]="explicit" if selected is not None else "highlighted_default" if any(samples[sid].get("highlight") for sid in cases) else "all_typed_default"
            c.text(8.45,9,"c  CDS inclusion by assigned R-M type",11,bold=True)
            c.text(8.45,16,"R(n): incomplete retained context; M/S(n): orphan · NA: unassessed · +: lower bound · -: not applicable",7.8)
            for ci,sid in enumerate(cases):
                case=samples[sid];rr=[r for r in rows if r["sample_id"]==sid]
                x=8.45+ci*(casew+gap)
                c.rect(x,23,casew,strainh,"#F3EFE4" if case.get("cohort") in ("external","external_validation") else "#F8E9EF")
                c.text(x+casew/2,23+strainh/2,labels[sid],9,bold=True,ha="center",owner=[x,23,casew,strainh])
                rolesw=sum(widths[f] for f in ("R","M","S","A"));groupy=23+strainh
                c.rect(x,groupy,casew-rolesw,4,"#F1F4F5");c.rect(x+casew-rolesw,groupy,rolesw,4,"#E4E8EB")
                c.text(x+(casew-rolesw)/2,groupy+2,"Shared-pool CDSs",8,bold=True,ha="center")
                c.text(x+casew-rolesw/2,groupy+2,"DNMB components",8,bold=True,ha="center")
                positions=[];cursor=x
                for field in fields:
                    cw=widths[field];positions.append((field,cursor,cw));c.rect(cursor,groupy+4,cw,6,"#ECEFF1" if field in ("R","M","S","A") else "#F4F7F6")
                    c.text(cursor+cw/2,groupy+7,alias(field,config),8,bold=True,ha="center");cursor+=cw
                for i,r in enumerate(rr):
                    y=bodytop+i*6;is_total=str(r.get("rm_type")).lower() in ("all","total")
                    for field,left,cw in positions:
                        fill="#DFE4E8" if is_total else "#F1F3F4" if field in ("R","M","S","A") else "#FAFBFB" if i%2 else "white"
                        c.rect(left,y,cw,6,fill)
                        text,value,state,subset=_type_display(r,field)
                        if status_kind(state)=="not_applicable":c.rect(left,y,cw,6,"#CCD3D8")
                        c.text(left+cw/2 if field=="Type" else left+cw-1.8,y+3,text,9.4,bold=is_total,ha="center" if field=="Type" else "right",owner=[left,y,cw,6],role="type_label" if field=="Type" else "type_numeric",ref={"sample_id":sid,"type":r["rm_type"],"field":field,"value":value,"included_subset":subset})
                        c.ledger["cells"].append({"kind":"type_cell","sample_id":sid,"rm_type":r["rm_type"],"field":field,
                            "value":value,"state":state,"included_subset":subset,"display":text,"bounds_mm":[left,y,cw,6],"numeric_font_pt":9.4})
                    c.line(x,y+6,x+casew,y+6)
            c.text(8.45,height-8,"Unique CDS totals are not sums of role columns. Context annotations do not establish biological activity.",7)
            pages.append(c)
    return pages

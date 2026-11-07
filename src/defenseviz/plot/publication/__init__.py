"""Editable publication figures over audited saved-result analysis outputs.

This package never classifies proteins, joins raw records or re-runs a detector.
The legacy :mod:`defenseviz.plot.renderer` API remains unchanged.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

from .drawing import Canvas, Fonts, STROKE_FLOOR_PT, status_kind, alias, source_color
from .atlas import atlas_pages, type_pages
from .manuscript_atlas import manuscript_atlas_pages
from .comparison import comparison_pages

__all__ = ["render_publication", "validate_publication_model"]

DEFAULTS = {"sections": ["atlas", "type_benchmarks", "comparison"],
    "font_family": "Arial", "font_fallbacks": ["Liberation Sans", "DejaVu Sans"],
    "body_font_size": 7.8, "atlas_page_size": [506.,365.],
    "max_rows_per_page": 60, "max_families_per_page": 64,
    "source_order": ["DefenseFinder","PADLOC","DefensePredictor"],
    "category_order": ["R-M","Modification","Nuclease","Cas","CRISPR-Cas","Abortive infection","Other","Other mechanisms","Unresolved"],
    "source_colors": {}, "category_colors": {}, "aliases": {}, "atlas_link_colors": {},
    "atlas_section_titles": {},
    "atlas_style": "manuscript", "atlas_title": None,
    "atlas_source_order": ["PADLOC","DefenseFinder","DefensePredictor"],
    "manuscript_workflow": {},
    "highlight_label": None, "sample_marker": None,
    "selection_label": "R subunits + nuclease effectors", "overlap_mode": "auto",
    "type_samples": None, "preview_png": False,
    "page_composition": "separate", "rm_fraction_scope": "assigned_I_IV", "composite_min_scale": .6}


def _integer(value, field):
    if value is not None and (isinstance(value,bool) or not isinstance(value,int) or value<0):
        raise ValueError(f"{field} must be a nonnegative integer or explicit null")


def validate_publication_model(result):
    if not isinstance(result,dict):raise ValueError("Publication result must be an analysis mapping")
    samples=[s.get("sample_id") for s in result.get("samples",[])]
    if any(not isinstance(s,str) or not s for s in samples) or len(samples)!=len(set(samples)):
        raise ValueError("Publication samples need unique nonempty sample_id values")
    sources=[s.get("id") for s in result.get("sources",[])]
    if any(not isinstance(s,str) or not s for s in sources) or len(sources)!=len(set(sources)):
        raise ValueError("Publication sources need unique nonempty id values")
    for scope,rows in [("categories",result.get("inventory_rows",{}).get("categories",[])),
                       ("families",result.get("inventory_rows",{}).get("families",[])),
                       ("cells",result.get("cells",[]))]:
        seen=set()
        for row in rows:
            key=(row.get("sample_id"),row.get("family_id",row.get("id"))) if scope=="cells" else row.get("id")
            if key in seen:raise ValueError(f"Duplicate {scope} identity: {key!r}")
            seen.add(key)
            if scope=="cells" and row.get("sample_id") not in samples:raise ValueError("Unknown cell sample")
            for field in ("N","G","associated_CDS","effector_CDS","motif_CDS"):_integer(row.get(field),field)
            n,a,e=row.get("N"),row.get("associated_CDS"),row.get("effector_CDS")
            if all(v is not None for v in (n,a,e)) and not 0<=e<=a<=n:
                raise ValueError("Publication counts require effector <= associated <= row union")
            for source,v in row.get("source_counts",{}).items():
                _integer(v,"source count")
                if source not in sources:raise ValueError(f"Unknown inventory source {source!r}")
                if v is not None and n is not None and v>n:raise ValueError("Source inclusion exceeds row union")
    for endpoint in ("overlap_regions","priority_overlap_regions"):
        overlap=result.get(endpoint,{})
        order=overlap.get("source_order",[]);seen=set();count=0
        if len(order)!=len(set(order)) or set(order)-set(sources):raise ValueError("Overlap needs unique declared source identities")
        for row in overlap.get("regions",[]):
            mask=row.get("mask");n=row.get("count")
            if not isinstance(mask,str) or len(mask)!=len(order) or set(mask)-{"0","1"}:
                raise ValueError("Invalid overlap membership mask")
            if mask in seen:raise ValueError("Duplicate overlap region")
            seen.add(mask);_integer(n,"overlap count")
            if n is None:raise ValueError("Overlap regions need known observed counts")
            count+=n
        if "union_count" in overlap and count!=overlap["union_count"]:
            raise ValueError("Overlap region sum differs from declared union count")
    if result.get("overlap_regions",{}).get("source_order",[]) != result.get("priority_overlap_regions",{}).get("source_order",[]):
        raise ValueError("Before/after overlap must preserve the same original source order")
    before=result.get("overlap_regions",{});after=result.get("priority_overlap_regions",{})
    if "sample_ids" in before and "sample_ids" in after and set(before["sample_ids"])!=set(after["sample_ids"]):
        raise ValueError("Before/after overlap must use the same sample scope")
    raw={r["mask"]:r["count"] for r in before.get("regions",[])}
    for r in after.get("regions",[]):
        if "1" in r["mask"] and r["count"]>raw.get(r["mask"],0):
            raise ValueError("Prioritized membership exceeds its original source intersection")
    seen_types=set()
    for row in result.get("type_tables",[]):
        if row.get("sample_id") not in samples:raise ValueError("Unknown type benchmark sample")
        key=(row["sample_id"],row.get("rm_type"))
        if key in seen_types:raise ValueError("Duplicate sample/type benchmark row")
        seen_types.add(key)
        for f in ("unique_CDS","pool_CDS","R","M","S","A","incomplete_R_context","orphan_M","orphan_S"):_integer(row.get(f),f)
        for source,n in row.get("source_counts",{}).items():
            _integer(n,"type source count")
            if source not in sources:raise ValueError("Unknown type source identity")
            if n is not None and row.get("pool_CDS") is not None and n>row["pool_CDS"]:raise ValueError("Type source count exceeds the shared typed pool")
        for field,total in (("incomplete_R_context","R"),("orphan_M","M"),("orphan_S","S")):
            if row.get(field) is not None and row.get(total) is not None and row[field]>row[total]:
                raise ValueError("Included type/context subset exceeds its role count")
    return {"samples":len(samples),"sources":len(sources),"activity_inferred":False}


def _config(title, options):
    if not isinstance(title,str) or not title.strip():raise ValueError("Publication title must be nonempty text")
    unknown=set(options)-set(DEFAULTS)
    if unknown:raise ValueError("Unknown publication layout options: "+", ".join(sorted(unknown)))
    config={**DEFAULTS,**options,"title":title}
    if config["page_composition"] not in ("separate","letter_composite"):raise ValueError("Unknown page_composition")
    if config["rm_fraction_scope"] not in ("assigned_I_IV","retained_I_III"):raise ValueError("Unknown rm_fraction_scope")
    scale=config["composite_min_scale"]
    if isinstance(scale,bool) or not isinstance(scale,(int,float)) or not math.isfinite(scale) or not 0<scale<=1:raise ValueError("composite_min_scale must be in (0,1]")
    if config["page_composition"]=="letter_composite" and config["atlas_style"]!="manuscript":raise ValueError("letter_composite requires manuscript atlas_style")
    if not isinstance(config["sections"],(list,tuple)) or set(config["sections"])-{"atlas","type_benchmarks","comparison"}:
        raise ValueError("Unknown publication section")
    if len(set(config["sections"]))!=len(config["sections"]):raise ValueError("Duplicate publication section")
    for field in ("max_rows_per_page","max_families_per_page"):
        if type(config[field]) is not int or config[field]<1:raise ValueError(field+" must be positive")
    if isinstance(config["body_font_size"],bool) or not isinstance(config["body_font_size"],(int,float)) or not math.isfinite(config["body_font_size"]) or config["body_font_size"]<7:
        raise ValueError("body_font_size must be finite and at least7 pt; paginate instead of shrinking")
    size=config["atlas_page_size"]
    if not isinstance(size,(list,tuple)) or len(size)!=2 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in size) or size[0]<270 or size[1]<145:
        raise ValueError("atlas_page_size needs finite dimensions at least270×145mm; smaller requests cannot contain the fixed-font workflow and evidence columns")
    if config["overlap_mode"] not in ("auto","intersection"):raise ValueError("overlap_mode must be auto or intersection")
    if config["atlas_style"] not in ("manuscript","compact"):raise ValueError("atlas_style must be manuscript or compact")
    if config["atlas_style"]=="manuscript" and (size[0]<447 or size[1]<200):
        raise ValueError("manuscript atlas_page_size needs at least447×200mm to preserve the five-card template; use compact for smaller pages")
    if config["atlas_title"] is not None and (not isinstance(config["atlas_title"],str) or not config["atlas_title"].strip()):raise ValueError("atlas_title must be nonempty text or null")
    if not isinstance(config["manuscript_workflow"],dict):raise ValueError("manuscript_workflow must be a mapping")
    if config["highlight_label"] is not None and (not isinstance(config["highlight_label"],str) or not config["highlight_label"].strip()):
        raise ValueError("highlight_label must be nonempty text or null")
    marker=config["sample_marker"]
    if marker is not None:
        if not isinstance(marker,dict) or set(marker)!={"field","label","color"} or any(not isinstance(v,str) or not v.strip() for v in marker.values()):
            raise ValueError("sample_marker needs nonempty field, label and color strings")
        from matplotlib.colors import to_rgba
        try:to_rgba(marker["color"])
        except ValueError as exc:raise ValueError("sample_marker color must be a valid color") from exc
    for field in ("source_colors","category_colors","aliases","atlas_section_titles","atlas_link_colors"):
        if not isinstance(config[field],dict):raise ValueError(field+" must be a mapping")
    for field in ("source_order","atlas_source_order","category_order","font_fallbacks"):
        if not isinstance(config[field],(list,tuple)) or any(not isinstance(v,str) or not v for v in config[field]) or len(set(config[field]))!=len(config[field]):raise ValueError(field+" needs unique nonempty strings")
    if config["type_samples"] is not None and (not isinstance(config["type_samples"],(list,tuple)) or any(not isinstance(s,str) for s in config["type_samples"])):raise ValueError("type_samples must be a sample-id list or null")
    if config["type_samples"] is not None and len(set(config["type_samples"])) != len(config["type_samples"]):raise ValueError("type_samples must contain distinct sample identities")
    if any(not isinstance(v,str) or not v for v in config["aliases"].values()):raise ValueError("Source aliases need nonempty strings")
    return config


def _size_note(config,physical):
    """Report measured physical sizes; never claim a print font floor."""
    def page(m):
        text=f"page {m['page']} at scale {m['export_scale']:.4f}, smallest text {m['min_text_pt']:.2f} pt" if m["min_text_pt"] is not None else f"page {m['page']} at scale {m['export_scale']:.4f}"
        return text+(f", thinnest line {m['min_stroke_pt']:.2f} pt" if m["min_stroke_pt"] is not None else "")
    measured="; ".join(page(m) for m in physical)
    if config["page_composition"]!="letter_composite":
        return "Measured physical sizes: "+measured+".\n\n"
    return ("Letter composition: atlas and type-table pages are scaled uniformly to fit Letter, with exported lines kept at "
            f"{STROKE_FLOOR_PT} pt or thicker. Measured sizes: "+measured+". Text below about 6 pt is small for print at Letter size; "
            "the default separate layout keeps native sizes on larger pages.\n\n")


def render_publication(result:dict, output_dir:Path, *, title="DefenseViz", external_type_cases=None, **layout_options):
    """Render arbitrary saved-result cohorts; return every generated artifact path.

    Every emitted element carries measured physical geometry and its source key.
    Data-dependent continuation pages preserve counts. Separate pages keep native
    fonts; optional Letter composition records its explicit uniform export scale.
    The parent pipeline owns transactionality; existing outputs are refused.
    """
    validate_publication_model(result)
    config=_config(title,layout_options)
    labels=[alias(s["id"],config) for s in result.get("sources",[])]
    if len(set(labels))!=len(labels):raise ValueError("Source aliases collide; supply distinct aliases")
    from .display_scope import type_case_projection
    type_result,type_selection,external_bindings=type_case_projection(result,external_type_cases,config["type_samples"],validate_publication_model,Path(output_dir).parent)
    type_config={**config,"type_samples":type_selection}
    model_json=json.dumps(result,sort_keys=True,ensure_ascii=False,allow_nan=False)
    model_hash=hashlib.sha256(model_json.encode("utf-8")).hexdigest()
    output_dir=Path(output_dir)
    if any(output_dir.glob("publication*")):raise FileExistsError("Publication output exists; use a fresh directory")
    fonts=Fonts(config["font_family"],config["font_fallbacks"])
    index=0
    def factory(width,height,section):
        nonlocal index
        index+=1
        return Canvas(width,height,fonts,section,index)
    pages=[]
    if config["page_composition"]=="letter_composite":
        from .composite import composite_pages
        pages=composite_pages(result,type_result,fonts,config,type_config,factory)
    for section in ([] if config["page_composition"]=="letter_composite" else config["sections"]):
        if section=="atlas":pages.extend((manuscript_atlas_pages if config["atlas_style"]=="manuscript" else atlas_pages)(result,fonts,config,factory))
        elif section=="type_benchmarks":pages.extend(type_pages(type_result,fonts,type_config,factory))
        else:pages.extend(comparison_pages(result,fonts,config,factory))
    if not pages:
        c=factory(215.9,279.4,"empty");c.text(8.45,12,title,12,bold=True)
        c.text(8.45,24,"No requested publication endpoints are available.",8);pages=[c]
    layouts=[];physical=[]
    for i,c in enumerate(pages,1):
        c.ledger["page_id"]=i
        layouts.append(c.validate())
        c.ledger["physical_metrics"]=c.physical_metrics()
        physical.append(dict(page=i,section=c.ledger.get("section"),**c.ledger["physical_metrics"]))
        if c.ledger.get("section") in ("atlas_composite","type_benchmark"):
            c.ledger["external_type_bindings"]=external_bindings
    external_labels=[s.get("label",s["sample_id"]) for s in type_result.get("samples",[])
                     if s.get("cohort")=="external" and any(s["sample_id"] in b["sample_ids"] for b in external_bindings)]
    output_dir.mkdir(parents=True,exist_ok=True)
    pdf=output_dir/"publication.pdf";outputs=[pdf]
    with matplotlib.rc_context({"svg.fonttype":"none","svg.hashsalt":"DefenseVizPublication",
                               "pdf.fonttype":42,"axes.unicode_minus":False}):
        with PdfPages(pdf,metadata={"Title":title,"Author":"DefenseViz","CreationDate":None,"ModDate":None}) as writer:
            for i,c in enumerate(pages,1):
                writer.savefig(c.fig)
                svg=output_dir/f"publication_page{i:03d}.svg"
                c.fig.savefig(svg,format="svg",metadata={"Title":title,"Date":None,"Creator":"DefenseViz"});outputs.append(svg)
                if config["preview_png"]:
                    png=output_dir/f"publication_page{i:03d}.png";c.fig.savefig(png,dpi=150);outputs.append(png)
                c.fig.clear()
    ledger=output_dir/"publication_layout.json"
    ledger.write_text(json.dumps({"schema_version":1,"units":"millimeters","bbox_convention":"x,y,width,height",
        "pages":layouts,"page_count":len(layouts),"input_model_sha256":model_hash,"config":config,
        "primary_input_model_sha256":model_hash,"external_type_bindings":external_bindings,
        "rm_fraction_scope":config["rm_fraction_scope"],"font_requested":fonts.preferred,"font_resolved":fonts.base,"font_files":list(fonts.used.values()),
        "source_labels":[{"id":s["id"],"label":s.get("label",s["id"]),"alias":alias(s["id"],config),"color":source_color(s["id"],config)} for s in result.get("sources",[])],
        "activity_inferred":False,"operation":"saved_result_publication_render"},indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    outputs.append(ledger)
    caption=output_dir/"publication_caption.md"
    coverage=[]
    for endpoint in ("overlap_regions","priority_overlap_regions"):
        item=result.get(endpoint,{})
        unassessed=item.get('coverage_status')=='unassessed'
        union_label='unavailable' if unassessed else f"{item.get('union_count',0):,}"
        outside_label='unavailable' if unassessed else str(item.get('outside_source_memberships',0))
        coverage.append(f"{endpoint}: observed union {union_label}; "
            f"coverage {item.get('coverage_status','unavailable').replace('unassessed','unavailable')}; "
            f"{len(item.get('sample_ids',[]))} included and {len(item.get('excluded_sample_ids',[]))} excluded samples. "
            f"Scope: {item.get('scope','not supplied')}. "
            f"Outside original source memberships: {outside_label}. "
            f"Priority CDSs outside the shared sample scope: {item.get('excluded_selected_count','unavailable')}.")
    caption.write_text("# DefenseViz publication figures\n\nDNMB annotations and source-bound evidence are summarized. "
        "Counts describe recorded annotations and support within the declared input scope, not exhaustive biological totals. Missing measurements use an em dash; pipeline coverage and unresolved candidates remain in the audit tables. Source agreement and reference homology do not establish biological activity.\n\n"
        "Atlas: category-colored disks encode binned CDS membership. A purple outline marks recorded putative nuclease-associated support; inner purple area is the supported effector/total fraction. An absent support mark does not establish absence of function. Pink rows are explicitly selected sample highlights. "
        "Effector CDSs include individually supported annotations and required cores of complete reviewed complexes; individual catalysis is not assigned to every core.\n\n"
        "R-M role rings mark support for at least one CDS; context rings mean every displayed group has supplied R support. A reference dot is separately curated context evidence. M(S)-only composition does not establish orphanhood. Intersection and assigned fractions retain separate numerator and denominator provenance.\n\n"
        "Type tables: Pool counts the CDSs placed in each R-M type by DNMB or, without a DNMB type, by a native R-M call. DF, PAD and DP count pool CDSs that each tool called in any defense family, so they show CDS inclusion rather than agreement on R-M type (for example, PADLOC PDC-S or DMS_other calls are counted). DNMB counts pool CDSs assigned to that type by DNMB. "
        "R(n) and M/S(n) give R CDSs with incomplete context and included orphan M/S CDSs within the count. (—) after a positive count means context or orphan status was not assessed; a count without parentheses has none recorded among assessed CDSs, and Type IV R context lies outside the Type I–III context rule. A bare em dash marks an unavailable or inapplicable value; gray cells are structurally inapplicable.\n\n"
        "Comparison: original source memberships are restricted to the declared priority set. Venn areas are schematic; intersection tables preserve exact memberships when Venn labels cannot fit. "
        "Category/family rows overlap. Equal bar tracks use each row union; associated totals and effector counts are separate columns.\n\n"
        +"Atlas R/total fraction scope: "+("strict retained Type I–III CDSs; reviewed Type IV excluded" if config["rm_fraction_scope"]=="retained_I_III" else "assigned Type I–IV CDSs, including supported reviewed Type IV components")+".\n\n"
        +(("External type cases (tan headers): "+"; ".join(external_labels)+". They are analyzed as separate projects and are not included in the atlas, overlap, family or workbook totals.\n\n") if external_bindings else "")
        +_size_note(config,physical)
        +"Priority definition: "+str(result.get("endpoint_definitions",{}).get("priority","declared analysis priority set; activity untested")).replace("unassessed","not evaluated")+".\n\n"
        +"\n\n".join(coverage)+"\n\n"
        +"Source labels: "+"; ".join(alias(s["id"],config)+": "+s.get("label",s["id"]) for s in result.get("sources",[]))+".\n",encoding="utf-8")
    outputs.append(caption)
    manifest=output_dir/"publication_manifest.json"
    manifest.write_text(json.dumps({"schema_version":1,"input_model_sha256":model_hash,"page_count":len(layouts),
        "outputs":[{"path":p.name,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in outputs],
        "resolved_font":fonts.base,"font_fallback":fonts.base!=fonts.preferred,"matplotlib_version":matplotlib.__version__,
        "external_type_bindings":external_bindings,"page_composition":config["page_composition"],"rm_fraction_scope":config["rm_fraction_scope"],
        "physical_metrics":physical,"letter_export_stroke_floor_pt":STROKE_FLOOR_PT if config["page_composition"]=="letter_composite" else None,
        "no_biological_inference":True},indent=2)+"\n",encoding="utf-8")
    outputs.append(manifest)
    return outputs

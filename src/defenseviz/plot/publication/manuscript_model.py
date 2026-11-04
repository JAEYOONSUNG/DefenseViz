"""Source-bound display projection for the detailed saved-results atlas.

This module adds no scientific assignments. Composition counts describe retained
members of explicit saved contexts; badges require separately bound display
evidence. They are not inferred from an effector flag or a reference hit.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import re


_TYPES = {"I": "Type I", "II": "Type II", "III": "Type III", "IV": "Type IV",
          "Type I": "Type I", "Type II": "Type II", "Type III": "Type III", "Type IV": "Type IV",
          "IIG": "Type II", "Type IIG": "Type II"}
_POSITIVE = {"supported", "accepted", "supported_annotation", "supported_putative_effector_annotation"}
_DISPLAY_STATES = {"supported", "not_supported", "unassessed"}
_CLASSES = ("I_MS", "I_RMS", "III_M", "III_MR")


def _key(row):
    k = row.get("sample_id"), row.get("gene_id")
    if any(not isinstance(x, str) or not x for x in k):
        raise ValueError("Atlas record requires exact sample_id and gene_id")
    return k


def _keys(values):
    return [list(k) for k in sorted(values)]


def _state(values):
    values = list(values)
    if values and all(x == "complete" for x in values):
        return "complete"
    return "partial" if any(x in {"complete", "partial"} for x in values) else "unassessed"


def _intersection_state(values):
    values = list(values)
    return _state(values) if values and all(x in {"complete", "partial"} for x in values) else "unassessed"


def _refs(row):
    refs = row.get("source_refs")
    if refs is None:
        refs = [row["source_ref"]] if row.get("source_ref") else []
    if not isinstance(refs, list) or not refs:
        raise ValueError("Atlas display evidence requires exact source references")
    for ref in refs:
        if (not isinstance(ref, dict) or not isinstance(ref.get("path"), str) or not ref["path"]
                or not re.fullmatch(r"[0-9a-f]{64}", ref.get("sha256", ""))):
            raise ValueError("Invalid atlas evidence source reference")
    return deepcopy(refs)


def _badge(states, mode):
    states = list(states)
    known = [x for x in states if x != "unassessed"]
    status = "complete" if states and len(known) == len(states) else "partial" if known else "unassessed"
    if mode == "all":
        value = False if "not_supported" in states else True if states and all(x == "supported" for x in states) else None
    else:
        value = True if "supported" in states else False if states and all(x == "not_supported" for x in states) else None
    return value, status


def _merge_display_state(previous, incoming, where):
    """Unassessed contributes no claim; incompatible known claims fail closed."""
    if incoming == "unassessed":
        return previous
    if previous not in {"unassessed", incoming}:
        raise ValueError("Contradictory atlas " + where + " evidence")
    return incoming


def _declared_channels(bundle):
    """Native final membership and reference matches are different endpoints."""
    native = {c["caller"] for c in bundle.get("calls", [])
              if c.get("evidence_kind", "final_member") == "final_member"}
    reference = {c["caller"] for c in bundle.get("calls", []) if c.get("evidence_kind") == "reference_hit"}
    coverage = bundle.get("coverage", [])
    native |= {c["caller"] for c in coverage if c.get("endpoint") == "native_calls"}
    reference |= {c["caller"] for c in coverage if c.get("endpoint") == "reference_matches"}
    # The legacy omitted endpoint means a final-call export except when the
    # provider is explicitly reference-only. Profile-only coverage is separate.
    native |= {c["caller"] for c in coverage if not c.get("endpoint") and
               (c["caller"] not in reference or c["caller"] in native)}
    return native, reference


def _channel_coverage(bundle, result, native, reference):
    fallback = {(c["sample_id"],c["caller"]):c["status"]
                for c in result.get("availability",{}).get("coverage",[])}
    channels = {"native_calls": {}, "reference_matches": {}}
    mixed = native & reference
    for endpoint, providers in [("native_calls", native), ("reference_matches", reference)]:
        channels[endpoint] = {k:v for k,v in fallback.items() if k[1] in providers and k[1] not in mixed}
    for c in bundle.get("coverage", []):
        endpoint = c.get("endpoint")
        if not endpoint:
            endpoint = "reference_matches" if c["caller"] in reference-native else "native_calls" if c["caller"] not in mixed else None
        if endpoint in channels:
            channels[endpoint][(c["sample_id"],c["caller"])] = c["status"]
    for c in bundle.get("calls", []):
        if c.get("evidence_kind", "final_member") == "final_member" and c.get("selected", True):
            channels["native_calls"].setdefault((c["sample_id"],c["caller"]), "partial")
    return channels["native_calls"], channels["reference_matches"]


def _retained_role_status(sample_id, typ, retained, assignments, saved_status):
    """Assess the displayed retained endpoint, not the broader native pool."""
    state = saved_status.get(sample_id, "unassessed")
    rows = [assignments[k] for k in retained if k[0] == sample_id]
    unknown = any(not any(t in _TYPES for t, role in rows_for_gene) or
                  (any(t == typ for t, role in rows_for_gene) and
                   not any(t == typ and role in {"R","M","S","RM","Aux","A","C","F"}
                           for t, role in rows_for_gene)) for rows_for_gene in rows)
    return "partial" if state == "complete" and unknown else state


def _taxonomy(samples, config):
    declared = config.get("taxonomy_legend", {})
    if not isinstance(declared, dict):
        raise ValueError("Atlas taxonomy_legend must be an object")
    result = {}
    for field, plural in [("genus","genera"),("species","species")]:
        items = declared.get(plural)
        if items is None:
            items=[];seen=set()
            for sample in samples:
                key=sample.get(field)
                if key and key not in seen:
                    display=sample.get("display",{})
                    items.append(dict(key=key,label=display.get(field+"_label",key),color=display.get(field+"_color")))
                    seen.add(key)
        if (not isinstance(items,list) or any(not isinstance(x,dict) or not x.get("key") or not x.get("label") for x in items)
                or len({x["key"]for x in items})!=len(items)):
            raise ValueError("Invalid atlas taxonomy legend records")
        result[plural]=[dict(deepcopy(x),genomes=sum(s.get(field)==x["key"]for s in samples))for x in items]
    return result


def _fallback(result):
    """Keep a usable, honest projection when raw saved input was not supplied."""
    samples = deepcopy(result.get("samples", []))
    native = [s["id"] for s in result.get("sources", []) if s.get("kind") in {"native_final_member", "mixed"}]
    references = [s["id"] for s in result.get("sources", []) if s.get("kind") in {"reference_annotation", "mixed"}]
    mixed = {s["id"] for s in result.get("sources", []) if s.get("kind") == "mixed"}
    source_order = [s for s in ("PADLOC", "DefenseFinder", "DefensePredictor") if s in native]
    source_order += [s for s in native if s not in source_order]
    cells = defaultdict(dict)
    for c in result.get("cells", []):
        cells[c["sample_id"]][c["family_id"]] = deepcopy(c)
    summaries = {r["sample_id"]: r for r in result.get("sample_summary", [])}
    rows = []
    for sample in samples:
        sid = sample["sample_id"]; s = summaries.get(sid, {})
        retained_total = s.get("retained_CDS")
        retained_scope = "all_saved_retained"
        retained_total_status = s.get("current_retained_counts_status",s.get("retained_counts_status",s.get("retained_ratio_status","unassessed")))
        if retained_total is None:
            retained_total = s.get("retained_RM_CDS")
            retained_scope = "saved_retained_RM_only" if retained_total is not None else "unavailable"
            retained_total_status = s.get("current_retained_RM_counts_status",s.get("retained_RM_counts_status",s.get("retained_ratio_status","unassessed")))
        if retained_total is None:
            retained_total_status = "unassessed"
        native_export_sources = s.get("rm_native_intersection_sources")
        compatible_native = (not mixed and isinstance(native_export_sources, list)
                             and len(native_export_sources) == len(source_order)
                             and set(native_export_sources) == set(source_order))
        numerator = s.get("rm_native_intersection_CDS") if compatible_native else None
        numerator_status = (s.get("rm_native_intersection_status", "unassessed")
                            if numerator is not None else "unassessed")
        pool_status = s.get("rm_pool_status", "unassessed")
        fraction_status = (_state([numerator_status, pool_status])
                           if numerator_status != "unassessed" else "unassessed")
        rows.append(dict(sample, display=deepcopy(sample.get("display", {})),
            rm_original=dict(source_counts={p:None if p in mixed else s.get("rm_source_counts", {}).get(p) for p in source_order},
                source_status={p:"unassessed" if p in mixed else s.get("rm_source_count_status", {}).get(p, "unassessed") for p in source_order},
                reference_counts={p:None if p in mixed else s.get("rm_source_counts", {}).get(p) for p in references},
                reference_status={p:"unassessed" if p in mixed else s.get("rm_source_count_status", {}).get(p,"unassessed") for p in references},input_pool_CDS=s.get("rm_pool_CDS", 0),
                pool_status=pool_status, intersection_CDS=numerator,
                intersection_count_status=numerator_status, intersection_status=fraction_status, intersection_sources=source_order),
            rm_roles={}, compositions={name:dict(count=0,status="unassessed",R_support_all=None,R_support_status="unassessed",
                reference_any=None,reference_status="unassessed",group_ids=[],member_keys=[]) for name in _CLASSES},
            retained=dict(R_CDS=s.get("retained_R_CDS"),total_CDS=retained_total,total_scope=retained_scope,
                R_status=s.get("current_retained_R_counts_status",s.get("retained_R_counts_status",s.get("retained_ratio_status","unassessed"))) if s.get("retained_R_CDS") is not None else "unassessed",
                total_status=retained_total_status,
                status=s.get("current_retained_ratio_status",s.get("retained_ratio_status","unassessed")) if retained_total is not None and s.get("retained_R_CDS") is not None else "unassessed",
                scope=s.get("current_retention_scope",retained_scope)),
            noncanonical=dict(source_counts={p:None if p in mixed else s.get("noncanonical_source_counts",{}).get(p) for p in source_order},
                source_status={p:"unassessed" if p in mixed else s.get("noncanonical_source_count_status",{}).get(p,"unassessed") for p in source_order},
                pool_CDS=s.get("noncanonical_CDS",0),pool_status="unassessed",masks=[],intersection_CDS=None,intersection_status="unassessed"),
            cells=cells[sid]))
    families = [deepcopy(f) for f in result.get("inventory_rows",{}).get("families",[])
                if f.get("category") != "R-M" and (f.get("N",0)>0 or f.get("counts_status") != "complete")]
    return dict(schema_version="manuscript-atlas-v1",samples=samples,taxonomy_legend=_taxonomy(samples,{}),source_order=source_order,reference_sources=references,rows=rows,
        families=families,group_records=[],provenance=dict(operation="result_only_display_projection",scientific_assignments_changed=False),
        warnings=["Original source memberships and all saved contexts require enrich_manuscript_model(bundle, result); unavailable channels remain unassessed."])


def build_manuscript_atlas_model(result):
    """Return the independent display model; never mutate the analysis result."""
    if not isinstance(result, dict):
        raise ValueError("Atlas result must be an object")
    model = result.get("manuscript_atlas")
    if model is None:
        return _fallback(result)
    if not isinstance(model, dict) or model.get("schema_version") != "manuscript-atlas-v1":
        raise ValueError("Invalid manuscript atlas model")
    return deepcopy(model)


def enrich_manuscript_model(bundle, result):
    """Attach a generic exact-set projection to a copied saved analysis result.

    Optional annotation records with evidence_kind rm_role_display or
    rm_context_display use a non-admitting reviewed_RM_component envelope:
    status=unassessed, include_in_union=false. Their display evidence states are
    explicit and independent of the analysis annotation status.
    """
    if not isinstance(bundle, dict) or not isinstance(result, dict):
        raise ValueError("Atlas enrichment requires canonical bundle and result objects")
    samples = deepcopy(result.get("samples", [])); sample_ids = [s["sample_id"] for s in samples]
    if set(sample_ids) != {s["sample_id"] for s in bundle.get("samples", [])}:
        raise ValueError("Atlas bundle/result sample identities differ")
    genes = {_key(g):g for g in bundle.get("genes", [])}
    ledger = {_key(g):g for g in result.get("gene_ledger", [])}; universe = set(ledger)
    for k, g in ledger.items():
        if k not in genes or (g.get("protein_sha256") and genes[k].get("protein_sha256") != g["protein_sha256"]):
            raise ValueError("Atlas bundle/result gene binding differs")
    config = (bundle.get("options") or {}).get("manuscript_atlas", {})
    if not isinstance(config, dict):
        raise ValueError("options.manuscript_atlas must be an object")
    sources = [s["id"] for s in result.get("sources", [])]
    source_sets = {s["source_id"]:{tuple(k) for k in s["member_keys"]} for s in result.get("source_sets", [])}
    native_sources, references = _declared_channels(bundle)
    default_order = [s for s in ("PADLOC", "DefenseFinder", "DefensePredictor") if s in native_sources]
    default_order += [s for s in sources if s in native_sources and s not in default_order]
    order = config.get("source_order", default_order)
    if (not isinstance(order,list) or len(order)!=len(set(order)) or any(s not in native_sources for s in order)):
        raise ValueError("Atlas source_order requires unique original native sources")
    reference_order = [s for s in sources if s in references]
    coverage, reference_coverage = _channel_coverage(bundle, result, native_sources, references)
    saved = {_key(r):r for r in bundle.get("saved_decisions", [])}
    retained = {k for k,g in ledger.items() if g.get("retained")}
    raw_rm = defaultdict(set); reference_rm = defaultdict(set)
    observed = defaultdict(set); native_sets = defaultdict(set); generic = defaultdict(set)
    for c in bundle.get("calls", []):
        if c.get("evidence_kind", "final_member") == "hmm_hit" or not c.get("selected",True):
            continue
        k = _key(c); provider = c["caller"]; observed[provider].add(k)
        kind = c.get("evidence_kind", "final_member")
        if kind == "final_member":
            native_sets[provider].add(k)
        if kind == "reference_hit":
            # The separate REBASE channel is original reference coverage;
            # an untyped reference is not silently reclassified as an R-M call.
            reference_rm[provider].add(k)
        elif c.get("category") == "R-M" or (c.get("rm_type") in _TYPES):
            raw_rm[provider].add(k)
        if kind == "final_member" and c.get("family") == "General defense":
            generic[provider].add(k)
    if any(observed[s] != source_sets.get(s,set()) for s in sources):
        raise ValueError("Atlas original source memberships disagree with analysis")
    scope = result.get("rm_candidate_scope", {})
    candidate_scope = {tuple(k) for k in scope.get("typed_member_keys", scope.get("member_keys", saved))
                       if (saved.get(tuple(k), {}).get("rm_type") or saved.get(tuple(k), {}).get("family")) in _TYPES}
    if not candidate_scope <= set(genes):
        raise ValueError("Atlas current R-M candidate scope contains unknown identities")
    raw_pool = candidate_scope | set().union(*raw_rm.values(), *reference_rm.values())
    for provider, keys in generic.items():
        raw_rm[provider] |= keys & raw_pool
    current_rm = set().union(*(set(map(tuple,r["member_keys"])) for r in result.get("inventory_rows",{}).get("categories",[]) if r.get("category") == "R-M"))
    noncanonical = universe - current_rm
    # Annotations already adjudicated by analyze_saved supply roles; the display
    # never upgrades a native family label into a retained/reviewed component.
    assignments = defaultdict(set)
    for k,g in ledger.items():
        records = ([saved[k]] if k in retained and k in saved else []) + [a for a in g.get("annotations",[]) if a.get("status") in _POSITIVE and a.get("assignment_kind") == "reviewed_RM_component"]
        for a in records:
            typ = _TYPES.get(a.get("rm_type") or a.get("family"))
            role = {"R_component":"R","A":"Aux"}.get(a.get("role"),a.get("role"))
            if typ:
                assignments[k].add((typ,role))
    context_assessments = {}
    for assessment in result.get("context_assessments", []):
        ids = assessment.get("group_ids", [assessment.get("group_id")])
        if not isinstance(ids, list) or not ids or len(set(ids)) != len(ids):
            raise ValueError("Atlas context assessment requires unique group_ids")
        for gid in ids:
            key = assessment.get("sample_id"), gid
            if key in context_assessments:
                raise ValueError("Duplicate atlas current context assessment")
            context_assessments[key] = assessment
    group_coverage = config.get("group_coverage",{})
    if not isinstance(group_coverage,dict) or any(k not in sample_ids or v not in {"complete","partial","unassessed"} for k,v in group_coverage.items()):
        raise ValueError("Invalid atlas saved-group coverage")
    current_review_members = defaultdict(set)
    for g in bundle.get("groups", []):
        members = {(g["sample_id"], m["gene_id"]) for m in g.get("members", [])}
        if (g.get("provider") == "Review" and g.get("membership_kind") == "saved_context"
                and members and members <= retained):
            current_review_members[g["sample_id"], _TYPES.get(g.get("rm_type") or g.get("family"))] |= members
    unresolved_group_types = set()
    for r in (bundle.get("suite") or {}).get("unresolved_contexts", []):
        sid, typ = r["sample_id"], _TYPES.get(r.get("family"))
        member = sid, r.get("gene_id")
        if member not in genes or (member in retained and not
                (group_coverage.get(sid) == "complete" and member in current_review_members[sid,typ])):
            unresolved_group_types.add((sid,typ))
    group_index = {}; group_rows = []; dedup = {}
    skipped_review_groups = []
    for g in bundle.get("groups", []):
        if g.get("membership_kind") != "saved_context" or g.get("provider") not in {"DNMB", "Review"}:
            continue
        sid = g["sample_id"]; gid = g["group_id"]; gkey = sid,gid
        if gkey in group_index:
            raise ValueError("Duplicate atlas saved context identity")
        original = {(sid,m["gene_id"]) for m in g.get("members",[])}
        if g.get("provider") == "Review" and (not original or not original <= retained):
            # The retained R-M panel cannot turn a general reviewed integrated
            # context into a DNMB-retained context by silently dropping members.
            skipped_review_groups.append(gkey)
            continue
        members = {(sid,m["gene_id"]) for m in g.get("members", [])
                   if m.get("retained") is True or ("retained" not in m and (sid,m["gene_id"]) in retained)}
        if not members <= retained:
            raise ValueError("Atlas saved context retained member conflicts with saved decision")
        typ = _TYPES.get(g.get("rm_type") or g.get("family"))
        assessment = context_assessments.get(gkey)
        if assessment is not None:
            if (assessment.get("provider") != g["provider"] or assessment.get("run_id") != g.get("run_id")
                    or set(map(tuple, assessment.get("original_member_keys", []))) != original
                    or set(map(tuple, assessment.get("member_keys", []))) != members):
                raise ValueError("Atlas current context assessment differs from exact saved group")
            context_state = assessment.get("context_state", "unassessed")
        else:
            context_state = g.get("context_state", "unassessed")
        current_composition = not (assessment and assessment.get("assessment_scope") == "pre_retention_context_not_reassessed")
        if (not current_composition and members and not
                (group_coverage.get(sid) == "complete" and members <= current_review_members[sid,typ])):
            # Keep old boundaries in the audit. Current source-bound review
            # contexts establish their own scope; this does not merge groups
            # or assert completeness of a displaced historical boundary.
            unresolved_group_types.add((sid,typ))
        unknown_member_role = typ is None or any(
            not any(t == typ and role in {"R","M","S","RM","Aux","C","F"}
                    for t,role in assignments[k]) for k in members)
        roles = set().union(*({"R","M"} if role == "RM" else {role}
                             for k in members for t,role in assignments[k] if t == typ))
        settled_class = ((typ == "Type I" and {"R","M","S"} <= roles) or
                         (typ == "Type III" and {"R","M"} <= roles))
        unknown_composition = unknown_member_role and not settled_class
        if current_composition and unknown_composition and members:
            unresolved_group_types.add((sid,typ))
        cls = None if unknown_composition else (("I_RMS" if "R" in roles else "I_MS") if typ == "Type I" and {"M","S"} <= roles else ("III_MR" if "R" in roles else "III_M") if typ == "Type III" and "M" in roles else None)
        row = dict(sample_id=sid,group_id=gid,group_ids=[gid],provider=g["provider"],run_id=g.get("run_id"),family=typ,
            composition_class=cls,member_keys=_keys(members),original_member_keys=_keys(original),
            current_composition=bool(current_composition),
            composition_status="unassessed" if unknown_composition or not current_composition else "complete",
            member_role_status="partial" if unknown_member_role else "complete",
            context_state=context_state,original_context_state=g.get("context_state","unassessed"),
            context_origin="reviewed_context" if g.get("provider") == "Review" else "native_saved_context",
            origin=g.get("origin"),rm_R_support="unassessed",rm_reference_link="unassessed",
            source_refs=deepcopy(g.get("source_refs",[])),evidence_refs=[])
        # Reciprocal exported partner rows can describe the same exact context.
        # Require known run and identical ORIGINAL members; do not merge merely
        # adjacent or overlapping groups, or groups sharing only retained genes.
        ident = (sid,g["provider"],g.get("run_id"),typ,tuple(sorted(original)))
        if g.get("run_id") and original and ident in dedup:
            previous = dedup[ident]
            if row["context_state"] != previous["context_state"] or row["member_keys"] != previous["member_keys"]:
                raise ValueError("Contradictory exact saved-context assessments")
            row = previous; row["group_ids"].append(gid)
        else:
            group_rows.append(row); dedup[ident] = row
        group_index[gkey] = row
    role_evidence = {}; display_refs = []
    for a in bundle.get("annotations", []):
        kind = a.get("evidence_kind")
        if kind not in {"rm_role_display","rm_context_display"}:
            continue
        if a.get("assignment_kind") != "reviewed_RM_component" or a.get("status") != "unassessed" or a.get("include_in_union") is not False:
            raise ValueError("Atlas display evidence requires a non-admitting annotation envelope")
        k = _key(a); sha = a.get("protein_sha256")
        if k not in genes or not isinstance(sha,str) or not re.fullmatch(r"[0-9a-f]{64}",sha) or genes[k].get("protein_sha256") != sha:
            raise ValueError("Atlas display evidence protein binding mismatch")
        refs = _refs(a); display_refs.extend(refs)
        names = ["rm_R_support"] if kind == "rm_role_display" else ["rm_R_support","rm_reference_link"]
        for name in names:
            if a.get(name,"unassessed") not in _DISPLAY_STATES:
                raise ValueError("Invalid atlas display evidence state")
        if kind == "rm_role_display":
            if k not in retained or not any(role in {"R","RM"} for _,role in assignments[k]):
                raise ValueError("R display support must bind a retained R/RM component")
            state = a.get("rm_R_support","unassessed")
            role_evidence[k] = _merge_display_state(role_evidence.get(k,"unassessed"), state, "R support")
        else:
            gkey = k[0],a.get("group_id")
            if gkey not in group_index or list(k) not in group_index[gkey]["member_keys"]:
                raise ValueError("Atlas group evidence must bind an exact retained group member")
            row = group_index[gkey]
            if not row["current_composition"] and any(a.get(n) == "supported" for n in names):
                raise ValueError("Pre-retention context cannot receive a current support/reference badge")
            if row["composition_class"] not in {"I_RMS","III_MR"} and any(a.get(n)=="supported" for n in names):
                raise ValueError("Complete R-M support/reference badge cannot extend to M/MS-only composition")
            if a.get("rm_reference_link") == "supported" and row["context_state"] not in {"complete","met","met_in_saved_context"}:
                raise ValueError("Curated reference badge requires an explicitly complete saved context")
            for name in names:
                state = a.get(name,"unassessed")
                row[name] = _merge_display_state(row[name], state, "group display")
            row["evidence_refs"].extend(refs)
    summaries = {s["sample_id"]:s for s in result.get("sample_summary",[])}
    availability = result.get("availability",{})
    saved_status = {s["sample_id"]:s["status"] for s in availability.get(
        "current_saved_decisions_by_sample", availability.get("saved_decisions_by_sample",[]))}
    saved_status.update({sid:s["current_retained_counts_status"] for sid,s in summaries.items()
                         if "current_retained_counts_status" in s})
    type_rows = {(s["sample_id"],s["rm_type"]):s for s in result.get("type_tables",[])}
    cells = defaultdict(dict)
    for c in result.get("cells",[]):
        cells[c["sample_id"]][c["family_id"]] = deepcopy(c)
    absence_families = config.get("absence_families",[])
    known_families = {f["id"]for f in result.get("inventory_rows",{}).get("families",[])}
    if not isinstance(absence_families,list) or len(set(absence_families))!=len(absence_families) or any(f not in known_families for f in absence_families):
        raise ValueError("Atlas absence_families must name explicit catalog families")
    rows=[]
    for sample in samples:
        sid=sample["sample_id"]; summary=summaries[sid]
        source_status={p:coverage.get((sid,p),"missing") for p in order}
        source_status={p:_state([v]) for p,v in source_status.items()}
        intersection = set.intersection(*(raw_rm[p] for p in order)) if order else set()
        numerator=sum(k[0]==sid for k in intersection)
        count_status=_intersection_state(source_status.values()); pool_status=summary.get("rm_pool_status","unassessed")
        native_export_sources = summary.get("rm_native_intersection_sources")
        if (not set(order) & references and isinstance(native_export_sources,list)
                and len(native_export_sources) == len(order) and set(native_export_sources) == set(order)):
            if (summary.get("rm_native_intersection_CDS") != numerator or
                    summary.get("rm_native_intersection_status") != count_status):
                raise ValueError("Atlas native intersection disagrees with analysis native endpoint")
        fraction_status=_state([count_status,pool_status]) if count_status!="unassessed" else "unassessed"
        rm_original=dict(source_counts={p:sum(k[0]==sid for k in raw_rm[p]) for p in order},source_status=source_status,
            reference_counts={p:sum(k[0]==sid for k in reference_rm[p]) for p in reference_order},
            reference_status={p:_state([reference_coverage.get((sid,p),"missing")]) for p in reference_order},
            input_pool_CDS=sum(k[0]==sid for k in raw_pool),pool_status=pool_status,intersection_CDS=numerator,
            intersection_count_status=count_status,intersection_status=fraction_status,intersection_sources=list(order))
        rm_roles={}
        for code,names in [("I",["R","M","S"]),("II",["R","M"]),("III",["R","M"]),("IV",["R","A"])]:
            typ="Type "+code; table=type_rows.get((sid,typ),{})
            for name in names:
                role="Aux" if name=="A" else name
                ks={k for k,aa in assignments.items() if k[0]==sid and (code=="IV" or k in retained) and any(t==typ and (v==role or v=="RM" and role in {"R","M"}) for t,v in aa)}
                status=(table.get("role_states",{}).get(name,"unassessed") if code=="IV" else
                        _retained_role_status(sid,typ,retained,assignments,saved_status))
                states=[role_evidence.get(k,"unassessed") for k in ks]
                evidence_status=_badge(states,"all")[1] if name=="R" and code!="IV" else "not_applicable"
                if name=="R" and code!="IV" and not ks and status=="complete":evidence_status="complete"
                rm_roles[code+"_"+name]=dict(count=len(ks),status=status,member_keys=_keys(ks),support_count=sum(x=="supported" for x in states),support_status=evidence_status,
                    status_scope="saved_retained_or_reviewed_TypeIV_roles" if code=="IV" else "saved_retained_roles")
        comps={}
        for cls in _CLASSES:
            gs=[g for g in group_rows if g["sample_id"]==sid and g["composition_class"]==cls and g["current_composition"]]
            status=group_coverage.get(sid,"partial" if gs else "unassessed")
            typ = "Type I" if cls.startswith("I_") else "Type III"
            if ((sid,typ) in unresolved_group_types or (sid,None) in unresolved_group_types) and status == "complete":
                status="partial"
            if saved_status.get(sid) in {"partial","unassessed"} and status == "complete":
                status="partial"
            supported,support_state=_badge([g["rm_R_support"]for g in gs],"all")
            linked,reference_state=_badge([g["rm_reference_link"]for g in gs],"any")
            comps[cls]=dict(count=len(gs),status=status,R_support_all=supported,R_support_status=support_state,
                reference_any=linked,reference_status=reference_state,group_ids=[g["group_id"]for g in gs],
                member_keys=_keys(set().union(*(set(map(tuple,g["member_keys"]))for g in gs))))
        nc={k for k in noncanonical if k[0]==sid}; masks=defaultdict(set)
        for k in nc:masks[''.join('1'if k in native_sets[p]else'0'for p in order)].add(k)
        mask_status=_intersection_state(source_status.values())
        mask_ids=[format(i,f"0{len(order)}b")for i in range(2**len(order))] if order else [""]
        nc_status=_state(source_status.values())
        if nc and nc_status=="unassessed":nc_status="partial"
        noncanonical_row=dict(source_counts={p:len(nc&native_sets[p])for p in order},source_status=source_status,
            pool_CDS=len(nc),pool_status=nc_status,masks=[dict(mask=m,count=len(masks[m]),status=mask_status,member_keys=_keys(masks[m]))for m in mask_ids],
            intersection_CDS=len(masks['1'*len(order)]) if order else 0,intersection_status=mask_status)
        for cell in cells[sid].values():
            members=set(map(tuple,cell["member_keys"]))
            cell["exclusive_source_counts"]={p:len(members & masks[''.join('1'if q==p else'0'for q in order)])for p in order}
            cell["exclusive_source_status"]=mask_status
        display=deepcopy(sample.get("display",{}))
        display["absence_family_ids"]=[f for f in absence_families if f in cells[sid] and cells[sid][f]["N"]==0 and cells[sid][f].get("counts_status")=="complete"]
        display["absence_label"]=config.get("absence_label","Observed absence")
        rows.append(dict(sample,display=display,rm_original=rm_original,rm_roles=rm_roles,compositions=comps,
            retained=dict(R_CDS=summary["retained_R_CDS"],total_CDS=summary["retained_CDS"],
                R_status=summary.get("current_retained_R_counts_status",summary.get("retained_R_counts_status","unassessed")),
                total_status=summary.get("current_retained_counts_status",summary.get("retained_counts_status","unassessed")),
                status=summary.get("current_retained_ratio_status",summary.get("retained_ratio_status","unassessed")),
                scope=summary.get("current_retention_scope","native saved candidate assessment")),
            noncanonical=noncanonical_row,cells=cells[sid]))
    families=[deepcopy(f)for f in result.get("inventory_rows",{}).get("families",[])if f.get("category")!="R-M" and (f.get("N",0)>0 or f.get("counts_status")!="complete")]
    model=dict(schema_version="manuscript-atlas-v1",samples=samples,taxonomy_legend=_taxonomy(samples,config),absence_label=config.get("absence_label","saved family"),source_order=list(order),reference_sources=reference_order,
        rows=rows,families=families,group_records=group_rows,
        provenance=dict(operation="exact_source_display_projection",scientific_assignments_changed=False,
            source_metadata=deepcopy(result.get("provenance",{}).get("source_metadata",[])),display_evidence_refs=display_refs),
        warnings=["R-M support and curated reference badges are independent display evidence, not effector or activity assignments.",
                  "Saved composition groups are direct explicit contexts; no proximity pairing or transitive expansion."] +
                 (["Current saved contexts with unresolved type or member role are preserved without guessed classes; affected composition counts are partial."] if unresolved_group_types else []) +
                 (["Reviewed contexts containing nonretained members are excluded from the retained R-M context panel; their integrated-system evidence is preserved."] if skipped_review_groups else []))
    enriched=deepcopy(result);enriched["manuscript_atlas"]=model
    return enriched

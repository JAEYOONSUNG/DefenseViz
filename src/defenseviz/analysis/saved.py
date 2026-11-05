"""Exact aggregation of saved DNMBsuite results, independent of legacy scoring.

No detector, sequence search, score predicate, proximity grouping, or activity
inference occurs here. Optional reviews and native component transfers retain
explicit evidence classes and source bindings.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
from importlib.resources import files
import json
import re

SUPPORTED = {"supported", "accepted", "supported_annotation", "supported_putative_effector_annotation"}
EVIDENCE_STATES = SUPPORTED | {"unassessed", "unsupported", "not_supported", "rejected", "incomplete", "conflicting", "out_of_scope", "not_applicable"}
DECISIONS = {"retained", "not_retained", "unassessed", "out_of_scope"}
COVERAGE = {"complete", "partial", "missing", "not_run"}
ROLES = {"R", "M", "S", "RM", "Aux", "A", "C", "F", "unknown", "Unknown"}
CONTEXT_MET = {"met", "complete", "met_in_saved_context"}
CONTEXT_MISSING = {"missing_retained_partner_in_saved_context", "missing"}


def _text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f"{name} must be nonempty text")
    return v


def _rows(v, name):
    if v is None:
        return []
    if not isinstance(v, (list, tuple)) or any(not isinstance(r, dict) for r in v):
        raise ValueError(f"{name} must be a list of objects")
    return v


def _sha(v, where, required=False):
    if v is None and not required:
        return None
    if not isinstance(v, str) or not re.fullmatch(r"[0-9a-f]{64}", v):
        raise ValueError(f"{where} must be lowercase SHA256")
    return v


def _refs(row, where, required=False):
    refs = row.get("source_refs")
    if refs is None:
        ref = row.get("source_ref")
        refs = [ref] if ref is not None else []
    _rows(refs, where + ".source_refs")
    if required and not refs:
        raise ValueError(f"{where} requires exact source_refs")
    for r in refs:
        _text(r.get("path"), where + ".source_ref.path")
        _sha(r.get("sha256"), where + ".source_ref.sha256", True)
    return deepcopy(refs)


def _keys(values):
    return [list(k) for k in sorted(values)]


def _role(v):
    return "Aux" if v == "A" else v


def _type(v):
    # This is an explicit nomenclature map, not a family-name classifier.
    # In particular, non-RM sentinels must never enter a hidden fifth type pool.
    return {"I": "Type I", "II": "Type II", "III": "Type III", "IV": "Type IV",
            "Type I": "Type I", "Type II": "Type II", "Type III": "Type III", "Type IV": "Type IV",
            "IIG": "Type II", "Type IIG": "Type II"}.get(v, "Unresolved")


def _status(states):
    values = list(states)
    if values and all(v == "complete" for v in values):
        return "complete"
    if any(v in {"complete", "partial"} for v in values):
        return "partial"
    return "unassessed"


def _intersection_status(states):
    values = list(states)
    return _status(values) if values and all(v in {"complete", "partial"} for v in values) else "unassessed"


def _registry():
    raw = files("defenseviz").joinpath("data/components.json").read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def analyze_saved(data):
    """Return source-traceable sets and plot tables from canonical saved inputs.

    Count fields are observed numbers; their companion status must be used to
    distinguish complete counts, lower bounds, and unassessed endpoints.
    """
    if not isinstance(data, dict):
        raise ValueError("saved data must be an object")
    d = deepcopy(data)
    samples = list(_rows(d.get("samples"), "samples"))
    sample_ids = []
    for r in samples:
        sid = _text(r.get("sample_id"), "samples.sample_id")
        _text(r.get("label"), "samples.label")
        if sid in sample_ids:
            raise ValueError(f"Duplicate sample_id {sid!r}")
        if "highlight" in r and type(r["highlight"]) is not bool:
            raise ValueError("sample highlight must be boolean")
        sample_ids.append(sid)
    sample_set = set(sample_ids)
    tables = {n: list(_rows(d.get(n), n)) for n in (
        "genes", "calls", "coverage", "saved_decisions", "annotations",
        "effector_evidence", "groups", "complexes", "source_metadata")}
    options = d.get("options") or {}
    if not isinstance(options, dict):
        raise ValueError("options must be an object")
    warnings = [("Stored R-M evidence reevaluated with the selected DNMB policy; no detector run."
                 if d.get("rm_evaluation") else "Saved-result import and aggregation; no detector or scoring run."),
                "Annotation evidence does not establish activity. Missing evidence is unassessed.",
                "Family/category rows can overlap; sums are not a unique CDS union."]

    def key(r, where):
        sid = _text(r.get("sample_id"), where + ".sample_id")
        if sid not in sample_set:
            raise ValueError(f"{where}: unknown sample_id {sid!r}")
        return sid, _text(r.get("gene_id"), where + ".gene_id")

    genes = {}
    for i, r in enumerate(tables["genes"]):
        k = key(r, f"genes[{i}]")
        if k in genes:
            raise ValueError(f"Duplicate gene identity {k!r}")
        _sha(r.get("protein_sha256"), "gene protein_sha256")
        if any(f in r for f in ("start", "end")):
            if any(type(r.get(f)) is not int for f in ("start", "end")) or not 1 <= r["start"] <= r["end"]:
                raise ValueError("Gene coordinates require 1-based inclusive start <= end")
            _text(r.get("contig"), "gene contig")
            if r.get("strand") not in {"+", "-"}:
                raise ValueError("Gene strand must be '+' or '-'")
        genes[k] = r
    explicit_genes = bool(genes)

    def bind(r, where, review=False):
        k = key(r, where)
        sha = _sha(r.get("protein_sha256"), where + ".protein_sha256", review)
        if k not in genes:
            if explicit_genes or review:
                raise ValueError(f"{where}: unknown gene {k!r}")
            genes[k] = {"sample_id": k[0], "gene_id": k[1], "protein_sha256": sha}
        current = genes[k].get("protein_sha256")
        if current and sha and current != sha:
            raise ValueError(f"{where}: protein hash mismatch for {k!r}")
        if review and not current:
            raise ValueError(f"{where}: review requires current gene protein hash")
        _refs(r, where, review)
        return k

    calls = []
    seen_calls = set()
    sources = []
    source_kind = defaultdict(set)
    source_members = defaultdict(set)
    native = set()
    call_by_gene = defaultdict(list)
    for i, r in enumerate(tables["calls"]):
        k = bind(r, f"calls[{i}]")
        provider = _text(r.get("caller"), "call caller")
        kind = r.get("evidence_kind", "final_member")
        if kind not in {"final_member", "reference_hit", "hmm_hit"}:
            raise ValueError("Unsupported saved call evidence_kind")
        if "selected" in r and type(r["selected"]) is not bool:
            raise ValueError("call selected must be boolean")
        if r.get("role") is not None and r["role"] not in ROLES:
            raise ValueError("Unsupported call role")
        ident = (k, provider, r.get("run_id"), r.get("system_id"), r.get("family"),
                 r.get("native_model"), kind, r.get("reference_id"))
        if ident in seen_calls:
            raise ValueError(f"Duplicate native observation {ident!r}")
        seen_calls.add(ident)
        r.update(evidence_kind=kind)
        calls.append(r); call_by_gene[k].append(r)
        if provider not in sources:
            sources.append(provider)
        if kind != "hmm_hit":
            if provider not in sources:
                sources.append(provider)
            source_kind[provider].add(kind)
            if r.get("selected", True):
                source_members[provider].add(k)
                if kind == "final_member":
                    native.add(k)
    cov = {}
    decision_coverage = {}
    for r in tables["coverage"]:
        sid = _text(r.get("sample_id"), "coverage sample_id")
        if sid not in sample_set:
            raise ValueError("Unknown coverage sample")
        provider = _text(r.get("caller"), "coverage caller")
        if r.get("status") not in COVERAGE:
            raise ValueError("Invalid coverage status")
        if r.get("endpoint") in {"saved_retention", "saved_decisions"}:
            if sid in decision_coverage and decision_coverage[sid] != r["status"]:
                raise ValueError("Conflicting duplicate saved-decision coverage")
            decision_coverage[sid] = r["status"]
            continue
        if r.get("endpoint") in {"hmm_hits", "profile_hits", "hmm_search", "raw_hmm"}:
            continue
        if r.get("endpoint") == "reference_matches":
            source_kind[provider].add("reference_hit")
        elif r.get("endpoint") == "native_calls" or not source_kind[provider]:
            source_kind[provider].add("final_member")
        ck = sid, provider
        if ck in cov and cov[ck] != r["status"]:
            raise ValueError("Conflicting duplicate source coverage; select one active run")
        cov[ck] = r["status"]
        if provider not in sources:
            sources.append(provider)
    for r in calls:
        if r["evidence_kind"] == "final_member" and r.get("selected", True):
            state = cov.get((r["sample_id"], r["caller"]), "missing")
            if state in {"missing", "not_run"}:
                # Unlike a claimed explicit missing artifact, omitted coverage
                # can accompany observed records; preserve a lower bound.
                if (r["sample_id"], r["caller"]) in cov:
                    raise ValueError("Native final member conflicts with missing/not_run coverage")
                cov[(r["sample_id"], r["caller"])] = "partial"
    order = options.get("source_order", sources)
    if not isinstance(order, list) or len(set(order)) != len(order) or any(x not in sources for x in order):
        raise ValueError("options.source_order must contain unique known sources")
    sources = list(order) + [s for s in sources if s not in order]

    saved = {}; retained = set(); decision_assessed = set()
    for i, r in enumerate(tables["saved_decisions"]):
        k = bind(r, f"saved_decisions[{i}]")
        if k in saved:
            raise ValueError("Duplicate saved decision: choose one active decision run")
        if r.get("state") not in DECISIONS:
            raise ValueError("Invalid saved decision state")
        if r.get("role") is not None and r["role"] not in ROLES:
            raise ValueError("Unsupported saved decision role")
        if "integrated_R_component_support" in r and type(r["integrated_R_component_support"]) is not bool:
            raise ValueError("integrated_R_component_support must be an explicit boolean")
        saved[k] = r
        if r["state"] != "unassessed":
            decision_assessed.add(k)
        if r["state"] == "retained":
            retained.add(k)
    capabilities = (d.get("suite") or {}).get("capabilities", [])
    capabilities = capabilities if isinstance(capabilities, list) else []
    empty_retention_assessed = {r.get("sample_id") for r in capabilities if r.get("saved_retention") is True and r.get("native_candidates") is True and r.get("candidate_count") == 0}
    empty_retention_assessed |= {sid for sid,state in decision_coverage.items() if state == "complete" and not any(k[0] == sid for k in saved)}
    saved_states = []
    for sid in sample_ids:
        rows = [r for k,r in saved.items() if k[0] == sid]
        known = sum(r["state"] != "unassessed" for r in rows)
        state = "complete" if (rows and known == len(rows)) or sid in empty_retention_assessed else "partial" if known else "unassessed"
        # Observed decisions cannot override an explicitly incomplete export.
        if sid in decision_coverage and decision_coverage[sid] != "complete":
            state = "partial" if known else "unassessed"
        saved_states.append(dict(sample_id=sid,status=state,record_count=len(rows),states=dict(Counter(r["state"] for r in rows))))
    saved_status_by_sample = {r["sample_id"]:r["status"] for r in saved_states}
    saved_availability = _status(r["status"] for r in saved_states)
    # A declared current selection is a different assessment universe from the
    # broader native candidate audit. Completeness here describes that exact
    # observed selection, never the exhaustiveness of discovery.
    current_scopes = {}
    for scope in (d.get("rm_evaluation") or {}).get("candidate_scope", []):
        if "gene_ids" not in scope:
            continue
        sid, ids = scope.get("sample_id"), scope["gene_ids"]
        if (sid not in sample_set or sid in current_scopes or not isinstance(ids, list)
                or any(not isinstance(g, str) or not g for g in ids) or len(set(ids)) != len(ids)):
            raise ValueError("Current candidate scope requires unique known samples and gene IDs")
        members = {(sid, g) for g in ids}
        if not members <= set(genes):
            raise ValueError("Current candidate scope contains unknown genes")
        current_scopes[sid] = members
    current_saved_states = []
    for sid in sample_ids:
        members = current_scopes.get(sid, {k for k in saved if k[0] == sid})
        known = members & decision_assessed
        state = ("complete" if known == members else "partial" if known else "unassessed")
        if sid not in current_scopes:
            state = saved_status_by_sample[sid]
        current_saved_states.append(dict(sample_id=sid, status=state, record_count=len(members),
            member_keys=_keys(members), unassessed_member_keys=_keys(members-known),
            scope="explicit current observed candidate selection" if sid in current_scopes else "native saved candidate assessment",
            exhaustive_discovery_claim=False))
    annotation_by_gene = defaultdict(list); reviewed = set(); blocked_effector = set()
    annotation_ids = set()
    for i, r in enumerate(tables["annotations"]):
        k = bind(r, f"annotations[{i}]", True)
        _text(r.get("rule_id"), "annotation rule_id")
        if r.get("assignment_kind") not in {"reviewed_RM_component", "noncatalytic_component"}:
            raise ValueError("Invalid annotation assignment_kind")
        if r.get("status") not in EVIDENCE_STATES:
            raise ValueError("Invalid annotation status")
        ident = (k, r["assignment_kind"], r["rule_id"], r.get("family"), r.get("role"))
        if ident in annotation_ids:
            raise ValueError("Duplicate reviewed annotation")
        annotation_ids.add(ident)
        if r.get("role") is not None and r["role"] not in ROLES | {"R_component"}:
            raise ValueError("Unsupported reviewed role")
        if "include_in_union" in r and type(r["include_in_union"]) is not bool:
            raise ValueError("include_in_union must be boolean")
        if r["assignment_kind"] == "noncatalytic_component" and r.get("include_in_union") is True:
            raise ValueError("Noncatalytic component adjudication is a veto, not an included membership annotation")
        annotation_by_gene[k].append(r)
        if r["status"] in SUPPORTED:
            if r["assignment_kind"] == "noncatalytic_component":
                blocked_effector.add(k)
            elif r.get("include_in_union", False):
                reviewed.add(k)
    for k, arr in annotation_by_gene.items():
        assignments_now = {(_type(r.get("rm_type") or r.get("family")), _role(r.get("role")))
                           for r in arr if r["status"] in SUPPORTED and r["assignment_kind"] == "reviewed_RM_component"}
        if len(assignments_now) > 1:
            raise ValueError("Contradictory reviewed role/type assignments require explicit resolution")
    U = native | retained | reviewed

    # Explicit evidence is separate from native model functional transfer.
    individual = set(); site = set(); other_role = set(); role_assessed = set()
    evidence_by_gene = defaultdict(list); site_grades = defaultdict(list)
    seen_evidence = set()
    for i, r in enumerate(tables["effector_evidence"]):
        k = bind(r, f"effector_evidence[{i}]", True)
        _text(r.get("rule_id"), "effector rule_id")
        if r.get("assignment_kind") not in {"individual_nuclease_effector", "nuclease_other_role", "uncertain_domain_function"}:
            raise ValueError("Invalid effector assignment_kind")
        if r.get("status") not in EVIDENCE_STATES:
            raise ValueError("Invalid effector evidence status")
        ident = (k, r["assignment_kind"], r["rule_id"])
        if ident in seen_evidence:
            raise ValueError("Duplicate effector evidence")
        seen_evidence.add(ident)
        if k not in U:
            raise ValueError("Effector evidence outside integrated universe; supply explicit membership annotation")
        evidence_by_gene[k].append(r)
        if r["status"] in SUPPORTED:
            role_assessed.add(k)
            if r["assignment_kind"] == "individual_nuclease_effector":
                if k in blocked_effector:
                    raise ValueError("Conflicting noncatalytic and individual-effector reviews require explicit resolution")
                individual.add(k)
            elif r["assignment_kind"] == "nuclease_other_role":
                other_role.add(k)
            if r.get("motif_state") in SUPPORTED:
                if not r.get("motif_grade"):
                    raise ValueError("Supported motif/site requires original motif_grade")
                site_grades[k].append(r["motif_grade"])
                if r["assignment_kind"] == "individual_nuclease_effector":
                    site.add(k)
    registry, registry_sha = _registry()
    automatic = []
    component_display = {}
    component_display_metadata = {}
    native_component = defaultdict(lambda: defaultdict(set))
    component_observations = []
    for r in calls:
        if r["evidence_kind"] != "final_member" or not r.get("selected", True) or r.get("origin") == "supplement":
            continue
        k = key(r, "native component")
        model = r.get("native_model")
        if not model or not _refs(r, "native component"):
            continue
        for rule in registry["component_rules"]:
            if r["caller"] != rule["provider"] or model not in rule["exact_models"]:
                continue
            native_component[k][rule["component"]].add(rule["rule_id"])
            if rule.get("display_family"):
                component_display[k] = (rule["display_family"], rule["display_category"])
                component_display_metadata[rule["display_family"]] = rule.get("display_metadata",{})
            component_observations.append((k, r["caller"], r.get("run_id"), r.get("system_id"), rule["component"]))
            if k in blocked_effector:
                continue
            assignment = rule["assignment_kind"]
            item = dict(sample_id=k[0], gene_id=k[1], rule_id=rule["rule_id"],
                        assignment_kind=assignment, status="supported",
                        evidence_grade="native_final_selected_specific_component_annotation",
                        source_refs=_refs(r, "native component"), native_model=model,
                        substrate=rule.get("substrate", "unresolved"),
                        primary_sources=rule["primary_sources"], motif_state="unassessed",
                        activity="unassessed", automatic=True)
            automatic.append(item); evidence_by_gene[k].append(item)
            if assignment == "individual_nuclease_effector":
                individual.add(k); role_assessed.add(k)
            elif assignment == "nuclease_other_role":
                other_role.add(k); role_assessed.add(k)
    # Final source-derived systems may be supplied separately. Derive them only
    # from explicit original IDs when canonical groups were omitted.
    groups = list(tables["groups"])
    if not groups:
        gs = {}
        for r in calls:
            if r["evidence_kind"] != "final_member" or not r.get("selected", True) or r.get("system_id") is None:
                continue
            raw_key = [r["sample_id"], r["caller"], r.get("run_id"), r.get("replicon_id"), r["system_id"]]
            gid = json.dumps(raw_key, ensure_ascii=False, separators=(",", ":"))
            g = gs.setdefault(gid, dict(group_id=gid, sample_id=r["sample_id"], provider=r["caller"], run_id=r.get("run_id"), membership_kind="native_system", native_system_id=r["system_id"], system_id=r["system_id"], members=[], source_refs=[]))
            k = key(r, "call group")
            if not any(m["gene_id"] == k[1] for m in g["members"]):
                g["members"].append(dict(gene_id=k[1], protein_sha256=genes[k].get("protein_sha256")))
            for ref in _refs(r, "call group"):
                if ref not in g["source_refs"]:
                    g["source_refs"].append(ref)
        groups = list(gs.values())
    group_index = {}; group_original = {}; group_current = {}
    for i, g in enumerate(groups):
        gid = _text(g.get("group_id"), f"groups[{i}].group_id")
        sid = _text(g.get("sample_id"), "group sample_id")
        group_key = (sid, gid)
        if sid not in sample_set or group_key in group_index:
            raise ValueError("Duplicate group_id or unknown group sample")
        if g.get("membership_kind") not in {"native_system", "saved_context"}:
            raise ValueError("Invalid group membership_kind")
        _text(g.get("provider"), "group provider")
        if g.get("origin") == "supplement" and g["membership_kind"] == "native_system":
            raise ValueError("Supplement cannot claim an original native system")
        review_context = g["membership_kind"] == "saved_context" and g["provider"] == "Review"
        group_refs = _refs(g, "group", review_context)
        native_identity_verified = g["membership_kind"] == "native_system" and bool(g.get("run_id")) and bool(g.get("native_system_id", g.get("system_id")))
        members = set(); current = set()
        for member in _rows(g.get("members"), "group members"):
            if "sample_id" in member and member["sample_id"] != sid:
                raise ValueError("Cross-sample group member")
            m = dict(member, sample_id=sid)
            if review_context:
                if "retained" in member:
                    raise ValueError("Review group cannot claim member retention")
                # The exact group artifact binds the listed protein identities.
                m["source_refs"] = group_refs
            k = bind(m, "group member", review_context)
            if k in members:
                raise ValueError("Duplicate member in one group")
            members.add(k)
            if g["membership_kind"] == "native_system":
                candidates = [c for c in call_by_gene[k] if c["caller"] == g["provider"] and c["evidence_kind"] == "final_member"
                              and c.get("selected", True) and c.get("origin") != "supplement"]
                if not candidates:
                    raise ValueError("Native group member lacks an original final call from its provider")
                native_id = g.get("native_system_id", g.get("system_id"))
                def compatible(c):
                    return all(a is None or b is None or a == b for a,b in
                               ((g.get("run_id"),c.get("run_id")), (native_id,c.get("system_id"))))
                if not any(compatible(c) for c in candidates):
                    raise ValueError("Native group member conflicts with its exact provider/run/system")
                native_identity_verified = native_identity_verified and any(
                    c.get("run_id") == g.get("run_id") and c.get("system_id") == native_id
                    for c in candidates)
                current.add(k)
            elif review_context:
                if k in U:
                    current.add(k)
            elif member.get("retained") is True or ("retained" not in member and k in retained):
                if k not in retained:
                    raise ValueError("Saved context retained member conflicts with saved decision")
                current.add(k)
        if g.get("origin") == "supplement" and g["membership_kind"] == "native_system":
            raise ValueError("Supplement cannot claim an original native system")
        g["native_identity_binding"] = "verified" if native_identity_verified and members else "unassessed"
        g["membership_scope"] = ("original_native_members" if g["membership_kind"] == "native_system" else
                                 "reviewed_current_integrated_members" if review_context else "retained_saved_context_members")
        group_index[group_key] = g; group_original[group_key] = members; group_current[group_key] = current
    complexes = list(tables["complexes"])
    # Versioned exact component requirements operate only inside the same native
    # group. They never infer a new boundary or transfer catalysis to all members.
    for gid, g in group_index.items():
        if (g["membership_kind"] != "native_system" or g.get("origin") == "supplement"
                or g.get("replicon_binding") == "unassessed" or g["native_identity_binding"] != "verified"):
            continue
        for rule in registry["complex_rules"]:
            native_id = g.get("native_system_id", g.get("system_id"))
            if native_id is None:
                continue
            found = {c: {k for k,provider,run,system,component in component_observations
                         if component == c and k in group_current[gid] and provider == g["provider"]
                         and run == g.get("run_id") and system == native_id}
                     for c in rule["required_components"]}
            if not all(len(v) == 1 for v in found.values()):
                continue
            chosen = set.union(*found.values())
            if any(len({c for c, ks in found.items() if k in ks}) > 1 for k in chosen):
                continue
            refs = _refs(g, "native complex")
            if not refs:
                continue
            complexes.append(dict(complex_id="registry:" + rule["rule_id"] + ":" + g["group_id"],
                group_id=g["group_id"], sample_id=g["sample_id"], status="complete",
                core_gene_ids=[k[1] for k in sorted(chosen)], source_refs=refs,
                rule_id=rule["rule_id"], primary_sources=rule["primary_sources"], automatic=True))
    cores = set(); module_groups = set(); seen_complexes = set(); complex_out = []
    for i, c in enumerate(complexes):
        cid = _text(c.get("complex_id"), "complex_id")
        cid_key = (c.get("sample_id"), cid)
        if cid_key in seen_complexes:
            raise ValueError("Duplicate complex_id")
        seen_complexes.add(cid_key)
        if c.get("status") not in {"complete", "incomplete", "unassessed"}:
            raise ValueError("Invalid complex status")
        _text(c.get("rule_id"), "complex rule_id")
        _refs(c, "complex", True)
        gid = (c.get("sample_id"), c.get("group_id"))
        if gid not in group_index or group_index[gid]["sample_id"] != c.get("sample_id"):
            raise ValueError("Complex does not bind its exact original group")
        ids = c.get("core_gene_ids")
        if not isinstance(ids, list) or not ids or len(ids) != len(set(ids)) or any(not isinstance(x, str) or not x for x in ids):
            raise ValueError("Complex requires unique explicit core_gene_ids")
        required = {(c["sample_id"], x) for x in ids}
        if "core_protein_sha256" in c:
            hashes = c["core_protein_sha256"]
            if not isinstance(hashes,dict) or set(hashes) != set(ids):
                raise ValueError("Complex core hash map must bind every core identity")
            for k in required:
                _sha(hashes[k[1]],"complex core protein_sha256",True)
                if genes.get(k,{}).get("protein_sha256") != hashes[k[1]]:
                    raise ValueError("Complex core protein hash mismatch")
        if not required <= group_original[gid] or not required <= U:
            raise ValueError("Complex core outside exact group or integrated universe")
        if c["status"] == "complete":
            if not required <= group_current[gid]:
                raise ValueError("Complete complex includes a nonretained saved-context member")
            cores |= required; module_groups.add(gid)
        complex_out.append(dict(c, individual_catalysis="unassigned_for_complex_core"))
    E = individual | cores
    site &= E
    evaluation = d.get("rm_evaluation") or {}
    prior_unassessed = {(r["sample_id"], r["gene_id"]) for r in evaluation.get("native_saved_decisions", [])
                        if r.get("state") == "unassessed"}
    replaced_samples = ({s["sample_id"] for s in evaluation.get("candidate_scope", [])}
                        if evaluation.get("decision_mode") == "replace_scoped" else set())
    pre_retention_groups = {gid for gid, g in group_index.items() if
        g.get("provider") == "DNMB" and g.get("membership_kind") == "saved_context" and bool(evaluation)
        and (g["sample_id"] in replaced_samples or bool(group_original[gid] & prior_unassessed))}
    selected_groups = []; G = set()
    for gid, g in group_index.items():
        if gid in pre_retention_groups and gid not in module_groups:
            # Old candidate context is retained in the audit, but cannot expand
            # a current nuclease-associated set without an explicit new review.
            continue
        own = group_current[gid] & individual
        if not own and gid not in module_groups:
            continue
        included = group_current[gid] & U
        G |= included
        selected_groups.append(dict(group_id=g["group_id"], sample_id=g["sample_id"], provider=g["provider"],
            membership_kind=g["membership_kind"], member_keys=_keys(included), seed_keys=_keys(own),
            membership_scope=g["membership_scope"], native_identity_binding=g["native_identity_binding"], origin=g.get("origin"),
            complete_reviewed_complex=gid in module_groups,
            outside_union_member_keys=_keys(group_current[gid] - U), source_refs=_refs(g, "group")))
    A = G | E; ungrouped = E - G

    # Roles are assignments, not evidence of complete systems. One fused RM is
    # one unique CDS and two explicit component memberships when rendered.
    assignments = defaultdict(set); role_records = defaultdict(list)
    for k, r in saved.items():
        if k in retained:
            typ = _type(r.get("rm_type") or r.get("family"))
            assignments[k].add((typ, _role(r.get("role"))))
            role_records[k].append(r)
    for k, arr in annotation_by_gene.items():
        for r in arr:
            if r["status"] in SUPPORTED and r["assignment_kind"] == "reviewed_RM_component" and k in U:
                role = "R" if r.get("role") == "R_component" else _role(r.get("role"))
                assignments[k].add((_type(r.get("rm_type") or r.get("family")), role))
                role_records[k].append(r)
    # A source-bound review may fill a missing role/type, but conflicting known
    # assignments require explicit resolution rather than a pooled fused role.
    for k, values in assignments.items():
        known_types = {t for t,role in values if t != "Unresolved"}
        known_roles = {role for t,role in values if role not in {None,"unknown","Unknown"}}
        if len(known_types) > 1 or len(known_roles) > 1:
            raise ValueError("Contradictory saved/reviewed role or type assignments require explicit resolution")
    R = {k for k, a in assignments.items() if any(role in {"R", "RM"} for _, role in a)}
    final_rm = {k for k, a in assignments.items() if any(t in {"Type I", "Type II", "Type III", "Type IV"} for t, _ in a)}
    # Annotation-based triage and the narrower nuclease-route support predicate
    # are distinct endpoints. An unmet route is not a negative R annotation.
    unmet_R = {k for k in R if saved.get(k, {}).get("integrated_R_component_support") is False}
    policy_supported_R = {k for k in R if saved.get(k, {}).get("integrated_R_component_support") is True}
    independently_reviewed_R = {k for k in R if any(
        r["status"] in SUPPORTED and r["assignment_kind"] == "reviewed_RM_component"
        and _role(r.get("role")) in {"R", "RM"} for r in annotation_by_gene[k])}
    support_qualified_R = policy_supported_R | independently_reviewed_R
    support_qualified_P = E | support_qualified_R
    P = E | R

    # A saved context describes its original membership. If filtering removes
    # members, its original positive assessment cannot silently describe a new
    # retained context. Explicit review contexts already bind the current union.
    # Preserve both states rather than rewriting the imported source artifact.
    context_assessments = []
    current_context_state = {}
    unique_contexts = {}
    for gid, g in group_index.items():
        if g["membership_kind"] != "saved_context":
            continue
        typ = _type(g.get("rm_type") or g.get("family"))
        original_state = g.get("context_state")
        state = original_state
        removed = group_original[gid] - group_current[gid]
        pre_retention = gid in pre_retention_groups
        requirements = g.get("required_roles")
        if requirements is not None:
            if (not isinstance(requirements, (list, tuple)) or not requirements
                    or any(r not in ROLES - {"unknown", "Unknown"} for r in requirements)):
                raise ValueError("Saved context required_roles must name known roles")
            requirements = set().union(*({"R", "M"} if r == "RM" else {_role(r)} for r in requirements))
        else:
            requirements = {"Type I": {"R", "M", "S"}, "Type III": {"R", "M"}}.get(typ)
        if pre_retention:
            state = "unassessed"
        elif removed and original_state in CONTEXT_MET and g.get("provider") != "Review":
            # Type IV has no inferred methylation-partner requirement.
            if typ != "Type IV" or requirements is not None:
                current_roles = set()
                unknown = False
                for k in group_current[gid]:
                    own = {role for family, role in assignments[k] if family == typ
                           and role in ROLES - {"unknown", "Unknown"}}
                    unknown |= not own
                    current_roles |= set().union(*({"R", "M"} if r == "RM" else {r} for r in own))
                if requirements is not None and requirements <= current_roles:
                    pass
                elif requirements is not None and not unknown:
                    state = "missing_retained_partner_in_saved_context"
                else:
                    state = "unassessed"
        current_context_state[gid] = state
        record = dict(sample_id=g["sample_id"], group_id=g["group_id"],
            group_ids=[g["group_id"]], provider=g["provider"], run_id=g.get("run_id"),
            family=typ, original_context_state=original_state, context_state=state,
            required_roles=sorted(requirements) if requirements is not None else None,
            original_member_keys=_keys(group_original[gid]), member_keys=_keys(group_current[gid]),
            removed_member_keys=_keys(removed), source_refs=_refs(g, "saved context"),
            assessment_scope=("pre_retention_context_not_reassessed" if pre_retention else
                "explicit_reviewed_current_context" if g.get("provider") == "Review" else "current_retained_context"))
        # Reciprocal exported rows describe one context only when their source
        # run, original boundary and type are identical. Different boundaries
        # remain distinct even if filtering leaves the same retained members.
        ident = (g["sample_id"], g["provider"], g.get("run_id"), typ,
                 tuple(sorted(group_original[gid]))) if g.get("run_id") and group_original[gid] else gid
        if ident in unique_contexts:
            previous = unique_contexts[ident]
            if any(previous[field] != record[field] for field in
                   ("original_context_state", "context_state", "member_keys", "required_roles")):
                raise ValueError("Contradictory exact saved-context assessments")
            previous["group_ids"].append(g["group_id"])
            for ref in record["source_refs"]:
                if ref not in previous["source_refs"]:
                    previous["source_refs"].append(ref)
        else:
            unique_contexts[ident] = record
            context_assessments.append(record)
    # Explicit family catalog and identity mappings are presentation membership.
    # They must not establish any of E, R, site, or activity above.
    catalog = list(_rows(options.get("family_catalog"), "family_catalog"))
    catalog_index = {}; families = []; pools = defaultdict(set)
    def family_add(fid, category, label=None, metadata=None):
        if fid not in catalog_index:
            defaults = registry.get("display_catalog_defaults", {}).get(fid, {})
            display_label = defaults.get("label",fid) if label in {None,fid} else label
            catalog_index[fid] = dict(id=fid, label=display_label, category=category or "Unresolved",
                                     **{k:v for k,v in defaults.items() if k != "label"})
            catalog_index[fid].update(metadata or {})
            families.append(fid)
        elif category and catalog_index[fid]["category"] != category:
            raise ValueError(f"Conflicting display categories for family {fid!r}")
    for r in catalog:
        fid = _text(r.get("id") or r.get("family_id"), "catalog family id")
        if fid in catalog_index:
            raise ValueError("Duplicate family catalog ID")
        family_add(fid, r.get("category") or r.get("category_id") or "Unresolved", r.get("label"),
                   {k:deepcopy(v) for k,v in r.items() if k not in {"id","family_id","category","category_id","label"}})
    memberships = options.get("family_memberships")
    if memberships is not None:
        seen_m = set()
        for m in _rows(memberships, "family_memberships"):
            k = bind(m, "family membership", True)
            if k not in U:
                raise ValueError("Display family membership outside integrated universe")
            fid = _text(m.get("family_id"), "membership family_id")
            if (k, fid) in seen_m:
                raise ValueError("Duplicate family membership")
            seen_m.add((k, fid))
            family_add(fid, m.get("category") or (catalog_index.get(fid) or {}).get("category") or "Unresolved")
            pools[fid].add(k)
        if set.union(set(), *pools.values()) != U:
            raise ValueError("Explicit family memberships must cover every integrated CDS")
    else:
        display_raw = json.loads(files("defenseviz").joinpath("data/families.json").read_text())
        aliases = {(r.get("caller"), r["family"]): (r.get("display_family",r["family"]),r["category"]) for r in display_raw}
        aliases.update({(r["provider"],r["native_family"]):(r["family"],r["category"]) for r in registry.get("display_aliases",[])})
        generic = set()
        for r in calls:
            k = key(r, "family call")
            if r["evidence_kind"] != "final_member" or not r.get("selected", True) or k not in U:
                continue
            fid = _text(r.get("family") or "Unresolved", "call family")
            category = r.get("category") or "Unresolved"
            if fid in {"General defense","Unresolved","Unknown"}:
                generic.add(k); continue
            mapped = None if r.get("category_source") == "family_map" else aliases.get((r["caller"],fid),aliases.get((None,fid)))
            if mapped:
                fid,category = mapped
            if category == "R-M" or _type(r.get("rm_type")) in {"Type I","Type II","Type III","Type IV"}:
                fid,category = "R-M","R-M"
            elif k in component_display:
                fid,category = component_display[k]
            meta = component_display_metadata.get(fid) or next((a.get("display_metadata",{}) for a in registry.get("display_aliases",[]) if a["family"] == fid),{})
            family_add(fid,category,metadata=meta);pools[fid].add(k)
        for k in U:
            if k in final_rm:
                family_add("R-M", "R-M"); pools["R-M"].add(k)
            if not any(k in st for st in pools.values()):
                fid = "General defense only" if k in generic else "Unresolved"
                family_add(fid, "Unresolved"); pools[fid].add(k)
    # All configured families, including observed-zero catalog entries, remain in
    # scientific output. Renderer may hide only explicitly known empty entries.
    category_pools = defaultdict(set)
    for fid in families:
        category_pools[catalog_index[fid]["category"]] |= pools[fid]
    native_sources = [s for s in sources if "final_member" in source_kind[s] or any(r["caller"] == s and r["evidence_kind"] == "hmm_hit" for r in calls)]
    # Inventory and functional coverage are scoped independently to each sample.
    # An explicit declaration describes review extent, never biological absence.
    axes = ("effector", "motif", "associated")
    assessment = options.get("assessment_coverage", {})
    scoped_assessment = options.get("assessment_coverage_by_sample", {})
    def validate_assessment(value, where):
        if not isinstance(value, dict) or set(value) - set(axes):
            raise ValueError(f"{where} must contain only effector, motif and associated")
        if any(not isinstance(v, str) or v not in {"complete", "partial", "unassessed"} for v in value.values()):
            raise ValueError(f"Invalid {where}")
    validate_assessment(assessment, "assessment_coverage")
    if not isinstance(scoped_assessment, dict) or set(scoped_assessment) - sample_set:
        raise ValueError("assessment_coverage_by_sample contains unknown sample IDs or is not an object")
    for sid, declaration in scoped_assessment.items():
        validate_assessment(declaration, f"assessment_coverage_by_sample[{sid!r}]")
    inventory_by_sample = {}
    assessment_by_sample = {}
    assessment_provenance = {}
    for sid in sample_ids:
        inventory_by_sample[sid] = _status(cov.get((sid, s), "missing") for s in native_sources)
        if not native_sources and any(k[0] == sid for k in retained | reviewed):
            inventory_by_sample[sid] = "partial"
        default = {"effector": "partial" if any(k[0] == sid for k in E | other_role) else "unassessed",
                   "motif": "partial" if any(k[0] == sid for k in site) else "unassessed",
                   "associated": "partial" if any(k[0] == sid for k in A) else "unassessed"}
        local = scoped_assessment.get(sid, {})
        assessment_by_sample[sid] = {}
        assessment_provenance[sid] = {}
        for axis in sorted(axes):
            if axis in local:
                value, origin = local[axis], "sample_declaration"
                option_path = f"options.assessment_coverage_by_sample[{sid!r}].{axis}"
            elif axis in assessment:
                value, origin = assessment[axis], "explicit_global_fallback"
                option_path = f"options.assessment_coverage.{axis}"
            else:
                value, origin, option_path = default[axis], "local_observed_evidence_default", None
            declarations = [m for m in tables["source_metadata"]
                            if m.get("source_kind") == "assessment_coverage_declaration"
                            and m.get("sample_id") == sid and m.get("declaration_scope") ==
                            ("sample" if origin == "sample_declaration" else "global")
                            and m.get("axes", {}).get(axis) == value] if option_path else []
            assessment_by_sample[sid][axis] = value
            assessment_provenance[sid][axis] = dict(origin=origin, option_path=option_path,
                source_refs=[ref for m in declarations for ref in _refs(m, "assessment coverage provenance")])
    inventory_state = _status(inventory_by_sample.values())
    role_global = _status(v["effector"] for v in assessment_by_sample.values())
    motif_global = _status(v["motif"] for v in assessment_by_sample.values())
    associated_state = _status(v["associated"] for v in assessment_by_sample.values())
    def aggregate(st, sid_scope=None):
        ids = sid_scope if sid_scope is not None else sample_ids
        sc = {s:len(st & source_members[s]) for s in sources}
        ss = {s:_status(cov.get((sid,s), "missing") for sid in ids) for s in sources}
        local_inventory = _status(inventory_by_sample[sid] for sid in ids)
        metric_status = {axis:_status(assessment_by_sample[sid][axis] for sid in ids) for axis in axes}
        return dict(member_keys=_keys(st), N=len(st), G=len({k[0] for k in st}), genome_count=len({k[0] for k in st}),
            source_counts=sc, source_count_status=ss, associated_CDS=len(st&A), effector_CDS=len(st&E), motif_CDS=len(st&site),
            individual_effector_CDS=len(st&individual), complex_core_CDS=len(st&cores), ungrouped_effector_CDS=len(st&ungrouped),
            counts_status=local_inventory, metric_status=dict(inventory=local_inventory, **metric_status))
    family_rows = []; cells = []
    for fid in families:
        meta = catalog_index[fid]
        z = dict(meta, **aggregate(pools[fid])); z["visible"] = bool(pools[fid]) or inventory_state != "complete"
        family_rows.append(z)
        for sid in sample_ids:
            st = {k for k in pools[fid] if k[0] == sid}
            cells.append(dict(sample_id=sid, family_id=fid, category=meta["category"], **aggregate(st, [sid])))
    category_rows = [dict(id=c, label=c, category=c, **aggregate(st)) for c,st in category_pools.items()]
    # Scope the overlap to samples with every configured source complete. Raw
    # observed sets remain separately available; incomplete absence is not a bit0.
    requested = options.get("comparison_sample_ids", sample_ids)
    if not isinstance(requested, list) or len(requested) != len(set(requested)) or not set(requested) <= sample_set:
        raise ValueError("Invalid comparison_sample_ids")
    common = {sid for sid in requested if sources and all(cov.get((sid,s)) == "complete" for s in sources)}
    def overlap(selected=None):
        scope = {k for k in (set.union(set(), *(source_members[s] for s in sources)) if selected is None else selected) if k[0] in common}
        buckets = defaultdict(set)
        for k in scope:
            buckets[''.join('1' if k in source_members[s] else '0' for s in sources)].add(k)
        mask_ids = [format(i, f"0{len(sources)}b") for i in range(1,2**len(sources))] if 0 < len(sources) <= 8 else sorted(buckets)
        if selected is not None and '0'*len(sources) in buckets:
            mask_ids = ['0'*len(sources)] + mask_ids
        return dict(source_order=sources, regions=[dict(mask=m,count=len(buckets[m]),member_keys=_keys(buckets[m])) for m in mask_ids],
            union_count=len(scope), source_sizes={s:len(scope&source_members[s]) for s in sources},
            sample_ids=[s for s in sample_ids if s in common], excluded_sample_ids=[s for s in requested if s not in common],
            excluded_selected_count=sum(k[0] not in common for k in selected) if selected is not None else None,
            coverage_status="complete" if common == set(requested) and requested else "partial" if common else "unassessed",
            outside_source_memberships=len(buckets['0'*len(sources)]) if sources else len(scope),
            scope="common samples with complete configured source endpoints; source agreement is not validation")
    # Case/type roles use unique source intersections, never pooled role sums.
    type_tables = []; components = []
    native_types = defaultdict(set)
    for call in calls:
        if call["evidence_kind"] == "final_member" and call.get("selected", True):
            typ = _type(call.get("rm_type"))
            if typ != "Unresolved" or call.get("category") == "R-M":
                native_types[key(call, "native type")].add(typ)
    type_pool = set(native_types) | final_rm
    type_bucket = {}
    type_disagreement = {}
    for k in type_pool:
        current = {t for t,role in assignments[k] if t in {"Type I","Type II","Type III","Type IV"}}
        observed = native_types[k] - {"Unresolved"}
        type_bucket[k] = next(iter(current)) if len(current) == 1 else next(iter(observed)) if not current and len(observed) == 1 else "Unresolved"
        type_disagreement[k] = len(observed) > 1 or bool(current and observed and current != observed)
    types = [t for t in ["Type I","Type II","Type III","Type IV","Unresolved"] if t in set(type_bucket.values())]
    types += ["Total"] if type_pool else []
    for sid in sample_ids:
        for typ in types:
            pool = {k for k in type_pool if k[0] == sid and (typ == "Total" or type_bucket[k] == typ)}
            st = pool & final_rm
            roles = {r:{k for k in st if any((typ == "Total" or type_bucket[k] == typ) and (v == r or v == "RM" and r in {"R","M"}) for t,v in assignments[k])} for r in ("R","M","S","Aux")}
            disposition_known = pool & (decision_assessed | set(role_records))
            roles_known = all(any(v in {"R","M","S","RM","Aux","A","C","F"} for t,v in assignments[k]) for k in st)
            row_assessment = ("complete" if roles_known and (disposition_known == pool if pool else saved_status_by_sample[sid] == "complete")
                              else "partial" if disposition_known else "unassessed")
            counts_status = ("complete" if disposition_known == pool and saved_status_by_sample[sid] == "complete"
                             else "partial" if disposition_known else "unassessed")
            pool_status = _status(cov.get((sid, s), "missing") for s in native_sources)
            if pool and pool_status == "unassessed":
                pool_status = "partial"
            # The shared-pool table conditions on these explicitly observed,
            # assigned-type members. Unassessed candidates outside this pool
            # still affect the broader retention endpoint, not this subset.
            within_pool_counts_status = ("complete" if disposition_known == pool and pool_status == "complete"
                                         else "partial" if disposition_known or not pool and pool_status == "partial"
                                         else "unassessed")
            within_pool_role_status = (within_pool_counts_status if roles_known else
                                       "partial" if st else "unassessed")
            applicable_roles = {"Type I":{"R","M","S"},"Type II":{"R","M"},"Type III":{"R","M"},"Type IV":{"R","A"}}.get(typ,{"R","M","S","A"})
            role_states = {r:("not_applicable" if r not in applicable_roles and not roles["Aux" if r == "A" else r] else row_assessment) for r in ("R","M","S","A")}
            within_pool_role_states = {r:("not_applicable" if role_states[r] == "not_applicable" else within_pool_role_status) for r in role_states}
            context_values = {}; orphan_values = defaultdict(list)
            for k in st:
                states = {r.get("context_state") for r in role_records[k] if (typ == "Total" or _type(r.get("rm_type") or r.get("family")) == typ) and r.get("context_state")}
                for gid,g in group_index.items():
                    if g["membership_kind"] == "saved_context" and k in group_current[gid] and _type(g.get("rm_type") or g.get("family")) in {t for t,_ in assignments[k]}:
                        state = current_context_state.get(gid)
                        if state in CONTEXT_MET:
                            states.add("met")
                        elif state in CONTEXT_MISSING:
                            states.add("missing")
                if states & CONTEXT_MET and states & CONTEXT_MISSING:
                    raise ValueError("Conflicting context assessments")
                if k in roles["R"]:
                    context_values[k] = "met" if states & CONTEXT_MET else "missing" if states & CONTEXT_MISSING else "unassessed"
                for r in role_records[k]:
                    if typ != "Total" and _type(r.get("rm_type") or r.get("family")) != typ:
                        continue
                    o = r.get("orphan")
                    if o is not None and type(o) is not bool:
                        raise ValueError("orphan must be explicit boolean or absent")
                    if o is not None:
                        for name in ("M","S"):
                            if k in roles[name]: orphan_values[name].append((k,o))
            if typ == "Total":
                context_values = {k:v for k,v in context_values.items() if type_bucket[k] in {"Type I","Type II","Type III"}}
            applicable = typ in {"Type I","Type II","Type III","Total"}
            context_denominator = len(context_values)
            assessed = [v for v in context_values.values() if v != "unassessed"] if applicable else []
            context_state = ("complete" if applicable and len(assessed) == context_denominator and row_assessment == "complete"
                             else "partial" if applicable and (assessed or not context_values and row_assessment == "partial") else "unassessed")
            count_missing = sum(v == "missing" for v in context_values.values()) if applicable and assessed else (0 if applicable and not context_values and row_assessment != "unassessed" else None)
            orphan_counts = {}; orphan_state = {}
            for name in ("M","S"):
                vals = defaultdict(set)
                for k,o in orphan_values[name]: vals[k].add(o)
                if any(len(v)>1 for v in vals.values()): raise ValueError("Conflicting orphan states")
                orphan_counts[name] = sum(True in v for v in vals.values()) if vals or not roles[name] and row_assessment != "unassessed" else None
                orphan_state[name] = ("complete" if len(vals) == len(roles[name]) and row_assessment == "complete"
                                      else "partial" if vals or not roles[name] and row_assessment == "partial" else "unassessed")
                if role_states[name] == "not_applicable":
                    orphan_counts[name] = None; orphan_state[name] = "not_applicable"
            z = aggregate(pool,[sid])
            type_tables.append(dict(sample_id=sid,rm_type=typ,is_total=typ=="Total",unique_CDS=len(st),member_keys=_keys(st),pool_CDS=len(pool),pool_member_keys=_keys(pool),type_disagreement_CDS=sum(type_disagreement[k] for k in pool),
                counts_status=counts_status,pool_status=pool_status,
                within_pool_counts_status=within_pool_counts_status,within_pool_role_states=within_pool_role_states,
                count_scope="observed assigned-type pool; broader candidate retention coverage remains separate",
                source_counts=z["source_counts"],source_count_status=z["source_count_status"],
                R=len(roles["R"]),M=len(roles["M"]),S=len(roles["S"]),A=len(roles["Aux"]),
                role_states=role_states,
                context_applicability="saved_I_III_context" if applicable else "not_assessed_by_I_III_rule",
                context_state=context_state,context_R_denominator=context_denominator if applicable else None,incomplete_R_context=count_missing,orphan_M=orphan_counts["M"],orphan_S=orphan_counts["S"],orphan_state=orphan_state,
                complete_contexts=sum(g["sample_id"] == sid and (typ == "Total" or g["family"] == typ) and g.get("context_state") in CONTEXT_MET for g in context_assessments) if any(g["sample_id"] == sid and g.get("context_state") is not None for g in context_assessments) else None,
                activity="unassessed"))
            for role,ks in roles.items():
                components.append(dict(sample_id=sid,rm_type=typ,role=role,count=len(ks),member_keys=_keys(ks),scope="saved retained or explicit reviewed annotation roles"))
    ledger = []
    for k in sorted(U):
        r = dict(sample_id=k[0],gene_id=k[1],protein_sha256=genes[k].get("protein_sha256"),
            native_member=k in native,retained=k in retained,reviewed_annotation_member=k in reviewed,
            final_RM=k in final_rm,restriction_role=k in R,support_qualified_restriction_role=k in support_qualified_R,
            priority_support_qualified=k in support_qualified_P,
            policy_R_component_support=saved.get(k, {}).get("integrated_R_component_support"),individual_effector=k in individual,
            complex_core=k in cores,effector=k in E,motif_site_supported=k in site,motif_grades=site_grades[k],
            associated=k in A,selected_group_member=k in G,ungrouped_effector=k in ungrouped,priority=k in P,
            source_membership={s:k in source_members[s] for s in sources},
            family_ids=[f for f in families if k in pools[f]],
            saved_decision=saved.get(k),annotations=annotation_by_gene[k],effector_evidence=evidence_by_gene[k],
            activity="unassessed",individual_catalysis="unassigned_complex_core" if k in cores-individual else "not_experimentally_assessed")
        if k in P:
            r["selection_route"] = "R_and_effector" if k in R&E else "R_only" if k in R else "effector_only"
            r["R_policy_support_state"] = ("supported" if k in policy_supported_R else
                "policy_nuclease_route_not_met" if k in unmet_R else "unassessed") if k in R else "not_applicable"
        ledger.append(r)
    source_out = [dict(id=s,label=s,kind="hmm_evidence" if not source_kind[s] else "mixed" if len(source_kind[s])>1 else "reference_annotation" if source_kind[s]=={"reference_hit"} else "native_final_member",coverage_status=_status(cov.get((sid,s),"missing") for sid in sample_ids)) for s in sources]
    global_summary = dict(integrated_CDS=len(U),native_union_CDS=len(native),retained_CDS=len(retained),reviewed_CDS=len(reviewed),
        DNMB_only_CDS=len(U-native),retained_only_CDS=len(retained-native),reviewed_only_CDS=len(reviewed-(native|retained)),
        final_RM_CDS=len(final_rm),restriction_role_CDS=len(R),individual_effector_CDS=len(individual),complex_core_CDS=len(cores),
        effector_CDS=len(E),motif_CDS=len(site),associated_CDS=len(A),system_member_CDS=len(G),ungrouped_effector_CDS=len(ungrouped),
        priority_CDS=len(P),R_and_effector_CDS=len(R&E),R_only_CDS=len(R-E),
        priority_support_qualified_CDS=len(support_qualified_P),policy_nuclease_route_unmet_R_CDS=len(unmet_R),
        R_only_policy_unmet_CDS=len(unmet_R-E),
        source_union_CDS=len(set.union(set(),*source_members.values())),
        metric_status=dict(inventory=inventory_state,associated=associated_state,effector=role_global,motif=motif_global))
    # Original R-M source membership is never backfilled from current DNMB roles.
    raw_rm_sources = defaultdict(set)
    native_rm_sources = defaultdict(set)
    for call in calls:
        if call["evidence_kind"] != "hmm_hit" and call.get("selected", True) and (call["evidence_kind"] == "reference_hit" or call.get("category") == "R-M" or _type(call.get("rm_type")) in {"Type I","Type II","Type III","Type IV"}):
            raw_rm_sources[call["caller"]].add(key(call,"original R-M source"))
            if call["evidence_kind"] == "final_member":
                native_rm_sources[call["caller"]].add(key(call,"native R-M source"))
    # Keep superseded/broader native candidates in the audit without silently
    # changing the explicitly supplied current R-M evaluation denominator.
    active_candidates = set(saved)
    for sid, scoped in current_scopes.items():
        active_candidates = {k for k in active_candidates if k[0] != sid} | scoped
    typed_candidates = {k for k in active_candidates
                        if _type(saved.get(k, {}).get("rm_type") or saved.get(k, {}).get("family"))
                        in {"Type I", "Type II", "Type III", "Type IV"}}
    raw_rm_source_members = set.union(set(), *raw_rm_sources.values())
    raw_rm_pool = raw_rm_source_members | typed_candidates
    generic_predictions = defaultdict(set)
    for call in calls:
        if call["evidence_kind"] == "final_member" and call.get("selected",True) and call.get("family") == "General defense":
            generic_predictions[call["caller"]].add(key(call,"general defense prediction"))
    for provider,members in generic_predictions.items():
        raw_rm_sources[provider] |= members & raw_rm_pool
        native_rm_sources[provider] |= members & raw_rm_pool
    # Keep the legacy all-configured-source endpoint, including references,
    # separate from the atlas's native-final endpoint. A mixed provider's
    # reference membership or coverage must not backfill its native channel.
    native_intersection_sources = [s for s in ("PADLOC", "DefenseFinder", "DefensePredictor")
                                   if "final_member" in source_kind[s]]
    native_intersection_sources += [s for s in sources if "final_member" in source_kind[s]
                                    and s not in native_intersection_sources]
    native_rm_coverage = {k:v for k,v in cov.items() if source_kind[k[1]] == {"final_member"}}
    for row in tables["coverage"]:
        provider = row["caller"]
        if row.get("endpoint") == "native_calls" or (not row.get("endpoint") and source_kind[provider] == {"final_member"}):
            native_rm_coverage[(row["sample_id"], provider)] = row["status"]
    for call in calls:
        if call["evidence_kind"] == "final_member" and call.get("selected", True):
            native_rm_coverage.setdefault((call["sample_id"], call["caller"]), "partial")
    all_rm_intersection = set.intersection(*(raw_rm_sources[s] for s in sources)) if sources else set()
    native_rm_intersection = (set.intersection(*(native_rm_sources[s] for s in native_intersection_sources))
                              if native_intersection_sources else set())
    current_rm_pool = category_pools.get("R-M",set())
    sample_summary = []
    for sid in sample_ids:
        st = {k for k in U if k[0] == sid}
        z = aggregate(st, [sid])
        retained_state = saved_status_by_sample[sid]
        retained_members = st & retained
        retained_rm_state = ("partial" if retained_state == "complete" and retained_members - final_rm else retained_state)
        retained_r_state = retained_rm_state
        if retained_r_state == "complete" and any(
                not any(role in {"R","M","S","RM","Aux","A","C","F"} for _,role in assignments[k])
                for k in retained_members & final_rm):
            retained_r_state = "partial"
        retention_assessment_members = current_scopes.get(sid, {k for k in saved if k[0] == sid}) | retained_members
        retention_assessed = retention_assessment_members & decision_assessed
        current_retained_state = ("complete" if retention_assessed == retention_assessment_members else
                                  "partial" if retention_assessed else "unassessed")
        if sid not in current_scopes:
            current_retained_state = saved_status_by_sample[sid]
        current_retained_rm_state = ("partial" if current_retained_state == "complete"
                                     and retained_members-final_rm else current_retained_state)
        current_retained_r_state = current_retained_rm_state
        if current_retained_r_state == "complete" and any(
                not any(role in {"R","M","S","RM","Aux","A","C","F"} for _,role in assignments[k])
                for k in retained_members & final_rm):
            current_retained_r_state = "partial"
        # These are original-source R-M endpoints, unlike the native-only
        # inventory/type pool. Reference coverage and saved candidate records
        # both participate, and assessment belongs to this sample alone.
        rm_source_states = [cov.get((sid, source), "missing") for source in sources]
        rm_pool_count = sum(k[0] == sid for k in raw_rm_pool)
        broad_pool_state = _status([*rm_source_states, saved_status_by_sample[sid]])
        if rm_pool_count and broad_pool_state == "unassessed":
            broad_pool_state = "partial"
        candidate_members = {k for k in active_candidates if k[0] == sid}
        # Retention can be unassessed while a bound Type-II candidate's pool
        # membership is known. An explicitly out-of-scope annotation supplies
        # no typed membership; this is not a biological-negative assertion.
        classified = {k for k in candidate_members if k in raw_rm_source_members or
                      (k in saved and (k in typed_candidates or saved[k]["state"] == "out_of_scope"))}
        membership_state = ("complete" if classified == candidate_members else
                            "partial" if classified else "unassessed")
        pool_sources = [s for s in sources if source_kind[s] & {"final_member", "reference_hit"}]
        membership_contribution = [membership_state] if candidate_members else []
        rm_pool_state = _status([*(cov.get((sid, source), "missing") for source in pool_sources), *membership_contribution])
        if rm_pool_count and rm_pool_state == "unassessed":
            rm_pool_state = "partial"
        rm_intersection_state = _intersection_status(rm_source_states)
        native_intersection_state = _intersection_status(
            native_rm_coverage.get((sid, source), "missing") for source in native_intersection_sources)
        sample_summary.append(dict(sample_id=sid, integrated_CDS=len(st), native_union_CDS=len(st&native),
            retained_CDS=len(st&retained), reviewed_CDS=len(st&reviewed), rm_CDS=len(st&final_rm),
            retained_counts_status=retained_state,retained_RM_counts_status=retained_rm_state,
            retained_R_counts_status=retained_r_state,retained_ratio_status=_status([retained_rm_state,retained_r_state]),
            current_retained_counts_status=current_retained_state,
            current_retained_RM_counts_status=current_retained_rm_state,
            current_retained_R_counts_status=current_retained_r_state,
            current_retained_ratio_status=_status([current_retained_state,current_retained_r_state]),
            current_retention_scope="declared current candidates plus explicitly retained active decisions outside that selection; not exhaustive discovery" if sid in current_scopes else "native saved candidate assessment",
            current_retention_assessment_member_keys=_keys(retention_assessment_members),
            current_retained_outside_candidate_member_keys=_keys(retained_members-current_scopes[sid]) if sid in current_scopes else [],
            associated_CDS=len(st&A), effector_CDS=len(st&E), motif_CDS=len(st&site), priority_CDS=len(st&P),
            source_counts={source:sum(k[0] == sid for k in source_members[source]) for source in sources},
            source_count_status=z["source_count_status"], metric_status=z["metric_status"],
            rm_source_counts={source:sum(k[0] == sid for k in raw_rm_sources[source]) for source in sources},
            rm_source_count_status=z["source_count_status"],rm_pool_CDS=rm_pool_count,
            rm_pool_status=rm_pool_state,rm_intersection_status=rm_intersection_state,
            rm_pool_broad_assessment_status=broad_pool_state,rm_candidate_membership_status=membership_state,
            rm_intersection_broad_assessment_status=(broad_pool_state if rm_source_states and
                all(state in {"complete", "partial"} for state in rm_source_states) else "unassessed"),
            rm_pool_count_scope="observed native/reference union and explicitly classified current candidates; not exhaustive discovery",
            rm_intersection_CDS=sum(k[0] == sid for k in all_rm_intersection),
            rm_intersection_sources=list(sources),
            rm_all_source_intersection_CDS=sum(k[0] == sid for k in all_rm_intersection),
            rm_all_source_intersection_status=rm_intersection_state,
            rm_native_intersection_CDS=sum(k[0] == sid for k in native_rm_intersection),
            rm_native_intersection_status=native_intersection_state,
            rm_native_intersection_sources=list(native_intersection_sources),
            retained_R_CDS=len(st & retained & R),retained_RM_CDS=len(st & retained & final_rm),
            noncanonical_CDS=len(st-current_rm_pool),
            noncanonical_source_counts={source:len(st & source_members[source] - current_rm_pool) for source in sources},
            noncanonical_source_count_status=z["source_count_status"],
            complete_saved_contexts=sum(g["sample_id"] == sid and g.get("context_state") in CONTEXT_MET for g in context_assessments) if any(g["sample_id"] == sid and g.get("context_state") is not None for g in context_assessments) else None,
            complete_context_state="observed_saved_assessment" if any(g["sample_id"] == sid and g["membership_kind"] == "saved_context" and g.get("context_state") is not None for g in group_index.values()) else "unassessed",
            rm_scope_definitions={"rm_pool":"original R-M source/reference union plus explicitly typed current DNMB candidates; untyped and broader historical nominees stay in audit", "intersection":"legacy all-configured-source intersection, including reference providers; rm_intersection_CDS and rm_all_source_intersection_CDS are aliases, with status from source coverage only", "native_intersection":"native final R-M members across rm_native_intersection_sources, including selected general predictions projected into the observed R-M pool; reference and HMM-only channels excluded; status from native endpoint coverage only", "general_prediction_RM":"saved general predictions intersect original R-M candidate pool; not a native RM type call", "noncanonical":"current integrated union minus current R-M display-category pool; source counts intersect this same pool"}))
    result = dict(schema_version="saved-results-v1",samples=samples,sources=source_out,sample_summary=sample_summary,
        selected_policy=deepcopy(evaluation.get("policy")) if evaluation else None,
        source_sets=[dict(source_id=s,endpoint=source_out[i]["kind"],member_keys=_keys(source_members[s]),count=len(source_members[s]),coverage_status=source_out[i]["coverage_status"]) for i,s in enumerate(sources)],
        global_summary=global_summary,overlap_regions=overlap(),priority_overlap_regions=overlap(P),
        inventory_rows=dict(categories=category_rows,families=family_rows),cells=cells,rm_components=components,type_tables=type_tables,
        gene_ledger=ledger,priority_ledger=[r for r in ledger if r["priority"]],groups_selected=selected_groups,complex_assignments=complex_out,
        context_assessments=context_assessments,
        rm_candidate_scope=dict(member_keys=_keys(active_candidates), historical_member_keys=_keys(saved),
            typed_member_keys=_keys(typed_candidates),untyped_member_keys=_keys(active_candidates-typed_candidates),
            unselected_historical_member_keys=_keys(set(saved)-active_candidates),
            selection="source_bound_evaluation_scope" if evaluation else "native_saved_candidates"),
        endpoint_definitions=dict(integrated="native final union OR saved retained OR explicitly included reviewed annotations",
            effector="individual putative nuclease-effector OR required complete reviewed restriction-complex core; not individual catalytic validation",
            associated="directly selected explicit group members within integrated union OR ungrouped effectors; no transitive expansion",
            priority="effector OR retained/reviewed restriction role; annotation-based triage, activity unassessed",
            priority_support_qualified="effector OR policy-supported R OR independently accepted reviewed R; narrower evidence tier, not experimental activity",
            motif="graded saved motif/reference-site support; not exhaustive or required for effector annotation"),
        warnings=warnings,availability=dict(coverage=[dict(sample_id=sid,caller=s,status=cov.get((sid,s),"missing")) for sid in sample_ids for s in sources],
            saved_decisions=saved_availability,saved_decisions_by_sample=saved_states,
            current_saved_decisions_by_sample=current_saved_states,
            effector=role_global,motif=motif_global,associated=associated_state,
            inventory_by_sample=inventory_by_sample, assessment_coverage_by_sample=assessment_by_sample,
            assessment_coverage_provenance_by_sample=assessment_provenance),
        workflow_capabilities=deepcopy((d.get("suite") or {}).get("capabilities",{})),
        provenance=dict(operation="stored_evidence_evaluation_and_aggregation" if evaluation else "saved_result_import_and_aggregation",detector_runs=0,
            scoring_runs=evaluation.get("scoring_runs", 0),
            component_registry_version=registry["version"],component_registry_sha256=registry_sha,
            automatic_component_assignments=automatic,source_metadata=tables["source_metadata"],input_sha256=d.get("input_sha256",{})))
    from .system_review import review_all_systems
    result["system_review"] = review_all_systems(d, result, rules=d.get("system_review_rules"))
    return result

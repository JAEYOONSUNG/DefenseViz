"""Generic integration of saved calls and DNMB evidence, never a detector run."""
from collections import defaultdict
from copy import deepcopy
from dataclasses import asdict
from .evidence import parse_evidence, parse_policy, evidence_dict, policy_metadata
from .policy import Evidence, Policy, Profile, evaluate

CATEGORIES = ("R-M", "Modification", "Nuclease", "CRISPR-Cas", "Abortive infection", "Other mechanisms", "Unresolved")
KINDS = {"final_member", "hmm_hit", "reference_hit"}
ROLES = {"R", "M", "S", "RM", "C", "Aux"}
RM_TYPES = {"Type I", "Type II", "Type III", "Type IV"}


def _string(value, context):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} must be a nonempty string")
    return value


def _records(value, context):
    if not isinstance(value, (list, tuple)) or any(not isinstance(r, dict) for r in value):
        raise ValueError(f"{context} must be a list of objects")
    return value


def integrate(samples, calls, evidence_records, *, coverage=None, genes=None, policy=None):
    """Integrate distinct identities using (sample_id, gene_id) tuple keys.

    Family cells count native final members. The integrated union additionally
    includes retained DNMB-only genes, explicitly reported separately. No pooled
    family/component counts are interpreted as complete biological systems.
    """
    p = parse_policy(policy)
    samples = deepcopy(_records(samples, "samples"))
    calls = deepcopy(_records(calls, "calls"))
    evidence_records = deepcopy(_records(evidence_records, "evidence_records"))
    sample_ids = []
    for i, s in enumerate(samples):
        sid = _string(s.get("sample_id"), f"samples[{i}].sample_id")
        _string(s.get("label"), f"samples[{i}].label")
        if sid in sample_ids:
            raise ValueError(f"Duplicate sample_id: {sid!r}")
        if "highlight" in s and type(s["highlight"]) is not bool:
            raise ValueError(f"Sample {sid!r}: highlight must be boolean")
        for k in ("genus", "species"):
            if k in s and s[k] is not None and not isinstance(s[k], str):
                raise ValueError(f"Sample {sid!r}: {k} must be a string")
        sample_ids.append(sid)
    valid_samples = set(sample_ids)
    callers = []
    warnings = ["Saved-call evaluation and DNMB evidence integration; no detector rerun.",
                "Historically calibrated evidence preset; no biological activity or independently validated performance is inferred.",
                "Family/category rows can overlap; row sums are not a CDS union. Native systems remain caller-specific."]
    def bound(row, where):
        sid = _string(row.get("sample_id"), where + ".sample_id")
        if sid not in valid_samples:
            raise ValueError(f"{where}: unknown sample_id {sid!r}")
        return sid
    def caller_name(name, where):
        _string(name, where)
        if name not in callers:
            callers.append(name)
        return name
    category_by_family = {}
    normalized = []
    exact_calls = set()
    system_family = {}
    for i, r in enumerate(calls):
        where = f"calls[{i}]"
        sid = bound(r, where)
        gid = _string(r.get("gene_id"), where + ".gene_id")
        caller = _string(r.get("caller"), where + ".caller")
        family = _string(r.get("family"), where + ".family")
        category = r.get("category") or "Unresolved"
        if category not in CATEGORIES:
            raise ValueError(f"{where}: unknown category {category!r}; map explicitly or use Unresolved")
        kind = r.get("evidence_kind", "final_member")
        if kind not in KINDS:
            raise ValueError(f"{where}: evidence_kind must be one of {sorted(KINDS)}")
        # Reference databases/providers are auxiliary evidence sources, not
        # independent native callers. An explicit coverage row can still opt a
        # provider into native caller scope without relying on its name.
        if kind != "reference_hit":
            caller_name(caller, where + ".caller")
        for k in ("system_id", "role", "rm_type", "source"):
            if r.get(k) is not None:
                _string(r[k], where + "." + k)
        if r.get("role") is not None and r["role"] not in ROLES:
            raise ValueError(f"{where}: unsupported role {r['role']!r}")
        r.update(category=category, evidence_kind=kind)
        for k in ("system_id", "role", "rm_type"):
            r.setdefault(k, None)
        if family in category_by_family and category_by_family[family] != category:
            raise ValueError(f"Family {family!r} has conflicting categories; use an explicit consistent map")
        category_by_family[family] = category
        identity = (sid, gid, caller, family, kind, r["system_id"], r["role"], r["rm_type"])
        if identity in exact_calls:
            raise ValueError(f"{where}: duplicate call identity {identity!r}")
        exact_calls.add(identity)
        if kind == "final_member" and r["system_id"] is not None:
            key = (sid, caller, r["system_id"])
            if key in system_family and system_family[key] != family:
                raise ValueError(f"{where}: native system_id assigned to conflicting families within sample/caller")
            system_family[key] = family
        normalized.append(r)
    cov = {}
    if coverage is None:
        cov = {(r["sample_id"], r["caller"]): "complete" for r in normalized if r["evidence_kind"] == "final_member"}
        warnings.append("Coverage was omitted: observed final-member sample/caller groups are treated as supplied complete native results; omitted groups and HMM/reference-only groups remain missing. Supply explicit coverage when this assumption is inappropriate.")
    else:
        for i, r in enumerate(_records(coverage, "coverage")):
            sid = bound(r, f"coverage[{i}]")
            caller = caller_name(r.get("caller"), f"coverage[{i}].caller")
            key = (sid, caller)
            if key in cov:
                raise ValueError(f"Duplicate coverage for sample/caller {key!r}")
            if r.get("status") not in {"complete", "missing", "not_run"}:
                raise ValueError(f"coverage[{i}]: invalid status")
            cov[key] = r["status"]
    for r in normalized:
        if r["evidence_kind"] == "final_member" and cov.get((r["sample_id"], r["caller"]), "missing") != "complete":
            raise ValueError(f"Call conflicts with incomplete coverage: {(r['sample_id'], r['caller'])!r}")
    coverage_out = [{"sample_id": sid, "caller": c, "status": cov.get((sid, c), "missing")} for sid in sample_ids for c in callers]
    # Validate geometry without inventing context: only a separately adjudicated
    # Evidence record may currently supply old_context/divergent_context.
    gene_records = deepcopy(_records(genes if genes is not None else [], "genes"))
    gene_keys = set()
    gene_sequence_hashes = {}
    for i, r in enumerate(gene_records):
        sid = bound(r, f"genes[{i}]")
        gid = _string(r.get("gene_id"), f"genes[{i}].gene_id")
        key = (sid, gid)
        if key in gene_keys:
            raise ValueError(f"Duplicate gene geometry: {key!r}")
        gene_keys.add(key)
        _string(r.get("contig"), f"genes[{i}].contig")
        if any(type(r.get(k)) is not int for k in ("start", "end")) or not 1 <= r["start"] <= r["end"]:
            raise ValueError(f"genes[{i}]: coordinates must be 1-based inclusive integers, start <= end")
        if r.get("strand") not in {"+", "-"}:
            raise ValueError(f"genes[{i}]: strand must be '+' or '-'")
        sha = r.get("sequence_sha256")
        if sha is not None and (not isinstance(sha, str) or len(sha) != 64 or any(c not in "0123456789abcdef" for c in sha)):
            raise ValueError(f"genes[{i}]: sequence_sha256 must be lowercase SHA-256")
        gene_sequence_hashes[key] = sha
    evidence_by_key = {}
    decisions_by_key = {}
    for i, r in enumerate(evidence_records):
        where = f"evidence_records[{i}]"
        sid = bound(r, where)
        gid = _string(r.get("gene_id"), where + ".gene_id")
        key = (sid, gid)
        if key in evidence_by_key:
            raise ValueError(f"Duplicate evidence for sample/gene {key!r}")
        if r.get("provenance") is not None and not isinstance(r["provenance"], dict):
            raise ValueError(f"{where}.provenance must be an object")
        evidence_sha = (r.get("provenance") or {}).get("sequence_sha256")
        if evidence_sha is not None:
            if not isinstance(evidence_sha, str) or len(evidence_sha) != 64 or any(c not in "0123456789abcdef" for c in evidence_sha):
                raise ValueError(f"{where}.provenance.sequence_sha256 must be lowercase SHA-256")
            geometry_sha = gene_sequence_hashes.get(key)
            if geometry_sha is not None and evidence_sha != geometry_sha:
                raise ValueError(f"{where}: sequence hash conflicts with gene geometry for {key!r}")
        # Incomplete records are explicitly unassessed rather than implicitly
        # assigning zero points to missing boolean assessments.
        if r.get("evidence") is None:
            evidence_by_key[key] = None
            decisions_by_key[key] = dict(sample_id=sid, gene_id=gid, status="unassessed", reason="policy_evidence_not_supplied", provenance=r.get("provenance", {}))
            continue
        try:
            e = parse_evidence(r["evidence"])
        except ValueError as exc:
            raise ValueError(f"{where}: {exc}") from exc
        evidence_by_key[key] = e
        if e.family not in RM_TYPES and e.family not in {"NA", "Unresolved", "Unknown"}:
            decisions_by_key[key] = dict(sample_id=sid, gene_id=gid, status="out_of_scope", evidence=evidence_dict(e), reason="DNMB preset is an R-M evidence policy")
            continue
        if e.role not in {"R", "M", "S", "RM"}:
            decisions_by_key[key] = dict(sample_id=sid, gene_id=gid, status="out_of_scope", evidence=evidence_dict(e), reason="unsupported R-M role")
            continue
        if not e.valid_sequence:
            decisions_by_key[key] = dict(sample_id=sid, gene_id=gid, status="unassessed", evidence=evidence_dict(e), reason="invalid_or_missing_sequence_cannot_supply_policy_support", provenance=r.get("provenance", {}))
            continue
        decision = evaluate(e, p)
        # Canonical JSON values without tuple/frozenset implementation details.
        decision = {k: list(v) if isinstance(v, tuple) else v for k, v in decision.items()}
        decisions_by_key[key] = dict(sample_id=sid, gene_id=gid, status="scored", evidence=evidence_dict(e), decision=decision, provenance=r.get("provenance", {}))
    native = [r for r in normalized if r["evidence_kind"] == "final_member"]
    evidence_calls = [r for r in normalized if r["evidence_kind"] != "final_member"]
    native_providers = {r["caller"] for r in native}
    evidence_sources = list(dict.fromkeys(r["caller"] for r in evidence_calls if r["caller"] not in native_providers))
    native_by_sample = defaultdict(set)
    caller_genes = defaultdict(set)
    family_genes = defaultdict(set)
    family_callers = defaultdict(set)
    family_systems = defaultdict(set)
    missing_system_ids = set()
    rm_pool = defaultdict(set)
    components = defaultdict(set)
    explicit_roles = defaultdict(set)
    families = []
    for r in native:
        sid, gid, caller, family = (r[k] for k in ("sample_id", "gene_id", "caller", "family"))
        native_by_sample[sid].add(gid)
        caller_genes[(sid, caller)].add(gid)
        family_genes[(sid, family)].add(gid)
        family_callers[(sid, family, caller)].add(gid)
        if family not in families:
            families.append(family)
        if r["system_id"] is not None:
            family_systems[(sid, family)].add((caller, r["system_id"]))
        else:
            missing_system_ids.add((sid, caller, family))
    for r in normalized:
        if r["evidence_kind"] != "hmm_hit" and r["category"] == "R-M":
            sid, gid = r["sample_id"], r["gene_id"]
            rm_pool[sid].add(gid)
            explicit_roles[(sid,gid)].add(r["role"])
            if r["role"] is not None:
                components[(sid, r["rm_type"] or "Unresolved", r["role"])].add(gid)
    for (sid, gid), e in evidence_by_key.items():
        if e is not None and e.family not in RM_TYPES and e.family not in {"NA", "Unresolved", "Unknown"}:
            continue
        rm_pool[sid].add(gid)
        if e is not None and e.role in ROLES:
            components[(sid, e.family if e.family in RM_TYPES else "Unresolved", e.role)].add(gid)
    for sid in sample_ids:
        for gid in sorted(rm_pool[sid]):
            roles = explicit_roles.get((sid,gid), set())
            control_only = bool(roles) and roles <= {"C", "Aux"}
            decisions_by_key.setdefault((sid, gid), dict(sample_id=sid, gene_id=gid,
                status="out_of_scope" if control_only else "unassessed",
                reason="explicit_control_or_auxiliary_role" if control_only else "policy_evidence_not_supplied"))
    cells = []
    for sid in sample_ids:
        for family in families:
            cells.append(dict(sample_id=sid, family=family, category=category_by_family[family],
                cds_count=len(family_genes[(sid, family)]),
                system_count=None if any((sid,c,family) in missing_system_ids for c in callers) or not any(cov.get((sid,c)) == "complete" for c in callers) else len(family_systems[(sid, family)]),
                known_native_system_count=len(family_systems[(sid, family)]),
                system_count_status="unassessed_missing_system_ids" if any((sid,c,family) in missing_system_ids for c in callers) else "observed" if any(cov.get((sid,c)) == "complete" for c in callers) else "unassessed_missing_coverage",
                count_status="observed" if any(cov.get((sid,c)) == "complete" for c in callers) else "unassessed",
                caller_counts={c: len(family_callers[(sid,family,c)]) if cov.get((sid,c)) == "complete" else None for c in callers}))
    summaries = []
    overlap = []
    for sid in sample_ids:
        scored = {gid for (s,gid), d in decisions_by_key.items() if s == sid and d["status"] == "scored"}
        retained = {gid for gid in scored if decisions_by_key[(sid,gid)]["decision"]["keep"]}
        native_set = native_by_sample[sid]
        summaries.append(dict(sample_id=sid, union_cds=len(native_set | retained), native_union_cds=len(native_set),
            retained_dnmb_only_cds=len(retained-native_set), rm_pool_cds=len(rm_pool[sid]), scored_cds=len(scored),
            retained_cds=len(retained), unassessed_cds=sum(decisions_by_key[(sid,g)]["status"] == "unassessed" for g in rm_pool[sid]),
            out_of_scope_cds=sum(decisions_by_key[(sid,g)]["status"] == "out_of_scope" for g in rm_pool[sid]),
            caller_counts={c: len(caller_genes[(sid,c)]) if cov.get((sid,c)) == "complete" else None for c in callers}))
        masks = defaultdict(list)
        for gid in sorted(native_set):
            masks[tuple(c for c in callers if gid in caller_genes[(sid,c)])].append(gid)
        for group, gids in masks.items():
            overlap.append(dict(sample_id=sid, callers=list(group), gene_ids=gids, count=len(gids),
                scope="observed_native_final_members", all_callers_complete=all(cov.get((sid,c)) == "complete" for c in callers)))
    rm_components = []
    for (sid, rm_type, role), gids in sorted(components.items(), key=lambda item: (sample_ids.index(item[0][0]), item[0][1], item[0][2])):
        supported = set()
        for gid in gids:
            d = decisions_by_key.get((sid,gid), {})
            e = evidence_by_key.get((sid,gid))
            if d.get("status") != "scored" or e is None or e.role != role or (e.family if e.family in RM_TYPES else "Unresolved") != rm_type:
                continue
            supported_by_policy = (d["decision"]["integrated_R_component_support"]
                                   if role in {"R", "RM"} else d["decision"]["keep"])
            if supported_by_policy:
                supported.add(gid)
        rm_components.append(dict(sample_id=sid, rm_type=rm_type, role=role, count=len(gids), supported_count=len(supported)))
    native_system_counts = []
    for sid in sample_ids:
        for caller in callers:
            for family in families:
                keys = {system_id for c,system_id in family_systems[(sid,family)] if c == caller}
                if cov.get((sid,caller)) == "complete":
                    missing = (sid,caller,family) in missing_system_ids
                    native_system_counts.append(dict(sample_id=sid, caller=caller, family=family,
                        count=None if missing else len(keys), known_count=len(keys),
                        count_status="unassessed_missing_system_ids" if missing else "observed"))
    return dict(samples=samples, callers=callers, families=[dict(family=f, category=category_by_family[f]) for f in families],
        cells=cells, cell_scope="native_final_members", sample_summary=summaries,
        decisions=[d for sid in sample_ids for (s,g),d in decisions_by_key.items() if s==sid], overlap=overlap,
        rm_components=rm_components, native_system_counts=native_system_counts, policy=asdict(p),
        policy_metadata=policy_metadata(p), warnings=warnings, coverage=coverage_out, evidence_calls=evidence_calls,
        evidence_sources=evidence_sources,
        native_calls=native, genes=gene_records, system_count_unit="sum_of_explicit_caller_specific_native_systems",
        rm_component_scope="explicit_role_candidates_in_RM_pool; not complete systems", activity="not_assessed")


__all__ = ["integrate", "Evidence", "Profile", "Policy", "evaluate", "parse_evidence", "parse_policy"]

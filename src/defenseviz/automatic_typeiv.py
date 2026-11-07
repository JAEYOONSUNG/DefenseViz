"""Source-bound Type IV component review, separate from strict DNMB scoring.

The packaged rules describe component architectures. They do not establish
activity, substrate specificity, catalytic residues, or a complete complex.
Only existing current R-M candidates are reviewed; no sequence search runs here.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
from importlib.resources import files
from itertools import product
import json
import re

ORIGIN = "automatic_typeiv_component_review"
RULE_VERSION = "source_bound_typeiv_components_v1"
SUPPORTED = {"supported", "accepted", "supported_annotation", "supported_putative_effector_annotation"}
HOLDS = ("alternative_annotation", "conflict", "partial", "hmm_role_ambiguous")


def _key(row):
    return row.get("sample_id"), row.get("gene_id")


def _refs(row):
    values = [row.get("source_ref"), *(row.get("source_refs") or [])]
    result = []
    for ref in values:
        if (isinstance(ref, dict) and ref.get("path")
                and re.fullmatch(r"[a-f0-9]{64}", str(ref.get("sha256", ""))) and ref not in result):
            result.append(deepcopy(ref))
    return result


def _all_refs(rows):
    result = []
    for row in rows:
        for ref in _refs(row):
            if ref not in result:
                result.append(ref)
    return result


def _bool(value):
    if type(value) is bool:
        return value
    return {"TRUE": True, "FALSE": False, "1": True, "0": False}.get(str(value).upper())


def _integer(value):
    if isinstance(value, bool) or not re.fullmatch(r"[0-9]+", str(value)):
        return None
    return int(value)


def _role(value):
    return {"R_component": "R", "A": "Aux"}.get(value, value)


def _type(value):
    return {"I": "Type I", "II": "Type II", "III": "Type III", "IV": "Type IV"}.get(value, value)


def _pinned(row, inventory):
    refs = _refs(row)
    return bool(refs) and all((row.get("sample_id"), ref["path"], ref["sha256"]) in inventory for ref in refs)


def _domains(gene, observations, inventory):
    """Use original accepted InterPro rows, never product/eggNOG descriptions."""
    accepted, errors = [], []
    for row in observations:
        if row.get("namespace") != "InterProScan" or row.get("evidence_type") != "domain_hit":
            continue
        native = (row.get("native_fields", {}) or {})
        start, end = _integer(native.get("start")), _integer(native.get("end"))
        length = _integer(native.get("protein_length"))
        valid = (row.get("sequence_binding") == "verified" and row.get("coordinate_binding") == "verified"
                 and row.get("protein_sha256") == gene.get("protein_sha256")
                 and native.get("protein_id") == gene.get("gene_id")
                 and str(native.get("sequence_md5", "")).lower() == gene.get("protein_md5")
                 and length == gene.get("protein_length") and length is not None
                 and start is not None and end is not None and 1 <= start <= end <= length
                 and _pinned(row, inventory))
        if not valid:
            errors.append("contradictory_raw_domain_binding")
            continue
        if native.get("status") != "T" or native.get("analysis") not in {"Pfam", "PANTHER"}:
            continue
        accepted.append(dict(analysis=native["analysis"], accession=native.get("signature_accession"),
            start=start, end=end, sequence_md5=native["sequence_md5"], protein_length=length,
            protein_sha256=gene["protein_sha256"], source_refs=_refs(row)))
    return accepted, sorted(set(errors))


def _architectures(domains, rules):
    selected, families = [], set()
    for rule in rules:
        options = [[d for d in domains if d["analysis"] == "Pfam" and d["accession"] == acc]
                   for acc in rule["all_pfams"]]
        combinations = [list(combo) for combo in product(*options)
                        if not rule["ordered"] or all(a["end"] < b["start"] for a, b in zip(combo, combo[1:]))]
        panther = [d for d in domains if d["analysis"] == "PANTHER"
                   and str(d["accession"]).split(":")[0] == rule.get("panther_family")]
        if not combinations or (rule.get("panther_family") and not panther) or rule["family"] in families:
            continue
        # Specific original architectures precede broader same-family routes.
        families.add(rule["family"])
        selected.append(dict(rule_id=rule["id"], component_family=rule["family"], role=_role(rule["role"]),
            domain_combinations=combinations, family_model_evidence=panther,
            requires_full_reference=rule.get("requires_full_reference", False),
            requires_native_typeiv=rule.get("requires_native_TypeIV", False),
            requires_native_domain_span=rule.get("requires_native_domain_span", False),
            required_partner=rule.get("required_partner"), source_urls=rule.get("source_urls", [])))
    return selected


def _native_typeiv(call):
    # Reuse the importer's exact named-model taxonomy, including versioned
    # aliases; no substring classifier of model names or protein descriptions.
    from .io.suite import RM, DF_RM_ALIASES, DF_NAMED_TYPE_IV_MODELS
    native = (call.get("native_fields", {}) or {})
    if call.get("caller") == "PADLOC":
        return RM.get(native.get("system")) == "Type IV"
    family = native.get("subtype") or call.get("native_family") or call.get("family")
    fqn = native.get("model_fqn")
    lineage = fqn.split("/") if isinstance(fqn, str) else []
    if lineage and family != lineage[-1]:
        return False
    named = DF_NAMED_TYPE_IV_MODELS.get(family)
    if named:
        return (not fqn or fqn == named[0]) and native.get("type") in {None, "RM", family, "Type IV", "RM_Type_IV", "RM_type_IV"}
    types = {value for value in (RM.get(family), DF_RM_ALIASES.get(family),
             RM.get(lineage[2]) if len(lineage) >= 3 and lineage[:2] == ["defense-finder-models", "RM"] else None,
             RM.get(native.get("type"))) if value}
    return types == {"Type IV"}


def _native_memberships(gene, calls, maps, inventory, lengths):
    result, problems = [], []
    for call in calls:
        caller = call.get("caller")
        if (caller not in {"DefenseFinder", "PADLOC"} or call.get("evidence_kind") != "final_member"
                or call.get("selected") is not True or call.get("rm_type") != "Type IV"):
            continue
        native = (call.get("native_fields", {}) or {})
        identity = native.get("hit_id" if caller == "DefenseFinder" else "target.name")
        matching = [m for m in maps if m.get("namespace") == caller.lower() and m.get("source_id") == identity]
        exact = (identity == gene["gene_id"] or
                 (matching and all(_key(m) == _key(gene) and _pinned(m, inventory) for m in matching)))
        # Native typing is explicit, not inferred from a generic HNH/Mrr label.
        native_type = _native_typeiv(call)
        valid = bool(exact and native_type and call.get("protein_sha256") == gene.get("protein_sha256") and _pinned(call, inventory))
        span = None
        if caller == "DefenseFinder":
            a, b, length = (_integer(native.get(x)) for x in ("hit_begin_match", "hit_end_match", "hit_seq_len"))
            valid = valid and length in lengths and a is not None and b is not None and 1 <= a <= b <= length
            if valid:
                span = [a, b]
        if not valid:
            problems.append("contradictory_native_typeiv_binding")
            continue
        result.append(dict(caller=caller, system_id=call.get("system_id"), native_id=identity,
                           protein_span=span, source_refs=_all_refs([call, *matching])))
    return result, sorted(set(problems))


def _normalization_binding(record, features, gene):
    """Respect the full evaluator's verified terminal-stop identity receipt.

    InterPro remains bound to the original source protein. Native HMM queries
    may omit its terminal stop, as already verified against the canonical
    translation, FASTA, and signed-by-hash stage manifest by full_evaluation.
    """
    binding = record.get("protein_binding")
    lengths = {gene.get("protein_length")}
    if not binding:
        return True, lengths
    sha = gene.get("protein_sha256")
    if binding.get("source_sha256") != sha:
        return False, lengths
    fields = {"source_sha256": "source_protein_sha256", "canonical_sha256": "canonical_protein_sha256", "hmm_sha256": "hmm_protein_sha256"}
    if any(not re.fullmatch(r"[a-f0-9]{64}", str(binding.get(k, ""))) for k in fields):
        return False, lengths
    if all(binding[k] == sha for k in fields):
        return True, lengths
    if (binding.get("normalization") != "terminal_stop_removal_only"
            or binding.get("source_receipt") not in {"exact_canonical_translation", "uppercase whitespace-free; terminal stop removed for HMM only"}):
        return False, lengths
    pins = {(r["path"], r["sha256"]) for r in _refs(record)}
    matching = [f for f in features if f.get("protein_sha256") == sha
                and any((r["path"], r["sha256"]) in pins for r in _refs(f))
                and all((f.get("native_fields") or {}).get(v) == binding[k] for k, v in fields.items())]
    if len(matching) != 1:
        return False, lengths
    clean_length = _integer(matching[0]["native_fields"].get("sequence_length_aa"))
    if clean_length is None or not 0 < clean_length <= gene["protein_length"]:
        return False, lengths
    lengths.add(clean_length)
    return True, lengths


def _reference(record, features):
    """Require the selected, source-bound curated Type IV R reference route."""
    evidence = (record.get("evidence", {}) or {})
    pins = {(x["path"], x["sha256"]) for x in _refs(record)}
    selected = [f for f in features if f.get("protein_sha256") == record.get("protein_sha256")
                and any((x["path"], x["sha256"]) in pins for x in _refs(f))]
    qualified = []
    for row in selected:
        n = (row.get("native_fields", {}) or {})
        if (n.get("curated_blast_family") == "Type IV" and n.get("curated_blast_role") == "R"
                and _bool(n.get("curated_blast_role_compatible")) is True
                and _bool(n.get("curated_blast_family_compatible")) is True):
            qualified.append(row)
    gate = bool((record.get("decision", {}) or {}).get("numeric_reference_pass") is True
                and evidence.get("curated_role_compatible") is True
                and evidence.get("curated_family_compatible") is True and len(qualified) == 1)
    return dict(qualified=gate, source_refs=_all_refs(qualified),
                subject=(qualified[0]["native_fields"].get("curated_blast_resolved_subject_id")
                         or qualified[0]["native_fields"].get("curated_blast_match")) if gate else None)


def _assignment_conflicts(gene, role, annotations, decisions):
    reasons = []
    for row in annotations:
        relevant = row.get("assignment_kind") == "reviewed_RM_component"
        if row.get("status") in SUPPORTED:
            if row.get("assignment_kind") == "noncatalytic_component" and role == "R":
                reasons.append("explicit_noncatalytic_review_conflict")
            if relevant:
                typ, old_role = _type(row.get("rm_type") or row.get("family")), _role(row.get("role"))
                if typ not in {None, "Unresolved", "Type IV"} or old_role not in {None, "unknown", "Unknown", role}:
                    reasons.append("explicit_review_assignment_conflict")
                if row.get("protein_sha256") != gene.get("protein_sha256"):
                    reasons.append("explicit_review_protein_binding_conflict")
        elif relevant and row.get("status") in {"rejected", "unsupported", "not_supported", "conflicting"}:
            reasons.append("explicit_review_veto")
    for row in decisions:
        if row.get("state") != "retained":
            continue
        typ, old_role = _type(row.get("rm_type") or row.get("family")), _role(row.get("role"))
        if typ not in {None, "Unresolved", "Type IV"} or old_role not in {None, "unknown", "Unknown", role}:
            reasons.append("retained_assignment_conflict")
        if row.get("protein_sha256") != gene.get("protein_sha256"):
            reasons.append("retained_protein_binding_conflict")
    return reasons


def apply_typeiv_component_review(data):
    """Return a copy with reviewed component annotations and an exact audit.

    This function does not modify saved decisions, scores, source calls, context
    groups, effector evidence, or motif evidence. Reapplication replaces only
    annotations it previously produced, and is deterministic and idempotent.
    """
    out = deepcopy(data)
    protocol_bytes = files("defenseviz").joinpath("data/typeiv_component_rules.json").read_bytes()
    protocol = json.loads(protocol_bytes)
    digest = hashlib.sha256(protocol_bytes).hexdigest()
    protocol_ref = dict(path="defenseviz:data/typeiv_component_rules.json", sha256=digest)
    out["annotations"] = [r for r in out.get("annotations", []) if r.get("origin") != ORIGIN]
    suite = (out.get("suite", {}) or {})
    inventory = {(r.get("sample_id"), r.get("path"), r.get("sha256")) for r in suite.get("sources", [])}
    indexes = {}
    for name, rows in (("genes", out.get("genes", [])), ("calls", out.get("calls", [])),
                       ("annotations", out["annotations"]), ("decisions", out.get("saved_decisions", [])),
                       ("profiles", suite.get("profiles", [])), ("features", suite.get("features", [])),
                       ("quarantine", suite.get("unmapped_records", []))):
        indexes[name] = defaultdict(list)
        for row in rows:
            indexes[name][_key(row)].append(row)
    maps = defaultdict(list)
    for row in suite.get("identity_maps", []):
        maps[row.get("sample_id")].append(row)
    candidates = defaultdict(list)
    for row in (out.get("rm_evaluation", {}) or {}).get("records", []):
        if _type((row.get("evidence") or {}).get("family", row.get("family"))) == "Type IV":
            candidates[_key(row)].append(row)
    records = []
    for key, current in sorted(candidates.items()):
        record = current[0]
        genes = indexes["genes"][key]
        gene = genes[0] if genes else {}
        evidence = (record.get("evidence", {}) or {})
        reasons = []
        sha, md5, length = gene.get("protein_sha256"), gene.get("protein_md5"), gene.get("protein_length")
        bound = (len(genes) == 1 and re.fullmatch(r"[a-f0-9]{64}", str(sha or ""))
                 and re.fullmatch(r"[a-f0-9]{32}", str(md5 or ""))
                 and type(length) is int and length > 0 and record.get("protein_sha256") == sha
                 and evidence.get("valid_sequence") is True and bool(_refs(record)) and bool(_refs(gene))
                 and record.get("assignment_status") in {"supported", "evaluated"})
        normalization_valid, native_lengths = _normalization_binding(record, indexes["features"][key], gene) if bound else (False, set())
        bound = bound and normalization_valid
        if record.get("protein_binding_status") not in {None, "verified"}:
            bound = False
        if not bound:
            reasons.append("current_protein_binding_unresolved")
        if len(current) != 1:
            reasons.append("ambiguous_current_typeiv_candidate")
        holds = {name: evidence.get(name) for name in HOLDS}
        if any(value is True for value in holds.values()):
            reasons.append("preserved_source_hold")
        if any(value is None for value in holds.values()):
            reasons.append("source_hold_assessment_incomplete")
        domains, errors = _domains(gene, indexes["profiles"][key], inventory) if bound else ([], [])
        reasons += errors
        if any(r.get("module") == "interproscan" or r.get("namespace") == "InterProScan" for r in indexes["quarantine"][key]):
            reasons.append("quarantined_raw_domain_evidence")
        rules = _architectures(domains, protocol["rules"])
        native, errors = _native_memberships(gene, indexes["calls"][key], maps[key[0]], inventory, native_lengths) if bound else ([], [])
        reasons += errors
        reference = _reference(record, indexes["features"][key])
        if not (reference["qualified"] or native):
            reasons.append("typeiv_corroboration_unresolved")
        if not rules:
            reasons.append("component_architecture_not_covered_or_incomplete")
        if len(rules) > 1:
            reasons.append("multiple_role_architectures_require_review")
        for match in rules:
            if match["requires_full_reference"] and not reference["qualified"]:
                reasons.append("architecture_route_needs_qualified_full_reference")
            if match["requires_native_typeiv"] and not native:
                reasons.append("hnh_route_needs_exact_native_typeiv_membership")
            if match["requires_native_domain_span"]:
                spans = [n["protein_span"] for n in native if n["caller"] == "DefenseFinder" and n["protein_span"]]
                if not any(a <= d["start"] <= d["end"] <= b for a, b in spans
                           for combo in match["domain_combinations"] for d in combo):
                    reasons.append("native_family_alignment_does_not_span_catalytic_domain")
        role = rules[0]["role"] if len(rules) == 1 else None
        if role:
            reasons += _assignment_conflicts(gene, role, indexes["annotations"][key], indexes["decisions"][key])
        reasons = sorted(set(reasons))
        supported = not reasons
        quarantine = any("conflict" in reason or "contradictory" in reason or "ambiguous" in reason or "quarantined" in reason for reason in reasons)
        sources = _all_refs([record, gene, *domains, *native, reference]) + [protocol_ref]
        result = dict(sample_id=key[0], gene_id=key[1], protein_sha256=record.get("protein_sha256"),
            status="supported_annotation" if supported else "quarantined" if quarantine else "unassessed",
            annotation_supported=supported, role=role if supported else None,
            rule_id=rules[0]["rule_id"] if supported else None, observed_architectures=rules,
            accepted_domain_count=len(domains), raw_domain_observation_count=sum(r.get("namespace") == "InterProScan" and r.get("evidence_type") == "domain_hit" for r in indexes["profiles"][key]),
            reference_gate=reference["qualified"], native_gate=bool(native), native_memberships=native,
            source_holds=holds, reasons=reasons, source_refs=sources,
            strict_policy_state=record.get("state"), strict_policy_score=(record.get("decision", {}) or {}).get("score_total"),
            activity="not_inferred", complete_complex="not_inferred", substrate_specificity="unresolved",
            catalytic_residue_validation="not_performed", motif_credit_added=False, effector_credit_added=False)
        records.append(result)
        if supported:
            out["annotations"].append(dict(sample_id=key[0], gene_id=key[1], protein_sha256=sha,
                origin=ORIGIN, rule_id=RULE_VERSION + ":" + result["rule_id"], protocol_sha256=digest,
                assignment_kind="reviewed_RM_component", family="Type IV", rm_type="Type IV", role=role,
                component_family=rules[0]["component_family"], status="supported_annotation", include_in_union=True,
                source_refs=sources, evidence_kind="source_bound_component_architecture",
                activity="not_inferred", complete_complex="not_inferred", substrate_specificity="unresolved"))
    out["automatic_typeiv_review"] = dict(schema_version=1, rule_version=RULE_VERSION, protocol_sha256=digest,
        protocol_source_ref=protocol_ref, scope=protocol["scope"], records=records,
        status_counts=dict(Counter(r["status"] for r in records)),
        role_counts=dict(Counter(r["role"] for r in records if r["annotation_supported"])),
        analysis_kind="fresh_evaluation_of_stored_exact_domain_and_current_rm_evidence",
        engine_runs=0, policy_changed=False, thresholds_changed=False, motif_credit_added=False,
        effector_credit_added=False, complete_complex_inferred=False)
    return out

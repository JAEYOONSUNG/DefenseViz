"""Lossless, additive review of saved defense detections and exact evidence.

This module does not change accepted annotations or E/A/P membership. A completed
review can legitimately have unassessed scientific axes. New rules describe a
bounded evidence contract, not a detector rerun or experimental activity.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
from importlib.resources import files
import json
from pathlib import Path


AXES = ("native_detection", "family_annotation", "component_role", "nuclease_domain",
        "motif_site", "system_architecture", "complex_core", "activity")
STATES = {"supported", "incomplete", "unassessed", "conflicting", "not_applicable"}
POSITIVE = {"supported", "accepted", "supported_annotation", "supported_putative_effector_annotation"}
ROUTE_ENDPOINTS = {
    "native_component": {"component_role", "family_annotation"},
    "qualified_domain": {"family_annotation", "nuclease_domain"},
    "qualified_site": {"motif_site"},
    "native_architecture": {"system_architecture", "complex_core"},
}
MATCH_FIELDS = {"provider", "exact_native_models", "exact_native_families", "evidence_kinds", "selected",
                "database_versions", "model_versions", "software_versions", "model_sha256", "namespace", "signature_accessions",
                "qualification_rules", "family_ids"}
REQUIRE_FIELDS = {"protein_binding", "source_binding", "same_native_unit", "all_components",
                  "qualified_domains", "qualified_sites"}
CONCLUSION_FIELDS = {"component", "role", "substrate", "family_id", "individual_catalysis", "excludes_roles"}


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _evidence_identity(value):
    """Relocating the same source file does not create new scientific evidence."""
    if isinstance(value, dict):
        return {k: _evidence_identity(v) for k, v in value.items() if k not in {"path", "resolved_path"}}
    if isinstance(value, list):
        return [_evidence_identity(v) for v in value]
    return value


def _id(prefix, value):
    return prefix + ":" + hashlib.sha256(_canonical(value).encode()).hexdigest()


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be nonempty text")
    return value


def _sha(value, name):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"{name} must be a lowercase SHA-256")
    return value


def _object(value, allowed, required, name):
    if not isinstance(value, dict) or set(value) - allowed or not required <= set(value):
        raise ValueError(f"{name}: unknown or missing fields; allowed {sorted(allowed)}, required {sorted(required)}")


def _strings(value, name):
    if not isinstance(value, list) or not value or any(not isinstance(v, str) or not v.strip() for v in value) or len(set(value)) != len(value):
        raise ValueError(f"{name} must be distinct nonempty strings")


def _values(value):
    return [value] if isinstance(value, str) else value


def _ref_schema(ref, name):
    if not isinstance(ref, dict) or not {"path", "sha256"} <= set(ref):
        raise ValueError(f"{name} requires path and sha256")
    _text(ref["path"], name + ".path")
    _sha(ref["sha256"], name + ".sha256")


def validate_system_review_rules(rules):
    """Validate finite, exact rules. No expression language or implicit wildcard."""
    _object(rules, {"schema_version", "registry_id", "version", "rules"},
            {"schema_version", "registry_id", "version", "rules"}, "system review registry")
    if type(rules["schema_version"]) is not int or rules["schema_version"] != 1:
        raise ValueError("System review rules require schema_version 1")
    _text(rules["registry_id"], "registry_id"); _text(rules["version"], "registry version")
    if not isinstance(rules["rules"], list):
        raise ValueError("System review rules must be a list")
    bundled = json.loads(files("defenseviz.data").joinpath("components.json").read_text())
    reserved = {r["rule_id"] for name in ("component_rules", "complex_rules") for r in bundled[name]}
    seen = set()
    required = {"rule_id", "version", "route", "endpoint", "primary_sources", "match", "requires", "conclusion"}
    for rule in rules["rules"]:
        _object(rule, required | {"architecture_contract", "defense_scope"}, required, "system review rule")
        rid = _text(rule["rule_id"], "rule_id")
        if rid in seen or rid in reserved:
            raise ValueError("Duplicate or reserved system review rule_id")
        seen.add(rid)
        _text(rule["version"], "rule version"); _strings(rule["primary_sources"], "primary_sources")
        route = rule["route"]
        if route not in ROUTE_ENDPOINTS or rule["endpoint"] not in ROUTE_ENDPOINTS[route]:
            raise ValueError("Unsupported system review route/endpoint combination")
        if rule.get("defense_scope", "unresolved") not in {"host_defense", "anti_defense", "unresolved"}:
            raise ValueError("Invalid rule defense_scope")
        match, req, conclusion = rule["match"], rule["requires"], rule["conclusion"]
        _object(match, MATCH_FIELDS, set(), "rule.match")
        _object(req, REQUIRE_FIELDS, {"protein_binding", "source_binding"}, "rule.requires")
        _object(conclusion, CONCLUSION_FIELDS, set(), "rule.conclusion")
        for field, value in match.items():
            if field in {"provider", "namespace"}:
                _text(value, "match." + field)
            elif field == "selected":
                if type(value) is not bool:
                    raise ValueError("match.selected must be boolean")
            else:
                values = _values(value) if field == "model_sha256" else value
                _strings(values, "match." + field)
                if field == "model_sha256":
                    for digest in values:
                        _sha(digest, "match.model_sha256")
        if set(match.get("evidence_kinds", [])) - {"final_member", "hmm_hit", "reference_hit"}:
            raise ValueError("Unsupported evidence kind in system review rule")
        if req["protein_binding"] is not True or req["source_binding"] is not True:
            raise ValueError("New review rules require protein_binding and source_binding")
        if "same_native_unit" in req and type(req["same_native_unit"]) is not bool:
            raise ValueError("same_native_unit must be boolean")
        for field, value in conclusion.items():
            if field == "excludes_roles":
                _strings(value, "conclusion.excludes_roles")
            else:
                _text(value, "conclusion." + field)
        if conclusion.get("substrate", "unknown") not in {"DNA", "RNA", "dual", "unknown", "unresolved"}:
            raise ValueError("Unsupported review substrate")
        if rule["endpoint"] == "component_role" and not conclusion.get("role"):
            raise ValueError("Component-role conclusion requires an explicit role")
        if rule["endpoint"] == "family_annotation" and not conclusion.get("family_id"):
            raise ValueError("Family annotation conclusion requires an explicit family_id")
        if route.startswith("native_"):
            if not {"provider", "evidence_kinds", "selected"} <= set(match):
                raise ValueError("Native rule requires exact provider, evidence_kinds and selected state")
            if match["selected"] is not True or match["evidence_kinds"] != ["final_member"]:
                raise ValueError("Native positive rules require selected final_member evidence")
            if not (match.get("database_versions") or match.get("model_versions") or match.get("model_sha256")):
                raise ValueError("Native rule requires exact database_versions, model_versions or model_sha256")
        if route == "native_component" and not match.get("exact_native_models"):
            raise ValueError("Native component rule requires exact_native_models")
        if route in {"qualified_domain", "qualified_site"} and not match.get("family_ids"):
            raise ValueError("Qualified rules require exact registry family_ids")
        if route in {"qualified_domain", "qualified_site"} and conclusion.get("family_id") and conclusion["family_id"] not in match["family_ids"]:
            raise ValueError("Qualified conclusion family_id must be one of its exact matched family_ids")
        route_match = ({"provider", "exact_native_models", "exact_native_families", "evidence_kinds", "selected", "database_versions", "model_versions", "software_versions", "model_sha256"}
                       if route.startswith("native_") else {"family_ids", "namespace", "signature_accessions", "qualification_rules", "model_sha256", "database_versions"})
        if set(match) - route_match:
            raise ValueError("Match predicate is not implemented for this rule route")
        if route.startswith("native_") and ("qualified_domains" in req or "qualified_sites" in req):
            raise ValueError("Qualified domain/site predicates require a qualified rule route")
        if route == "qualified_domain" and "qualified_sites" in req:
            raise ValueError("Site requirements require the qualified_site route")
        if route == "native_architecture":
            if not match.get("exact_native_families") or req.get("same_native_unit") is not True:
                raise ValueError("Architecture requires exact_native_families and same_native_unit")
            contract = rule.get("architecture_contract")
            _object(contract, {"contract_id", "version", "source_ref"}, {"contract_id", "version", "source_ref"}, "architecture_contract")
            _text(contract["contract_id"], "architecture contract_id"); _text(contract["version"], "architecture version")
            _ref_schema(contract["source_ref"], "architecture source_ref")
            components = req.get("all_components")
            if not isinstance(components, list) or not components:
                raise ValueError("Architecture requires nonempty all_components")
            component_ids = set()
            for component in components:
                fields = {"component", "component_rule_ids", "minimum_members", "maximum_members"}
                _object(component, fields, fields, "architecture component")
                name = _text(component["component"], "component")
                if name in component_ids:
                    raise ValueError("Duplicate required component")
                component_ids.add(name)
                _strings(component["component_rule_ids"], "component_rule_ids")
                lo, hi = component["minimum_members"], component["maximum_members"]
                if type(lo) is not int or type(hi) is not int or not 1 <= lo <= hi:
                    raise ValueError("Component multiplicities require 1 <= minimum <= maximum")
            if rule["endpoint"] == "complex_core" and conclusion.get("individual_catalysis") != "unassigned":
                raise ValueError("Complex cores cannot assign individual catalysis")
        elif "architecture_contract" in rule or "all_components" in req or "same_native_unit" in req:
            raise ValueError("Architecture-only fields used outside native_architecture")
        # Domain/site predicates consume exact qualified registry evidence only.
        for field, allowed, required_fields in (
            ("qualified_domains", {"namespace", "signature_accession", "model_sha256", "database_version", "qualification_rule", "minimum_model_coverage"}, {"namespace", "signature_accession", "qualification_rule"}),
            ("qualified_sites", {"reference_sha256", "required_positions"}, {"reference_sha256", "required_positions"}),
        ):
            if field not in req:
                continue
            if not isinstance(req[field], list) or not req[field]:
                raise ValueError(field + " must be a nonempty list")
            for predicate in req[field]:
                _object(predicate, allowed, required_fields, field)
                for key, value in predicate.items():
                    if key.endswith("sha256"):
                        _sha(value, field + "." + key)
                    elif key == "minimum_model_coverage":
                        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
                            raise ValueError("Invalid minimum_model_coverage")
                    elif key == "required_positions":
                        if not isinstance(value, list) or not value or any(type(p) is not int or p < 1 for p in value) or len(set(value)) != len(value):
                            raise ValueError("required_positions must be distinct positive integers")
                    else:
                        _text(value, field + "." + key)
                if field == "qualified_domains" and not (predicate.get("model_sha256") or predicate.get("database_version")):
                    raise ValueError("Qualified domain requires model digest or database release")
    index = {r["rule_id"]: r for r in rules["rules"]}
    for rule in rules["rules"]:
        for component in rule["requires"].get("all_components", []):
            for rid in component["component_rule_ids"]:
                other = index.get(rid)
                if other is None or other["route"] != "native_component" or other["conclusion"].get("component") != component["component"]:
                    raise ValueError("Architecture references an unknown or incompatible component rule")
    return None


def _refs(record):
    refs = list(record.get("source_refs") or [])
    if record.get("source_ref") and record["source_ref"] not in refs:
        refs.append(record["source_ref"])
    return refs


def _axis(reason="no_applicable_rule", coverage="unassessed"):
    return dict(status="unassessed", reason_codes=[reason], supporting_evidence_ids=[],
                conflicting_evidence_ids=[], evaluated_rule_ids=[], coverage=coverage)


def _axes():
    return {name: _axis("activity_not_assessed" if name == "activity" else "no_applicable_rule") for name in AXES}


def _add_axis(axes, name, state, reason, evidence=None, rule_id=None):
    state = "conflicting" if state == "conflict" else state
    if state not in STATES:
        state = "unassessed"
    item = axes[name]
    if item["reason_codes"] == ["no_applicable_rule"]:
        item["reason_codes"] = []
    if reason and reason not in item["reason_codes"]:
        item["reason_codes"].append(reason)
    ranking = {"not_applicable": 0, "unassessed": 1, "incomplete": 2, "supported": 3, "conflicting": 4}
    if ranking[state] > ranking[item["status"]] or (state == "not_applicable" and item["status"] == "unassessed" and not evidence):
        item["status"] = state
    if evidence:
        key = "conflicting_evidence_ids" if state == "conflicting" else "supporting_evidence_ids"
        if evidence not in item[key]:
            item[key].append(evidence)
    if rule_id and rule_id not in item["evaluated_rule_ids"]:
        item["evaluated_rule_ids"].append(rule_id)


def _safe_state(values):
    values = list(values)
    if values and all(v == "complete" for v in values):
        return "complete"
    if any(v in {"complete", "partial"} for v in values):
        return "partial"
    return "unassessed"


class _Bindings:
    def __init__(self, data):
        self.inputs = data.get("input_sha256", {})
        self.verified = {}
        self.rules_base = Path(data.get("system_review_rules_base", "."))
        self.aliases = defaultdict(set)
        self.alias_hashes = defaultdict(set)
        for source in (data.get("suite") or {}).get("sources", []):
            if source.get("path") and source.get("sha256") and source.get("resolved_path"):
                self.aliases[(source["path"], source["sha256"])].add(source["resolved_path"])
                self.alias_hashes[source["path"]].add(source["sha256"])

    def refs(self, refs, *, require_imported=False):
        if not refs:
            return False
        all_available = True
        for ref in refs:
            _ref_schema(ref, "review source reference")
            path, expected = ref["path"], ref["sha256"]
            aliases = self.aliases.get((path, expected), set())
            if len(aliases) > 1:
                raise ValueError("Ambiguous exact source path/digest alias")
            if path in self.alias_hashes and not aliases:
                raise ValueError("Source alias hash differs from imported source")
            if aliases:
                path = next(iter(aliases))
            if require_imported and path not in self.inputs:
                all_available = False
                continue
            if path in self.inputs and self.inputs[path] != expected:
                raise ValueError("Review source hash conflicts with imported source")
            if path in self.verified:
                if self.verified[path] != expected:
                    raise ValueError("Conflicting review source digests")
                continue
            p = Path(path)
            if not p.is_file():
                if path in self.inputs:
                    raise ValueError("Pinned imported review source file unavailable")
                all_available = False
                continue
            actual = hashlib.sha256(p.read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError("System review source hash mismatch")
            self.verified[path] = actual
        return all_available


def _version(call, axis):
    source = call.get("source_version") or {}
    if not isinstance(source, dict):
        source = {}
    if axis == "database":
        return call.get("database_version") or source.get("db_version")
    if axis == "model":
        origin = call.get("model_version_source", "")
        fallback = origin.endswith("database_version") or origin.endswith("db_version")
        return (None if fallback else call.get("model_version")) or source.get("models_version")
    return call.get("software_version") or source.get("software_version")


def _defense_scope(call):
    """Read exact declared namespaces, never Anti_* or family-name substrings."""
    native = call.get("native_fields") or {}
    namespace = call.get("native_namespace", native.get("namespace"))
    fqn = native.get("model_fqn", call.get("model_fqn"))
    path = fqn.split("/") if isinstance(fqn, str) else []
    if call.get("caller") == "DefenseFinder":
        if namespace == "ADF" or len(path) > 1 and path[:2] == ["defense-finder-models", "ADF"]:
            return "anti_defense", "exact_native_ADF_namespace"
        if namespace in {"DefenseFinder", "CasFinder"} or len(path) > 1 and path[:2] in (["defense-finder-models", "DefenseFinder"], ["defense-finder-models", "RM"]) or path and path[0] == "CasFinder":
            return "host_defense", "exact_native_host_namespace"
    return "unresolved", "native_namespace_unresolved"


def _protein_binding(call, gene):
    current = gene.get("protein_sha256")
    native = call.get("native_fields") or {}
    producer = native.get("protein_sha256", native.get("sequence_sha256"))
    if producer is not None and current is not None and producer != current:
        raise ValueError("Producer protein hash differs from current gene")
    if producer and current:
        origin = "producer_reported"
    elif call.get("protein_binding_origin") == "importer_attached":
        origin = "importer_attached"
    elif call.get("protein_sha256"):
        origin = "canonical_supplied"
    else:
        origin = "canonical_gene_identity" if current else "unavailable"
    return dict(origin=origin, independent_producer_binding=origin == "producer_reported",
                current_protein_sha256=current, source_refs=deepcopy(_refs(call)))


def _native_scope(rule, call):
    m = rule["match"]
    return (call.get("caller") == m["provider"]
            and (not m.get("exact_native_models") or call.get("native_model") in m["exact_native_models"])
            and (not m.get("exact_native_families") or call.get("native_family", call.get("family")) in m["exact_native_families"]))


def _native_requirements(rule, call, gene, bindings):
    reasons = []
    m = rule["match"]
    if call.get("evidence_kind", "final_member") not in m["evidence_kinds"]:
        reasons.append("auxiliary_evidence_not_final_membership")
    if call.get("selected", True) != m["selected"]:
        reasons.append("unselected_detection_preserved")
    current = gene.get("protein_sha256")
    if not current:
        reasons.append("protein_binding_unavailable")
    else:
        _sha(current, "current protein_sha256")
        if call.get("protein_sha256") and call["protein_sha256"] != current:
            raise ValueError("System review protein hash mismatch")
    if not bindings.refs(_refs(call), require_imported=True):
        reasons.append("source_binding_unavailable")
    # Explicitly supplied constraints are all enforced. Either an exact database
    # release or model hash is mandatory; both supplied means both must agree.
    if m.get("database_versions"):
        v = _version(call, "database")
        if v is None:
            reasons.append("producer_version_unavailable")
        elif v not in m["database_versions"]:
            reasons.append("producer_version_mismatch")
    if m.get("model_versions"):
        v = _version(call, "model")
        if v is None:
            reasons.append("producer_version_unavailable")
        elif v not in m["model_versions"]:
            reasons.append("producer_version_mismatch")
    if m.get("software_versions"):
        v = _version(call, "software")
        if v is None:
            reasons.append("software_version_unavailable")
        elif v not in m["software_versions"]:
            reasons.append("software_version_mismatch")
    if m.get("model_sha256"):
        h = call.get("model_sha256")
        if h is None:
            reasons.append("model_digest_unavailable")
        else:
            _sha(h, "native call model_sha256")
            if h not in _values(m["model_sha256"]):
                reasons.append("model_digest_mismatch")
    scope, _ = _defense_scope(call)
    requested = rule.get("defense_scope", "unresolved")
    if scope == "anti_defense" and requested != "anti_defense":
        reasons.append("defense_scope_mismatch")
    elif scope == "host_defense" and requested == "anti_defense":
        reasons.append("defense_scope_mismatch")
    return reasons


def _verified_registry(data, genes, bindings):
    """Validate wrapper linkage and return exact candidate/entry/assessment rows."""
    pack = data.get("system_review_evidence")
    if pack is None:
        return []
    if not isinstance(pack, dict) or pack.get("schema_version") != 1 or pack.get("source_files_verified") is not True:
        raise ValueError("System review evidence must come from verified registry evaluation")
    _sha(pack.get("registry_sha256"), "system review registry_sha256")
    registry, candidates, assessments = pack.get("registry"), pack.get("candidates"), pack.get("assessments")
    if not isinstance(registry, dict) or not isinstance(candidates, list) or not isinstance(assessments, list):
        raise ValueError("Malformed verified system review evidence")
    hashes = pack.get("input_sha256")
    if not isinstance(hashes, dict) or not hashes or pack["registry_sha256"] not in hashes.values():
        raise ValueError("Verified registry lacks actual input hash bindings")
    if not bindings.refs([{"path": p, "sha256": h} for p, h in hashes.items()]):
        raise ValueError("Verified registry source file unavailable")
    # The adapter verifies full registry/raw-candidate files. Bind their content
    # here too, so replacing in-memory rows after evaluation cannot retain trust.
    if pack.get("registry_path"):
        if hashes.get(pack["registry_path"]) != pack["registry_sha256"]:
            raise ValueError("Registry path is not bound to exact registry hash")
        registry_files = [pack["registry_path"]]
    else:
        registry_files = [p for p, h in hashes.items() if h == pack["registry_sha256"]]
    if not any(json.loads(Path(p).read_text()) == registry for p in registry_files):
        raise ValueError("Verified registry content differs from pinned file")
    candidate_binding = pack.get("candidates_sha256")
    if pack.get("candidates_path"):
        if not candidate_binding or hashes.get(pack["candidates_path"]) != candidate_binding:
            raise ValueError("Candidates path is not bound to exact candidate hash")
        candidate_files = [pack["candidates_path"]]
    else:
        candidate_files = [p for p, h in hashes.items() if candidate_binding is None or h == candidate_binding]
    bound_candidates = False
    for path in candidate_files:
        try:
            if json.loads(Path(path).read_text()) == candidates:
                bound_candidates = True
                break
        except (UnicodeDecodeError, json.JSONDecodeError):
            pass
    if not bound_candidates:
        raise ValueError("Verified candidates differ from pinned source")
    by_candidate = {}
    for candidate in candidates:
        key = (candidate.get("sample_id"), candidate.get("query_id"))
        if key in by_candidate or key not in genes:
            raise ValueError("Duplicate or unknown verified review candidate")
        sha = candidate.get("sequence_sha256")
        if sha != genes[key].get("protein_sha256") or hashlib.sha256(candidate.get("sequence", "").encode()).hexdigest() != sha:
            raise ValueError("Verified review candidate protein mismatch")
        by_candidate[key] = candidate
    by_assessment = {}
    for assessment in assessments:
        key = (assessment.get("sample_id"), assessment.get("query_id"))
        if key in by_assessment or key not in by_candidate or assessment.get("source_files_verified") is not True:
            raise ValueError("Duplicate, unknown or unverified registry assessment")
        if assessment.get("sequence_sha256") != by_candidate[key]["sequence_sha256"]:
            raise ValueError("Registry assessment query mismatch")
        if any(hashes.get(p) != h for p, h in assessment.get("source_sha256", {}).items()):
            raise ValueError("Registry assessment sources differ from wrapper")
        by_assessment[key] = assessment
    if set(by_assessment) != set(by_candidate):
        raise ValueError("Registry assessment coverage differs from candidate list")
    from .registry import assess_domains, assess_registry
    entries = {entry["family_id"]: entry for entry in registry["families"]}
    if len(entries) != len(registry["families"]):
        raise ValueError("Duplicate registry family")
    out = []
    # Recompute all assertions with the public verifier. A serialized boolean or
    # injected supported assessment is never a trust boundary.
    registry_dir = Path(pack.get("registry_dir") or Path(registry_files[0]).parent)
    source_root = Path(pack.get("source_root") or registry_dir)
    recomputed = assess_registry(registry, candidates, registry_dir=registry_dir, source_root=source_root, verify_sources=True)
    if assessments != recomputed:
        raise ValueError("Registry assessments differ from exact supplied evidence evaluation")
    for assessment in recomputed:
        if any(hashes.get(p) != h for p, h in assessment.get("source_sha256", {}).items()):
            raise ValueError("Reverified registry source missing from wrapper hash bindings")
    for key, candidate in by_candidate.items():
        assessment = by_assessment[key]
        for family in assessment.get("families", []):
            entry = entries.get(family.get("family_id"))
            if entry is None:
                raise ValueError("Registry assessment family not registered")
            _, usable, _ = assess_domains(candidate, entry)
            out.append(dict(key=key, candidate=candidate, entry=entry, family=family, usable=usable,
                            source_sha256=assessment.get("source_sha256", {}), registry_sha256=pack["registry_sha256"], registry_dir=registry_dir))
    return out


def _qualified_match(rule, row):
    m, req = rule["match"], rule["requires"]
    if row["entry"]["family_id"] not in m["family_ids"]:
        return None
    if rule["conclusion"].get("family_id") and rule["conclusion"]["family_id"] != row["entry"]["family_id"]:
        return "unassessed", ["qualified_conclusion_family_mismatch"]
    usable = row["usable"]
    def hit_matches(u, predicate):
        hit, model = u["hit"], u["model"]
        namespace = predicate.get("namespace")
        accession = predicate.get("signature_accession")
        qrule = predicate.get("qualification_rule")
        coverage = ((hit["hmm_end"] - hit["hmm_start"] + 1) / hit["hmm_length"]
                    if hit.get("hmm_length") and hit.get("hmm_start") and hit.get("hmm_end") else hit.get("profile_coverage", 0))
        return ((namespace is None or model.get("source") == namespace)
                and (accession is None or model.get("accession", model.get("profile")) == accession)
                and (not predicate.get("model_sha256") or model.get("model_sha256") == predicate["model_sha256"])
                and (not predicate.get("database_version") or hit.get("database_version") == predicate["database_version"])
                and (qrule is None or hit.get("qualification", {}).get("rule") == qrule)
                and coverage >= predicate.get("minimum_model_coverage", 0))
    filtered = [u for u in usable if (not m.get("namespace") or u["model"].get("source") == m["namespace"])
                and (not m.get("signature_accessions") or u["model"].get("accession", u["model"].get("profile")) in m["signature_accessions"])
                and (not m.get("model_sha256") or u["model"].get("model_sha256") in _values(m["model_sha256"]))
                and (not m.get("qualification_rules") or u["hit"].get("qualification", {}).get("rule") in m["qualification_rules"])
                and (not m.get("database_versions") or u["hit"].get("database_version") in m["database_versions"])]
    reasons = []
    # A site rule carrying a domain-selector constraint cannot fall back to an
    # unrelated independently supported reference-site route.
    if not filtered and set(m) & {"namespace", "signature_accessions", "qualification_rules", "model_sha256", "database_versions"}:
        reasons.append("qualified_domain_required")
    if any(not any(hit_matches(u, p) for u in filtered) for p in req.get("qualified_domains", [])):
        reasons.append("qualified_domain_required")
    if req.get("qualified_domains"):
        filtered = [u for u in filtered if any(hit_matches(u, p) for p in req["qualified_domains"])]
    if rule["route"] == "qualified_domain":
        if not filtered:
            reasons.append("qualified_domain_required")
        if rule["endpoint"] == "nuclease_domain" and not any(u["model"].get("evidence_kind") == "nuclease_domain" and u["model"].get("role") == "nuclease" for u in filtered):
            reasons.append("qualified_nuclease_domain_required")
    from .registry import assess_domains, assess_reference
    narrowed = {**row["candidate"], "domain_hits": [u["hit"] for u in filtered]}
    family_assessment, narrowed_usable, _ = assess_domains(narrowed, row["entry"])
    domain_constrained = bool(set(m) & {"namespace", "signature_accessions", "qualification_rules", "model_sha256", "database_versions"} or req.get("qualified_domains"))
    references = []
    for reference in row["entry"].get("references", []):
        references.append(assess_reference(row["candidate"], row["entry"], reference, filtered, row["registry_dir"]))
        if not domain_constrained:
            references.append(assess_reference(row["candidate"], row["entry"], reference, filtered, row["registry_dir"], reference_route=True))
    for p in req.get("qualified_sites", []):
        registered = {r["accession"] for r in row["entry"].get("references", []) if r["sequence_sha256"] == p["reference_sha256"]}
        if not any(r.get("state") == "supported" and r.get("reference") in registered and set(p["required_positions"]) <= {s["reference_position"] for s in r.get("sites", []) if s.get("conserved")} for r in references):
            reasons.append("qualified_site_required")
    field = {"family_annotation": "family_evidence", "nuclease_domain": "domain_evidence", "motif_site": "catalytic_site_evidence"}[rule["endpoint"]]
    assessed = row["family"].get(field, {})
    if rule["endpoint"] == "family_annotation":
        assessed = family_assessment
    if rule["endpoint"] == "nuclease_domain":
        domain_supported = any(u["model"].get("evidence_kind") == "nuclease_domain" and u["model"].get("role") == "nuclease" for u in narrowed_usable)
        assessed = {"state": "supported" if domain_supported else "unassessed",
                    "reason": "qualified_narrowed_nuclease_domain" if domain_supported else "selected_model_does_not_support_nuclease_domain"}
    if rule["endpoint"] == "motif_site":
        status = ("supported" if any(r.get("state") == "supported" for r in references) else
                  "conflict" if any(r.get("state") == "conflict" for r in references) else
                  "incomplete" if any(r.get("state") == "incomplete" for r in references) else "unassessed")
        assessed = {"state": status, "reason": "qualified_exact_reference_site_contract" if status == "supported" else "qualified_site_contract_unmet"}
    status = "unassessed" if reasons else assessed.get("state", "unassessed")
    return status, reasons or [assessed.get("reason", "qualified_evidence_unassessed")]


def review_all_systems(data, result, rules=None):
    """Return every observed unit's review without modifying either input."""
    if rules is not None:
        validate_system_review_rules(rules)
    rule_list = rules["rules"] if rules else []
    samples = {s["sample_id"] for s in data.get("samples", [])}
    genes = {(g["sample_id"], g["gene_id"]): g for g in data.get("genes", [])}
    ledger = {(g["sample_id"], g["gene_id"]): g for g in result.get("gene_ledger", [])}
    bindings = _Bindings(data)
    if rules is not None and data.get("system_review_rule_source"):
        source = data["system_review_rule_source"]
        if not bindings.refs([source]) or json.loads(Path(source["path"]).read_text()) != rules:
            raise ValueError("System review rules differ from pinned rule source")
    units, call_index, detection_units, accounting, gene_units = {}, {}, {}, [], defaultdict(set)
    evidence, evaluations = [], []
    evidence_index = {}
    def append_evidence(record):
        eid = record["evidence_id"]
        if eid in evidence_index:
            if _evidence_identity(evidence_index[eid]) != _evidence_identity(record):
                raise ValueError("Conflicting records share an evidence identity")
            return
        evidence_index[eid] = record
        evidence.append(record)
    raw_rows = 0
    def ensure_unit(identity, kind, sid, provider=None, run=None, replicon=None, system=None):
        uid = _id("unit", identity)
        if uid not in units:
            units[uid] = dict(unit_id=uid, unit_kind=kind, sample_id=sid, provider=provider,
                run_id=run, replicon_id=replicon, native_system_id=system,
                native_families=[], native_models=[], native_defense_scopes=[], defense_scope="unresolved", member_keys=[], selected_member_keys=[],
                detection_ids=[], source_refs=[], context_group_keys=[],
                review_execution={"status": "completed"}, axes=_axes(), reason_codes=[],
                native_identity_binding="exact_canonical" if run is not None and system is not None and replicon is not None else "unassessed",
                native_identity_binding_scope="canonical_source_scoping_not_independent_producer_validation",
                run_binding_origins=[], replicon_binding_origins=[])
        return units[uid]
    for call in data.get("calls", []):
        key = (call["sample_id"], call["gene_id"])
        if key[0] not in samples:
            raise ValueError("Unknown sample in system review call")
        genes.setdefault(key, {"sample_id": key[0], "gene_id": key[1]})
        kind = call.get("evidence_kind", "final_member")
        provider, run, replicon, system = call.get("caller"), call.get("run_id"), call.get("replicon_id"), call.get("system_id")
        # Explicit gene contig can bind an absent provider replicon; it does not
        # merge a missing identity with another contig or invent coordinates.
        if replicon is None:
            replicon = genes[key].get("contig")
        if kind == "final_member" and system is not None:
            unit_kind = "native_system"
            run_scope = run if run is not None else ["unbound_run", sorted({ref.get("sha256") for ref in _refs(call) if ref.get("sha256")})]
            identity = [unit_kind, key[0], provider, run_scope, replicon, system]
        else:
            unit_kind = "ungrouped_detection" if kind == "final_member" else "auxiliary_detection"
            identity = [unit_kind, key[0], provider, run, replicon, key[1], kind]
        unit = ensure_unit(identity, unit_kind, key[0], provider, run, replicon, system)
        for destination, origin in (("run_binding_origins", call.get("run_id_kind", "canonical_supplied" if run is not None else "unavailable")),
                                    ("replicon_binding_origins", "native_call" if call.get("replicon_id") is not None else "gene_contig" if replicon is not None else "unavailable")):
            if origin not in unit[destination]:
                unit[destination].append(origin)
        detection_identity = [key, provider, run, replicon, system, kind, call.get("native_family", call.get("family")), call.get("family"), call.get("native_model"), call.get("reference_id"), call.get("selected", True)]
        did = _id("detection", detection_identity)
        if did in call_index:
            raise ValueError("Duplicate canonical detection in system review")
        call_index[did] = call
        detection_units[did] = unit["unit_id"]
        unit["detection_ids"].append(did)
        scope, scope_reason = _defense_scope(call)
        if scope not in unit["native_defense_scopes"]:
            unit["native_defense_scopes"].append(scope)
        known_scopes = set(unit["native_defense_scopes"]) - {"unresolved"}
        unit["defense_scope"] = next(iter(known_scopes)) if len(known_scopes) == 1 else "unresolved"
        if list(key) not in unit["member_keys"]:
            unit["member_keys"].append(list(key))
        if kind == "final_member" and call.get("selected", True) and list(key) not in unit["selected_member_keys"]:
            unit["selected_member_keys"].append(list(key))
        for dest, value in (("native_families", call.get("native_family", call.get("family"))), ("native_models", call.get("native_model"))):
            if value is not None and value not in unit[dest]:
                unit[dest].append(value)
        for ref in _refs(call):
            if ref not in unit["source_refs"]:
                unit["source_refs"].append(deepcopy(ref))
        count = len(call.get("native_observations") or [call.get("native_fields")])
        raw_rows += count
        accounting.append(dict(detection_id=did, unit_id=unit["unit_id"], sample_id=key[0], gene_id=key[1],
            provider=provider, evidence_kind=kind, selected=call.get("selected", True), raw_observation_count=count,
            defense_scope=scope, defense_scope_reason=scope_reason, protein_binding=_protein_binding(call, genes[key])))
        gene_units[key].add(unit["unit_id"])
        _add_axis(unit["axes"], "native_detection", "supported", "observed_saved_detection", did)
        if not call.get("selected", True):
            unit["reason_codes"].append("unselected_detection_preserved")
        if kind != "final_member":
            unit["reason_codes"].append("auxiliary_evidence_not_final_membership")
    stored_reviews = defaultdict(list)
    for table in ("saved_decisions", "annotations", "effector_evidence"):
        for record in data.get(table, []):
            key = (record["sample_id"], record["gene_id"])
            if key[0] not in samples:
                raise ValueError("Unknown sample in saved review record")
            genes.setdefault(key, {"sample_id": key[0], "gene_id": key[1], "protein_sha256":record.get("protein_sha256")})
            stored_reviews[key].append(dict(table=table, record=deepcopy(record)))
    target_keys = set(gene_units) | set(ledger) | set(stored_reviews)
    for key in sorted(target_keys - set(gene_units)):
        unit = ensure_unit(["review_only_gene", *key], "review_only_gene", key[0], "DNMB")
        unit["member_keys"] = [list(key)]
        gene_units[key].add(unit["unit_id"])
    # Context links are informative, never the unit discovery mechanism. Native
    # contexts require exact identity and member inclusion; no wildcard matches.
    for group in data.get("groups", []):
        members = {(group["sample_id"], m["gene_id"]) for m in group.get("members", [])}
        for uid in set.union(set(), *(gene_units.get(k, set()) for k in members)):
            unit = units[uid]
            if group.get("membership_kind") == "native_system":
                if (group.get("provider"), group.get("run_id"), group.get("native_system_id", group.get("system_id"))) != (unit["provider"], unit["run_id"], unit["native_system_id"]):
                    continue
                if not members <= {tuple(k) for k in unit["member_keys"]}:
                    continue
            unit["context_group_keys"].append([group["sample_id"], group["group_id"]])
    gene_reviews = {}
    roles, excluded_roles = defaultdict(set), defaultdict(set)
    role_sources, exclusion_sources = defaultdict(lambda: defaultdict(set)), defaultdict(lambda: defaultdict(set))
    for key in sorted(target_keys):
        existing = ledger.get(key, {})
        review = dict(sample_id=key[0], gene_id=key[1], protein_sha256=genes.get(key, existing).get("protein_sha256"),
            unit_ids=sorted(gene_units[key]), in_integrated_universe=key in ledger, axes=_axes(),
            review_execution={"status": "completed"}, preserved_assignment_ids=[], reason_codes=[],
            accepted_membership={name: bool(existing.get(name)) for name in ("effector", "associated", "priority", "motif_site_supported", "restriction_role", "individual_effector", "complex_core")})
        native_scopes = sorted({scope for uid in gene_units[key] for scope in units[uid]["native_defense_scopes"]})
        known_scopes = set(native_scopes) - {"unresolved"}
        review["native_defense_scopes"] = native_scopes
        review["defense_scope"] = next(iter(known_scopes)) if len(known_scopes) == 1 else "unresolved"
        gene_reviews[key] = review
        if stored_reviews.get(key):
            review["stored_review_records"] = stored_reviews[key]
        if any(units[u]["detection_ids"] for u in gene_units[key]):
            _add_axis(review["axes"], "native_detection", "supported", "observed_saved_detection")
        accepted = [("annotations", x) for x in existing.get("annotations", []) if x.get("status") in POSITIVE]
        accepted += [("effector_evidence", x) for x in existing.get("effector_evidence", []) if x.get("status") in POSITIVE]
        if isinstance(existing.get("saved_decision"), dict) and existing["saved_decision"].get("state") == "retained":
            accepted.append(("saved_decision", existing["saved_decision"]))
        for kind, record in accepted:
            eid = _id("existing", [key, kind, _evidence_identity(record)])
            if eid not in review["preserved_assignment_ids"]:
                review["preserved_assignment_ids"].append(eid)
            append_evidence(dict(evidence_id=eid, sample_id=key[0], gene_id=key[1], kind="preserved_accepted_assignment", record=deepcopy(record)))
            assignment = record.get("assignment_kind")
            if assignment == "noncatalytic_component":
                roles[key].add("noncatalytic")
                excluded_roles[key].add("nuclease")
                exclusion_sources[key]["nuclease"].add(eid)
            elif assignment == "individual_nuclease_effector":
                roles[key].add("nuclease")
                role_sources[key]["nuclease"].add(eid)
            elif record.get("role"):
                roles[key].add(record["role"])
                role_sources[key][record["role"]].add(eid)
            role_bearing = bool(record.get("role")) or assignment in {"individual_nuclease_effector", "nuclease_other_role", "noncatalytic_component"}
            if role_bearing:
                _add_axis(review["axes"], "component_role", "supported", "existing_supported_assignment_preserved", eid, record.get("rule_id"))
            if record.get("family"):
                _add_axis(review["axes"], "family_annotation", "supported", "existing_supported_assignment_preserved", eid, record.get("rule_id"))
            if record.get("motif_state") in POSITIVE:
                _add_axis(review["axes"], "motif_site", "supported", "existing_graded_motif_assignment_preserved", eid, record.get("rule_id"))
        # Flags preserve accepted endpoints even when legacy fixtures do not
        # provide all evidence fields. They cannot establish new stronger axes.
        for flag, axis in (("individual_effector", "component_role"), ("restriction_role", "component_role"), ("complex_core", "complex_core"), ("motif_site_supported", "motif_site")):
            if existing.get(flag):
                eid = _id("accepted_endpoint", [key, flag])
                append_evidence(dict(evidence_id=eid, sample_id=key[0], gene_id=key[1], kind="preserved_accepted_endpoint", endpoint=flag, source="analysis.gene_ledger", source_refs=deepcopy(_refs(existing.get("saved_decision") or {}))))
                _add_axis(review["axes"], axis, "supported", "existing_supported_assignment_preserved", eid)
        if existing.get("individual_effector"):
            roles[key].add("nuclease")
            role_sources[key]["nuclease"].add(_id("accepted_endpoint", [key, "individual_effector"]))
    # Index raw context by exact gene identity without copying large descriptions
    # or treating transferred text as qualified observed domains.
    for table in ("profiles", "annotations", "features", "motifs", "structures"):
        for row in (data.get("suite") or {}).get(table, []):
            key = (row.get("sample_id"), row.get("gene_id"))
            if key in gene_reviews:
                native = row.get("native_fields") or {}
                fields = ("protein_id", "sequence_md5", "protein_length", "analysis", "signature_accession", "signature_description",
                          "start", "end", "score", "status", "interpro_accession", "interpro_description", "go_annotations", "pathway_annotations",
                          "residue", "site_start", "site_end", "site_description", "seed_ortholog", "evalue", "eggNOG_OGs", "max_annot_lvl",
                          "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "PFAMs")
                summary = {field: native[field] for field in fields if field in native}
                gene_reviews[key].setdefault("context_evidence", []).append(dict(table=table, namespace=row.get("namespace"),
                    evidence_type=row.get("evidence_type"), source_refs=deepcopy(_refs(row)), native_summary=summary,
                    protein_sha256=row.get("protein_sha256"), sequence_binding=row.get("sequence_binding", "unassessed"),
                    coordinate_binding=row.get("coordinate_binding", "unassessed"), functional_inference="context_only_not_qualified_rule_evidence"))
    qualified = _verified_registry(data, genes, bindings)
    qualified_by_gene = defaultdict(list)
    for row in qualified:
        key = row["key"]
        if key not in gene_reviews:
            # Verified all-proteome evidence is retained as input but cannot
            # introduce a defense detection or accepted gene by itself.
            continue
        qualified_by_gene[key].append(row)
        eid = _id("registry", [key, row["registry_sha256"], row["entry"]["family_id"]])
        append_evidence(dict(evidence_id=eid, sample_id=key[0], gene_id=key[1], kind="verified_registry_assessment", family_id=row["entry"]["family_id"], registry_sha256=row["registry_sha256"], source_sha256=row["source_sha256"], assessment=deepcopy(row["family"])))
        for field, axis in (("family_evidence", "family_annotation"), ("domain_evidence", "nuclease_domain"), ("catalytic_site_evidence", "motif_site"), ("reference_site_evidence", "motif_site")):
            assessment = row["family"].get(field, {})
            _add_axis(gene_reviews[key]["axes"], axis, assessment.get("state", "unassessed"), assessment.get("reason", "verified_registry_unassessed"), eid)
    component_support = defaultdict(lambda: defaultdict(set))
    component_observed = defaultdict(lambda: defaultdict(set))
    def record_rule(rule, key, unit, status, reasons, source_ids, component_members=None):
        axis = rule["endpoint"]
        rid = rule["rule_id"]
        eid = _id("rule_evaluation", [rid, rule["version"], list(key) if key else None, unit["unit_id"] if unit else None, source_ids])
        item = dict(evaluation_id=eid, rule_id=rid, rule_version=rule["version"], route=rule["route"], endpoint=axis,
                    status="conflicting" if status == "conflict" else status, reason_codes=sorted(set(reasons)),
                    source_evidence_ids=source_ids, conclusion=deepcopy(rule["conclusion"]),
                    primary_sources=list(rule["primary_sources"]), defense_scope=rule.get("defense_scope", "unresolved"))
        if key:
            item.update(sample_id=key[0], gene_id=key[1])
        if unit:
            item.update(unit_id=unit["unit_id"], sample_id=unit["sample_id"])
        if component_members is not None:
            item["component_members"] = component_members
            item["architecture_contract"] = deepcopy(rule["architecture_contract"])
            item["scope"] = "registered_component_set_not_detector_grammar_replay"
        evaluations.append(item)
        target = gene_reviews[key]["axes"] if key else unit["axes"]
        for reason in reasons:
            _add_axis(target, axis, status, reason, eid, rid)
        if key and unit:
            # Only this exact unit's native rule receives this support. Intrinsic
            # gene evidence from another unit is exposed separately below.
            for reason in reasons:
                _add_axis(unit["axes"], axis, status, reason, eid, rid)
            if status == "supported":
                supported = unit["axes"][axis].setdefault("supported_member_keys", [])
                if list(key) not in supported:
                    supported.append(list(key))
        if key and status == "supported":
            if rule["conclusion"].get("role"):
                roles[key].add(rule["conclusion"]["role"])
                excluded_roles[key].update(rule["conclusion"].get("excludes_roles", []))
                role_sources[key][rule["conclusion"]["role"]].add(eid)
                for role in rule["conclusion"].get("excludes_roles", []):
                    exclusion_sources[key][role].add(eid)
                if rule["conclusion"]["role"] == "noncatalytic":
                    excluded_roles[key].add("nuclease")
                    exclusion_sources[key]["nuclease"].add(eid)
            if unit:
                component_support[unit["unit_id"]][rid].add(key)
        return item
    for did, call in call_index.items():
        key = (call["sample_id"], call["gene_id"])
        uid = detection_units[did]
        for rule in rule_list:
            if rule["route"] != "native_component" or not _native_scope(rule, call):
                continue
            component_observed[uid][rule["rule_id"]].add(key)
            reasons = _native_requirements(rule, call, genes[key], bindings)
            item = record_rule(rule, key, units[uid], "unassessed" if reasons else "supported", reasons or ["exact_registered_native_component"], [did])
            item["protein_binding"] = _protein_binding(call, genes[key])
            item["version_provenance"] = dict(model_version=_version(call,"model"), model_version_source=call.get("model_version_source", "canonical_model_version"),
                database_version=_version(call,"database"), database_version_source=call.get("database_version_source", "canonical_database_version"),
                software_version=_version(call,"software"), source_refs=deepcopy(_refs(call)))
    for key, rows in qualified_by_gene.items():
        for rule in rule_list:
            if rule["route"] not in {"qualified_domain", "qualified_site"}:
                continue
            for row in rows:
                match = _qualified_match(rule, row)
                if match is not None:
                    status, reasons = match
                    scopes = set(gene_reviews[key]["native_defense_scopes"]) - {"unresolved"}
                    requested_scope = rule.get("defense_scope", "unresolved")
                    if requested_scope == "host_defense" and scopes == {"anti_defense"} or requested_scope == "anti_defense" and scopes == {"host_defense"}:
                        status, reasons = "unassessed", ["defense_scope_mismatch"]
                    source_id = _id("registry", [key, row["registry_sha256"], row["entry"]["family_id"]])
                    record_rule(rule, key, None, status, reasons, [source_id])
    for unit in units.values():
        for rule in rule_list:
            if rule["route"] != "native_architecture" or unit["unit_kind"] != "native_system":
                continue
            eligible = [call_index[did] for did in unit["detection_ids"] if _native_scope(rule, call_index[did])]
            if not eligible:
                continue
            reasons = []
            if unit["native_identity_binding"] != "exact_canonical":
                reasons.append("native_unit_identity_unassessed")
            if set(unit["native_families"]) - set(rule["match"]["exact_native_families"]):
                reasons.append("native_family_contract_conflict")
            contract_ref = deepcopy(rule["architecture_contract"]["source_ref"])
            if not Path(contract_ref["path"]).is_absolute():
                contract_ref["path"] = str((bindings.rules_base / contract_ref["path"]).resolve())
            if not bindings.refs([contract_ref]):
                reasons.append("architecture_contract_source_unavailable")
            for call in eligible:
                reasons.extend(_native_requirements(rule, call, genes[(call["sample_id"], call["gene_id"])], bindings))
            components, chosen = [], []
            for requirement in rule["requires"]["all_components"]:
                members = set.union(set(), *(component_support[unit["unit_id"]].get(rid, set()) for rid in requirement["component_rule_ids"]))
                observed = set.union(set(), *(component_observed[unit["unit_id"]].get(rid, set()) for rid in requirement["component_rule_ids"]))
                components.append(dict(component=requirement["component"], member_keys=[list(k) for k in sorted(members)], observed_member_keys=[list(k) for k in sorted(observed)], minimum_members=requirement["minimum_members"], maximum_members=requirement["maximum_members"]))
                chosen.extend(members)
                if len(members) < requirement["minimum_members"]:
                    if observed - members:
                        reasons.append("required_component_unbound")
                    if len(observed) < requirement["minimum_members"]:
                        reasons.append("required_component_unobserved")
                if len(members) > requirement["maximum_members"]:
                    reasons.append("ambiguous_component_multiplicity")
            if len(chosen) != len(set(chosen)):
                reasons.append("one_gene_cannot_supply_distinct_components_without_fusion_contract")
            state = ("incomplete" if reasons and set(reasons) <= {"required_component_unobserved", "ambiguous_component_multiplicity", "one_gene_cannot_supply_distinct_components_without_fusion_contract"} else "unassessed") if reasons else "supported"
            item = record_rule(rule, None, unit, state, reasons or ["complete_registered_component_set"], unit["detection_ids"], components)
            if state == "supported" and rule["endpoint"] == "complex_core":
                for key in set(chosen):
                    _add_axis(gene_reviews[key]["axes"], "complex_core", "supported", "registered_complex_core_individual_catalysis_unassigned", item["evaluation_id"], rule["rule_id"])
    # Explicit exclusions are scoped claims. Different roles alone are not a
    # contradiction: fused and multi-domain proteins can have several functions.
    evaluations_by_gene, evaluations_by_unit = defaultdict(list), defaultdict(list)
    for evaluation in evaluations:
        if evaluation.get("gene_id") is not None:
            evaluations_by_gene[(evaluation["sample_id"], evaluation["gene_id"])].append(evaluation)
        if evaluation.get("unit_id") is not None:
            evaluations_by_unit[evaluation["unit_id"]].append(evaluation)
    for key, review in gene_reviews.items():
        conflict = roles[key] & excluded_roles[key]
        if conflict:
            _add_axis(review["axes"], "component_role", "conflicting", "conflicting_supported_assignments")
            review["axes"]["component_role"]["conflicting_evidence_ids"] = sorted(set.union(set(), *(role_sources[key][role] | exclusion_sources[key][role] for role in conflict)))
            review["conflicting_roles"] = sorted(conflict)
        review["proposed_or_preserved_roles"] = sorted(roles[key])
        review["registered_proposal_scopes"] = sorted({e["defense_scope"] for e in evaluations_by_gene[key] if e["status"] == "supported"})
    for unit in units.values():
        local = [e for e in evaluations_by_unit[unit["unit_id"]] if e["status"] == "supported" and e["endpoint"] == "component_role"]
        local_roles = {e["conclusion"]["role"] for e in local if e["conclusion"].get("role")}
        exclusions = {r for e in local for r in e["conclusion"].get("excludes_roles", [])}
        if "noncatalytic" in local_roles:
            exclusions.add("nuclease")
        conflicts = local_roles & exclusions
        if conflicts:
            _add_axis(unit["axes"], "component_role", "conflicting", "conflicting_supported_unit_assignments")
            unit["axes"]["component_role"]["conflicting_evidence_ids"] = [e["evaluation_id"] for e in local if e["conclusion"].get("role") in conflicts or set(e["conclusion"].get("excludes_roles", [])) & conflicts or e["conclusion"].get("role") == "noncatalytic"]
    # Gene intrinsic or preserved evidence is informative about a member but is
    # never a source/unit-specific claim. Complex completion in P1 cannot lend
    # its supported status to S1 merely because both contain the same CDS.
    for unit in units.values():
        unit["member_gene_evidence"] = [dict(member_key=member,
            scope="gene_intrinsic_or_other_unit_evidence_not_this_unit_contract",
            axes=deepcopy({axis:gene_reviews[tuple(member)]["axes"][axis] for axis in AXES if axis != "native_detection"}),
            unit_ids=gene_reviews[tuple(member)]["unit_ids"]) for member in unit["member_keys"]]
    coverage_rows = []
    options = data.get("options") or {}
    legacy = options.get("assessment_coverage", {})
    scoped = options.get("assessment_coverage_by_sample", {})
    allowed_coverage = {"complete", "partial", "unassessed"}
    if not isinstance(legacy, dict) or set(legacy) - {"effector", "associated", "motif"} or any(s not in allowed_coverage for s in legacy.values()):
        raise ValueError("Invalid legacy assessment_coverage declaration")
    if not isinstance(scoped, dict) or set(scoped) - samples:
        raise ValueError("Invalid sample assessment coverage identities")
    for declaration in scoped.values():
        if not isinstance(declaration, dict) or set(declaration) - {"effector", "associated", "motif"} or any(s not in allowed_coverage for s in declaration.values()):
            raise ValueError("Invalid scoped assessment coverage declaration")
    availability = result.get("availability") or {}
    resolved = availability.get("assessment_coverage_by_sample", {})
    resolved_provenance = availability.get("assessment_coverage_provenance_by_sample", {})
    if not isinstance(resolved, dict) or set(resolved) - samples:
        raise ValueError("Invalid resolved assessment coverage identities")
    for declaration in resolved.values():
        if not isinstance(declaration, dict) or set(declaration) - {"effector", "associated", "motif"} or any(s not in allowed_coverage for s in declaration.values()):
            raise ValueError("Invalid resolved assessment coverage declaration")
    if not isinstance(resolved_provenance, dict) or set(resolved_provenance) - samples:
        raise ValueError("Invalid resolved assessment coverage provenance")
    for sid in sorted(samples):
        sample_units = [u for u in units.values() if u["sample_id"] == sid]
        sample_genes = [r for k, r in gene_reviews.items() if k[0] == sid]
        source_rows = [r for r in data.get("coverage", []) if r.get("sample_id") == sid]
        observed = [r for k,r in ledger.items() if k[0] == sid]
        defaults = dict(effector="partial" if any(r.get("effector") or any(e.get("assignment_kind") == "nuclease_other_role" and e.get("status") in POSITIVE for e in r.get("effector_evidence", [])) for r in observed) else "unassessed",
                        associated="partial" if any(r.get("associated") for r in observed) else "unassessed",
                        motif="partial" if any(r.get("motif_site_supported") for r in observed) else "unassessed")
        declared = {axis: resolved.get(sid, {}).get(axis, scoped.get(sid, {}).get(axis, legacy.get(axis, defaults[axis]))) for axis in ("effector", "associated", "motif")}
        origins = {axis: deepcopy(resolved_provenance.get(sid, {}).get(axis, dict(
            origin="sample_declaration" if axis in scoped.get(sid,{}) else "explicit_global_fallback" if axis in legacy else "local_observed_evidence_default",
            option_path=f"options.assessment_coverage_by_sample[{sid!r}].{axis}" if axis in scoped.get(sid,{}) else f"options.assessment_coverage.{axis}" if axis in legacy else None,
            source_refs=[]))) for axis in declared}
        annotation_capabilities = [deepcopy(c) for c in (data.get("suite") or {}).get("capabilities", []) if c.get("sample_id") == sid and any(k.startswith(("raw_interpro", "raw_domain", "raw_eggnog")) for k in c)]
        row = dict(sample_id=sid, unit_count=len(sample_units), distinct_gene_count=len(sample_genes),
                   input_coverage=availability.get("inventory_by_sample", {}).get(sid, _safe_state(r.get("status") for r in source_rows)),
                   review_execution="completed", declared_assessment_coverage=declared,
                   coverage_origin={axis:origins[axis]["origin"] for axis in origins}, assessment_coverage_provenance=origins,
                   annotation_capabilities=annotation_capabilities,
                   axes={axis: dict(status_counts=dict(Counter(r["axes"][axis]["status"] for r in sample_genes)), assessed_gene_count=sum(r["axes"][axis]["status"] != "unassessed" for r in sample_genes)) for axis in AXES})
        coverage_rows.append(row)
        for review in sample_genes:
            for axis in AXES:
                coverage_axis = "motif" if axis == "motif_site" else "effector" if axis in {"component_role", "complex_core", "nuclease_domain"} else None
                review["axes"][axis]["coverage"] = declared[coverage_axis] if coverage_axis else row["input_coverage"] if axis == "native_detection" else "unassessed"
        for unit in sample_units:
            unit["axes"]["native_detection"]["coverage"] = _safe_state(r.get("status") for r in source_rows if r.get("caller") == unit["provider"])
            for item in unit["member_gene_evidence"]:
                item["axes"] = deepcopy({axis:gene_reviews[tuple(item["member_key"])]["axes"][axis] for axis in AXES if axis != "native_detection"})
    queue = []
    for kind, rows in (("unit", units.values()), ("gene", gene_reviews.values())):
        for row in rows:
            for axis, assessment in row["axes"].items():
                if axis in {"native_detection", "activity"} or assessment["status"] not in {"unassessed", "incomplete", "conflicting"}:
                    continue
                # Missing optional mechanistic axes are not demands that every
                # defense protein must contain a nuclease or catalytic motif.
                # Keep their unassessed scientific state, but queue only a
                # concrete applicable claim or unknown identity/mechanism.
                explicit_claim = bool(assessment["evaluated_rule_ids"]) or assessment["status"] in {"incomplete", "conflicting"}
                if kind == "unit" and axis in {"family_annotation", "component_role", "nuclease_domain", "motif_site"}:
                    continue  # Gene-level evidence owns these tasks.
                if kind == "gene" and axis == "system_architecture":
                    continue
                if not explicit_claim:
                    unknown_gene = kind == "gene" and axis in {"family_annotation", "component_role"}
                    unknown_system = kind == "unit" and row["unit_kind"] == "native_system" and axis == "system_architecture"
                    if not (unknown_gene or unknown_system):
                        continue
                queue_class = ("conflicting_claims" if assessment["status"] == "conflicting" else
                               "applicable_evidence_contract_unmet" if explicit_claim else "unregistered_mechanism")
                for reason in assessment["reason_codes"]:
                    target = row["unit_id"] if kind == "unit" else [row["sample_id"], row["gene_id"]]
                    requirement = ("registered_exact_rule" if reason == "no_applicable_rule" else
                                   "producer_database_version_or_model_hash" if "version" in reason or "digest" in reason else
                                   "qualified_site_map" if "site" in reason or axis == "motif_site" else
                                   "bound_native_component_members" if "component" in reason or axis in {"system_architecture", "complex_core"} else
                                   "source_bound_adjudication" if assessment["status"] == "conflicting" else "source_bound_axis_evidence")
                    queue.append(dict(queue_id=_id("queue", [kind, target, axis, reason]), target_kind=kind, target=target,
                        sample_id=row["sample_id"], axis=axis, status=assessment["status"], reason_code=reason,
                        queue_class=queue_class, required_for_claim=explicit_claim,
                        candidate_rule_ids=assessment["evaluated_rule_ids"], evidence_ids=assessment["supporting_evidence_ids"] + assessment["conflicting_evidence_ids"], required_evidence=[requirement]))
    queued_evaluations = {eid for item in queue for eid in item["evidence_ids"]}
    for evaluation in evaluations:
        if evaluation["status"] not in {"unassessed", "incomplete", "conflicting"} or evaluation["evaluation_id"] in queued_evaluations:
            continue
        # Support from a broad route does not erase the unmet requirements of an
        # explicitly requested narrower rule on the same gene/axis.
        for reason in evaluation["reason_codes"]:
            queue.append(dict(queue_id=_id("queue", ["rule_evaluation", evaluation["evaluation_id"], reason]),
                target_kind="rule_evaluation", target=evaluation["evaluation_id"], sample_id=evaluation["sample_id"],
                axis=evaluation["endpoint"], status=evaluation["status"], reason_code=reason,
                queue_class="conflicting_claims" if evaluation["status"] == "conflicting" else "applicable_evidence_contract_unmet",
                required_for_claim=True, candidate_rule_ids=[evaluation["rule_id"]],
                evidence_ids=[evaluation["evaluation_id"], *evaluation["source_evidence_ids"]], required_evidence=["rule_required_evidence"]))
    quarantine = deepcopy((data.get("suite") or {}).get("unmapped_records", []))
    expected_detections = len(data.get("calls", []))
    if len(accounting) != expected_detections or len({r["detection_id"] for r in accounting}) != expected_detections or not set(ledger) <= set(gene_reviews):
        raise ValueError("System review accounting lost original detections or accepted genes")
    totals = dict(canonical_detection_count=expected_detections, accounted_detection_count=len(accounting),
                  raw_observation_count=raw_rows, unit_count=len(units), distinct_reviewed_gene_count=len(gene_reviews),
                  accepted_gene_count=len(ledger), accepted_genes_accounted=len(set(ledger) & set(gene_reviews)),
                  saved_review_gene_count=len(stored_reviews), saved_review_genes_accounted=len(set(stored_reviews) & set(gene_reviews)),
                  quarantined_raw_record_count=len(quarantine), unaccounted_detection_count=0, unaccounted_accepted_gene_count=0,
                  unit_counts=dict(Counter(u["unit_kind"] for u in units.values())))
    verified_candidate_keys = {(c["sample_id"], c["query_id"]) for c in (data.get("system_review_evidence") or {}).get("candidates", [])}
    totals.update(verified_evidence_candidate_count=len(verified_candidate_keys),
                  verified_evidence_candidates_outside_detection_scope=[list(k) for k in sorted(verified_candidate_keys - target_keys)])
    return dict(schema_version=1, analysis_kind="evaluation_of_saved_system_evidence", review_execution={"status": "completed", "unit_count": len(units)},
                units=sorted(units.values(), key=lambda u: u["unit_id"]), gene_reviews=[gene_reviews[k] for k in sorted(gene_reviews)],
                evidence=evidence, rule_evaluations=evaluations, pending_queue=queue, detection_accounting=accounting,
                quarantine=quarantine, sample_coverage=coverage_rows, accounting=totals,
                endpoint_policy="review_proposals_only_no_accepted_membership_changes",
                provenance=dict(detector_runs=0, scoring_runs=0, input_sha256=deepcopy(data.get("input_sha256", {})),
                    verified_evidence_sha256=bindings.verified, rules_sha256=hashlib.sha256(_canonical(rules).encode()).hexdigest() if rules else None,
                    rules_source=deepcopy(data.get("system_review_rule_source")),
                    existing_component_registry_sha256=result.get("provenance", {}).get("component_registry_sha256")))

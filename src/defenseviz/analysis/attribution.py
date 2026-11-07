"""Source-bound family attribution, separate from saved R-M score decisions.

The input is a canonical, hash-verified suite import. Curation's explicit
unproven R-M association plus a complete native component contract can support
another system scope. Neither methylase similarity nor caller agreement is a
system attribution rule. Original arrays are immutable; projection is explicit.
No thresholds are fitted, sequence searches run, or activity inferred here.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
import json
import re


CANONICAL = {"Type I", "Type II", "Type III", "Type IV"}
UNPROVEN = {"DNA_MTase_RM_association_unproven", "RM_association_unproven"}
SUPPORTED_RM = {"classic_RM_supported", "fused_RM_supported", "modification_dependent_restriction_supported"}
CONFLICT_FIELDS = ("other_defense_conflict", "evidence_conflict", "role_conflict",
                   "unresolved_role_conflict", "family_conflict", "non_dna_mtase_conflict")
POSITIVE = {"supported", "accepted", "supported_annotation", "supported_putative_effector_annotation"}

# Producer-defined component requirements and maximum separation, not a new
# biological distance cutoff. A native final unit still must bind every member.
# DISARM DrmE is optional in the producer contract, so it is not newly required.
NATIVE_CONTEXT_CONTRACTS = (
    dict(rule_id="padloc_brex_type_III_scope_v1", provider="PADLOC", native_family="brex_type_III",
         family="BREX", category="Modification", target_component="PglXI",
         required_components=["BrxA", "BrxC", "BrxF", "PglXI", "BrxHII", "PglZ"],
         prohibited_components=[], maximum_separation=3,
         definition_sha256="438d123279bf6a16c6ce786473940bcd37536823a8ae606ab9201b4dffd81272"),
    dict(rule_id="padloc_disarm_type_II_scope_v1", provider="PADLOC", native_family="disarm_type_II",
         family="DISARM", category="Modification", target_component="DrmMII",
         required_components=["DrmA", "DrmB", "DrmC", "DrmMII"],
         prohibited_components=["DrmD", "DrmMI"], maximum_separation=1,
         definition_sha256="a13046d89b8c92332676237f826561d56adcebc243df33b57bbc859f41350a51"),
    dict(rule_id="padloc_hma_scope_v1", provider="PADLOC", native_family="hma",
         family="Hma", category="Nuclease", target_component="HmaB",
         required_components=["HmaA", "HmaB", "HmaC"], prohibited_components=[], maximum_separation=3,
         definition_sha256="0a9e474979443884edf3bee472de7a14812d6feedb431c92392d3bdab9b96d45"),
)


def _key(row):
    return row.get("sample_id"), row.get("gene_id")


def _refs(row):
    refs = list(row.get("source_refs") or [])
    if row.get("source_ref") and row["source_ref"] not in refs:
        refs.append(row["source_ref"])
    return refs


def _digest(value):
    return isinstance(value, str) and re.fullmatch("[0-9a-f]{64}", value) is not None


def _stable(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _true(value):
    return value is True or isinstance(value, str) and value.upper() == "TRUE"


class _Evidence:
    """Bind canonical records to the importer's exact file/hash inventory.

    Files are verified by suite intake, not reopened by this pure analysis.
    A digest found under a different unrelated path does not bind a source.
    """
    def __init__(self, data):
        self.hashes = data.get("input_sha256") or {}
        self.aliases = defaultdict(set)
        for row in (data.get("suite") or {}).get("sources", []):
            if row.get("path") and row.get("resolved_path"):
                self.aliases[(row["path"], row.get("sha256"))].add(row["resolved_path"])

    def refs(self, refs):
        if not refs:
            return False
        for ref in refs:
            path, sha = ref.get("path"), ref.get("sha256")
            if not path or not _digest(sha):
                return False
            if path in self.hashes:
                if self.hashes[path] != sha:
                    return False
                continue
            paths = self.aliases.get((path, sha), set())
            if len(paths) != 1 or self.hashes.get(next(iter(paths))) != sha:
                return False
        return True

    def bound(self, row, gene):
        return (gene is not None and _key(row) == _key(gene)
                and _digest(row.get("protein_sha256"))
                and row["protein_sha256"] == gene.get("protein_sha256")
                and self.refs(_refs(row)) and self.refs(_refs(gene)))


def _coordinates(gene):
    parts = gene.get("location_parts")
    if parts is not None:
        if not isinstance(parts, list) or len(parts) > 1:
            return False
        if parts:
            part = parts[0]
            if (not isinstance(part, (list, tuple)) or len(part) != 3
                    or any(type(value) is not int for value in part)
                    or list(part[:2]) != [gene.get("start"), gene.get("end")]
                    or part[2] != {"+":1, "-":-1}.get(gene.get("strand"))):
                return False
    return (gene.get("contig") and gene.get("strand") in {"+", "-"}
            and type(gene.get("start")) is int and type(gene.get("end")) is int
            and 1 <= gene["start"] <= gene["end"]
            and (parts is None or isinstance(parts, list) and len(parts) <= 1))


def _context(rule, unit, genes, evidence):
    reasons, members = [], []
    components = Counter()
    for call in unit:
        gene = genes.get(_key(call)); native = call.get("native_fields") or {}
        if not evidence.bound(call, gene):
            reasons.append("unbound_native_component")
            continue
        if not _coordinates(gene):
            reasons.append("unresolved_component_coordinates")
            continue
        if (native.get("system") != rule["native_family"]
                or native.get("hmm.name") != call.get("native_model")):
            reasons.append("native_component_metadata_conflict")
            continue
        if (native.get("seqid") != gene["contig"] or native.get("strand") != gene["strand"]
                or str(native.get("start")) != str(gene["start"])
                or str(native.get("end")) != str(gene["end"])):
            reasons.append("native_coordinate_binding_conflict")
            continue
        component = native.get("protein.name")
        components[component] += 1
        members.append(dict(sample_id=call["sample_id"], gene_id=call["gene_id"],
                            protein_sha256=call["protein_sha256"], component=component,
                            contig=gene["contig"], start=gene["start"], end=gene["end"],
                            source_refs=deepcopy(_refs(call))))
    if not set(rule["required_components"]) <= set(components):
        reasons.append("incomplete_native_component_set")
    if set(rule["prohibited_components"]) & set(components):
        reasons.append("prohibited_native_component")
    if any(components[c] != 1 for c in rule["required_components"] if c in components):
        reasons.append("ambiguous_native_component_multiplicity")
    if len({_key(row) for row in members}) != len(members):
        reasons.append("ambiguous_native_member_identity")
    if len({m["contig"] for m in members}) != 1:
        reasons.append("different_or_missing_contigs")
    # Read the native protein ordinal, which is the coordinate system used by
    # the PADLOC separation contract. Genomic spans remain independently bound.
    positions = []
    for call in unit:
        try:
            value = (call.get("native_fields") or {})["relative.position"]
            if isinstance(value, bool) or str(int(value)) != str(value) or int(value) < 1:
                raise ValueError
            positions.append(int(value))
        except (KeyError, TypeError, ValueError):
            reasons.append("unassessed_native_gene_separation")
    if len(positions) == len(unit):
        ordered = sorted(zip(positions, unit), key=lambda item: item[0])
        if len(set(positions)) != len(positions):
            reasons.append("ambiguous_native_gene_positions")
        # PADLOC counts intervening proteins: ordinal difference <= space+1.
        elif any(b-a > rule["maximum_separation"] + 1 for a,b in zip(sorted(positions), sorted(positions)[1:])):
            reasons.append("outside_producer_gene_separation")
        if all(genes.get(_key(c)) and _coordinates(genes[_key(c)]) for _,c in ordered):
            starts = [genes[_key(c)]["start"] for _,c in ordered]
            if starts != sorted(starts):
                reasons.append("native_position_coordinate_conflict")
    return dict(rule_id=rule["rule_id"], family=rule["family"], category=rule["category"],
                provider=rule["provider"], native_family=rule["native_family"],
                system_id=unit[0]["system_id"], run_id=unit[0]["run_id"],
                definition_source=dict(url="https://github.com/padlocbio/padloc-db/blob/v2.0.0/sys/" + rule["native_family"] + ".yaml",
                                       database_version="v2.0.0", sha256=rule["definition_sha256"]),
                required_components=rule["required_components"], maximum_separation=rule["maximum_separation"],
                target_component=rule["target_component"], members=members,
                status="supported" if not reasons else "unassessed", reason_codes=sorted(set(reasons)))


def _rm_protection(key, curation, decision, data, evidence, gene):
    reasons = []
    if curation.get("rm_association_class") in SUPPORTED_RM:
        reasons.append("explicit_canonical_RM_association")
    if decision.get("role") in {"R", "RM", "S"} or curation.get("final_role") in {"R", "RM", "S"}:
        reasons.append("independent_canonical_RM_role")
    if curation.get("final_family") in CANONICAL and curation["final_family"] != decision.get("rm_type", decision.get("family")):
        reasons.append("curation_saved_family_conflict")
    if any(curation.get(f) is True for f in CONFLICT_FIELDS):
        reasons.append("unresolved_curation_conflict")
    for row in [decision, decision.get("native_context") or {}, curation.get("native_fields") or {}]:
        if (row.get("context_state") in {"met", "complete", "met_in_saved_context"}
                or row.get("context_support") is True
                or any(_true(row.get(f)) for f in ("typei_operon_supported", "typeii_operon_supported", "typeiii_operon_supported"))):
            reasons.append("supported_canonical_RM_context")
    for group in data.get("groups", []):
        if (gene is not None and group.get("sample_id") == key[0]
                and group.get("rm_type", group.get("family")) in CANONICAL
                and group.get("context_state") in {"met", "complete", "met_in_saved_context"}
                and any(m.get("gene_id") == key[1] and m.get("protein_sha256") == gene.get("protein_sha256")
                        for m in group.get("members", []))
                and evidence.refs(_refs(group))):
            reasons.append("independent_saved_canonical_context")
    for row in data.get("annotations", []):
        if (_key(row) == key and row.get("status") in POSITIVE
                and row.get("assignment_kind") == "reviewed_RM_component"
                and evidence.bound(row, gene)):
            # A review explicitly describing an orphan M establishes chemistry,
            # not cognate canonical-system membership. Other reviews are kept.
            if not (row.get("orphan") is True and row.get("role") == "M"
                    and row.get("family", row.get("rm_type")) == decision.get("family", decision.get("rm_type"))):
                reasons.append("independent_reviewed_RM_assignment")
    for row in data.get("calls", []):
        if (_key(row) == key and row.get("selected", True) and evidence.bound(row, gene)
                and (row.get("rm_type") in CANONICAL or row.get("category") == "R-M")):
            if row.get("evidence_kind") == "final_member":
                reasons.append("independent_native_RM_membership")
            elif _true((row.get("native_fields") or {}).get("rebase_reference_is_gold_standard")):
                reasons.append("characterized_canonical_reference_conflict")
    return sorted(set(reasons))


def apply_family_attribution(data):
    """Return a deep copy with additive, source-bound ``family_attribution``.

    The canonical calls, decisions, annotations, and configured memberships are
    unchanged. Unsupported/malformed/ambiguous evidence produces no override.
    The three finite contracts are annotation rules, not exhaustive discovery.
    """
    if not isinstance(data, dict):
        raise ValueError("Family attribution requires a canonical data object")
    out = deepcopy(data)
    genes = {}
    for row in data.get("genes", []):
        if _key(row) in genes:
            raise ValueError("Duplicate gene identity in family attribution")
        genes[_key(row)] = row
    evidence = _Evidence(data)
    decisions = {}
    for row in data.get("saved_decisions", []):
        if _key(row) in decisions:
            raise ValueError("Duplicate saved decision in family attribution")
        decisions[_key(row)] = row
    curations = defaultdict(list)
    for row in data.get("curation_records", []):
        curations[_key(row)].append(row)
    units = defaultdict(list)
    for call in data.get("calls", []):
        if (call.get("evidence_kind") == "final_member" and call.get("selected", True)
                and call.get("system_id") and call.get("run_id")):
            units[(call.get("sample_id"), call.get("caller"), call["run_id"],
                   call["system_id"], call.get("native_family"))].append(call)
    contexts = defaultdict(list)
    for rule in NATIVE_CONTEXT_CONTRACTS:
        for (sid, provider, _, _, family), calls in units.items():
            if provider != rule["provider"] or family != rule["native_family"]:
                continue
            context = _context(rule, calls, genes, evidence)
            for call in calls:
                if (call.get("native_fields") or {}).get("protein.name") == rule["target_component"]:
                    contexts[_key(call)].append(context)
    rows = []
    for key, context_rows in sorted(contexts.items()):
        decision = decisions.get(key, {})
        if decision.get("state") != "retained" or decision.get("rm_type", decision.get("family")) not in CANONICAL:
            continue
        gene = genes.get(key)
        reasons = []
        curation_rows = curations.get(key, [])
        if len(curation_rows) != 1:
            reasons.append("missing_or_ambiguous_curation")
        curation = curation_rows[0] if len(curation_rows) == 1 else {}
        if (curation.get("schema_version") != 1 or curation.get("source_kind") != "DNMB_REBASEfinder_curation"
                or curation.get("usable_for_attribution") is not True
                or curation.get("binding_status") != "bound"
                or not evidence.bound(curation, gene)):
            reasons.append("unbound_curation")
        if not evidence.bound(decision, gene):
            reasons.append("unbound_saved_decision")
        native = curation.get("native_fields") or {}
        if native.get("query", native.get("locus_tag", native.get("gene_id"))) != key[1]:
            reasons.append("curation_native_query_conflict")
        if any(curation.get(f) != native.get(f) for f in ("final_family", "final_role", "rm_association_class")):
            reasons.append("curation_projection_conflict")
        if any((curation.get(f) is True) != _true(native.get(f)) for f in CONFLICT_FIELDS):
            reasons.append("curation_conflict_projection_mismatch")
        protections = _rm_protection(key, curation, decision, data, evidence, gene)
        supported = [c for c in context_rows if c["status"] == "supported"]
        families = {c["family"] for c in supported}
        if not supported:
            reasons.append("no_complete_supported_native_context")
        if len(families) > 1:
            protections.append("conflicting_noncanonical_families")
        if (curation.get("rm_association_class") not in UNPROVEN
                or curation.get("final_role") != "M" or decision.get("role") != "M"):
            reasons.append("canonical_scope_not_explicitly_unproven_methylase")
        action = ("retain_dual_or_conflicting" if supported and protections
                  else "assign_noncanonical_scope" if not reasons and not protections
                  else "retain_original_unassessed")
        target = supported[0] if action == "assign_noncanonical_scope" else None
        identity = dict(sample_id=key[0], gene_id=key[1], protein_sha256=(gene or {}).get("protein_sha256"))
        record = dict(**identity, schema_version=1, action=action,
                      status="supported" if target else "conflicting" if supported and protections else "unassessed",
                      reason_codes=sorted(set(reasons+protections)),
                      target_family=target["family"] if target else None,
                      target_category=target["category"] if target else None,
                      preserved_component_role="M" if target else decision.get("role"),
                      previous_assignment=deepcopy(decision), curation=deepcopy(curation),
                      native_contexts=deepcopy(context_rows), activity="unassessed",
                      score_changed=False, native_observations_changed=False)
        record["transition_id"] = "family-attribution:" + _stable(record)
        rows.append(record)
    out["family_attribution"] = rows
    out["family_attribution_summary"] = dict(schema_version=1, actions=dict(Counter(r["action"] for r in rows)),
        override_keys=[list(_key(r)) for r in rows if r["action"] == "assign_noncanonical_scope"],
        evidence_scope="source-bound curation and complete native component annotation; no activity inference",
        engine_runs=0, scoring_runs=0)
    return out


def project_family_attribution(data):
    """Return an effective aggregation view, recomputing every scope override.

    Stored transition flags alone cannot authorize an override. Original records
    are nested in the projection archive; calls and curation stay unchanged.
    A saved score remains unchanged while its derived R-M state is out_of_scope.
    """
    if data.get("family_attribution_projection"):
        raise ValueError("Family attribution projection requires the unprojected canonical bundle")
    out = apply_family_attribution(data)
    overrides = {_key(r): r for r in out["family_attribution"] if r["action"] == "assign_noncanonical_scope"}
    archive = dict(saved_decisions=[], annotations=[], family_memberships=[])
    for decision in out.get("saved_decisions", []):
        transition = overrides.get(_key(decision))
        if not transition:
            continue
        original = deepcopy(decision); archive["saved_decisions"].append(original)
        decision.update(state="out_of_scope", family=transition["target_family"], rm_type=None,
                        role=None, final_component=None, prior_saved_decision=original,
                        decision_reason="source_bound_noncanonical_system_scope",
                        family_attribution_transition_id=transition["transition_id"],
                        derived_assignment=True, activity="unassessed")
    for row in out.get("annotations", []):
        transition = overrides.get(_key(row))
        if transition and row.get("assignment_kind") == "reviewed_RM_component":
            archive["annotations"].append(deepcopy(row))
            row.update(status="out_of_scope", include_in_union=False,
                       family_attribution_transition_id=transition["transition_id"],
                       derived_assignment=True)
    options = out.get("options") or {}
    if options.get("family_memberships") is not None:
        catalog = {r.get("id", r.get("family_id")): r for r in options.get("family_catalog", [])}
        memberships = []
        for row in options["family_memberships"]:
            fid = row.get("family_id")
            category = row.get("category") or catalog.get(fid, {}).get("category")
            if _key(row) in overrides and (category == "R-M" or fid in CANONICAL | {"R-M"}):
                archive["family_memberships"].append(deepcopy(row))
            else:
                memberships.append(row)
        existing = {(_key(r), r.get("family_id")) for r in memberships}
        for key, t in overrides.items():
            if (key, t["target_family"]) not in existing:
                refs = _refs(t["curation"])
                memberships.append(dict(sample_id=key[0], gene_id=key[1], protein_sha256=t["protein_sha256"],
                    family_id=t["target_family"], category=t["target_category"], source_refs=deepcopy(refs),
                    family_attribution_transition_id=t["transition_id"], derived_assignment=True))
        options["family_memberships"] = memberships
    out["family_attribution_projection"] = dict(schema_version=1, operation="derived_effective_RM_scope",
        original_records=archive, override_keys=[list(k) for k in sorted(overrides)],
        original_scores_preserved=True, native_calls_preserved=True,
        original_saved_decisions_preserved_in="original_records.saved_decisions and prior_saved_decision")
    return out

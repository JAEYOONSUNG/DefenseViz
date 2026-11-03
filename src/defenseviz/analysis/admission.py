"""Source-bound reference-anchor admission without sample or locus whitelists."""
import hashlib
import json
from copy import deepcopy
from dataclasses import replace, asdict
from pathlib import Path
from .evidence import parse_evidence, parse_policy
from .policy import Profile, profile_pass
from . import registry as core


def entry_sha256(entry):
    """Digest of the exact registered family contract, with stable JSON encoding."""
    return hashlib.sha256(json.dumps(entry, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                                     allow_nan=False).encode("utf-8")).hexdigest()


def validate_admission_rules(registry, rules, *, registry_dir=Path("."), source_root=None):
    """Validate every rule before candidate filtering, even for an empty batch."""
    core.validate_registry(registry, registry_dir)
    if not isinstance(rules, list):
        raise ValueError("Admission rules must be a list")
    allowed = {"family_id", "RM_families", "roles", "entry_sha256", "model_accession", "model_sha256",
               "reference_accession", "required_anchor_subset", "source", "source_sha256"}
    entries = {x["family_id"]: x for x in registry["families"]}
    root = Path(source_root) if source_root is not None else Path(registry_dir)
    sources = {}
    for entry in registry["families"]:
        for ref in entry.get("references", []):
            if "sequence" not in ref:
                path = Path(ref["sequence_file"])
                if not path.is_absolute():
                    path = Path(registry_dir) / path
                path = path.resolve()
                sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    for rule in rules:
        if not isinstance(rule, dict) or set(rule) != allowed:
            raise ValueError("Admission rule needs the exact documented family/model/reference/source fields; identity whitelists are forbidden")
        for name in allowed - {"RM_families", "roles", "required_anchor_subset"}:
            if not isinstance(rule[name], str) or not rule[name]:
                raise ValueError(f"Admission {name} must be a nonempty string")
        for name in ("RM_families", "roles", "required_anchor_subset"):
            expected_type = int if name == "required_anchor_subset" else str
            if not isinstance(rule[name], list) or not rule[name] or any(type(v) is not expected_type for v in rule[name]) or len(rule[name]) != len(set(rule[name])):
                raise ValueError(f"Admission {name} must be a nonempty unique list")
        if any(x not in {"Type I", "Type II", "Type III"} for x in rule["RM_families"]):
            raise ValueError("Reference-anchor score admission supports explicitly registered Type I/II/III families")
        if any(x not in {"R", "RM"} for x in rule["roles"]):
            raise ValueError("Anchor score admission requires an R or RM role")
        if any(type(x) is not int or x < 1 for x in rule["required_anchor_subset"]):
            raise ValueError("Anchor positions must be positive 1-based integers")
        core._verify_file_binding(rule["source"], rule["source_sha256"], root, sources)
        entry = entries.get(rule["family_id"])
        if entry is None or entry_sha256(entry) != rule["entry_sha256"]:
            raise ValueError("Admission rule does not bind the exact registered family entry")
        core._sha256(rule["model_sha256"], "Admission model_sha256")
        ref = next((r for r in entry.get("references", []) if r["accession"] == rule["reference_accession"]), None)
        if ref is None:
            raise ValueError("Admission rule reference is absent from its family entry")
        required = {site["position"] for site in ref.get("catalytic_sites", []) if site.get("role") == "nuclease" and site.get("required_for_reference_subset", True) is not False}
        if set(rule["required_anchor_subset"]) != required:
            raise ValueError("Admission anchor subset must equal the registered required subset")
        if not any(m.get("accession") == rule["model_accession"] and m.get("model_sha256") == rule["model_sha256"] and m.get("evidence_kind") == "nuclease_domain" and m.get("role") == "nuclease" for m in entry.get("domain_models", [])):
            raise ValueError("Admission rule must bind an exact registered nuclease-domain model")
    return sources


def admit_reference_anchors(evidence, candidate, registry, rules, *, registry_dir=Path("."),
                            source_root=None, policy=None):
    """Return (Evidence, provenance) after verified, coherent-domain admission.

    Each rule explicitly binds an entry digest, model, reference and required
    anchor subset to compatible R-M types/roles. Rules are portable family-level
    contracts; unknown fields (including locus/sample selectors) are rejected.
    No score is assigned here. evaluate() applies the unchanged DNMB conditions.
    """
    e = parse_evidence(evidence)
    p = parse_policy(policy)
    assessed = core.assess_candidate(candidate, registry, registry_dir, source_root=source_root,
                                     verify_sources=True)
    sources = validate_admission_rules(registry, rules, registry_dir=registry_dir, source_root=source_root)
    sources.update(assessed["source_sha256"])
    entries = {x["family_id"]: x for x in registry["families"]}
    evaluations = {x["family_id"]: x for x in assessed["families"]}
    selected = []
    seen_domains = set()
    for rule in rules:
        entry = entries[rule["family_id"]]
        family = evaluations.get(rule["family_id"])
        if e.family not in rule["RM_families"] or e.role not in rule["roles"] or entry["substrate"] not in {"DNA", "dual"} or family is None:
            continue
        if any(family[k]["state"] != "supported" for k in ("family_evidence", "domain_evidence", "catalytic_site_evidence")):
            continue
        ref = next((r for r in entry.get("references", []) if r["accession"] == rule["reference_accession"]), None)
        if ref is None:
            raise ValueError("Admission rule reference is absent from its family entry")
        required = {s["position"] for s in ref.get("catalytic_sites", []) if s.get("role") == "nuclease" and s.get("required_for_reference_subset", True) is not False}
        if set(rule["required_anchor_subset"]) != required:
            raise ValueError("Admission anchor subset must equal the registered required subset")
        _, usable, _ = core.assess_domains(candidate, entry)
        for assessment in family["catalytic_site_evidence"].get("references", []):
            if assessment["state"] != "supported" or assessment.get("evidence_route") != "qualified_domain" or assessment["reference"] != ref["accession"]:
                continue
            sites = assessment["sites"]
            if {s["reference_position"] for s in sites} != required or not all(s["conserved"] is True for s in sites):
                continue
            for item in usable:
                d, model = item["hit"], item["model"]
                if model["evidence_kind"] != "nuclease_domain" or model["role"] != "nuclease":
                    continue
                if d.get("accession") != rule["model_accession"] or d.get("model_sha256") != rule["model_sha256"]:
                    continue
                if not all(d["qstart"] <= s["query_position"] <= d["qend"] for s in sites):
                    continue
                ev = d.get("independent_evalue")
                if not core.finite_number(ev) or ev < 0:
                    continue
                if not core.valid_interval(d.get("hmm_start"), d.get("hmm_end"), d.get("hmm_length")):
                    continue
                profile = Profile(True, True, ev, (d["qend"]-d["qstart"]+1)/len(candidate["sequence"]),
                                  (d["hmm_end"]-d["hmm_start"]+1)/d["hmm_length"])
                if not profile_pass(profile, p):
                    continue
                key = (rule["family_id"], d["accession"], d["model_sha256"], d["qstart"], d["qend"], d["hmm_start"], d["hmm_end"])
                if key in seen_domains:
                    continue
                seen_domains.add(key)
                selected.append(dict(profile=asdict(profile), family_id=rule["family_id"], entry_sha256=rule["entry_sha256"],
                    model_accession=d["accession"], model_sha256=d["model_sha256"], query_interval=[d["qstart"], d["qend"]],
                    reference_accession=ref["accession"], required_anchor_subset=sorted(required)))
    if selected:
        start = len(e.profiles)
        for i, source in enumerate(selected):
            source['profile_index'] = start+i
        e = replace(e, profiles=e.profiles+tuple(Profile(**r["profile"]) for r in selected),
                    reference_nuclease_anchor_support=True,
                    reference_nuclease_anchor_profile_indices=e.reference_nuclease_anchor_profile_indices+tuple(range(start, start+len(selected))))
    return e, dict(sample_id=candidate["sample_id"], gene_id=candidate["query_id"], sequence_sha256=candidate["sequence_sha256"],
                   sources=selected, source_sha256=sources, source_files_verified=True,
                   activity="not_assessed", score_changed=False)


def apply_reference_anchors(evidence_records, candidates, registry, admission_rules, *,
                            source_root=None, registry_dir=Path("."), policy=None):
    """Enrich exact sample/gene/sequence matches and return (records, audit).

    Records require provenance.sequence_sha256 to admit new evidence. Candidate
    gene_id is explicit; query_id is accepted as an exact alias only if the two
    fields agree. Identifiers are never split, normalized or joined by sequence
    alone. Candidates outside the evidence pool are reported, not injected.
    """
    from .evidence import evidence_dict
    if not isinstance(evidence_records, list) or not isinstance(candidates, list):
        raise ValueError("Evidence records and candidates must be lists")
    validated_sources = validate_admission_rules(registry, admission_rules, registry_dir=registry_dir, source_root=source_root)
    parse_policy(policy)
    rows = deepcopy(evidence_records)
    index = {}
    for i, r in enumerate(rows):
        if not isinstance(r, dict):
            raise ValueError("Evidence record must be an object")
        key = (r.get("sample_id"), r.get("gene_id"))
        if any(not isinstance(x, str) or not x for x in key):
            raise ValueError("Evidence records require explicit sample_id and gene_id")
        if key in index:
            raise ValueError(f"Duplicate evidence sample/gene identity: {key!r}")
        index[key] = i
    seen = set()
    audit = []
    all_sources = dict(validated_sources)
    for original in candidates:
        if not isinstance(original, dict):
            raise ValueError("Candidate must be an object")
        c = deepcopy(original)
        gene_id = c.get("gene_id", c.get("query_id"))
        if "query_id" in c and c["query_id"] != gene_id:
            raise ValueError("Candidate gene_id/query_id alias conflict")
        key = (c.get("sample_id"), gene_id)
        if any(not isinstance(x, str) or not x for x in key):
            raise ValueError("Candidates require explicit sample_id and gene_id/query_id")
        if key in seen:
            raise ValueError(f"Duplicate candidate sample/gene identity: {key!r}")
        seen.add(key)
        c["query_id"] = gene_id
        verified = core.assess_candidate(c, registry, registry_dir, source_root=source_root, verify_sources=True)
        all_sources.update(verified["source_sha256"])
        if key not in index:
            audit.append(dict(sample_id=key[0], gene_id=key[1], status="outside_evidence_pool"))
            continue
        row = rows[index[key]]
        if row.get("evidence") is None:
            audit.append(dict(sample_id=key[0], gene_id=key[1], status="unassessed_missing_policy_evidence"))
            continue
        provenance = row.get("provenance") or {}
        expected = provenance.get("sequence_sha256")
        core._sha256(expected, "evidence provenance.sequence_sha256 required for enrichment")
        if expected != c.get("sequence_sha256"):
            raise ValueError(f"Exact sample/gene sequence mismatch for {key!r}")
        enriched, binding = admit_reference_anchors(row["evidence"], c, registry, admission_rules,
            registry_dir=registry_dir, source_root=source_root, policy=policy)
        previous = provenance.get("reference_anchor_admission")
        if previous is not None:
            if not isinstance(previous, dict) or not isinstance(previous.get("sources"), list):
                raise ValueError("Existing anchor-admission provenance is malformed")
            old_e = parse_evidence(row["evidence"])
            old_sources = previous["sources"]
            def semantic_sources(items):
                return [{k:v for k,v in item.items() if k != "profile_index"} for item in items]
            if semantic_sources(old_sources) != semantic_sources(binding["sources"]):
                raise ValueError("Existing anchor admission differs from supplied evidence; use original unextended policy records")
            for item in old_sources:
                i = item.get("profile_index")
                if type(i) is not int or i not in old_e.reference_nuclease_anchor_profile_indices or not 0 <= i < len(old_e.profiles) or asdict(old_e.profiles[i]) != item.get("profile"):
                    raise ValueError("Existing anchor admission does not bind its actual linked profile")
            # Source verification still occurred. Reuse the actual existing
            # profile indices rather than appending identical profile evidence.
            enriched = old_e
            binding["sources"] = deepcopy(old_sources)
        row["evidence"] = evidence_dict(enriched)
        row["provenance"] = dict(provenance, reference_anchor_admission=binding)
        all_sources.update(binding["source_sha256"])
        audit.append(dict(sample_id=key[0], gene_id=key[1], status="admitted" if binding["sources"] else "no_qualified_anchor_admission",
                          admitted_domains=len(binding["sources"]), sequence_sha256=expected))
    return rows, dict(records=audit, source_sha256=all_sources,
        analysis_kind="source_bound_reference_anchor_evidence_admission_not_detector_run",
        activity="not_assessed", thresholds_changed=False)

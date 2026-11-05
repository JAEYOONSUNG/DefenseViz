"""Independent identity/evidence controls for the additive review engine."""
from copy import deepcopy
import hashlib
import json

import pytest

from defenseviz.analysis.registry import assess_registry
from defenseviz.analysis.system_review import review_all_systems, validate_system_review_rules


def sha(value):
    return hashlib.sha256(value.encode()).hexdigest()


def saved(tmp_path, *ids):
    source = tmp_path / "native.tsv"
    source.write_text("synthetic observations\n")
    ref = {"path": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "row": 2}
    genes = [dict(sample_id="s", gene_id=g, protein_sha256=sha("ADCEFGHIKLMN"), protein_length=12, contig="contig") for g in ids]
    calls = [dict(sample_id="s", gene_id=g, caller="Tool", native_family="Future", family="Future", native_model="F__"+g,
                  system_id="sys", run_id="run", replicon_id="contig", evidence_kind="final_member", selected=True,
                  database_version="1", source_ref=deepcopy(ref)) for g in ids]
    data = dict(samples=[dict(sample_id="s", label="Synthetic")], genes=genes, calls=calls, groups=[],
                coverage=[dict(sample_id="s", caller="Tool", status="complete")], input_sha256={str(source): ref["sha256"]})
    result = dict(gene_ledger=[dict(sample_id="s", gene_id=g, protein_sha256=sha("ADCEFGHIKLMN"), effector=False, associated=False, priority=False) for g in ids], provenance={})
    return data, result


def component(g="a", **changes):
    rule = dict(rule_id="rule_"+g, version="1", route="native_component", endpoint="component_role", primary_sources=["synthetic:test"],
                match=dict(provider="Tool", exact_native_models=["F__"+g], evidence_kinds=["final_member"], selected=True, database_versions=["1"]),
                requires=dict(protein_binding=True, source_binding=True), conclusion=dict(component=g, role="nuclease", substrate="DNA"))
    rule.update(changes)
    return rule


def registry(*rules):
    return dict(schema_version=1, registry_id="synthetic", version="1", rules=list(rules))


def gene(review, gid="a"):
    return next(r for r in review["gene_reviews"] if r["gene_id"] == gid)


def pack(tmp_path, data, *, broad=False, extra_other=False):
    source = tmp_path / "evidence.txt"
    source.write_text("synthetic qualified evidence\n")
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    sequence = "ADCEFGHIKLMN"
    model = dict(accession="TEST.1", source="SyntheticDB", model_sha256="1"*64, evidence_kind="nuclease_domain", role="nuclease", family_specific=not broad, min_model_coverage=.7)
    reference = dict(accession="REF", sequence=sequence, sequence_sha256=sha(sequence), numbering_origin="1-based full reference sequence",
                     independent_mapping_methods=["method_A", "method_B"], catalytic_sites=[dict(position=2, residue="D", role="nuclease", evidence_source="synthetic"), dict(position=4, residue="E", role="nuclease", evidence_source="synthetic")])
    entry = dict(family_id="FutureEvidence", component="effector", substrate="DNA", domain_models=[model], references=[reference])
    reg = dict(schema_version="1.0", families=[entry])
    hit = dict(accession="TEST.1", source="SyntheticDB", model_sha256="1"*64, query_sequence_sha256=sha(sequence),
               qstart=1, qend=12, hmm_start=1, hmm_end=12, hmm_length=12, target_length=12,
               qualification=dict(passed=True, rule="registered_qualification", source=str(source), source_sha256=source_sha))
    candidate = dict(candidate_id="candidate", sample_id="s", query_id="a", sequence=sequence, sequence_sha256=sha(sequence), domain_hits=[hit],
                     site_maps=[dict(reference_accession="REF", reference_sequence_sha256=sha(sequence), query_sequence_sha256=sha(sequence), method=m,
                                     positions={"2": 2, "4": 4}, source=str(source), source_sha256=source_sha) for m in ("method_A", "method_B")])
    if extra_other:
        other = dict(accession="OTHER.1", source="SyntheticDB", model_sha256="2"*64, evidence_kind="other_domain", role="motor")
        entry["domain_models"].append(other)
        candidate["domain_hits"].append({**deepcopy(hit), "accession":"OTHER.1", "model_sha256":"2"*64})
    reg_path, candidates_path = tmp_path / "registry.json", tmp_path / "candidates.json"
    reg_path.write_text(json.dumps(reg)); candidates_path.write_text(json.dumps([candidate]))
    assessments = assess_registry(reg, [candidate], registry_dir=tmp_path, source_root=tmp_path)
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, reg_path, candidates_path)}
    wrapper = dict(schema_version=1, registry=reg, candidates=[candidate], assessments=assessments, input_sha256=hashes,
                   registry_sha256=hashes[str(reg_path)], candidates_sha256=hashes[str(candidates_path)], source_files_verified=True,
                   registry_dir=str(tmp_path), source_root=str(tmp_path))
    data["system_review_evidence"] = wrapper
    return wrapper


def qualified_rule(route, endpoint, **match_extra):
    return dict(rule_id="qualified", version="1", route=route, endpoint=endpoint, primary_sources=["synthetic:test"],
                match={"family_ids":["FutureEvidence"], **match_extra}, requires=dict(protein_binding=True, source_binding=True),
                conclusion={"family_id":"FutureEvidence"} if endpoint == "family_annotation" else {"role":"nuclease", "substrate":"DNA"})


def test_all_calls_and_accepted_genes_without_partial_group_loss(tmp_path):
    data, result = saved(tmp_path, "a", "b")
    data["calls"][1]["system_id"] = "other"
    data["groups"] = [dict(sample_id="s", group_id="only_a", provider="Tool", run_id="run", system_id="sys", membership_kind="native_system", members=[{"gene_id":"a"}])]
    data["genes"].append(dict(sample_id="s", gene_id="retained", protein_sha256=sha("ADCEFGHIKLMN")))
    result["gene_ledger"].append(dict(sample_id="s", gene_id="retained", restriction_role=True))
    before = deepcopy((data, result))
    out = review_all_systems(data, result)
    assert (data, result) == before
    assert out["accounting"]["accounted_detection_count"] == 2
    assert out["accounting"]["accepted_genes_accounted"] == 3
    assert out["accounting"]["unit_counts"] == {"native_system":2, "review_only_gene":1}
    assert gene(out)["axes"]["native_detection"]["status"] == "supported"
    assert gene(out)["axes"]["component_role"]["status"] == "unassessed"
    assert out["review_execution"]["status"] == "completed"


def test_source_poor_legacy_input_completes_unassessed(tmp_path):
    data, result = saved(tmp_path, "a")
    data["calls"][0].pop("source_ref")
    data["genes"][0].pop("protein_sha256")
    out = review_all_systems(data, result, registry(component()))
    assert gene(out)["axes"]["component_role"]["status"] == "unassessed"
    reasons = out["rule_evaluations"][0]["reason_codes"]
    assert {"source_binding_unavailable", "protein_binding_unavailable"} <= set(reasons)


def test_raw_observation_rows_are_not_distinct_gene_counts(tmp_path):
    data, result = saved(tmp_path, "a")
    data["calls"][0]["native_observations"] = [{"row":2}, {"row":3}]
    out = review_all_systems(data, result)
    assert out["accounting"]["canonical_detection_count"] == 1
    assert out["accounting"]["raw_observation_count"] == 2
    assert out["accounting"]["distinct_reviewed_gene_count"] == 1


def test_qualified_registry_axes_do_not_promote_membership(tmp_path):
    data, result = saved(tmp_path, "a")
    pack(tmp_path, data)
    before = deepcopy((data, result))
    out = review_all_systems(data, result)
    assert (data, result) == before
    axes = gene(out)["axes"]
    assert axes["family_annotation"]["status"] == axes["nuclease_domain"]["status"] == axes["motif_site"]["status"] == "supported"
    assert axes["activity"]["status"] == "unassessed"
    assert gene(out)["accepted_membership"]["effector"] is False


def test_generic_domain_does_not_establish_family_or_reference_sites(tmp_path):
    data, result = saved(tmp_path, "a")
    pack(tmp_path, data, broad=True)
    axes = gene(review_all_systems(data, result))["axes"]
    assert axes["nuclease_domain"]["status"] == "supported"
    assert axes["family_annotation"]["status"] == axes["motif_site"]["status"] == "unassessed"


@pytest.mark.parametrize("change", ["assessment", "candidate", "source", "protein", "missing_assessment"])
def test_tampered_qualified_wrappers_fail_closed(tmp_path, change):
    data, result = saved(tmp_path, "a")
    wrapper = pack(tmp_path, data)
    if change == "assessment":
        wrapper["assessments"][0]["families"][0]["domain_evidence"]["reason"] = "forged"
    elif change == "candidate":
        wrapper["candidates"][0]["site_maps"][0]["positions"]["2"] = 3
    elif change == "source":
        (tmp_path / "evidence.txt").write_text("changed")
    elif change == "protein":
        data["genes"][0]["protein_sha256"] = "9"*64
    else:
        wrapper["assessments"] = []
    with pytest.raises(ValueError):
        review_all_systems(data, result)


def test_qualified_extension_exact_family_and_site_contract(tmp_path):
    data, result = saved(tmp_path, "a")
    pack(tmp_path, data)
    rule = qualified_rule("qualified_site", "motif_site", signature_accessions=["TEST.1"], namespace="SyntheticDB")
    rule["requires"]["qualified_sites"] = [dict(reference_sha256=sha("ADCEFGHIKLMN"), required_positions=[2,4])]
    out = review_all_systems(data, result, registry(rule))
    assert out["rule_evaluations"][0]["status"] == "supported"
    rule["match"]["family_ids"] = ["AnotherFamily"]
    assert review_all_systems(data, result, registry(rule))["rule_evaluations"] == []


@pytest.mark.parametrize("route,endpoint", [("qualified_site","motif_site"), ("qualified_domain","family_annotation"), ("qualified_domain","nuclease_domain")])
def test_other_domain_cannot_borrow_family_or_site_support(tmp_path, route, endpoint):
    data, result = saved(tmp_path, "a")
    pack(tmp_path, data, extra_other=True)
    rule = qualified_rule(route, endpoint, signature_accessions=["OTHER.1"], namespace="SyntheticDB")
    out = review_all_systems(data, result, registry(rule))
    assert out["rule_evaluations"][0]["status"] != "supported"


def test_new_rule_does_not_overwrite_existing_noncatalytic_veto(tmp_path):
    data, result = saved(tmp_path, "a")
    result["gene_ledger"][0]["annotations"] = [dict(status="supported", assignment_kind="noncatalytic_component", rule_id="reviewed_veto")]
    out = review_all_systems(data, result, registry(component()))
    axes = gene(out)["axes"]
    assert axes["component_role"]["status"] == "conflicting"
    assert len(axes["component_role"]["conflicting_evidence_ids"]) == 2
    assert gene(out)["accepted_membership"]["effector"] is False


def test_exact_digest_scalar_cannot_pass_by_substring(tmp_path):
    data, result = saved(tmp_path, "a")
    rule = component()
    rule["match"].pop("database_versions")
    rule["match"]["model_sha256"] = "1"*64
    data["calls"][0]["model_sha256"] = "1"
    with pytest.raises(ValueError, match="SHA-256"):
        review_all_systems(data, result, registry(rule))


@pytest.mark.parametrize("change", ["ignored_domain_on_native", "ignored_provider_on_qualified", "ignored_site_label"])
def test_unsupported_safeguards_are_rejected_not_ignored(change):
    if change == "ignored_domain_on_native":
        rule = component()
        rule["requires"]["qualified_domains"] = []
    elif change == "ignored_provider_on_qualified":
        rule = qualified_rule("qualified_domain", "nuclease_domain", provider="Tool")
    else:
        rule = qualified_rule("qualified_site", "motif_site")
        rule["requires"]["qualified_sites"] = [dict(reference_sha256="1"*64, required_positions=[2], site_set_id="unimplemented")]
    with pytest.raises(ValueError):
        validate_system_review_rules(registry(rule))


def test_scope_metadata_never_uses_external_flag(tmp_path):
    data, result = saved(tmp_path, "a")
    data["options"] = {"assessment_coverage":{"effector":"partial"}, "assessment_coverage_by_sample":{"s":{"effector":"complete"}}}
    out = review_all_systems(data, result)
    data["samples"][0]["cohort"] = "external"
    assert review_all_systems(data, result) == out
    assert gene(out)["axes"]["component_role"]["coverage"] == "complete"
    assert gene(out)["axes"]["component_role"]["status"] == "unassessed"


def test_manifest_relative_source_alias_requires_exact_path_and_digest(tmp_path):
    data, result = saved(tmp_path, "a")
    ref = data["calls"][0]["source_ref"]
    absolute = ref["path"]
    ref["path"] = "runs/arbitrary/native.tsv"
    data["suite"] = {"sources":[dict(path=ref["path"], sha256=ref["sha256"], resolved_path=absolute)]}
    assert gene(review_all_systems(data, result, registry(component())))["axes"]["component_role"]["status"] == "supported"
    ref["sha256"] = "9"*64
    with pytest.raises(ValueError, match="alias hash"):
        review_all_systems(data, result, registry(component()))


def test_unaccepted_saved_records_are_reviewed_without_becoming_members(tmp_path):
    data, result = saved(tmp_path, "a")
    data["saved_decisions"] = [dict(sample_id="s", gene_id="never_accepted", state="unassessed")]
    out = review_all_systems(data, result)
    row = gene(out, "never_accepted")
    assert row["in_integrated_universe"] is False
    assert row["stored_review_records"][0]["record"]["state"] == "unassessed"
    assert out["accounting"]["saved_review_gene_count"] == 1
    assert out["accounting"]["saved_review_genes_accounted"] == 1


def test_all_emitted_axis_evidence_ids_resolve(tmp_path):
    data, result = saved(tmp_path, "a")
    result["gene_ledger"][0].update(individual_effector=True, complex_core=True, motif_site_supported=True, restriction_role=True)
    out = review_all_systems(data, result, registry(component()))
    known = {r["evidence_id"] for r in out["evidence"]} | {r["evaluation_id"] for r in out["rule_evaluations"]} | {r["detection_id"] for r in out["detection_accounting"]}
    for row in out["units"] + out["gene_reviews"]:
        for axis in row["axes"].values():
            assert set(axis["supporting_evidence_ids"] + axis["conflicting_evidence_ids"]) <= known


def test_queue_does_not_require_nuclease_or_motif_for_known_motor(tmp_path):
    data, result = saved(tmp_path, "a")
    rule = component()
    rule["conclusion"]["role"] = "motor"
    out = review_all_systems(data, result, registry(rule))
    assert gene(out)["axes"]["motif_site"]["status"] == "unassessed"
    assert not any(q["axis"] in {"motif_site", "nuclease_domain", "complex_core"} for q in out["pending_queue"])
    assert not any(q["target_kind"] == "gene" and q["axis"] == "system_architecture" for q in out["pending_queue"])
    assert not any(q["target_kind"] == "unit" and q["axis"] == "component_role" for q in out["pending_queue"])


def test_queue_preserves_applicable_failed_site_claim(tmp_path):
    data, result = saved(tmp_path, "a")
    pack(tmp_path, data)
    rule = qualified_rule("qualified_site", "motif_site")
    rule["requires"]["qualified_sites"] = [dict(reference_sha256=sha("ADCEFGHIKLMN"), required_positions=[2,4,8])]
    out = review_all_systems(data, result, registry(rule))
    # Another valid route can support the general site axis while this narrower
    # explicit claim still needs additional evidence in its rule evaluation.
    assert out["rule_evaluations"][0]["status"] == "unassessed"
    pending = [q for q in out["pending_queue"] if q["target_kind"] == "rule_evaluation"]
    assert pending and all(q["required_for_claim"] for q in pending)
    assert pending[0]["queue_class"] == "applicable_evidence_contract_unmet"

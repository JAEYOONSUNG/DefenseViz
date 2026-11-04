"""Source identity and missing-axis boundaries for fresh stored-evidence evaluation."""
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json

import pytest

from defenseviz.analysis.evidence import evidence_dict, parse_evidence
from defenseviz.analysis.policy import Evidence, Policy, evaluate
from defenseviz.analysis.saved import analyze_saved
from defenseviz.rm_evaluation import apply_rm_evidence


def sha(value):
    return hashlib.sha256(value.encode()).hexdigest()


def evidence(**updates):
    value = dict(family="Type II", role="M", valid_sequence=True,
                 identity=None, query_coverage=None, reference_coverage=None,
                 evalue=None, role_compatible=False, family_compatible=None,
                 hmm_role_ambiguous=False, profiles=(), raw_nuclease_classes=frozenset(),
                 legacy_nuclease=False, legacy_motor=False, legacy_mtase_motif=False,
                 old_strong=True, old_moderate=False, old_profile=False,
                 annotation=False, old_motif=True, old_context=False,
                 strict_structure=False, alternative_annotation=False,
                 conflict=False, partial=False)
    value.update(updates)
    return evidence_dict(Evidence(**value))


def setup(tmp_path, count=1, *, sample="sample"):
    native = tmp_path / "native.tsv"
    native.write_text("gene\tstate\n" + "".join(f"gene{i}\tlegacy\n" for i in range(count)))
    source = tmp_path / "stored_evidence.json"
    source.write_text("{\"description\":\"synthetic complete measured evidence fixture\"}")
    source_ref = dict(path=source.name, sha256=sha(source.read_text()), record_index=0)
    native_ref = dict(path=str(native), sha256=sha(native.read_text()), row=2)
    genes = [dict(sample_id=sample, gene_id=f"gene{i}", protein_sha256=sha(f"protein{i}")) for i in range(count)]
    d = dict(samples=[dict(sample_id=sample, label="Synthetic sample")], genes=genes,
             calls=[dict(sample_id=sample, gene_id=g["gene_id"], caller="OriginalTool", family="Type II",
                         rm_type="Type II", category="R-M", evidence_kind="final_member", source_ref=native_ref)
                    for g in genes],
             coverage=[dict(sample_id=sample, caller="OriginalTool", status="complete", endpoint="native_calls")],
             saved_decisions=[dict(sample_id=sample, gene_id=g["gene_id"], state="unassessed",
                                   family="Type II", role="M", native_context="legacy_context", source_ref=native_ref)
                              for g in genes],
             annotations=[], effector_evidence=[], groups=[], complexes=[],
             inputs=[native], input_sha256={str(native):native_ref["sha256"]},
             input_roles={str(native):[dict(sample_id=sample, native_result=True, kind="saved_candidates", module="rebasefinder")]})
    doc = dict(schema_version=1, input_sha256={native.name:native_ref["sha256"]},
               selected_policy=asdict(Policy()),
               records=[dict(**g, evidence=evidence(), source_refs=[source_ref]) for g in genes],
               candidate_scope=[dict(sample_id=sample, status="complete",
                                     gene_ids=[g["gene_id"] for g in genes], source_refs=[source_ref])])
    path = tmp_path / "rm_evidence.json"
    return d, doc, path


def apply(d, doc, path):
    path.write_text(json.dumps(doc))
    return apply_rm_evidence(d, path)


def test_evaluates_exact_frozen_predicate_and_preserves_native_data(tmp_path):
    d, doc, path = setup(tmp_path, 2)
    doc["records"][1]["evidence"] = evidence(old_strong=False, old_motif=False)
    before = deepcopy(d)
    out = apply(d, doc, path)
    assert d == before
    assert out is not d
    for field in ("calls", "groups", "genes", "samples", "annotations", "effector_evidence"):
        assert out[field] == d[field]
    assert [r["state"] for r in out["saved_decisions"]] == ["retained", "not_retained"]
    for raw, measured in zip(doc["records"], out["rm_evaluation"]["records"]):
        expected = json.loads(json.dumps(evaluate(parse_evidence(raw["evidence"]), Policy())))
        assert measured["decision"] == expected
        assert measured["previous_saved_decision"]["state"] == "unassessed"
    meta = out["rm_evaluation"]
    assert meta["operation"] == "fresh_evaluation_of_stored_RM_evidence"
    assert meta["scoring_runs"] == 1 and meta["policy_evaluations"] == 2
    assert meta["detector_runs"] == meta["structure_prediction_runs"] == 0
    assert meta["thresholds_changed"] is False and meta["activity_inferred"] is False
    assert out["input_sha256"][str(path)] == sha(path.read_text())
    assert analyze_saved(out)["global_summary"]["retained_CDS"] == 1


def test_140_legacy_candidates_23_evaluated_are_partial_not_na_or_complete(tmp_path):
    d, doc, path = setup(tmp_path, 140)
    doc["records"] = doc["records"][:23]
    for r in doc["records"][11:]:
        r["evidence"] = evidence(old_strong=False, old_motif=False)
    doc["candidate_scope"][0].update(status="partial", gene_ids=[r["gene_id"] for r in doc["records"]])
    out = apply(d, doc, path)
    assert len(out["saved_decisions"]) == 140
    assert sum(r["state"] == "unassessed" for r in out["saved_decisions"]) == 117
    scope = out["rm_evaluation"]["candidate_scope"][0]
    assert scope["evaluated_count"] == 23 and scope["active_candidate_count"] == 140
    assert len(scope["outside_declared_scope_gene_ids"]) == 117
    result = analyze_saved(out)
    summary = result["sample_summary"][0]
    assert summary["retained_CDS"] == 11 and summary["retained_counts_status"] == "partial"
    assert result["availability"]["saved_decisions_by_sample"][0]["status"] == "partial"
    assert summary["source_counts"]["OriginalTool"] == 140


@pytest.mark.parametrize("field", ["old_context", "strict_structure", "divergent_context", "curated_role_compatible", "reference_nuclease_anchor_support"])
def test_absent_axes_never_default_to_false(tmp_path, field):
    d, doc, path = setup(tmp_path)
    del doc["records"][0]["evidence"][field]
    with pytest.raises(ValueError, match="full Evidence"):
        apply(d, doc, path)


@pytest.mark.parametrize("field,value", [("old_context", None), ("strict_structure", "false"), ("identity", 1.2)])
def test_missing_or_malformed_assessments_rejected(tmp_path, field, value):
    d, doc, path = setup(tmp_path)
    doc["records"][0]["evidence"][field] = value
    with pytest.raises(ValueError):
        apply(d, doc, path)


def test_changed_threshold_or_incomplete_policy_rejected(tmp_path):
    d, doc, path = setup(tmp_path)
    doc["selected_policy"]["r_identity_min"] = .7
    with pytest.raises(ValueError, match="unchanged selected policy"):
        apply(d, doc, path)
    del doc["selected_policy"]["r_identity_min"]
    with pytest.raises(ValueError, match="full policy"):
        apply(d, doc, path)


@pytest.mark.parametrize("mutation,message", [
    ("protein", "protein hash"), ("source", "source hash"),
    ("native", "native artifact binding"), ("native_changed", "native artifact changed"),
    ("no_native_role", "native result artifact"), ("sample", "known sample_id"),
    ("gene", "candidate_scope"), ("duplicate", "duplicate evaluated"),
    ("old_duplicate", "Duplicate saved decision"), ("unknown", "unknown fields")])
def test_identity_and_source_boundary_rejections(tmp_path, mutation, message):
    d, doc, path = setup(tmp_path)
    if mutation == "protein": doc["records"][0]["protein_sha256"] = "a"*64
    elif mutation == "source": doc["records"][0]["source_refs"][0]["sha256"] = "a"*64
    elif mutation == "native": doc["input_sha256"]["native.tsv"] = "a"*64
    elif mutation == "native_changed": (tmp_path/"native.tsv").write_text("changed")
    elif mutation == "no_native_role": d["input_roles"] = {}
    elif mutation == "sample": doc["records"][0]["sample_id"] = "unknown"
    elif mutation == "gene": doc["records"][0]["gene_id"] = "unknown"
    elif mutation == "duplicate": doc["records"].append(deepcopy(doc["records"][0]))
    elif mutation == "old_duplicate": d["saved_decisions"].append(deepcopy(d["saved_decisions"][0]))
    elif mutation == "unknown": doc["records"][0]["activity"] = True
    with pytest.raises(ValueError, match=message):
        apply(d, doc, path)


def test_each_sample_needs_own_imported_native_binding(tmp_path):
    d, doc, path = setup(tmp_path)
    d["samples"].append(dict(sample_id="external"))
    d["genes"].append(dict(sample_id="external", gene_id="x", protein_sha256=sha("x")))
    doc["candidate_scope"].append(dict(sample_id="external", status="complete", gene_ids=["x"], source_refs=doc["candidate_scope"][0]["source_refs"]))
    with pytest.raises(ValueError, match="native-bound"):
        apply(d, doc, path)


def test_complete_scope_cannot_hide_omitted_or_invalid_legacy_candidate(tmp_path):
    d, doc, path = setup(tmp_path, 2)
    doc["records"].pop()
    with pytest.raises(ValueError, match="unassessed native candidates"):
        apply(d, doc, path)
    doc["candidate_scope"][0]["gene_ids"].pop()
    with pytest.raises(ValueError, match="omits"):
        apply(d, doc, path)


def test_genuine_complete_zero_is_distinct_from_missing(tmp_path):
    d, doc, path = setup(tmp_path, 0)
    assert analyze_saved(d)["availability"]["saved_decisions"] == "unassessed"
    out = apply(d, doc, path)
    result = analyze_saved(out)
    assert result["availability"]["saved_decisions"] == "complete"
    assert result["sample_summary"][0]["retained_CDS"] == 0
    assert out["rm_evaluation"]["scoring_runs"] == 0


@pytest.mark.parametrize("prior", [dict(state="not_retained"), dict(state="retained", role="R"),
                                   dict(state="retained", family="Type I"), dict(state="retained", family=None, rm_type="Type I")])
def test_assessed_conflicts_require_explicit_current_source_selection(tmp_path, prior):
    d, doc, path = setup(tmp_path)
    d["saved_decisions"][0].update(prior)
    with pytest.raises(ValueError, match="conflicting assessed"):
        apply(d, doc, path)


def test_same_assessed_state_preserved_in_provenance_and_reapplication_rejected(tmp_path):
    d, doc, path = setup(tmp_path)
    d["saved_decisions"][0]["state"] = "retained"
    out = apply(d, doc, path)
    assert out["rm_evaluation"]["native_saved_decisions"] == d["saved_decisions"]
    with pytest.raises(ValueError, match="already applied"):
        apply_rm_evidence(out, path)


def test_unassessed_state_cannot_mask_stale_saved_protein_binding(tmp_path):
    d, doc, path = setup(tmp_path)
    d["saved_decisions"][0]["protein_sha256"] = "a"*64
    with pytest.raises(ValueError, match="protein hash contradicts"):
        apply(d, doc, path)


@pytest.mark.parametrize("update,state", [(dict(role="Aux"),"out_of_scope"), (dict(role="helicase"),"out_of_scope"),
                                          (dict(valid_sequence=False),"unassessed")])
def test_out_of_scope_and_invalid_sequences_never_scored(tmp_path, update, state):
    d, doc, path = setup(tmp_path)
    doc["candidate_scope"][0]["status"] = "partial"
    doc["records"][0]["evidence"] = evidence(**update)
    out = apply(d, doc, path)
    assert out["saved_decisions"][0]["state"] == state
    assert out["rm_evaluation"]["policy_evaluations"] == 0
    assert out["rm_evaluation"]["scoring_runs"] == 0
    assert out["rm_evaluation"]["records"][0]["decision"] is None
    assert analyze_saved(out)["global_summary"]["retained_CDS"] == 0


def test_old_context_is_not_current_completeness_and_groups_untouched(tmp_path):
    d, doc, path = setup(tmp_path)
    d["groups"] = [dict(group_id="raw-context", sample_id="sample", provider="DNMB", membership_kind="saved_context",
                        members=[dict(gene_id="gene0", retained=True)], status="complete")]
    d["saved_decisions"][0].update(context_state="complete", orphan=False)
    doc["records"][0]["evidence"]["old_context"] = True
    out = apply(d, doc, path)
    assert out["groups"] == d["groups"]
    assert "context_state" not in out["saved_decisions"][0]
    assert "orphan" not in out["saved_decisions"][0]
    assert out["saved_decisions"][0]["native_context"] == "legacy_context"


def test_explicit_current_context_requires_separate_hashed_evidence(tmp_path):
    d, doc, path = setup(tmp_path)
    record = doc["records"][0]
    record.update(context_state="missing_retained_partner_in_saved_context", orphan=True)
    with pytest.raises(ValueError, match="source_refs"):
        apply(d, doc, path)
    record["context_source_refs"] = deepcopy(record["source_refs"])
    out = apply(d, doc, path)
    assert out["saved_decisions"][0]["orphan"] is True
    assert out["saved_decisions"][0]["context_state"] == record["context_state"]
    record["context_source_refs"][0]["sha256"] = "b"*64
    with pytest.raises(ValueError, match="source hash"):
        apply(d, doc, path)


def test_identifiers_and_record_order_do_not_change_policy_outputs(tmp_path):
    d, doc, path = setup(tmp_path, 3)
    doc["records"][1]["evidence"] = evidence(old_strong=False)
    original = apply(d, doc, path)
    for table in (d["genes"], d["calls"], d["saved_decisions"], doc["records"]):
        for row in table:
            row["gene_id"] = "opaque::" + row["gene_id"]
    doc["candidate_scope"][0]["gene_ids"] = ["opaque::" + x for x in doc["candidate_scope"][0]["gene_ids"]]
    doc["records"].reverse()
    permuted = apply(d, doc, path)
    a = {r["gene_id"]: r["decision"] for r in original["rm_evaluation"]["records"]}
    b = {r["gene_id"].removeprefix("opaque::"): r["decision"] for r in permuted["rm_evaluation"]["records"]}
    assert a == b


@pytest.mark.parametrize("kind", ["unmapped", "unbound_native_candidate", "cached_live_candidate_count_mismatch",
                                  "sentinel_live_retained_count_mismatch", "sentinel_live_row_count_mismatch"])
def test_complete_scope_cannot_override_native_candidate_integrity_conflicts(tmp_path, kind):
    d, doc, path = setup(tmp_path)
    issue = dict(sample_id="sample", module="rebasefinder", source_ref=d["saved_decisions"][0]["source_ref"])
    if kind == "unmapped":
        d["suite"] = dict(unmapped_records=[dict(**issue, assessment_status="unassessed_unbound_identity")])
    else:
        d["suite"] = dict(conflicts=[dict(**issue, code=kind, observed=1, expected=2)])
    with pytest.raises(ValueError, match="native candidate identity/count integrity"):
        apply(d, doc, path)
    doc["candidate_scope"][0]["status"] = "partial"
    out = apply(d, doc, path)
    assert out["rm_evaluation"]["candidate_scope"][0]["native_candidate_integrity_issues"]
    assert analyze_saved(out)["sample_summary"][0]["retained_counts_status"] == "partial"


def test_unrelated_sample_or_annotation_gap_does_not_block_candidate_scope(tmp_path):
    d, doc, path = setup(tmp_path)
    d["suite"] = dict(conflicts=[dict(sample_id="other", code="unbound_native_candidate")],
                      unmapped_records=[dict(sample_id="sample", module="eggnog", source_ref=dict(path="orthology.tsv")),
                                        dict(sample_id="sample", module="rebasefinder", source_ref=dict(path="structure.tsv"))])
    out = apply(d, doc, path)
    assert out["rm_evaluation"]["candidate_scope"][0]["coverage_status"] == "complete"
    assert out["rm_evaluation"]["candidate_scope"][0]["native_candidate_integrity_issues"] == []


def test_complete_scope_cannot_add_non_native_candidates(tmp_path):
    d, doc, path = setup(tmp_path, 2)
    d["saved_decisions"].pop()
    with pytest.raises(ValueError, match="outside the original native candidate"):
        apply(d, doc, path)
    doc["candidate_scope"][0]["status"] = "partial"
    out = apply(d, doc, path)
    assert out["rm_evaluation"]["candidate_scope"][0]["outside_native_candidate_gene_ids"] == ["gene1"]


def test_optional_expected_decision_checks_exact_historical_parity(tmp_path):
    d, doc, path = setup(tmp_path)
    doc["records"][0]["expected_decision"] = dict(keep=True, score_total=55, source_refs=deepcopy(doc["records"][0]["source_refs"]))
    out = apply(d, doc, path)
    record = out["rm_evaluation"]["records"][0]
    assert record["parity_status"] == "verified"
    assert record["expected_decision"]["source_refs"][0]["path"] == str(tmp_path/"stored_evidence.json")
    assert out["rm_evaluation"]["expected_decisions_verified"] == 1
    doc["records"][0]["expected_decision"]["score_total"] = 54
    with pytest.raises(ValueError, match="parity mismatch"):
        apply(d, doc, path)
    doc["records"][0]["expected_decision"].update(score_total=55, keep=False)
    with pytest.raises(ValueError, match="parity mismatch"):
        apply(d, doc, path)


@pytest.mark.parametrize("mutation", ["keep", "score", "source", "unscored"])
def test_expected_decision_cannot_assert_unverified_parity(tmp_path, mutation):
    d, doc, path = setup(tmp_path)
    record = doc["records"][0]
    record["expected_decision"] = dict(keep=True, score_total=55, source_refs=deepcopy(record["source_refs"]))
    if mutation == "keep": record["expected_decision"]["keep"] = "true"
    elif mutation == "score": record["expected_decision"]["score_total"] = True
    elif mutation == "source": record["expected_decision"]["source_refs"][0]["sha256"] = "a"*64
    else: record["evidence"]["valid_sequence"] = False
    with pytest.raises(ValueError):
        apply(d, doc, path)


def test_native_retention_coverage_and_row_references_are_preserved(tmp_path):
    d, doc, path = setup(tmp_path)
    native_coverage = dict(sample_id="sample", caller="DNMB", endpoint="saved_retention", status="partial",
                           source_ref=deepcopy(d["saved_decisions"][0]["source_ref"]))
    d["coverage"].append(native_coverage)
    before = deepcopy(d)
    out = apply(d, doc, path)
    assert d == before
    assert out["rm_evaluation"]["native_saved_retention_coverage"] == [native_coverage]
    assert out["saved_decisions"][0]["native_source_ref"] == d["saved_decisions"][0]["source_ref"]
    assert out["saved_decisions"][0]["source_ref"]["path"] == str(path)
    assert [r["status"] for r in out["coverage"] if r.get("endpoint") == "saved_retention"] == ["complete"]


def test_identical_hash_paralogs_remain_distinct_exact_gene_identities(tmp_path):
    d, doc, path = setup(tmp_path, 2)
    d["genes"][1]["protein_sha256"] = d["genes"][0]["protein_sha256"]
    doc["records"][1]["protein_sha256"] = doc["records"][0]["protein_sha256"]
    doc["records"][1]["evidence"] = evidence(old_strong=False, old_motif=False)
    d["saved_decisions"][1]["source_ref"] = {**d["saved_decisions"][1]["source_ref"], "row":3}
    out = apply(d, doc, path)
    rows = {r["gene_id"]:r for r in out["saved_decisions"]}
    assert rows["gene0"]["state"] == "retained" and rows["gene1"]["state"] == "not_retained"
    assert rows["gene0"]["native_source_ref"]["row"] == 2
    assert rows["gene1"]["native_source_ref"]["row"] == 3
    assert len(out["rm_evaluation"]["records"]) == 2


def test_retained_R_role_is_separate_from_policy_R_component_support(tmp_path):
    d, doc, path = setup(tmp_path)
    doc["records"][0]["evidence"] = evidence(role="R")
    out = apply(d, doc, path)
    row = out["saved_decisions"][0]
    assert row["state"] == "retained" and row["role"] == "R"
    assert row["integrated_R_component_support"] is False
    assert row["score_total"] == 55


def test_explicit_replace_scoped_resets_omitted_old_retention_and_preserves_originals(tmp_path):
    d, doc, path = setup(tmp_path, 3)
    for row in d["saved_decisions"]:
        row.update(state="retained", score_total=100, context_state="complete", orphan=False)
    doc["decision_mode"] = "replace_scoped"
    doc["records"] = doc["records"][:2]
    doc["records"][0]["evidence"] = evidence(old_strong=False, old_motif=False)
    doc["records"][0]["expected_decision"] = dict(keep=False, score_total=0, source_refs=deepcopy(doc["records"][0]["source_refs"]))
    doc["candidate_scope"][0].update(status="partial", gene_ids=["gene0", "gene1"])
    before = deepcopy(d)
    out = apply(d, doc, path)
    assert d == before
    assert out["rm_evaluation"]["native_saved_decisions"] == before["saved_decisions"]
    assert [r["state"] for r in out["saved_decisions"]] == ["not_retained", "retained", "unassessed"]
    assert all("context_state" not in r and "orphan" not in r for r in out["saved_decisions"])
    assert out["saved_decisions"][2]["score_total"] is None
    transitions = out["rm_evaluation"]["decision_transitions"]
    assert all(r["previous_state"] == "retained" for r in transitions)
    assert transitions[2]["current_state"] == "unassessed" and transitions[2]["evaluated"] is False
    assert analyze_saved(out)["sample_summary"][0]["retained_CDS"] == 1


def test_replace_scoped_does_not_change_other_samples_or_license_complete_omissions(tmp_path):
    d, doc, path = setup(tmp_path, 2)
    d["samples"].append(dict(sample_id="outside", label="Untouched"))
    d["genes"].append(dict(sample_id="outside", gene_id="gene0", protein_sha256=sha("outside")))
    other = dict(sample_id="outside", gene_id="gene0", family="Type II", role="M", state="retained")
    d["saved_decisions"].append(other)
    doc["decision_mode"] = "replace_scoped"
    out = apply(d, doc, path)
    assert out["saved_decisions"][-1] == other
    doc["records"].pop()
    with pytest.raises(ValueError, match="unassessed native candidates"):
        apply(d, doc, path)


def test_unknown_decision_mode_is_not_an_implicit_override(tmp_path):
    d, doc, path = setup(tmp_path)
    doc["decision_mode"] = "override"
    with pytest.raises(ValueError, match="decision_mode"):
        apply(d, doc, path)

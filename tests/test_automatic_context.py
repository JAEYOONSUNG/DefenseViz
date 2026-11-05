"""Current retained context uses exact original boundaries and known dispositions."""
from copy import deepcopy
import hashlib

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.automatic_context import apply_current_context_review


def fixture(tmp_path, typ="Type I", roles=("R", "M", "S")):
    source = tmp_path / "contexts.tsv"
    source.write_text("synthetic exact source boundary and evaluated records\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    ref = dict(path=str(source), sha256=digest, row=2)
    genes = [dict(sample_id="sample", gene_id=f"g{i}", protein_sha256=hashlib.sha256(f"protein{i}".encode()).hexdigest())
             for i in range(len(roles))]
    decisions = [dict(**gene, state="retained", family=typ, rm_type=typ, role=role,
                      source_refs=[ref], integrated_R_component_support=role in {"R", "RM"})
                 for gene, role in zip(genes, roles)]
    group = dict(group_id="original", sample_id="sample", provider="DNMB", run_id=digest,
                 membership_kind="saved_context", family=typ, context_state="complete_mrs",
                 members=[dict(gene_id=g["gene_id"], protein_sha256=g["protein_sha256"]) for g in genes], source_refs=[ref])
    return dict(samples=[dict(sample_id="sample", label="Synthetic")], genes=genes,
                calls=[], saved_decisions=decisions, groups=[group],
                coverage=[dict(sample_id="sample", caller="DNMB", endpoint="saved_retention", status="complete")],
                input_sha256={str(source):digest}, rm_evaluation=dict(
                    evaluation_id="evaluation", decision_mode="evaluate_unassessed", records=deepcopy(decisions),
                    native_saved_decisions=[dict(**g, state="unassessed") for g in genes], candidate_scope=[]))


def generated(out):
    return [g for g in out["groups"] if g.get("origin") == "automatic_current_context_review"]


def test_complete_current_type_i_is_composed_without_policy_or_native_changes(tmp_path):
    data = fixture(tmp_path)
    before = deepcopy(data)
    out = apply_current_context_review(data)
    assert data == before
    assert out["groups"][:1] == data["groups"]
    assert len(generated(out)) == 1
    group = generated(out)[0]
    assert group["provider"] == "Review" and group["required_roles"] == ["M", "R", "S"]
    assert group["context_state"] == "met"
    for old, new in zip(data["saved_decisions"], out["saved_decisions"]):
        assert new["state"] == old["state"]
        assert new["integrated_R_component_support"] == old["integrated_R_component_support"]
        assert new["context_state"] == "met"
        assert new["context_source_refs"]
        assert "activity" not in new
    assert out["saved_decisions"][1]["orphan"] is False
    assert out["rm_evaluation"]["current_context_review"]["activity_inferred"] is False
    result = analyze_saved(out)
    current = [r for r in result["context_assessments"] if r["provider"] == "Review"]
    assert len(current) == 1 and current[0]["context_state"] == "met"
    assert result["global_summary"]["retained_CDS"] == 3


def test_filter_removes_partner_then_context_is_missing_not_complete(tmp_path):
    data = fixture(tmp_path)
    data["saved_decisions"][1]["state"] = "not_retained"
    out = apply_current_context_review(data)
    group = generated(out)[0]
    assert group["context_state"] == "missing_retained_partner_in_saved_context"
    assert {m["gene_id"] for m in group["members"]} == {"g0", "g2"}
    assert len(out["rm_evaluation"]["current_context_review"]["records"][0]["original_dispositions"]) == 3
    result = analyze_saved(out)
    current = next(r for r in result["context_assessments"] if r["provider"] == "Review")
    assert current["context_state"] == "missing_retained_partner_in_saved_context"


def test_unassessed_partner_prevents_negative_context_and_orphan(tmp_path):
    data = fixture(tmp_path)
    data["saved_decisions"][0]["state"] = "unassessed"
    out = apply_current_context_review(data)
    assert generated(out)[0]["context_state"] == "unassessed"
    assert all("orphan" not in d and "context_state" not in d for d in out["saved_decisions"])


def test_filtering_known_original_R_marks_only_scoped_retained_orphan(tmp_path):
    data = fixture(tmp_path)
    data["saved_decisions"][0]["state"] = "not_retained"
    out = apply_current_context_review(data)
    for decision in out["saved_decisions"][1:]:
        assert decision["orphan"] is True
        assert decision["orphan_scope"] == "within_exact_retained_saved_context_not_genome_wide_absence"
    assert "orphan" not in out["saved_decisions"][0]


def test_missing_boundary_does_not_fabricate_singletons_or_orphans(tmp_path):
    data = fixture(tmp_path)
    data["groups"] = []
    out = apply_current_context_review(data)
    assert generated(out) == []
    assert all("orphan" not in d and "context_state" not in d for d in out["saved_decisions"])


@pytest.mark.parametrize("mutation", ["run", "source", "protein", "decision_source"])
def test_unbound_boundary_or_partner_evidence_stays_unknown(tmp_path, mutation):
    data = fixture(tmp_path)
    if mutation == "run": data["groups"][0].pop("run_id")
    elif mutation == "source": data["input_sha256"] = {}
    elif mutation == "protein": data["groups"][0]["members"][1]["protein_sha256"] = "f"*64
    else: data["saved_decisions"][1]["source_refs"] = []
    out = apply_current_context_review(data)
    if mutation == "decision_source":
        assert generated(out)[0]["context_state"] == "unassessed"
    else:
        assert generated(out) == []
    assert all("orphan" not in d for d in out["saved_decisions"])


def test_singleton_M_is_not_claimed_orphan_by_unobserved_global_R_absence(tmp_path):
    data = fixture(tmp_path, typ="Type III", roles=("M",))
    out = apply_current_context_review(data)
    assert generated(out)[0]["context_state"] == "missing_retained_partner_in_saved_context"
    assert "orphan" not in out["saved_decisions"][0]


@pytest.mark.parametrize("typ,roles", [("Type III", ("R", "M")), ("Type I", ("RM", "S"))])
def test_type_specific_contract_and_explicit_fusion_roles(tmp_path, typ, roles):
    out = apply_current_context_review(fixture(tmp_path, typ, roles))
    assert generated(out)[0]["context_state"] == "met"


def test_type_ii_requires_explicit_contract(tmp_path):
    data = fixture(tmp_path, typ="Type II", roles=("R", "M"))
    assert generated(apply_current_context_review(data))[0]["context_state"] == "unassessed"
    data["groups"][0]["required_roles"] = ["R", "M"]
    assert generated(apply_current_context_review(data))[0]["context_state"] == "met"


def test_type_iv_never_infers_methylation_requirement(tmp_path):
    data = fixture(tmp_path, typ="Type IV", roles=("R",))
    out = apply_current_context_review(data)
    assert generated(out)[0]["context_state"] == "unassessed"
    assert "required_roles" not in generated(out)[0]
    data["groups"][0]["required_roles"] = ["R"]
    out = apply_current_context_review(data)
    assert generated(out)[0]["context_state"] == "met"
    assert generated(out)[0]["required_roles"] == ["R"]


def test_reciprocal_rows_same_boundary_are_not_double_counted(tmp_path):
    data = fixture(tmp_path)
    reciprocal = deepcopy(data["groups"][0]); reciprocal["group_id"] = "reciprocal"
    reciprocal["members"].reverse()
    reciprocal["source_refs"][0]["row"] = 3
    data["groups"].append(reciprocal)
    out = apply_current_context_review(data)
    assert len(generated(out)) == 1
    assert generated(out)[0]["original_group_ids"] == ["original", "reciprocal"]
    assert len(generated(out)[0]["source_refs"]) == 2
    assert out["groups"][:2] == data["groups"]


def test_explicit_existing_context_is_not_overwritten(tmp_path):
    data = fixture(tmp_path)
    data["saved_decisions"][1].update(context_state="met", orphan=False, context_source_refs=[{"proof": "explicit"}])
    out = apply_current_context_review(data)
    assert out["saved_decisions"][1] == data["saved_decisions"][1]


def test_idempotent_but_stale_decisions_cannot_reuse_review(tmp_path):
    data = fixture(tmp_path)
    once = apply_current_context_review(data)
    assert apply_current_context_review(once) == once
    once["saved_decisions"][1]["state"] = "not_retained"
    with pytest.raises(ValueError, match="stale"):
        apply_current_context_review(once)


def test_unchanged_previously_assessed_boundary_does_not_duplicate_current_review(tmp_path):
    data = fixture(tmp_path)
    data["rm_evaluation"]["records"] = []
    out = apply_current_context_review(data)
    assert generated(out) == []


def test_conflicting_overlapping_exact_boundaries_stay_unknown(tmp_path):
    data = fixture(tmp_path)
    alternative = deepcopy(data["groups"][0])
    alternative["group_id"] = "alternative_without_S"
    alternative["members"] = alternative["members"][:2]
    data["groups"].append(alternative)
    out = apply_current_context_review(data)
    assert len(generated(out)) == 2
    assert {g["context_state"] for g in generated(out)} == {"unassessed"}
    records = out["rm_evaluation"]["current_context_review"]["records"]
    assert {r["evaluated_context_state"] for r in records} == {"met", "missing_retained_partner_in_saved_context"}
    assert all("competing_exact_saved_contexts" in r["reason_codes"] for r in records)
    analyze_saved(out)  # Competing boundaries cannot inject contradictory query claims.


def test_unknown_or_conflicting_partner_type_does_not_supply_missing_role(tmp_path):
    data = fixture(tmp_path)
    data["saved_decisions"][1].update(family="Type II", rm_type="Type II")
    out = apply_current_context_review(data)
    assert generated(out)[0]["context_state"] == "unassessed"
    assert "original_boundary_role_or_type_unresolved" in out["rm_evaluation"]["current_context_review"]["records"][0]["reason_codes"]


def test_explicit_registered_derived_evidence_can_bind_current_context(tmp_path):
    data = fixture(tmp_path)
    artifact = tmp_path / "fresh_alignment.tsv"
    artifact.write_text("fresh independently recorded alignment\n")
    digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
    data["derived_artifact_sha256"] = {str(artifact): digest}
    data["saved_decisions"][1]["source_refs"].append(dict(path=str(artifact), sha256=digest))
    out = apply_current_context_review(data)
    assert generated(out)[0]["context_state"] == "met"
    assert str(artifact) not in data["input_sha256"]
    artifact.write_text("modified after context review\n")
    with pytest.raises(ValueError, match="artifact hash mismatch"):
        apply_current_context_review(out)


def test_wrong_derived_artifact_digest_fails_closed(tmp_path):
    data = fixture(tmp_path)
    artifact = tmp_path / "fresh_alignment.tsv"
    artifact.write_text("changed or untrusted derived evidence\n")
    data["derived_artifact_sha256"] = {str(artifact): "f" * 64}
    with pytest.raises(ValueError, match="artifact hash mismatch"):
        apply_current_context_review(data)


def loaded_context(tmp_path):
    from defenseviz.io.suite import load_suite
    from tests.io.test_suite import native, write_table
    root = native(tmp_path)
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "final_family", "final_role", "curation_keep", "typeiii_context_status", "typeiii_context_partners"],
                [["g:a", "Type III", "R", "NA", "complete_mr", "g_b:M:annotation:+1g:8bp"],
                 ["g_b", "Type III", "M", "NA", "complete_mr", "g:a:R:annotation:-1g:8bp"]])
    data = load_suite(root)
    original = deepcopy(data["saved_decisions"])
    for row in data["saved_decisions"]:
        row.update(state="retained", family="Type III", rm_type="Type III",
                   role="R" if row["gene_id"] == "g:a" else "M")
        row.pop("context_state", None)
    data["rm_evaluation"] = dict(evaluation_id="current", records=deepcopy(data["saved_decisions"]),
                                 native_saved_decisions=original, candidate_scope=[])
    return data


def test_actual_suite_relative_source_refs_bind_via_exact_source_registry(tmp_path):
    from pathlib import Path
    data = loaded_context(tmp_path)
    original_group = next(g for g in data["groups"] if g["provider"] == "DNMB")
    assert not Path(original_group["source_refs"][0]["path"]).is_absolute()
    assert all(Path(p).is_absolute() for p in data["input_sha256"])
    out = apply_current_context_review(data)
    assert len(generated(out)) == 1
    assert generated(out)[0]["context_state"] == "met"
    assert {m["gene_id"] for m in generated(out)[0]["members"]} == {"g:a", "g_b"}
    assert any(r["provider"] == "Review" and r["context_state"] == "met" for r in analyze_saved(out)["context_assessments"])


def test_suite_source_basename_cannot_replace_exact_resolved_binding(tmp_path):
    data = loaded_context(tmp_path)
    for source in data["suite"]["sources"]:
        if source["module"] == "rebasefinder":
            source["resolved_path"] = "/different-root/" + source["resolved_path"].split("/")[-1]
    out = apply_current_context_review(data)
    assert generated(out) == []

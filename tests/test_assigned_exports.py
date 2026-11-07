"""Assigned components, policy retention, and effector evidence stay distinct."""
from copy import deepcopy
import csv
import hashlib

from openpyxl import load_workbook
import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.plot.publication.manuscript_model import enrich_manuscript_model
from defenseviz.suite_exports import export_suite_workbook


REF = {"path": "review.tsv", "sha256": "a" * 64, "row": 2}


def bundle():
    genes = [{"sample_id": "s", "gene_id": gid,
              "protein_sha256": hashlib.sha256(gid.encode()).hexdigest()}
             for gid in ["methylase", "restriction", "motor", "unresolved"]]
    decisions = [dict(sample_id="s", gene_id=g["gene_id"],
                      family="Type II" if g["gene_id"] == "methylase" else "Type IV",
                      role="M" if g["gene_id"] == "methylase" else "R",
                      state="retained" if g["gene_id"] == "methylase" else "not_retained",
                      score_total=55 if g["gene_id"] == "methylase" else 35, source_ref=REF)
                 for g in genes]
    annotations = [dict(**g, family="Type IV", role=role, status="supported",
                        assignment_kind="reviewed_RM_component", include_in_union=True,
                        rule_id="synthetic_bound_architecture", source_refs=[REF])
                   for g, role in [(genes[1], "R"), (genes[2], "Aux")]]
    return dict(samples=[dict(sample_id="s", label="Synthetic")], genes=genes,
                saved_decisions=decisions, annotations=annotations,
                calls=[dict(sample_id="s", gene_id="restriction", caller="DefenseFinder",
                            family="RM_Type_IV", rm_type="Type IV", role="R", category="R-M",
                            evidence_kind="final_member", source_ref=REF)],
                coverage=[dict(sample_id="s", caller="DefenseFinder", status="complete")])


def read_tsv(path):
    with path.open(newline="") as stream:
        return {r["gene_id"]: r for r in csv.DictReader(stream, delimiter="\t")}


def test_reviewed_typeiv_reaches_assigned_export_without_changing_policy_or_effector(tmp_path):
    data = bundle()
    original = deepcopy(data)
    result = analyze_saved(data)
    assert data == original
    summary = result["sample_summary"][0]
    assert (summary["assigned_R_CDS"], summary["assigned_RM_CDS"]) == (1, 3)
    assert (summary["retained_R_CDS"], summary["retained_CDS"]) == (0, 1)
    assert summary["source_counts"]["DefenseFinder"] == 1
    assert summary["effector_CDS"] == summary["motif_CDS"] == 0
    assert summary["priority_CDS"] == 1
    row = enrich_manuscript_model(data, result)["manuscript_atlas"]["rows"][0]
    assert (row["retained"]["R_CDS"], row["retained"]["total_CDS"]) == (1, 3)

    paths = export_suite_workbook(data, result, tmp_path)
    assigned = read_tsv(paths["assigned"])
    assert set(assigned) == {"methylase", "restriction", "motor"}
    assert set(read_tsv(paths["retained"])) == {"methylase"}
    assert set(read_tsv(paths["priority"])) == {"restriction"}
    assert assigned["motor"]["assigned_rm_role"] == "Aux"
    assert assigned["motor"]["R_subunit"] == "False"
    assert assigned["motor"]["effector_CDS"] == "False"
    assert assigned["restriction"]["retention_state"] == "not_retained"
    assert assigned["restriction"]["rm_score"] == "35"
    assert assigned["restriction"]["assigned_rm_type"] == "Type IV"
    wb = load_workbook(paths["workbook"])
    assert wb["R-M assigned"].max_row == 4
    assert wb["R-M retained"].max_row == 2
    assert wb["R-M evaluation"].max_row == 5


def test_fused_roles_count_one_cds_and_unknown_role_remains_assigned(tmp_path):
    data = bundle()
    data["saved_decisions"][0]["role"] = "RM"
    data["annotations"][1]["role"] = None
    result = analyze_saved(data)
    summary = result["sample_summary"][0]
    assert (summary["assigned_R_CDS"], summary["assigned_RM_CDS"]) == (2, 3)
    assert summary["assigned_counts_status"] == "partial"
    assigned = read_tsv(export_suite_workbook(data, result, tmp_path)["assigned"])
    assert len(assigned) == 3
    assert assigned["motor"]["assigned_rm_type"] == "Type IV"
    assert assigned["motor"]["assigned_rm_role"] == ""


def test_inconsistent_final_assignment_ledger_fails_without_writing(tmp_path):
    data = bundle()
    result = analyze_saved(data)
    result["gene_ledger"][0]["rm_assignments"] = []
    with pytest.raises(ValueError, match="Assigned R-M component and CDS ledgers disagree"):
        export_suite_workbook(data, result, tmp_path)
    assert not (tmp_path / "tables").exists()


def test_legacy_export_does_not_admit_stale_reviewed_annotation(tmp_path):
    data = bundle()
    result = analyze_saved(data)
    for row in result["gene_ledger"]:
        row.pop("rm_assignments")
        row.pop("final_RM")
    result.pop("rm_components")
    data["annotations"][0]["protein_sha256"] = "0" * 64
    assigned = read_tsv(export_suite_workbook(data, result, tmp_path)["assigned"])
    assert set(assigned) == {"methylase", "motor"}


def test_untyped_restriction_role_does_not_exceed_assigned_denominator():
    data = bundle()
    data["annotations"] = []
    data["saved_decisions"][0].update(family="Unresolved", role="R")
    result = analyze_saved(data)
    summary = result["sample_summary"][0]
    assert (summary["assigned_R_CDS"], summary["assigned_RM_CDS"]) == (0, 0)
    assert summary["retained_R_CDS"] == summary["priority_CDS"] == 1
    row = enrich_manuscript_model(data, result)["manuscript_atlas"]["rows"][0]
    assert (row["retained"]["R_CDS"], row["retained"]["total_CDS"]) == (0, 0)
    assert row["retained"]["status"] == "partial"


@pytest.mark.parametrize("typ", ["Type I", "Type II", "Type III", "Type IV"])
def test_reviewed_component_projection_and_equivalent_role_alias(typ):
    data = bundle()
    data["annotations"] = data["annotations"][:1]
    data["annotations"][0].update(family=typ, role="R_component")
    data["annotations"].append(dict(data["annotations"][0], rule_id="second_bound_source", role="R"))
    result = analyze_saved(data)
    row = enrich_manuscript_model(data, result)["manuscript_atlas"]["rows"][0]
    role = row["rm_roles"][typ.replace("Type ", "") + "_R"]
    assert role["count"] == 1
    assert result["sample_summary"][0]["assigned_R_CDS"] == 1
    assert result["sample_summary"][0]["retained_R_CDS"] == 0
    assert result["sample_summary"][0]["effector_CDS"] == 0
    assert row["compositions"]["I_RMS"]["count"] == 0

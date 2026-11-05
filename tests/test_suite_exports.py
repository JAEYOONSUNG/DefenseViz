from copy import deepcopy
import csv
import builtins

import pytest

from defenseviz.suite_exports import export_suite_workbook


def fixture():
    genes = [{"sample_id": s, "gene_id": g, "protein_sha256": p, "protein_length": 120, "start": 10, "end": 370}
             for s, g, p in [("A", "same", "a"), ("B", "same", "b"), ("A", "eff", "e"), ("A", "other", "o"), ("A", "rejected", "x")]]
    retained = [{"sample_id": s, "gene_id": "same", "state": "retained", "rm_type": "Type I", "role": "R",
                 "score_total": 55, "integrated_R_component_support": s == "A"} for s in ["A", "B"]]
    data = {"samples": [{"sample_id": "A", "label": "Sample alpha"}, {"sample_id": "B", "label": "Sample beta"}],
            "genes": genes, "saved_decisions": retained + [{"sample_id": "A", "gene_id": "rejected", "state": "not_retained", "role": "M", "score_total": 25}],
            "calls": [{"sample_id": "A", "gene_id": "same", "caller": "DF", "rm_type": "Type III", "role": "M", "selected": True},
                      {"sample_id": "A", "gene_id": "same", "caller": "PAD", "rm_type": "Type I", "role": "R", "selected": True}],
            "annotations": [{"sample_id": "A", "gene_id": "same", "protein_sha256": "wrong", "product": "wrong protein"}],
            "suite": {"annotations": [
                {"sample_id": "A", "gene_id": "same", "native_fields": {"product": "restriction protein", "gene": "hsdR", "protein_id": "pA", "PFAM": "PF00001"},
                 "source_ref": {"path": "alpha/input.tsv", "sha256": "hash", "row": 3}},
                {"sample_id": "B", "gene_id": "same", "native_fields": {"product": "beta restriction protein", "protein_id": "pB"}},
                {"sample_id": "A", "gene_id": "other", "native_fields": {"product": "=1+1"}},
            ]}, "effector_evidence": [{"sample_id": "A", "gene_id": "same", "rule_id": "qualified_R"}]}
    ledger = []
    for gene in genes[:-1]:
        s, g = gene["sample_id"], gene["gene_id"]
        ledger.append({"sample_id": s, "gene_id": g, "retained": g == "same", "priority": g != "other", "effector": s == "A" and g != "other",
                       "associated": g != "other", "restriction_role": g == "same", "priority_support_qualified": s == "A" and g != "other",
                       "source_membership": {"DF": g == "same", "PAD": s == "A"}, "family_ids": ["R-M" if g == "same" else "Example"], "activity": "unassessed"})
    result = {"gene_ledger": ledger, "priority_ledger": [r for r in ledger if r["priority"]],
              "automatic_system_review": {"gene_assessments": [dict(sample_id=r["sample_id"], gene_id=r["gene_id"], review_execution="completed",
                 endpoints={"effector": {"status": "supported" if r["effector"] else "unresolved"}, "associated": {"status": "supported" if r["associated"] else "unresolved"}}) for r in ledger]},
              "inventory_rows": {"families": [{"id": "R-M", "N": 2, "G": 2, "effector_CDS": 1, "associated_CDS": 2, "source_counts": {"DF": 2, "PAD": 1}}],
                                 "categories": [{"id": "R-M", "label": "R-M", "member_keys": [["A", "same"], ["B", "same"]], "N": 2, "G": 2}]},
              "system_review": {"units": [{"sample_id": "A", "provider": "DF", "unit_id": "u", "native_families": ["R-M"],
                  "member_keys": [["A", "same"], ["A", "other"]], "selected_member_keys": [["A", "same"]], "review_execution": {"status": "completed"}}]},
              "type_tables": [{"sample_id": "A", "rm_type": "Type I", "unique_CDS": 1, "R": 1, "M": 0, "S": 0, "A": 0}]}
    return data, result


def tsv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sheet_records(ws):
    values = list(ws.values)
    return [dict(zip(values[0], row)) for row in values[1:]]


def test_exact_gene_keys_and_separate_selection_scopes(tmp_path):
    data, result = fixture()
    paths = export_suite_workbook(data, result, tmp_path)
    assert set(paths) == {"workbook", "classified", "retained", "priority"}
    classified = tsv(paths["classified"])
    assert len(classified) == 4
    assert {(r["sample_id"], r["gene_id"]) for r in tsv(paths["retained"])} == {("A", "same"), ("B", "same")}
    assert {(r["sample_id"], r["gene_id"]) for r in tsv(paths["priority"])} == {("A", "same"), ("B", "same"), ("A", "eff")}
    a = next(r for r in classified if r["sample_id"] == "A" and r["gene_id"] == "same")
    b = next(r for r in classified if r["sample_id"] == "B")
    assert a["product"] == "restriction protein" and b["product"] == "beta restriction protein"
    assert a["excluded_protein_mismatch_records"] == "1"
    assert a["native_types"] == "DF: Type III; PAD: Type I"
    assert a["native_roles"] == "DF: M; PAD: R"
    assert a["rm_type"] == "Type I" and a["rm_role"] == "R"
    assert a["domain_annotations"] == "PF00001"
    assert "alpha/input.tsv" in a["source_refs"]
    assert b["effector_CDS"] == "False" and b["priority_CDS"] == "True"
    assert b["activity_assessment"] == "unassessed"


def test_workbook_is_readable_numeric_and_preserves_literal_annotation(tmp_path):
    from openpyxl import load_workbook
    data, result = fixture()
    paths = export_suite_workbook(data, result, tmp_path)
    wb = load_workbook(paths["workbook"])
    assert wb.sheetnames == ["Overview", "Classified CDSs", "R-M retained", "Priority CDSs", "R-M evaluation", "Systems", "Families", "Categories", "Type summary"]
    assert wb["Classified CDSs"].freeze_panes == "D2"
    assert wb["Classified CDSs"].auto_filter.ref.endswith("5")
    assert wb["R-M evaluation"].max_row == 4  # Includes rejected nominee outside integrated CDSs.
    overview = {r["metric"]: r["value"] for r in sheet_records(wb["Overview"])}
    assert overview["Classified CDSs"] == 4 and overview["R-M retained CDSs"] == 2 and overview["Priority CDSs"] == 3
    rows = sheet_records(wb["Classified CDSs"])
    assert all(isinstance(r["protein length"], int) for r in rows)
    assert all(isinstance(r["retained RM"], bool) for r in rows)
    assert "source refs" not in rows[0]
    assert next(r for r in rows if r["gene id"] == "other")["product"] == "=1+1"
    for ws in wb:
        for row in ws:
            for cell in row:
                assert cell.font.name == "Arial"
                assert cell.data_type != "f"
                if isinstance(cell.value, (int, float)) and not isinstance(cell.value, bool):
                    assert cell.alignment.horizontal == "right"
    system = sheet_records(wb["Systems"])[0]
    assert system["member CDSs"] == 2 and system["selected member CDSs"] == 1


@pytest.mark.parametrize("kind", ["missing_gene", "duplicate_gene", "retention_disagreement", "priority_disagreement"])
def test_inconsistent_identity_or_selection_fails_before_writing(tmp_path, kind):
    data, result = fixture()
    if kind == "missing_gene":
        data["genes"] = data["genes"][1:]
    elif kind == "duplicate_gene":
        data["genes"].append(deepcopy(data["genes"][0]))
    elif kind == "retention_disagreement":
        result["gene_ledger"][0]["retained"] = False
    else:
        result["priority_ledger"] = []
    with pytest.raises(ValueError):
        export_suite_workbook(data, result, tmp_path)
    assert not (tmp_path / "tables").exists()


def test_empty_inventory_has_explicit_headers_and_zero_overview(tmp_path):
    from openpyxl import load_workbook
    paths = export_suite_workbook({"genes": [], "saved_decisions": []}, {"gene_ledger": [], "priority_ledger": []}, tmp_path)
    wb = load_workbook(paths["workbook"])
    assert wb["Classified CDSs"].max_row == 1
    assert tsv(paths["classified"]) == []
    assert sheet_records(wb["Overview"])[1]["value"] == 0


@pytest.mark.parametrize("text", ["+enzyme annotation", "-enzyme annotation", "@enzyme annotation"])
def test_annotation_prefixes_remain_literal_excel_strings(tmp_path, text):
    from openpyxl import load_workbook
    data, result = fixture()
    data["suite"]["annotations"][-1]["native_fields"]["product"] = text
    wb = load_workbook(export_suite_workbook(data, result, tmp_path)["workbook"])
    values = [cell for row in wb["Classified CDSs"] for cell in row if cell.value == text]
    assert len(values) == 1 and values[0].data_type == "s"


def test_existing_deliverables_are_not_overwritten(tmp_path):
    data, result = fixture()
    paths = export_suite_workbook(data, result, tmp_path)
    before = paths["workbook"].read_bytes()
    with pytest.raises(FileExistsError):
        export_suite_workbook(data, result, tmp_path)
    assert paths["workbook"].read_bytes() == before


def test_missing_excel_dependency_is_actionable(tmp_path, monkeypatch):
    original = builtins.__import__
    def blocked(name, *args, **kwargs):
        if name == "openpyxl":
            raise ImportError("missing")
        return original(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(RuntimeError, match=r"pip install 'defenseviz\[suite\]'"):
        export_suite_workbook({}, {}, tmp_path)

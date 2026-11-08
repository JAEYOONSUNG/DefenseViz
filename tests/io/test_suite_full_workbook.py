"""Complete native workbooks are authoritative annotation inputs, not saved calls."""
import csv
import hashlib
import json
from xml.etree import ElementTree
from zipfile import ZipFile

import openpyxl
import pytest

from defenseviz.io.suite import load_suite
from defenseviz.native_pipeline import _prepare_sample


HEADERS = ["locus_tag", "translation", "product", "EggNOG_PFAMs",
           "Signature accession_Pfam", "dbAPIS_description", "note", "number", "flag"]


def native(tmp_path, name="isolate"):
    root = tmp_path / name
    module = root / "dnmb_module_rebasefinder"
    module.mkdir(parents=True)
    reduced = module / "rebasefinder_input.tsv"
    reduced.write_text("locus_tag\ttranslation\tproduct\ng1\tMAAA\treduced annotation\n")
    return root, reduced


def workbook(root, rows=None, name="isolate_total.xlsx", sheet="1.GenBank_table"):
    book = openpyxl.Workbook()
    table = book.active
    table.title = sheet
    table.append(HEADERS)
    rows = rows or [["g1", "MAAA", "complete annotation", "Domain_A", "PF12345",
                     "optional native annotation", "literal\tnote\nwith \"quotes\"", 1.25, True],
                    ["g2", "MCCC", "second protein", "Domain_B", "PF54321", "", "", 7, False],
                    ["rna1", "", "untranslated feature", "", "", "", "", "", ""]]
    for row in rows:
        table.append(row)
    book.create_sheet("Other_native_output").append(["not", "gene", "annotations"])
    path = root / name
    book.save(path)
    book.close()
    return path, [dict(zip(HEADERS, row)) for row in rows]


def annotations(data):
    return [row for row in data["suite"]["annotations"] if row.get("namespace") is None]


def declare(root, genes):
    path = root.parent / "suite.json"
    path.write_text(json.dumps({"schema_version": 1, "samples": [
        {"sample_id": root.name, "path": root.name, "genes": str(genes.relative_to(root.parent))}]}))
    return path


def test_parent_selects_complete_workbook_before_valid_reduced_input(tmp_path):
    root, reduced = native(tmp_path)
    path, expected = workbook(root)
    data = load_suite(tmp_path, strict_inputs=True)
    rows = annotations(data)
    assert [row["native_fields"] for row in rows] == expected
    assert all(list(row["native_fields"]) == HEADERS for row in rows)
    assert [row["source_ref"]["row"] for row in rows] == [2, 3, 4]
    assert all(row["source_ref"]["sheet"] == "1.GenBank_table" for row in rows)
    assert all(row["source_ref"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest() for row in rows)
    assert str(reduced.resolve()) not in data["input_sha256"]
    selected = next(s for s in data["suite"]["sources"] if s["role"] == "gene_annotations")
    assert selected["selection_reason"] == "complete_native_workbook"
    assert selected["selected"] and selected["row_count"] == 3
    assert {g["gene_id"] for g in data["genes"]} == {"g1", "g2", "rna1"}


def test_unselected_malformed_reduced_input_does_not_fail_strict_workbook(tmp_path):
    root, reduced = native(tmp_path)
    reduced.write_text("locus_tag\ttranslation\tproduct\ng1\tMAAA\tunescaped\ncontinuation\n")
    workbook(root)
    data = load_suite(tmp_path, strict_inputs=True)
    assert len(annotations(data)) == 3
    assert str(reduced.resolve()) not in data["input_sha256"]
    assert not data["suite"]["conflicts"]


def test_explicit_genes_override_both_workbook_and_automatic_ambiguity(tmp_path):
    root, reduced = native(tmp_path)
    workbook(root)
    workbook(root, name="other_total.xlsx")
    data = load_suite(tmp_path, manifest=declare(root, reduced), strict_inputs=True)
    assert [r["native_fields"]["product"] for r in annotations(data)] == ["reduced annotation"]
    assert not any(p.endswith("_total.xlsx") for p in data["input_sha256"])
    selected = next(s for s in data["suite"]["sources"] if s["role"] == "gene_annotations")
    assert selected["selection_reason"] == "explicit_manifest_genes"


def test_explicit_native_workbook_disambiguates_standard_multi_sheet_outputs(tmp_path):
    root, _ = native(tmp_path)
    workbook(root)
    chosen, expected = workbook(root, name="chosen_total.xlsx", rows=[
        ["chosen", "MYYY", "chosen native source", "Domain_Z", "PF67890", "", "", "", ""]])
    data = load_suite(tmp_path, manifest=declare(root, chosen), strict_inputs=True)
    assert [r["native_fields"] for r in annotations(data)] == expected
    assert annotations(data)[0]["source_ref"]["sheet"] == "1.GenBank_table"
    assert annotations(data)[0]["source_ref"]["sha256"] == hashlib.sha256(chosen.read_bytes()).hexdigest()
    assert not any(p.endswith("isolate_total.xlsx") for p in data["input_sha256"])


def test_explicit_generic_single_sheet_workbook_keeps_its_declared_source(tmp_path):
    root, _ = native(tmp_path)
    workbook(root)
    book = openpyxl.Workbook()
    book.active.append(["locus_tag", "translation", "product"])
    book.active.append(["selected", "MFFF", "explicit single sheet"])
    chosen = root / "custom.xlsx"
    book.save(chosen)
    book.close()
    data = load_suite(tmp_path, manifest=declare(root, chosen))
    assert [r["gene_id"] for r in annotations(data)] == ["selected"]


def test_explicit_malformed_tsv_does_not_fall_back_to_workbook(tmp_path):
    root, reduced = native(tmp_path)
    workbook(root)
    reduced.write_text("locus_tag\ttranslation\ng1\tMAAA\textra\n")
    with pytest.raises(ValueError, match="row/header column count differs"):
        load_suite(tmp_path, manifest=declare(root, reduced))


def test_multiple_automatic_workbooks_fail_instead_of_narrowing_to_tsv(tmp_path):
    root, _ = native(tmp_path)
    workbook(root)
    workbook(root, name="other_total.xlsx")
    with pytest.raises(ValueError, match="multiple total workbooks"):
        load_suite(tmp_path)


@pytest.mark.parametrize("problem", ["missing_sheet", "duplicate_gene", "unsafe_formula"])
def test_invalid_selected_workbook_never_falls_back_to_reduced_tsv(tmp_path, problem):
    root, _ = native(tmp_path)
    rows = [["g1", "MAAA", "protein", "", "", "", "", "", ""]]
    if problem == "duplicate_gene":
        rows.append(rows[0][:])
    if problem == "unsafe_formula":
        rows[0][2] = "=1+1"
    workbook(root, rows, sheet="wrong" if problem == "missing_sheet" else "1.GenBank_table")
    with pytest.raises(ValueError, match={"missing_sheet": "worksheet", "duplicate_gene": "duplicate gene annotation", "unsafe_formula": "formulas are not imported"}[problem]):
        load_suite(tmp_path)


def test_tsv_only_import_preserves_existing_contract(tmp_path):
    root, reduced = native(tmp_path)
    data = load_suite(tmp_path, strict_inputs=True)
    assert len(annotations(data)) == 1
    assert data["input_sha256"][str(reduced)] == hashlib.sha256(reduced.read_bytes()).hexdigest()
    reduced.write_text("locus_tag\ttranslation\ng1\tMAAA\textra\n")
    with pytest.raises(ValueError, match="row/header column count differs"):
        load_suite(tmp_path)


def test_large_literal_cells_and_hyperlink_labels_reach_full_preparation(tmp_path):
    root, _ = native(tmp_path)
    path, expected = workbook(root)
    # openxlsx can emit cells larger than Excel's editor limit. Patch only the
    # synthetic XML fixture so openpyxl's writer cannot truncate the test value.
    literal = "A" * 150_000 + '\tquoted "text"\nsecond line'
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(path) as z:
        contents = {n: z.read(n) for n in z.namelist()}
    xml = ElementTree.fromstring(contents["xl/worksheets/sheet1.xml"])
    xml.find(".//x:c[@r='G2']/x:is/x:t", namespace).text = literal
    # A literal-only hyperlink exercises presentation labels without evaluation.
    cell = xml.find(".//x:c[@r='C2']", namespace)
    cell.clear()
    cell.set("r", "C2")
    ElementTree.SubElement(cell, "{" + namespace["x"] + "}f").text = 'HYPERLINK("https://example.org/item","literal label")'
    contents["xl/worksheets/sheet1.xml"] = ElementTree.tostring(xml)
    with ZipFile(path, "w") as z:
        for name, content in contents.items():
            z.writestr(name, content)
    expected[0]["note"], expected[0]["product"] = literal, "literal label"
    data = load_suite(tmp_path, strict_inputs=True)
    assert [r["native_fields"] for r in annotations(data)] == expected
    transform = next(t for t in data["suite"]["table_transforms"] if t["operation"] == "literal_hyperlink_label")
    assert transform["formula_evaluated"] is False and transform["label"] == "literal label"
    # Discovery has already run; the preparation helper only needs the pinned
    # original-file presence. No native code or detector runs in this test.
    genome = root / "source.gbff"
    genome.write_text("synthetic source bytes")
    out = tmp_path / "prepared"
    out.mkdir()
    assert _prepare_sample(data, root.name, root, out) == genome
    previous = csv.field_size_limit()
    try:
        csv.field_size_limit(1_000_000)
        with (out / "annotations.tsv").open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            actual = list(reader)
            assert reader.fieldnames == HEADERS
    finally:
        csv.field_size_limit(previous)
    assert actual == [{key: str(value) for key, value in row.items()} for row in expected[:2]]
    binding = json.loads((out / "protein_binding.json").read_text())
    assert binding["source_protein_sha256"] == {"g1": hashlib.sha256(b"MAAA").hexdigest(), "g2": hashlib.sha256(b"MCCC").hexdigest()}


def test_workbook_protein_mutation_conflicts_with_exact_fasta_binding(tmp_path):
    root, _ = native(tmp_path)
    workbook(root)
    (root / "proteins.faa").write_text(">g1\nMAWA\n>g2\nMCCC\n")
    with pytest.raises(ValueError, match="conflicting protein sequence"):
        load_suite(tmp_path, strict_inputs=True)

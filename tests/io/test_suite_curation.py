"""Upstream curation is an exact observation, never a manufactured assignment."""
from copy import deepcopy
import hashlib

import pytest

from defenseviz.io import suite
from defenseviz.io.suite import load_suite
from defenseviz.io.suite_curation import native_curation_record
from tests.io.test_suite import native, write_table
from tests.io.test_suite_review_scope import rewrite


def curation(root, **updates):
    path = root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv"

    def mutate(headers, rows):
        rows[0].update(updates)
        return list(dict.fromkeys([*headers, *updates])), rows

    rewrite(path, mutate)
    return path


def test_native_curation_preserves_exact_scope_flags_sources_and_calls(tmp_path, monkeypatch):
    root = native(tmp_path)
    path = curation(root, rm_association_class="DNA_MTase_RM_association_unproven",
                    annotation_class="RM_specific", curation_decision="retain_high_confidence",
                    curation_reason="near_full_length_REBASE_match", other_defense_confirmed="FALSE",
                    other_defense_conflict="TRUE", evidence_conflict="FALSE", evidence_axes="reference+annotation",
                    independent_evidence_axes="2")
    data = load_suite(root)
    record = data["curation_records"][0]
    assert len(data["curation_records"]) == 2
    assert record["source_kind"] == "DNMB_REBASEfinder_curation"
    assert record["decision_namespace"] == "DNMB_native_curation"
    assert record["downstream_policy_reevaluated"] is False
    assert record["rm_association_class"] == "DNA_MTase_RM_association_unproven"
    assert record["other_defense_confirmed"] is False
    assert record["other_defense_conflict"] is True
    assert record["independent_evidence_axes"] == 2
    assert record["curation_keep"] is True and record["curation_score"] == 7
    assert record["binding_status"] == "bound" and record["usable_for_attribution"]
    assert record["protein_sha256"] == hashlib.sha256(b"MAAA").hexdigest()
    assert record["identity_binding"]["protein_binding_origin"] == "importer_attached"
    assert record["identity_binding"]["gene_source_refs"]
    assert record["source_ref"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert record["source_ref"]["row"] == 2
    raw = next(x for x in data["suite"]["features"] if x["gene_id"] == "g:a")
    assert record["native_fields"] == raw["native_fields"]
    assert record["source_ref"] == raw["source_ref"]
    coverage = next(c for c in data["suite"]["capabilities"] if "native_curation_status" in c)
    assert coverage["native_curation_status"] == "partial"
    assert coverage["native_curation_artifact_read_complete"] is True
    assert coverage["native_curation_artifact_scope"] == "all_evidence_export"
    assert (coverage["native_curation_source_rows"], coverage["native_curation_observed_rows"],
            coverage["native_curation_bound_rows"], coverage["native_curation_missing_rows"]) == (3, 2, 2, 1)
    assert coverage["native_curation_review_exhaustive"] is None
    monkeypatch.setattr(suite, "native_curation_record", lambda *args: None)
    without_projection = load_suite(root)
    for key in ("saved_decisions", "calls", "groups", "annotations", "effector_evidence"):
        assert data[key] == without_projection[key]


@pytest.mark.parametrize("updates,code", [
    ({"protein_sha256": "0" * 64}, "curation_sequence_digest_mismatch"),
    ({"protein_sha256": "not-a-hash"}, "invalid_curation_sequence_digest"),
    ({"translation": "MAAT"}, "curation_sequence_digest_mismatch"),
    ({"protein_length": "5"}, "curation_gene_metadata_mismatch"),
    ({"start": "2"}, "curation_gene_metadata_mismatch"),
    ({"start": "0"}, "invalid_curation_coordinates"),
    ({"end": "inf"}, "invalid_curation_coordinates"),
    ({"direction": "?"}, "invalid_curation_coordinates"),
    ({"locus_tag": "g_b"}, "curation_identity_field_conflict"),
    ({"other_defense_confirmed": "maybe"}, "invalid_curation_boolean"),
    ({"independent_evidence_axes": "1.5"}, "invalid_curation_number"),
])
def test_conflicting_curation_is_quarantined_without_rewriting_native_decision(tmp_path, updates, code):
    root = native(tmp_path)
    curation(root, **updates)
    data = load_suite(root)
    record = data["curation_records"][0]
    assert record["binding_status"] == "conflicting" and not record["usable_for_attribution"]
    assert code in {x["code"] for x in record["validation_issues"]}
    coverage = next(c for c in data["suite"]["capabilities"] if "native_curation_status" in c)
    assert coverage["native_curation_conflicting_rows"] == 1
    assert data["saved_decisions"][0]["state"] == "retained"
    assert data["saved_decisions"][0]["score_total"] == 7
    assert any(x["code"] == "invalid_native_curation_binding" and x["action"] == "raw_curation_quarantined"
               for x in data["suite"]["conflicts"])
    with pytest.raises(ValueError, match="invalid_native_curation_binding"):
        load_suite(root, strict_inputs=True)


def test_exact_source_sequence_digest_does_not_claim_catalysis(tmp_path):
    root = native(tmp_path)
    curation(root, protein_sha256=hashlib.sha256(b"MAAA").hexdigest(), translation="ma aa",
             rm_association_class="other_defense_system", other_defense_confirmed="TRUE")
    data = load_suite(root)
    record = data["curation_records"][0]
    assert record["identity_binding"]["protein_binding_origin"] == "native_digest_verified"
    assert record["usable_for_attribution"]
    assert data["annotations"] == data["effector_evidence"] == []


def test_missing_sequence_does_not_become_a_bound_review(tmp_path):
    root = native(tmp_path)
    path = root / "dnmb_module_rebasefinder/rebasefinder_input.tsv"
    rewrite(path, lambda headers, rows: (headers, [{**r, "translation": ""} for r in rows]))
    data = load_suite(root)
    assert all(r["binding_status"] == "unassessed" and not r["usable_for_attribution"]
               for r in data["curation_records"])
    assert data["saved_decisions"][0]["state"] == "retained"


def test_unbound_curation_keeps_exact_native_query_and_raw_row(tmp_path):
    root = native(tmp_path)
    curation(root, query="unmapped:locus")
    data = load_suite(root)
    record = data["curation_records"][0]
    assert record["gene_id"] == "unmapped:locus" and record["binding_status"] == "unassessed"
    assert not record["usable_for_attribution"] and record["protein_sha256"] is None
    assert any(r["native_fields"]["query"] == "unmapped:locus" for r in data["suite"]["unmapped_records"])


@pytest.mark.parametrize("filename,sheet", [
    ("R-M_REBASE_analysis.xlsx", "All_Augmented_Evidence"),
    ("DNMB_REBASEfinder_augmented_hits.xlsx", "All_evidence"),
])
def test_actual_upstream_workbook_curation_schema(tmp_path, filename, sheet):
    import openpyxl
    root = native(tmp_path)
    (root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv").unlink()
    workbook = openpyxl.Workbook()
    page = workbook.active
    page.title = sheet
    page.append(["locus_tag", "curation_keep", "curation_tier", "curation_decision", "final_family",
                 "final_role", "rm_association_class", "other_defense_confirmed", "contig", "start", "end", "direction"])
    page.append(["g:a", False, "other_defense", "separate_other_defense", "other_defense", "other_defense",
                 "other_defense_system", True, "chr A", 1, 12, "+"])
    path = root / "dnmb_module_rebasefinder" / filename
    workbook.save(path)
    data = load_suite(root)
    record, = data["curation_records"]
    assert record["curation_keep"] is False and record["other_defense_confirmed"] is True
    assert record["usable_for_attribution"]
    assert record["source_ref"]["sheet"] == sheet
    assert record["source_ref"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert record["final_family"] == "other_defense"  # No invented BREX/DISARM subtype.


def test_subset_workbook_cannot_stand_in_for_all_evidence(tmp_path):
    import openpyxl
    root = native(tmp_path)
    (root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv").unlink()
    workbook = openpyxl.Workbook()
    workbook.active.title = "Other_defense"
    workbook.active.append(["query", "curation_keep", "curation_tier"])
    workbook.active.append(["g:a", False, "other_defense"])
    workbook.save(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.xlsx")
    assert load_suite(root)["curation_records"] == []


def test_product_or_family_text_alone_is_not_curation(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "product", "final_family", "final_role"], [["g:a", "BREX methylase", "BREX", "M"]])
    data = load_suite(root)
    assert data["curation_records"] == []
    assert data["annotations"] == data["effector_evidence"] == []
    coverage = next(c for c in data["suite"]["capabilities"] if "native_curation_status" in c)
    assert coverage["native_curation_status"] == "not_exported"
    assert coverage["native_curation_missing_rows"] == 1


@pytest.mark.parametrize("has_export,expected", [(False, "missing"), (True, "empty_export")])
def test_empty_export_and_missing_source_are_distinct_without_claiming_review_completion(tmp_path, has_export, expected):
    root = native(tmp_path)
    path = root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv"
    path.unlink()
    if has_export:
        write_table(path, ["query", "curation_keep", "curation_tier", "rm_association_class"])
    data = load_suite(root)
    assert data["curation_records"] == []
    coverage = next(c for c in data["suite"]["capabilities"] if "native_curation_status" in c)
    assert coverage["native_curation_status"] == expected
    assert coverage["native_curation_artifact_read_complete"] is has_export
    assert coverage["native_curation_review_exhaustive"] is None


def test_projection_does_not_mutate_inputs_and_rejects_unpinned_source():
    native_row = {"query": "x", "curation_keep": False, "curation_tier": "review"}
    gene = {"sample_id": "s", "gene_id": "x", "protein_sha256": hashlib.sha256(b"MAAA").hexdigest(),
            "source_refs": [{"path": "input.tsv", "sha256": "a" * 64, "row": 2}]}
    original_row, original_gene = deepcopy(native_row), deepcopy(gene)
    record = native_curation_record("s", "x", native_row, {"path": "calls.tsv", "row": 2}, gene)
    assert not record["usable_for_attribution"]
    assert record["validation_issues"][0]["code"] == "invalid_curation_source_binding"
    assert native_row == original_row and gene == original_gene

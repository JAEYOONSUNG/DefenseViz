"""Reference-observation coverage is distinct from unresolved candidate identity."""
import csv
import hashlib
import json

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.io.suite import SuiteReader, load_suite
from defenseviz.rm_evaluation import _candidate_integrity_issues


def table(path, headers, rows=()):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="," if path.suffix == ".csv" else "\t")
        writer.writerow(headers)
        writer.writerows(rows)


def sample(tmp_path, name="sample", rows=None, headers=None):
    root = tmp_path / name
    table(root / "dnmb_module_rebasefinder/rebasefinder_input.tsv",
          ["locus_tag", "translation", "contig"], [["geneA", "MAAA", "chr"]])
    table(root / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv",
          ["sys_id", "type", "hit_id"], [["system1", "Wadjet", "geneA"]])
    table(root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv",
          ["system.number", "seqid", "system", "target.name"], [[1, "chr", "Wadjet", "geneA"]])
    table(root / "dnmb_module_defensepredictor/defense_predictor_output.csv",
          ["locus_tag", "mean_log_odds", "predicted_general_defense"], [["geneA", 5, "True"]])
    headers = headers or ["query", "curation_keep", "family_id", "hit_label", "blast_identity"]
    rows = rows if rows is not None else [["geneA", "FALSE", "Type II", "referenceA", "0.8"]]
    table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv", headers, rows)
    return root


def reference_status(data, name="sample"):
    return next(r["status"] for r in data["coverage"] if r["caller"] == "REBASE" and r["sample_id"] == name
                and r["endpoint"] == "reference_matches")


def detail(data, name="sample"):
    return next(r for r in data["suite"]["coverage"] if r["module"] == "rebasefinder" and r["sample_id"] == name)


def test_unbound_no_hit_candidate_keeps_reference_endpoint_complete_and_candidate_partial(tmp_path):
    root = sample(tmp_path, rows=[["geneA", "FALSE", "Type II", "referenceA", "0.8"],
                                 ["missing-protein", "FALSE", "NA", "NA", "NA"]])
    source = root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv"
    original = source.read_bytes()
    data = load_suite(root)
    assert reference_status(data) == "complete"
    assert detail(data)["candidate_binding_status"] == "partial"
    assert detail(data)["unbound_candidate_row_count"] == 1
    assert detail(data)["unbound_reference_row_count"] == 0
    assert {c["gene_id"] for c in data["calls"] if c["caller"] == "REBASE"} == {"geneA"}
    assert len(data["saved_decisions"]) == 1
    assert data["suite"]["unmapped_records"][0]["has_coherent_reference_pair"] is False
    assert _candidate_integrity_issues(data, "sample", {})  # Complete retention override remains blocked.
    assert source.read_bytes() == original
    assert data["input_sha256"][str(source)] == hashlib.sha256(original).hexdigest()
    with pytest.raises(ValueError, match="unbound_native_candidate"):
        load_suite(root, strict_inputs=True)


@pytest.mark.parametrize("subject_key,identity_key", [("raw_hit_label", "raw_blast_identity"),
    ("blast_match", "blast_identity"), ("hit_label", "blast_identity"),
    ("curated_blast_match", "curated_blast_identity")])
def test_unbound_actual_reference_match_keeps_endpoint_partial(tmp_path, subject_key, identity_key):
    root = sample(tmp_path, headers=["query", "curation_keep", subject_key, identity_key],
                  rows=[["missing-protein", "FALSE", "referenceA", "0.8"]])
    data = load_suite(root)
    assert reference_status(data) == "partial"
    assert detail(data)["unbound_reference_row_count"] == 1
    assert data["suite"]["unmapped_records"][0]["has_coherent_reference_pair"] is True
    assert not data["suite"]["reference_matches"]
    assert not any(c["caller"] == "REBASE" for c in data["calls"])


def test_reference_pair_never_borrows_identity_from_another_origin(tmp_path):
    root = sample(tmp_path, headers=["query", "raw_hit_label", "raw_blast_identity", "curated_blast_match", "curated_blast_identity"],
                  rows=[["missing-protein", "raw-reference", "NA", "NA", "0.8"]])
    data = load_suite(root)
    assert reference_status(data) == "partial"
    assert detail(data)["unbound_reference_row_count"] == 1
    assert detail(data)["incomplete_reference_row_count"] == 1
    assert not data["suite"]["reference_matches"]


@pytest.mark.parametrize("gene", ["geneA", "missing-protein"])
@pytest.mark.parametrize("subject_key,identity_key", [("raw_hit_label", "raw_blast_identity"),
    ("blast_match", "blast_identity"), ("hit_label", "blast_identity"),
    ("curated_blast_match", "curated_blast_identity")])
def test_incomplete_reference_observation_is_partial_when_bound_or_unbound(tmp_path, gene, subject_key, identity_key):
    root = sample(tmp_path, headers=["query", subject_key, identity_key], rows=[[gene, "referenceA", "NA"]])
    data = load_suite(root)
    assert reference_status(data) == "partial"
    assert detail(data)["incomplete_reference_row_count"] == 1
    assert not data["suite"]["reference_matches"]


@pytest.mark.parametrize("gene", ["geneA", "missing-protein"])
def test_explicit_hmm_only_generic_label_does_not_invent_reference_match(tmp_path, gene):
    root = sample(tmp_path, headers=["query", "source", "hit_label", "blast_identity"],
                  rows=[[gene, "hmm", "Synthetic_HMM_profile", "NA"]])
    data = load_suite(root)
    assert reference_status(data) == "complete"
    assert detail(data)["incomplete_reference_row_count"] == 0
    assert detail(data)["unbound_reference_row_count"] == 0
    assert not data["suite"]["reference_matches"]


@pytest.mark.parametrize("mode,label", [("operon_context", "typeI_context:R:sequence_motif"),
    ("operon_context", "typeII_context:M:annotation"), ("operon_context", "typeIII_context:R:sequence_motif"),
    ("annotation_motif_candidate", "typeIV_context:Mrr_or_modified_DNA_REase")])
def test_native_context_candidate_marker_is_preserved_without_reference_membership(tmp_path, mode, label):
    root = sample(tmp_path, headers=["query", "source", "evidence_mode", "hit_label", "blast_match", "blast_identity"],
                  rows=[["geneA", "rebasefinder", mode, label, label, "NA"]])
    data = load_suite(root)
    assert reference_status(data) == "complete"
    assert detail(data)["incomplete_reference_row_count"] == 0
    assert not data["suite"]["reference_matches"]
    assert data["suite"]["features"][0]["native_fields"]["hit_label"] == label


@pytest.mark.parametrize("source,mode,label,identity", [
    ("rebasefinder", "operon_context", "Actual.reference", "NA"),
    ("rebasefinder", "blast", "typeI_context:R:sequence_motif", "NA"),
    ("other_producer", "operon_context", "typeI_context:R:sequence_motif", "NA"),
    ("rebasefinder", "operon_context", "typeI_context:R:sequence_motif", "0.8"),
])
def test_context_mode_does_not_hide_actual_or_unsupported_reference_observations(tmp_path, source, mode, label, identity):
    root = sample(tmp_path, headers=["query", "source", "evidence_mode", "hit_label", "blast_identity"],
                  rows=[["geneA", source, mode, label, identity]])
    data = load_suite(root)
    assert reference_status(data) == "partial"
    assert detail(data)["incomplete_reference_row_count"] == 1
    assert not data["suite"]["reference_matches"]


def test_bound_observation_uses_same_raw_first_pair_selection(tmp_path):
    root = sample(tmp_path, headers=["query", "raw_hit_label", "raw_blast_identity", "curated_blast_match", "curated_blast_identity"],
                  rows=[["geneA", "raw-reference", "0.4", "curated-reference", "0.9"]])
    data = load_suite(root)
    observation = data["suite"]["reference_matches"][0]
    assert (observation["reference_id"], observation["identity"], observation["match_origin"]) == ("raw-reference", 0.4, "original_raw")


@pytest.mark.parametrize("rows", [[], [["missing-protein", "FALSE", "NA", "NA", "NA"]]])
def test_reference_schema_with_no_observed_matches_establishes_zero(tmp_path, rows):
    data = load_suite(sample(tmp_path, rows=rows))
    assert reference_status(data) == "complete"
    assert not data["suite"]["reference_matches"]
    assert detail(data)["reference_pair_schema"] is True


@pytest.mark.parametrize("rows", [[], [["geneA", "FALSE", "nuclease_profile", "1e-8"]]])
def test_hmm_only_export_cannot_establish_reference_match_absence(tmp_path, rows):
    data = load_suite(sample(tmp_path, headers=["query", "curation_keep", "hmm_hit", "hmm_evalue"], rows=rows))
    assert reference_status(data) == "partial"
    assert detail(data)["reference_pair_schema"] is False
    assert not data["suite"]["reference_matches"]


def test_declared_partial_reference_run_is_not_upgraded_by_no_hit_candidates(tmp_path):
    sample(tmp_path, rows=[["missing-protein", "FALSE", "NA", "NA", "NA"]])
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"schema_version": 1, "samples": [{"sample_id": "sample", "path": "sample",
        "modules": {"rebasefinder": {"status": "partial"}}}]}))
    data = load_suite(tmp_path, manifest=manifest)
    assert reference_status(data) == "partial"


def test_failed_native_run_stays_partial_even_when_present_match_rows_bind(tmp_path, monkeypatch):
    root = sample(tmp_path)
    monkeypatch.setattr(SuiteReader, "metadata", lambda *args: {"cache": {"modules": {"REBASEfinder": {"ok": False}}}})
    data = load_suite(root)
    assert reference_status(data) == "partial"
    assert detail(data)["candidate_binding_status"] == "complete"


@pytest.mark.parametrize("unbound_reference,expected_samples", [(False, ["one", "two"]), (True, ["one"])])
def test_four_set_common_scope_uses_reference_completeness_not_candidate_binding(tmp_path, unbound_reference, expected_samples):
    sample(tmp_path, "one")
    missing = ["missing-protein", "FALSE", "Type II", "referenceB", "0.6"] if unbound_reference else ["missing-protein", "FALSE", "NA", "NA", "NA"]
    sample(tmp_path, "two", rows=[["geneA", "FALSE", "Type II", "referenceA", "0.8"], missing])
    result = analyze_saved(load_suite(tmp_path))
    overlap = result["overlap_regions"]
    assert set(overlap["source_order"]) == {"PADLOC", "DefenseFinder", "DefensePredictor", "REBASE"}
    assert overlap["sample_ids"] == expected_samples
    assert overlap["union_count"] == len(expected_samples)
    assert overlap["source_sizes"] == {s: len(expected_samples) for s in overlap["source_order"]}
    assert next(r for r in overlap["regions"] if r["mask"] == "1111")["count"] == len(expected_samples)

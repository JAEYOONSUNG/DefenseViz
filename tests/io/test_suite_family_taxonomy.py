"""Provider taxonomy normalizes family display without borrowing activity roles."""
import csv
import json

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.io.suite import SuiteReader, load_suite


def imported(tmp_path, provider, family, native=None, mapping=None, kind="final_member"):
    reader = SuiteReader(tmp_path, None, None, False)
    reader.gene("sample", "geneA")
    source = {"path": "native.tsv", "sha256": "0" * 64, "row": 2}
    call = reader.call("sample", "geneA", provider, family, "systemA", source,
                       native=native, mapping=mapping, evidence_kind=kind, model="native_component")
    return reader, call


@pytest.mark.parametrize("family", ["CBASS_I", "CBASS_II", "CBASS_III_HNH", "CBASS_IV_TIR", "CBASS_II_Phospholipase"])
@pytest.mark.parametrize("fqn", [False, True])
def test_df_cbass_parent_and_bounded_native_subtypes_normalize(tmp_path, family, fqn):
    native = {"type": "CBASS", "subtype": family}
    if fqn:
        native["model_fqn"] = "defense-finder-models/DefenseFinder/CBASS/" + family
    reader, call = imported(tmp_path, "DefenseFinder", family, native)
    assert (call["family"], call["category"]) == ("CBASS", "Abortive infection")
    assert call["native_family"] == family and call["native_fields"] == native
    assert call["family_mapping_rule"] == "defensefinder_cbass_family_display_v1"
    assert call["role"] is call["rm_type"] is None
    assert reader.data["saved_decisions"] == reader.data["effector_evidence"] == []


@pytest.mark.parametrize("family", ["CAS_Class1-Subtype-IV-A", "CAS_Class1-Subtype-I-F3", "CAS_Class2-Subtype-VI", "CAS_Cluster"])
@pytest.mark.parametrize("fqn", [False, True])
def test_df_casfinder_native_namespace_and_legacy_grammar_normalize(tmp_path, family, fqn):
    native = {"type": "NA", "subtype": family}
    if fqn:
        native["model_fqn"] = "CasFinder/" + family
    _, call = imported(tmp_path, "DefenseFinder", family, native)
    assert (call["family"], call["category"]) == ("CRISPR-Cas", "CRISPR-Cas")
    assert call["native_family"] == family and call["role"] is None


@pytest.mark.parametrize("family,display,category", [
    ("cas_type_IV-A", "CRISPR-Cas", "CRISPR-Cas"),
    ("cas_type_I-F3", "CRISPR-Cas", "CRISPR-Cas"),
    ("cas_type_V-K", "CRISPR-Cas", "CRISPR-Cas"),
    ("cas_type_arrays", "CRISPR-Cas", "CRISPR-Cas"),
    ("cbass_type_III", "CBASS", "Abortive infection"),
    ("cbass_type_IIs", "CBASS", "Abortive infection"),
    ("cbass_type_IV", "CBASS", "Abortive infection"),
    ("PDC-M69", "PDC-M (pooled)", "Unresolved"),
    ("PDC-M01", "PDC-M (pooled)", "Unresolved"),
    ("PDC-S73", "PDC-S (pooled)", "Unresolved"),
])
def test_padloc_native_family_grammar_is_not_limited_to_one_cohort(tmp_path, family, display, category):
    native = {"system": family, "hmm.name": "Synthetic_profile"}
    _, call = imported(tmp_path, "PADLOC", family, native)
    assert (call["family"], call["category"]) == (display, category)
    assert call["native_family"] == family and call["native_fields"] == native
    assert call["system_id"] == "systemA" and call["native_model"] == "native_component"
    assert call["role"] is call["rm_type"] is None


@pytest.mark.parametrize("family,display,category", [
    ("pAgo_LongA", "pAgo LongA", "Nuclease"),
    ("pAgo_LongB_SIRim", "pAgo LongB-SIRim", "Abortive infection"),
])
def test_exact_pago_native_subtype_is_not_generic_argonaute(tmp_path, family, display, category):
    _, call = imported(tmp_path, "DefenseFinder", family, {
        "type": "pAgo", "subtype": family, "model_fqn": "defense-finder-models/DefenseFinder/pAgo/" + family})
    assert (call["family"], call["category"]) == (display, category)
    assert call["role"] is None


def test_padloc_solo_argonaute_does_not_claim_long_a_or_nuclease_role(tmp_path):
    reader, call = imported(tmp_path, "PADLOC", "argonaute_solo", {"system": "argonaute_solo", "protein.name": "pAgo"})
    assert (call["family"], call["category"]) == ("Argonaute", "Unresolved")
    assert call["native_family"] == "argonaute_solo" and call["role"] is None
    assert not reader.data["effector_evidence"]


@pytest.mark.parametrize("provider,family,native", [
    ("PADLOC", "cas_type_I-C_like", {}),
    ("PADLOC", "cas_type_VII-A", {}),
    ("PADLOC", "cbass_type_V", {}),
    ("PADLOC", "PDC-M69_like", {}),
    ("PADLOC", "PDC-MX", {}),
    ("PADLOC", "argonaute_type_I", {}),
    ("PADLOC", "argonaute_solo_like", {}),
    ("PADLOC", "CBASS_II", {}),
    ("DefenseFinder", "CBASS_V", {"type": "CBASS"}),
    ("DefenseFinder", "CBASS_II_unknown", {"type": "CBASS"}),
    ("DefenseFinder", "CBASS_II", {"type": "CBASS", "model_fqn": "custom/CBASS_II"}),
    ("DefenseFinder", "CAS_Class9-Subtype-IV-A", {}),
    ("DefenseFinder", "pAgo_LongA_like", {"type": "pAgo"}),
    ("DefenseFinder", "pAgo_LongB", {"type": "pAgo"}),
    ("DefenseFinder", "PDC-M69", {}),
    ("REBASE", "cas_type_IV-A", {}),
    ("DefensePredictor", "CBASS_II", {}),
    ("PADLOC", "SspABCD", {}),
    ("PADLOC", "SspBCD", {}),
])
def test_unknown_lookalikes_and_other_endpoints_are_not_inferred(tmp_path, provider, family, native):
    _, call = imported(tmp_path, provider, family, native)
    assert call["family"] == family and call["category"] == "Unresolved"
    assert "family_mapping_rule" not in call


def test_hmm_profile_is_not_a_native_system_taxonomy_call(tmp_path):
    _, call = imported(tmp_path, "DefenseFinder", "CBASS_II", {"type": "CBASS"}, kind="hmm_hit")
    assert call["family"] == "CBASS_II" and "family_mapping_rule" not in call
    assert call["selected"] is False


@pytest.mark.parametrize("key", [(None, "CAS_Class1-Subtype-IV-A"), ("DefenseFinder", "CAS_Class1-Subtype-IV-A")])
def test_explicit_user_family_mapping_wins_over_native_display_default(tmp_path, key):
    mapping = {key: dict(category="Other mechanisms", display_family="User grouping", rm_type=None,
                        priority=1, source="family_map", caller=key[0])}
    _, call = imported(tmp_path, "DefenseFinder", "CAS_Class1-Subtype-IV-A", {
        "subtype": "CAS_Class1-Subtype-IV-A", "model_fqn": "CasFinder/CAS_Class1-Subtype-IV-A"}, mapping=mapping)
    assert (call["family"], call["category"]) == ("User grouping", "Other mechanisms")
    assert call["native_family"] == "CAS_Class1-Subtype-IV-A" and call["category_source"] == "family_map"
    assert "family_mapping_rule" not in call


def test_normalizing_display_never_merges_distinct_native_call_identities(tmp_path):
    reader, _ = imported(tmp_path, "DefenseFinder", "CBASS_I", {"type": "CBASS"})
    reader.call("sample", "geneA", "DefenseFinder", "CBASS_II", "systemA",
                {"path": "native.tsv", "sha256": "0" * 64, "row": 3}, native={"type": "CBASS"}, model="native_component")
    assert len(reader.data["calls"]) == 2
    assert {c["family"] for c in reader.data["calls"]} == {"CBASS"}
    assert {c["native_family"] for c in reader.data["calls"]} == {"CBASS_I", "CBASS_II"}


def test_suite_import_unifies_known_source_family_display_and_keeps_systems(tmp_path):
    root = tmp_path / "sample"
    def table(relative, header, rows):
        p = root / relative
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="," if p.suffix == ".csv" else "\t")
            writer.writerow(header); writer.writerows(rows)
    table("dnmb_module_rebasefinder/rebasefinder_input.tsv", ["locus_tag", "translation"], [["geneA", "MAAA"]])
    table("dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv", ["hit_id", "sys_id", "type", "subtype", "model_fqn"],
          [["geneA", "df_system", "NA", "CAS_Class1-Subtype-IV-A", "CasFinder/CAS_Class1-Subtype-IV-A"]])
    table("dnmb_module_padloc/padloc_query_proteins_padloc.csv", ["system.number", "seqid", "system", "target.name"],
          [[1, "chr", "cas_type_IV-A", "geneA"]])
    data = load_suite(root)
    result = analyze_saved(data)
    family = next(f for f in result["inventory_rows"]["families"] if f["id"] == "CRISPR-Cas")
    assert family["N"] == 1 and family["source_counts"]["DefenseFinder"] == family["source_counts"]["PADLOC"] == 1
    assert len(data["groups"]) == 2
    assert {g["family"] for g in data["groups"]} == {"CAS_Class1-Subtype-IV-A", "cas_type_IV-A"}
    assert result["global_summary"]["effector_CDS"] == result["global_summary"]["retained_CDS"] == 0


def test_category_only_user_mapping_survives_import_and_analysis_catalog_aliases(tmp_path):
    root = tmp_path / "sample"
    genes = root / "dnmb_module_rebasefinder/rebasefinder_input.tsv"
    genes.parent.mkdir(parents=True)
    genes.write_text("locus_tag\ttranslation\ngeneA\tMAAA\n")
    calls = root / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv"
    calls.parent.mkdir()
    calls.write_text("hit_id\tsys_id\ttype\tsubtype\tmodel_fqn\tgene_name\n"
                     "geneA\tnative-system\tCBASS\tCBASS_II_HNH\t"
                     "defense-finder-models/DefenseFinder/CBASS/CBASS_II_HNH\tnative-component\n")
    baseline = load_suite(root)
    override = tmp_path / "families.tsv"
    override.write_text("family\tcaller\tcategory\nCBASS_II_HNH\tDefenseFinder\tOther mechanisms\n")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"schema_version": 1, "family_map": "families.tsv",
        "samples": [{"sample_id": "sample", "path": "sample"}]}))
    data = load_suite(root, manifest=manifest)
    before, after = baseline["calls"][0], data["calls"][0]
    assert (before["family"], before["category"]) == ("CBASS", "Abortive infection")
    assert (after["family"], after["category"], after["category_source"]) == (
        "CBASS_II_HNH", "Other mechanisms", "family_map")
    for field in ("sample_id", "gene_id", "caller", "native_family", "native_fields", "system_id",
                  "run_id", "protein_sha256", "native_model", "role", "rm_type", "selected"):
        assert before[field] == after[field]
    assert before["source_ref"]["sha256"] == after["source_ref"]["sha256"]
    assert data["saved_decisions"] == baseline["saved_decisions"] == []
    assert data["effector_evidence"] == baseline["effector_evidence"] == []
    result = analyze_saved(data)
    family = next(r for r in result["inventory_rows"]["families"] if r["id"] == "CBASS_II_HNH")
    assert family["category"] == "Other mechanisms" and family["N"] == 1
    assert not any(r["id"] == "CBASS" and r["N"] for r in result["inventory_rows"]["families"])
    assert result["global_summary"]["effector_CDS"] == result["global_summary"]["retained_CDS"] == 0

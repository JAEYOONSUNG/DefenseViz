"""Named native model taxonomy does not confer retained roles or effectors."""
import csv
import hashlib

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.io.suite import SuiteReader, load_suite


MODELS = [
    ("McrBC", "defense-finder-models/RM/McrBC/McrBC"),
    ("GmrSD_RM_Type_IV", "defense-finder-models/RM/GmrSD_RM_Type_IV/GmrSD_RM_Type_IV"),
    ("EcoKMcrA", "defense-finder-models/DefenseFinder/EcoKMcrA/EcoKMcrA"),
]


def fixture(tmp_path, rows):
    root = tmp_path / "isolate"
    genes = root / "dnmb_module_rebasefinder/rebasefinder_input.tsv"
    genes.parent.mkdir(parents=True)
    genes.write_text("locus_tag\ttranslation\ngeneA\tMAAA\n")
    calls = root / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv"
    calls.parent.mkdir()
    with calls.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["sys_id", "hit_id", "type", "subtype", "model_fqn", "gene_name"])
        writer.writerows(rows)
    return root, calls


@pytest.mark.parametrize("family,fqn", MODELS)
@pytest.mark.parametrize("with_fqn", [False, True])
def test_exact_named_native_models_are_type_iv_without_retention_or_activity(tmp_path, family, fqn, with_fqn):
    root, source = fixture(tmp_path, [["systemA", "geneA", family, family, fqn if with_fqn else "", "Synthetic_component"]])
    original = source.read_bytes()
    data = load_suite(root)
    call = next(c for c in data["calls"] if c["caller"] == "DefenseFinder")
    assert call["rm_type"] == "Type IV" and call["category"] == "R-M"
    assert call["native_family"] == family and call["native_type"] == family
    assert call["role"] is None
    assert call["mapping_rule"] == "defensefinder_exact_named_type_iv_model_v1"
    assert call["mapping_reference"].startswith("https://github.com/mdmparis/defense-finder-models/blob/")
    assert call["source_ref"]["sha256"] == hashlib.sha256(original).hexdigest()
    assert source.read_bytes() == original
    if family == "EcoKMcrA":
        assert "https://doi.org/10.1093/nar/gky731" in call["mapping_references"]
    result = analyze_saved(data)
    row = next(r for r in result["type_tables"] if r["rm_type"] == "Type IV")
    assert row["pool_CDS"] == 1 and row["unique_CDS"] == 0
    assert row["R"] == row["M"] == row["S"] == row["A"] == 0
    assert row["activity"] == "unassessed"
    assert result["global_summary"]["retained_CDS"] == result["global_summary"]["effector_CDS"] == 0
    assert not data["saved_decisions"] and not data["effector_evidence"]


def test_named_type_iv_keeps_other_native_system_membership(tmp_path):
    root, _ = fixture(tmp_path, [
        ["systemA", "geneA", "McrBC", "McrBC", "defense-finder-models/RM/McrBC/McrBC", "RM_Type_IV__Type_IV_21"],
        ["systemB", "geneA", "CoCoNut", "CoCoNut_I-A", "defense-finder-models/DefenseFinder/CoCoNut/CoCoNut_I-A", "CoCoNut-I-A__I-A_CnuB"],
    ])
    data = load_suite(root)
    calls = [c for c in data["calls"] if c["caller"] == "DefenseFinder"]
    assert len(calls) == 2
    assert {c["native_family"] for c in calls} == {"McrBC", "CoCoNut_I-A"}
    assert next(c for c in calls if c["native_family"] == "CoCoNut_I-A")["rm_type"] is None
    assert len([g for g in data["groups"] if g["provider"] == "DefenseFinder"]) == 2
    result = analyze_saved(data)
    assert result["global_summary"]["native_union_CDS"] == 1
    assert result["global_summary"]["effector_CDS"] == 0


@pytest.mark.parametrize("family,fqn", MODELS)
def test_same_model_basename_in_unrecognized_namespace_is_not_promoted(tmp_path, family, fqn):
    root, _ = fixture(tmp_path, [["systemA", "geneA", family, family, "custom/" + family, "Synthetic_component"]])
    call = load_suite(root)["calls"][0]
    assert call["rm_type"] is None
    assert call.get("mapping_rule") != "defensefinder_exact_named_type_iv_model_v1"


def test_unknown_named_rm_model_stays_type_unresolved(tmp_path):
    root, _ = fixture(tmp_path, [["systemA", "geneA", "RM", "McrBC_like_new", "defense-finder-models/RM/McrBC_like_new/McrBC_like_new", "Synthetic_component"]])
    data = load_suite(root)
    call = data["calls"][0]
    assert call["category"] == "R-M" and call["rm_type"] is None
    assert next(r for r in analyze_saved(data)["type_tables"] if r["rm_type"] == "Unresolved")["pool_CDS"] == 1


@pytest.mark.parametrize("family,fqn", MODELS)
def test_named_model_cannot_override_conflicting_native_type(tmp_path, family, fqn):
    root, _ = fixture(tmp_path, [["systemA", "geneA", "RM_Type_I", family, fqn, "Synthetic_component"]])
    with pytest.raises(ValueError, match="contradicts native taxonomy"):
        load_suite(root)


@pytest.mark.parametrize("provider,kind", [("PADLOC", "final_member"), ("REBASE", "reference_hit"), ("DefenseFinder", "hmm_hit")])
def test_named_df_alias_does_not_transfer_to_other_endpoints(tmp_path, provider, kind):
    reader = SuiteReader(tmp_path, None, None, False)
    reader.gene("sample", "geneA")
    call = reader.call("sample", "geneA", provider, "EcoKMcrA", None,
        {"path": "synthetic.tsv", "sha256": "0" * 64}, native={"type": "EcoKMcrA"}, evidence_kind=kind)
    assert call["rm_type"] is None and call["role"] is None
    assert call.get("mapping_rule") != "defensefinder_exact_named_type_iv_model_v1"

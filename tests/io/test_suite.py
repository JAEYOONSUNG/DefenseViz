"""Native, synthetic suite contracts: no manuscript identifiers or data."""
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from defenseviz.io.suite import load_suite


def write_table(path, columns, rows=()):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="," if path.suffix == ".csv" else "\t")
        writer.writerow(columns)
        writer.writerows(rows)
    return path


def native(tmp_path, name="isolate α"):
    root = tmp_path / name
    write_table(root / "dnmb_module_rebasefinder/rebasefinder_input.tsv",
                ["locus_tag", "contig", "start", "end", "direction", "translation", "product"],
                [["g:a", "chr A", 1, 12, "+", "MAAA", "unknown"],
                 ["g_b", "chr A", 21, 32, "-", "MBBB", "unknown"],
                 ["g_c", "plasmid", 1, 12, "+", "MCCC", "unknown"]])
    write_table(root / "dnmb_module_defensefinder/defensefinder_id_map.tsv",
                ["defensefinder_id", "locus_tag"], [["native1", "g:a"], ["native2", "g_b"], ["native3", "g_c"]])
    write_table(root / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv",
                ["hit_id", "sys_id", "type", "subtype", "gene_name", "hit_status"],
                [["native1", "sys A", "Wadjet", "Wadjet_I", "Wadjet__JetD_I", "mandatory"],
                 ["native2", "sys A", "Wadjet", "Wadjet_I", "Wadjet__JetC_I", "mandatory"]])
    write_table(root / "dnmb_module_defensefinder/defensefinder_systems.tsv",
                ["sys_id", "type", "subtype", "protein_in_syst", "genes_count"],
                [["sys A", "Wadjet", "Wadjet_I", "native1,native2", 2]])
    write_table(root / "dnmb_module_defensefinder/defensefinder_hmmer.tsv",
                ["hit_id", "gene_name", "i_eval", "hit_score", "hit_profile_cov", "hit_seq_cov"],
                [["native3", "Other_Profile", "1e-8", 40, 0.9, 0.8]])
    write_table(root / "dnmb_module_padloc/padloc_query_proteins_id_map.tsv",
                ["query_id", "locus_tag", "contig_id", "contig"], [["protein__safeA", "g:a", "chr A", "chr A"], ["protein__safeB", "g_b", "chr A", "chr A"], ["proteinC", "g_c", "plasmid", "plasmid"]])
    write_table(root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv",
                ["system.number", "seqid", "system", "target.name", "hmm.name", "protein.name"],
                [[1, "chr A", "RM_Type_II", "protein__safeA", "restriction_hmm", "Restriction"],
                 [1, "plasmid", "RM_Type_II", "proteinC", "restriction_hmm_2", "Restriction"],
                 [2, "chr A", "RM_Type_II", "protein__safeA", "restriction_hmm_3", "Restriction"]])
    write_table(root / "dnmb_module_defensepredictor/defense_predictor_query_id_map.tsv",
                ["query_id", "locus_tag"], [["protein__safeA", "g:a"], ["protein__safeB", "g_b"], ["proteinC", "g_c"]])
    write_table(root / "dnmb_module_defensepredictor/defense_predictor_output.csv",
                ["locus_tag", "mean_log_odds", "sd_log_odds"], [["g:a", 8, 0.5], ["g_b", -3, 0.7]])
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "curation_keep", "final_family", "final_role", "curation_tier", "curation_score",
                 "family_id", "hit_label", "blast_identity", "typeii_context_status", "typeii_context_partners"],
                [["g:a", "TRUE", "Type II", "R", "high", 7, "Type II", "R.Ref", 0.8, "complete_mr", "g_b:M:annotation:+1g:8bp"],
                 ["g_b", "FALSE", "Type II", "M", "review", 2, "Type II", "NA", "NA", "NA", "NA"],
                 ["g_c", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]])
    return root


def test_native_identity_evidence_and_completion(tmp_path):
    root = native(tmp_path)
    data = load_suite(root)
    assert data["samples"][0]["sample_id"] == "isolate α"
    calls = data["calls"]
    assert sum(c["caller"] == "PADLOC" for c in calls) == 3
    assert len({c["gene_id"] for c in calls if c["caller"] == "PADLOC"}) == 2
    assert len({c["system_id"] for c in calls if c["caller"] == "PADLOC"}) == 3
    df = [c for c in calls if c["caller"] == "DefenseFinder" and c["evidence_kind"] == "final_member"]
    assert {c["gene_id"] for c in df} == {"g:a", "g_b"}
    assert df[0]["native_model"] == "Wadjet__JetD_I"
    assert df[0]["native_fields"]["hit_status"] == "mandatory"
    assert any(c["gene_id"] == "g_c" and c["evidence_kind"] == "hmm_hit" for c in calls)
    assert [c["gene_id"] for c in calls if c["caller"] == "REBASE"] == ["g:a"]
    assert [d["state"] for d in data["saved_decisions"]] == ["retained", "not_retained", "unassessed"]
    assert data["annotations"] == data["effector_evidence"] == data["complexes"] == []
    assert data["saved_decisions"][0]["protein_sha256"] == hashlib.sha256(b"MAAA").hexdigest()
    dp = next(c for c in data["suite"]["coverage"] if c["module"] == "defensepredictor")
    assert (dp["mapped_gene_count"], dp["scored_gene_count"], dp["missing_gene_count"]) == (3, 2, 1)
    assert len(data["suite"]["scores"]) == 2
    assert all(c["status"] == "complete" for c in data["coverage"] if c["endpoint"] != "hmm_hits")
    groups = [g for g in data["groups"] if g["provider"] == "DefenseFinder"]
    assert groups[0]["native_system_id"] == "sys A"
    assert groups[0]["run_id"] == df[0]["run_id"]
    for path in data["inputs"]:
        assert data["input_sha256"][str(path)] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_missing_header_zero_and_explicit_disabled_are_distinct(tmp_path):
    root = tmp_path / "empty"
    write_table(root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv",
                ["system.number", "seqid", "system", "target.name"])
    config = tmp_path / "manifest.json"
    config.write_text(json.dumps({"schema_version": 1, "samples": [{"sample_id": "z", "path": "empty",
                      "modules": {"defensepredictor": {"status": "not_run"}}}]}))
    data = load_suite(tmp_path, manifest=config)
    assert {x["caller"]: x["status"] for x in data["coverage"] if x["endpoint"] != "hmm_hits"} == {
        "DefenseFinder": "missing", "PADLOC": "complete", "DefensePredictor": "not_run", "REBASE": "missing"}
    assert not data["calls"]


def test_actual_completed_zero_rds_without_artifact(tmp_path):
    rscript = shutil.which("Rscript")
    if not rscript:
        pytest.skip("Rscript unavailable")
    root = tmp_path / "native"
    md = root / "dnmb_module_padloc"
    md.mkdir(parents=True)
    target = md / ".dnmb_padloc_complete"
    subprocess.run([rscript, "--vanilla", "-e", "args<-commandArgs(TRUE);saveRDS(list(alias='PADLOC',n_hits=0L,elapsed_s=1,written_at='date',sentinel_version=1L),args[1])", str(target)], check=True)
    data = load_suite(root, rscript=rscript)
    assert next(c for c in data["coverage"] if c["caller"] == "PADLOC")["status"] == "complete"
    assert not data["calls"]
    assert data["suite"]["commands"][0]["exit_code"] == 0


def test_live_members_survive_stale_rds_and_strict_rejects(tmp_path, monkeypatch):
    root = native(tmp_path)
    (root / ".dnmb_module_stage_cache.rds").write_bytes(b"fake RDS fixture")
    payload = {"reader_schema_version": 1, "sentinels": {}, "cache": {"modules": {
        "REBASEfinder": {"ok": True, "table_counts": {"all_hits": 99, "hits": 8}}}}}
    def fake(argv, **kwargs):
        assert argv[1] == "--vanilla" and argv[2].endswith("read_suite_metadata.R")
        assert kwargs.get("shell") is None
        return subprocess.CompletedProcess(argv, 0, json.dumps(payload), "")
    monkeypatch.setattr(subprocess, "run", fake)
    data = load_suite(root, rscript="fixed-R")
    assert len(data["saved_decisions"]) == 3
    assert data["suite"]["conflicts"][0]["code"] == "cached_live_candidate_count_mismatch"
    with pytest.raises(ValueError, match="strict input consistency"):
        load_suite(root, rscript="fixed-R", strict_inputs=True)


def test_mapping_is_exact_not_label_parsing(tmp_path):
    root = native(tmp_path)
    path = root / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv"
    path.write_text(path.read_text().replace("native1", "native1__g:a"))
    with pytest.raises(ValueError, match="ID map has no entry"):
        load_suite(root)


@pytest.mark.parametrize("rows", [[["x", "g:a"], ["x", "g_b"]], [["x", "g:a"], ["y", "g:a"]]])
def test_non_bijective_map_fails(tmp_path, rows):
    root = native(tmp_path)
    write_table(root / "dnmb_module_defensefinder/defensefinder_id_map.tsv", ["defensefinder_id", "locus_tag"], rows)
    with pytest.raises(ValueError, match="one-to-one"):
        load_suite(root)


def test_coordinate_and_sequence_conflicts_fail(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_defensefinder/defensefinder_id_map.tsv", ["defensefinder_id", "locus_tag", "start"], [["native1", "g:a", 2]])
    with pytest.raises(ValueError, match="conflicting start"):
        load_suite(root)
    root = native(tmp_path)
    (root / "proteins.faa").write_text(">g:a\nMAAT\n")
    with pytest.raises(ValueError, match="conflicting protein"):
        load_suite(root)


def test_legacy_workbook_is_unassessed_and_broken_dimension_is_ignored(tmp_path):
    import openpyxl
    import zipfile
    root = native(tmp_path)
    (root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv").unlink()
    path = root / "dnmb_module_rebasefinder/R-M_REBASE_analysis.xlsx"
    w = openpyxl.Workbook()
    s = w.active
    s.title = "RM_Comprehensive"
    s.append(["locus_tag", "passed_blast_filter", "rm_type", "subunit", "blast_match", "blast_identity"])
    s.append(["g:a", True, "Type II", "R", "R.Ref", 0.8])
    s.append(["g_b", False, "Type II", "M", None, None])
    w.save(path)
    with zipfile.ZipFile(path) as z:
        contents = {name: z.read(name) for name in z.namelist()}
    contents["xl/worksheets/sheet1.xml"] = contents["xl/worksheets/sheet1.xml"].replace(b'ref="A1:F3"', b'ref="A1"')
    with zipfile.ZipFile(path, "w") as z:
        for name, content in contents.items():
            z.writestr(name, content)
    data = load_suite(root)
    assert len(data["saved_decisions"]) == 2
    assert {d["state"] for d in data["saved_decisions"]} == {"unassessed"}
    assert all(d["family"] is None for d in data["saved_decisions"])


def test_reference_recognition_cannot_promote_candidate(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "rec_seq", "reference_rec_seq", "hit_label", "family_id"], [["g:a", "ACTG", "ACTG", "Annotation", "Type II"]])
    data = load_suite(root)
    assert data["saved_decisions"][0]["state"] == "unassessed"
    assert not [c for c in data["calls"] if c["caller"] == "REBASE"]


def test_live_tsv_preferred_to_redundant_workbook(tmp_path):
    root = native(tmp_path)
    (root / "dnmb_module_rebasefinder/R-M_REBASE_analysis.xlsx").write_bytes(b"unread redundant workbook")
    data = load_suite(root)
    assert len(data["saved_decisions"]) == 3
    assert not any(str(p).endswith("R-M_REBASE_analysis.xlsx") for p in data["inputs"])


def test_system_membership_disagreement_reported(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_defensefinder/defensefinder_systems.tsv", ["sys_id", "protein_in_syst"], [["sys A", "native1"]])
    data = load_suite(root)
    assert any(c["code"] == "genes_systems_membership_mismatch" for c in data["suite"]["conflicts"])
    assert len([c for c in data["calls"] if c["caller"] == "DefenseFinder" and c["evidence_kind"] == "final_member"]) == 2


def test_missing_map_allowed_only_when_native_ids_exactly_bind_known_genes(tmp_path):
    root = native(tmp_path)
    (root / "dnmb_module_defensefinder/defensefinder_id_map.tsv").unlink()
    with pytest.raises(ValueError, match="exact gene binding"):
        load_suite(root)


def test_manifest_arbitrary_cohort_order_and_explicit_threshold(tmp_path):
    native(tmp_path, "folder-one")
    native(tmp_path, "folder-two")
    config = tmp_path / "manifest.json"
    config.write_text(json.dumps({"schema_version": 1, "samples": [
        {"sample_id": "sample/Z", "path": "folder-two", "label": "Second", "cohort": "external", "highlight": True,
         "modules": {"defensepredictor": {"score_threshold": 9}}},
        {"sample_id": "sample/A", "path": "folder-one", "label": "First"}]}))
    data = load_suite(tmp_path, manifest=config)
    assert [s["sample_id"] for s in data["samples"]] == ["sample/Z", "sample/A"]
    assert data["samples"][0]["cohort"] == "external"
    assert not [c for c in data["calls"] if c["sample_id"] == "sample/Z" and c["caller"] == "DefensePredictor"]
    assert len(data["genes"]) == 6


@pytest.mark.parametrize("score", ["nan", "inf", "bad", ""])
def test_invalid_dp_scores_are_errors(tmp_path, score):
    root = native(tmp_path)
    write_table(root / "dnmb_module_defensepredictor/defense_predictor_output.csv", ["locus_tag", "mean_log_odds"], [["g:a", score]])
    with pytest.raises(ValueError, match="score|numeric"):
        load_suite(root)


def test_native_domain_and_motif_rows_preserved_without_rescore(tmp_path):
    root = native(tmp_path)
    path = root / "dnmb_module_interproscan/interproscan_results.tsv"
    path.parent.mkdir()
    path.write_text("g:a\tmd5\t4\tPfam\tPFexample\tNuclease\t1\t4\t1e-5\tT\tdate\t-\t-\t-\t-\n")
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_motif_hits_raw.tsv", ["locus_tag", "motif", "start_aa", "end_aa"], [["g:a", "pattern", 1, 3]])
    data = load_suite(root)
    assert any(p["namespace"] == "InterProScan" and p["native_fields"]["signature_accession"] == "PFexample" for p in data["suite"]["profiles"])
    assert data["suite"]["motifs"][0]["native_fields"]["start_aa"] == "1"
    assert data["effector_evidence"] == []
    assert data["saved_decisions"][0]["score_total"] == 7


def test_unknown_manifest_and_ambiguous_sources_fail(tmp_path):
    root = native(tmp_path)
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps({"schema_version": 1, "samples": [{"sample_id": "s", "path": root.name}], "typo": True}))
    with pytest.raises(ValueError, match="unknown suite manifest"):
        load_suite(tmp_path, manifest=path)
    (root / "one.faa").write_text(">g:a\nMAAA\n")
    (root / "two.faa").write_text(">g:a\nMAAA\n")
    with pytest.raises(ValueError, match="ambiguous root protein"):
        load_suite(root)


def test_reference_selected_is_separate_from_native_or_saved_retention(tmp_path):
    from defenseviz.analysis.saved import analyze_saved
    root = native(tmp_path)
    write_table(root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv",
                ["system.number", "seqid", "system", "target.name", "hmm.name"], [[1, "chr A", "RM_Type_II", "protein__safeA", "restriction_hmm"]])
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "curation_keep", "final_family", "final_role", "family_id", "hit_label", "blast_identity"],
                [["g_c", "FALSE", "Type II", "M", "Type II", "M.Reference", 0.7]])
    data = load_suite(root)
    reference = next(c for c in data["calls"] if c["caller"] == "REBASE")
    assert reference["selected"] is True and reference["evidence_kind"] == "reference_hit"
    assert reference["reference_id"] == "M.Reference"
    assert reference["category"] == "R-M" and reference["rm_type"] == "Type II"
    result = analyze_saved(data)
    assert next(x for x in result["source_sets"] if x["source_id"] == "REBASE")["count"] == 1
    assert not any(x["gene_id"] == "g_c" for x in result["gene_ledger"])
    # Direct source membership must not turn a not-retained reference candidate
    # into a native final call or an explicitly retained DNMB decision.
    assert not [c for c in data["calls"] if c["gene_id"] == "g_c" and c["evidence_kind"] == "final_member"]
    assert data["saved_decisions"][0]["state"] == "not_retained"


@pytest.mark.parametrize("parent,subtype,fqn,category,rm_type", [
    ("RM", "RM_Type_IV_1", "defense-finder-models/RM/RM/RM_Type_IV_1", "R-M", "Type IV"),
    ("RM", "Unspecified_RM_model", "", "R-M", None),
    ("RM_Type_II", "DNA_met_RM_Type_II_PDDEXK", "defense-finder-models/RM/RM_Type_II/DNA_met_RM_Type_II_PDDEXK", "R-M", "Type II"),
    ("NA", "N6_RM_Type_I", "defense-finder-models/RM/RM_Type_I/N6_RM_Type_I", "R-M", "Type I"),
    ("NA", "N6_N4_RM_Type_III", "defense-finder-models/RM/RM_Type_III/N6_N4_RM_Type_III", "R-M", "Type III"),
    ("RM_Type_IIG", "N6_RM_Type_IIG", "defense-finder-models/RM/RM_Type_IIG/N6_RM_Type_IIG", "R-M", "Type II"),
    ("NA", "CAS_Class1-Subtype-IV-A", "CasFinder/CAS_Class1-Subtype-IV-A", "CRISPR-Cas", None),
    ("NA", "CAS_Class1-Subtype-I-G", "CasFinder/CAS_Class1-Subtype-I-G", "CRISPR-Cas", None),
    ("CAS", "Unspecified_CAS_model", "", "CRISPR-Cas", None),
    ("Other", "RM_Type_IV_1", "Other/RM_Type_IV_1", "Unresolved", None),
    ("Other", "CAS_Class1-Subtype-IV-A", "Other/CAS_Class1-Subtype-IV-A", "Unresolved", None),
    ("Anti_RM", "unrelated_RM_Type_IV_like", "Other/Anti_RM/unrelated_RM_Type_IV_like", "Unresolved", None),
    ("CAS_like", "unrelated", "Other/CasFinder/unrelated", "Unresolved", None),
])
def test_defensefinder_native_taxonomy_uses_exact_parent_and_namespace(tmp_path, parent, subtype, fqn, category, rm_type):
    root = native(tmp_path)
    columns = ["hit_id", "sys_id", "type", "subtype", "model_fqn"]
    values = ["native1", "sys A", parent, subtype, fqn]
    write_table(root / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv", columns, [values])
    write_table(root / "dnmb_module_defensefinder/defensefinder_systems.tsv", ["sys_id", "protein_in_syst"], [["sys A", "native1"]])
    data = load_suite(root)
    call = next(c for c in data["calls"] if c["caller"] == "DefenseFinder" and c["evidence_kind"] == "final_member")
    assert (call["category"], call["rm_type"]) == (category, rm_type)
    assert call["native_family"] == subtype and call["native_type"] == parent
    assert call["native_fields"] == dict(zip(columns, values))
    assert call["role"] is None and data["effector_evidence"] == []
    if subtype == "RM_Type_IV_1" and parent == "RM":
        assert call["mapping_rule"] == "defensefinder_2.0.2_exact_rm_subtype_v1"
        assert "/2.0.2/" in call["mapping_reference"]
        assert call["mapping_inputs"]["model_fqn"] == fqn
        assert call["model_version"] is None  # Mapping reference never fabricates producer metadata.


@pytest.mark.parametrize("parent,subtype,fqn", [
    ("RM", "RM_Type_IV_1", "defense-finder-models/RM/RM/RM_Type_IV"),
    ("RM_Type_I", "N6_N4_RM_Type_III", "defense-finder-models/RM/RM_Type_III/N6_N4_RM_Type_III"),
    ("RM", "CAS_Class1-Subtype-IV-A", "CasFinder/CAS_Class1-Subtype-IV-A"),
])
def test_defensefinder_rejects_contradictory_native_taxonomy(tmp_path, parent, subtype, fqn):
    root = native(tmp_path)
    write_table(root / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv",
                ["hit_id", "sys_id", "type", "subtype", "model_fqn"], [["native1", "sys A", parent, subtype, fqn]])
    with pytest.raises(ValueError, match="taxonomy|model_fqn"):
        load_suite(root)


@pytest.mark.parametrize("family", ["RM_Type_IIG", "RM_type_IIG", "Type IIG"])
def test_native_iig_is_type_ii_without_assigning_a_role(tmp_path, family):
    root = native(tmp_path)
    write_table(root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv",
                ["system.number", "seqid", "system", "target.name"], [[1, "chr A", family, "protein__safeA"]])
    data = load_suite(root)
    call = next(c for c in data["calls"] if c["caller"] == "PADLOC")
    assert (call["category"], call["rm_type"], call["role"]) == ("R-M", "Type II", None)
    assert call["native_family"] == family


def test_literal_negative_strand_and_source_binding_roles(tmp_path):
    data = load_suite(native(tmp_path))
    assert next(g for g in data["genes"] if g["gene_id"] == "g_b")["strand"] == "-"
    roles = data["input_roles"]
    assert all(not r["native_result"] for p, rs in roles.items() if p.endswith("id_map.tsv") for r in rs)
    assert all(r["native_result"] for p, rs in roles.items() if p.endswith("best_solution_genes.tsv") for r in rs)


def test_duplicate_candidate_rows_preserve_observations_not_duplicate_genes(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "family_id", "hit_label", "blast_identity"],
                [["g:a", "Type II", "R.Ref", 0.8], ["g:a", "Type II", "R.Ref", 0.7]])
    data = load_suite(root)
    assert len(data["suite"]["reference_matches"]) == 2
    assert len(data["saved_decisions"]) == 1
    reference = [c for c in data["calls"] if c["caller"] == "REBASE"]
    assert len(reference) == 1 and len(reference[0]["source_refs"]) == 2


def test_unbound_candidate_is_quarantined_and_partial(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "curation_keep", "final_family", "final_role"], [["not an identified gene", "NA", "NA", "NA"]])
    data = load_suite(root)
    assert not data["saved_decisions"]
    assert len(data["suite"]["unmapped_records"]) == 1
    assert next(c for c in data["coverage"] if c["caller"] == "REBASE")["status"] == "partial"
    assert not any(g["gene_id"] == "not an identified gene" for g in data["genes"])
    with pytest.raises(ValueError, match="unbound_native_candidate"):
        load_suite(root, strict_inputs=True)


def test_legacy_unquoted_newline_falls_back_to_native_workbook(tmp_path):
    import openpyxl
    root = native(tmp_path)
    path = root / "dnmb_module_rebasefinder/rebasefinder_input.tsv"
    path.write_text("locus_tag\ttranslation\tproduct\ng:a\tMAAA\tunquoted\nnote\n")
    w = openpyxl.Workbook()
    s = w.active
    s.title = "1.GenBank_table"
    s.append(["locus_tag", "translation", "product"])
    s.append(["g:a", "MAAA", "unquoted\nnote"])
    s.append(["g_b", "MBBB", "other"])
    s.append(["g_c", "MCCC", "other"])
    w.save(root / "example_total.xlsx")
    data = load_suite(root)
    assert any(c["code"] == "malformed_gene_tsv_fallback" for c in data["suite"]["conflicts"])
    assert len(data["genes"]) == 3
    assert any(not s["selected"] for s in data["suite"]["sources"] if s["path"].endswith("rebasefinder_input.tsv"))


def test_retention_zero_sentinel_does_not_prove_reference_zero(tmp_path, monkeypatch):
    root = tmp_path / "empty"
    (root / "dnmb_module_rebasefinder").mkdir(parents=True)
    (root / ".dnmb_module_stage_cache.rds").write_bytes(b"fixture")
    payload = {"reader_schema_version": 1, "sentinels": {"REBASEfinder": {"alias": "REBASEfinder", "sentinel_version": 1, "n_hits": 0}},
               "cache": {"modules": {"REBASEfinder": {"ok": True, "table_counts": {"hits": 0, "all_hits": 9}}}}}
    monkeypatch.setattr(subprocess, "run", lambda argv, **kw: subprocess.CompletedProcess(argv, 0, json.dumps(payload), ""))
    data = load_suite(root, rscript="fixed-R")
    assert next(c for c in data["coverage"] if c["caller"] == "REBASE")["status"] == "partial"


def test_equal_sequences_do_not_merge_loci_or_samples(tmp_path):
    root = native(tmp_path)
    p = root / "dnmb_module_rebasefinder/rebasefinder_input.tsv"
    p.write_text(p.read_text().replace("MBBB", "MAAA"))
    data = load_suite(root)
    a, b = [next(g for g in data["genes"] if g["gene_id"] == key) for key in ["g:a", "g_b"]]
    assert a["protein_sha256"] == b["protein_sha256"]
    assert a["gene_id"] != b["gene_id"] and len(data["genes"]) == 3


def test_reference_subject_and_measurements_stay_coherent(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "blast_match", "blast_identity", "curated_blast_match", "curated_blast_identity", "family_id"],
                [["g:a", "Original", 0, "Curated_without_identity", "NA", "Type II"],
                 ["g_b", "Original_without_identity", "NA", "Curated", 0.7, "Type II"],
                 ["g_c", "NA", "NA", "Candidate_without_identity", "NA", "Type II"]])
    data = load_suite(root)
    references = {c["gene_id"]: c for c in data["calls"] if c["caller"] == "REBASE"}
    assert references["g:a"]["reference_id"] == "Original"
    assert references["g:a"]["match_identity"] == 0
    assert references["g:a"]["match_origin"] == "legacy_blast"
    assert references["g_b"]["reference_id"] == "Curated"
    assert references["g_b"]["match_identity"] == 0.7
    assert "g_c" not in references
    assert {c["reference_scope"] for c in references.values()} == {"exported_reference_observations"}


def test_bundled_resources_are_software_provenance_not_native_inputs(tmp_path):
    data = load_suite(native(tmp_path))
    assert not any(p.name == "families.json" for p in data["inputs"])
    assert any(s["resource"] == "package:defenseviz/data/families.json" and len(s["sha256"]) == 64 for s in data["suite"]["software_resources"])


def eggnog_export(root, text):
    path = root / "dnmb_module_eggnog/eggnog_out.emapper.annotations"
    path.parent.mkdir(exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def test_eggnog_preserves_header_driven_annotations_without_role_inference(tmp_path):
    root = native(tmp_path)
    path = eggnog_export(root, "## emapper test-version\n#query\tDescription\tseed_ortholog\tevalue\tPFAMs\n"
                         "g:a\tPutative nuclease\tortholog-A\t0\tPFexample\n"
                         "g_c\tUnknown\tortholog-C\t1e-10\t-\n## 2 queries processed\n")
    data = load_suite(root)
    rows = [r for r in data["suite"]["annotations"] if r.get("namespace") == "EggNOG"]
    assert [r["gene_id"] for r in rows] == ["g:a", "g_c"]
    assert rows[0]["native_fields"]["evalue"] == "0"
    assert rows[0]["native_fields"]["Description"] == "Putative nuclease"
    assert rows[0]["source_ref"]["row"] == 3
    assert rows[0]["source_ref"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    source = next(s for s in data["suite"]["sources"] if s["module"] == "eggnog")
    assert source["row_count"] == 2 and len(source["comments"]) == 2
    assert data["input_roles"][str(path.resolve())] == [dict(sample_id=root.name, module="eggnog", kind="functional_annotations", native_result=True)]
    assert data["annotations"] == data["effector_evidence"] == []


def test_eggnog_unbound_query_is_quarantined_without_prefix_guessing(tmp_path):
    root = native(tmp_path)
    eggnog_export(root, "#query\tseed_ortholog\nprotein__safeA\tortholog\ng:a\tortholog\n")
    data = load_suite(root)
    rows = [r for r in data["suite"]["annotations"] if r.get("namespace") == "EggNOG"]
    assert [r["gene_id"] for r in rows] == ["g:a"]
    assert any(r["module"] == "eggnog" and r["native_fields"]["#query"] == "protein__safeA" for r in data["suite"]["unmapped_records"])
    capability = next(c for c in data["suite"]["capabilities"] if "raw_eggnog_rows" in c)
    assert capability["raw_eggnog_status"] == "partial" and capability["raw_eggnog_unbound_rows"] == 1
    with pytest.raises(ValueError, match="unbound_eggnog_query"):
        load_suite(root, strict_inputs=True)


@pytest.mark.parametrize("values", [("first", "second"), ("second", "first")])
def test_eggnog_conflicting_queries_preserve_all_rows_and_report_conflict(tmp_path, values):
    root = native(tmp_path)
    eggnog_export(root, "#query\tseed_ortholog\n" + "".join(f"g:a\t{value}\n" for value in values))
    data = load_suite(root)
    rows = [r for r in data["suite"]["annotations"] if r.get("namespace") == "EggNOG"]
    assert {r["native_fields"]["seed_ortholog"] for r in rows} == set(values)
    assert any(c["code"] == "conflicting_eggnog_query" for c in data["suite"]["conflicts"])
    assert data["effector_evidence"] == []
    with pytest.raises(ValueError, match="conflicting_eggnog_query"):
        load_suite(root, strict_inputs=True)


@pytest.mark.parametrize("text", ["g:a\tortholog\n", "#query\tseed_ortholog\ng:a\n",
                                  "#query\tseed_ortholog\tseed_ortholog\n", "#query\tDescription\n"])
def test_eggnog_rejects_ambiguous_or_malformed_headers(tmp_path, text):
    root = native(tmp_path)
    eggnog_export(root, text)
    with pytest.raises(ValueError):
        load_suite(root)


def test_eggnog_complete_empty_and_missing_are_distinct(tmp_path):
    root = native(tmp_path)
    missing = load_suite(root)
    assert any(c.get("raw_eggnog_status") == "missing" for c in missing["suite"]["capabilities"])
    eggnog_export(root, "## emapper version\n#query\tseed_ortholog\n## 0 queries processed\n")
    data = load_suite(root)
    capability = next(c for c in data["suite"]["capabilities"] if "raw_eggnog_rows" in c)
    assert capability["raw_eggnog_status"] == "complete" and capability["raw_eggnog_rows"] == 0


def test_native_replicon_binding_rejects_contradictory_gene_unit(tmp_path):
    root = native(tmp_path)
    path = root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv"
    path.write_text(path.read_text().replace('1,plasmid,RM_Type_II,proteinC', '1,plasmid,RM_Type_II,protein__safeA'))
    with pytest.raises(ValueError, match="replicon contradicts"):
        load_suite(root)


def test_unknown_replicon_alias_is_not_guessed_equal_or_different(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_padloc/padloc_query_proteins_id_map.tsv", ["query_id", "locus_tag"], [["protein__safeA", "g:a"]])
    write_table(root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv", ["system.number", "seqid", "system", "target.name"], [[1, "unmapped_alias", "Novel", "protein__safeA"]])
    data = load_suite(root)
    group = next(g for g in data["groups"] if g["provider"] == "PADLOC")
    assert group["replicon_binding"] == "unassessed"


def test_declared_sample_unit_cannot_be_overridden_by_native_row(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_defensefinder/defensefinder_id_map.tsv", ["sample_id", "defensefinder_id", "locus_tag"], [["other_sample", "native1", "g:a"]])
    with pytest.raises(ValueError, match="sample unit"):
        load_suite(root)


def test_native_call_coordinates_must_match_bound_gene(tmp_path):
    root = native(tmp_path)
    write_table(root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv", ["system.number", "seqid", "system", "target.name", "start", "end", "strand"],
                [[1, "chr A", "Novel", "protein__safeA", 22, 33, "+"]])
    with pytest.raises(ValueError, match="conflicting start"):
        load_suite(root)


@pytest.mark.parametrize("states", [("complete", "missing"), ("missing", "complete")])
def test_duplicate_candidate_context_conflicts_never_use_row_order(tmp_path, states):
    root = native(tmp_path)
    write_table(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
                ["query", "curation_keep", "final_family", "final_role", "typeii_context_status"],
                [["g:a", "TRUE", "Type II", "R", state] for state in states])
    with pytest.raises(ValueError, match="contradictory saved context payloads"):
        load_suite(root)

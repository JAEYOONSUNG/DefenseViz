"""Automatic cohort discovery uses native boundaries and source metadata."""
import hashlib
import json
from pathlib import Path

import pytest

from defenseviz.io import suite
from defenseviz.io.suite_discovery import discover_suite


def native(path, assembly=None, organism="Example bacterium", strain="S-1"):
    module = path / "dnmb_module_rebasefinder"
    module.mkdir(parents=True)
    (module / "rebasefinder_input.tsv").write_text(
        "locus_tag\tcontig\tstart\tend\tdirection\ttranslation\tproduct\n"
        "g1\tctg\t1\t12\t+\tMAAA\tunknown\n", encoding="utf-8")
    if assembly:
        (path / "source.gbff").write_text(genbank(assembly, organism, strain), encoding="utf-8")
    return path


def genbank(assembly, organism, strain):
    return ("LOCUS       EXAMPLE                  12 bp    DNA     linear   BCT 01-JAN-2000\n"
            "DEFINITION  Deliberately unrelated display title.\n"
            "ACCESSION   EXAMPLE\n"
            "DBLINK      BioProject: PRJNA000001\n"
            f"            Assembly: {assembly}\n"
            f"  ORGANISM  {organism}\n"
            "            Bacteria; Example.\n"
            "FEATURES             Location/Qualifiers\n"
            "     source          1..12\n"
            f'                     /organism="{organism}"\n'
            f'                     /strain="{strain}"\n'
            "     CDS             1..12\n"
            '                     /strain="not a source strain"\n'
            '                     /locus_tag="g1"\n'
            '                     /translation="MAAA"\n'
            "ORIGIN\n        1 atggctgctgct\n//\n")


def specs(root):
    return [row["spec"] for row in discover_suite(root)["samples"]]


def test_nested_flat_and_single_native_roots_with_source_identity(tmp_path):
    first = native(tmp_path / "group" / "deep" / "folder notes", "GCF_000000001.9")
    second = native(tmp_path / "flat arbitrary")
    found = specs(tmp_path)
    assert {row["sample_id"] for row in found} == {"GCF_000000001.9", second.name}
    assert specs(first) == [next(row for row in found if row["sample_id"].startswith("GCF_"))]
    assert {row["path"] for row in found} == {str(first), str(second)}
    # Deterministic independent of filesystem enumeration order.
    assert found == specs(tmp_path)


def test_accepted_native_root_stops_module_and_plasmid_descent(tmp_path):
    sample = native(tmp_path / "sample")
    native(sample / "plasmid" / "run", "GCF_000000002.1")
    native(sample / "dnmb_module_rebasefinder" / "nested", "GCF_000000003.1")
    assert [row["sample_id"] for row in specs(tmp_path)] == ["sample"]


def test_hidden_subtrees_are_excluded(tmp_path):
    native(tmp_path / "kept")
    ignored = native(tmp_path / ".hidden" / "deep" / "ignored")
    assert [row["sample_id"] for row in specs(tmp_path)] == ["kept"]
    # An explicitly selected native root is still a legitimate input.
    assert specs(ignored)[0]["sample_id"] == "ignored"


@pytest.mark.parametrize("name", ["output", "Outputs", "staging", "receipts", "reviews", "artifacts"])
def test_ordinary_output_folder_names_do_not_hide_native_samples(tmp_path, name):
    sample = native(tmp_path / name / "deep" / "sample")
    assert specs(tmp_path)[0]["path"] == str(sample)
    direct = native(tmp_path / "other" / name)
    assert len(specs(tmp_path)) == 2
    assert specs(direct)[0]["sample_id"] == name


def test_identified_generated_output_is_pruned_but_native_root_takes_precedence(tmp_path):
    native(tmp_path / "kept")
    generated = tmp_path / "any name"
    native(generated / "copied input")
    (generated / "analysis.json").write_text("{}", encoding="utf-8")
    receipt = generated / "receipt.json"
    receipt.write_text('{"tool":"OtherTool"}', encoding="utf-8")
    assert len(specs(tmp_path)) == 2
    receipt.write_text('{"tool":"DefenseViz"}', encoding="utf-8")
    assert [row["sample_id"] for row in specs(tmp_path)] == ["kept"]
    (generated / "dnmb_module_padloc").mkdir()
    assert {row["sample_id"] for row in specs(tmp_path)} == {"kept", "any name"}


def test_symlink_aliases_and_ancestor_cycles_do_not_duplicate_samples(tmp_path):
    sample = native(tmp_path / "group" / "sample", "GCF_000000004.1")
    try:
        (tmp_path / "alias").symlink_to(sample, target_is_directory=True)
        (tmp_path / "group" / "cycle").symlink_to(tmp_path, target_is_directory=True)
    except OSError:
        pytest.skip("Directory symlinks unavailable for this user/platform")
    inventory = discover_suite(tmp_path)
    assert len(inventory["samples"]) == 1
    assert inventory["samples"][0]["directory"] == str(sample)
    assert any(row["reason"] == "duplicate_or_cycle" for row in inventory["skipped"])


def test_symlinks_cannot_promote_nested_plasmids_or_hidden_inputs(tmp_path):
    sample = native(tmp_path / "sample")
    plasmid = native(sample / "plasmid" / "nested")
    hidden = native(tmp_path / ".hidden" / "input")
    try:
        (tmp_path / "a-plasmid-alias").symlink_to(plasmid, target_is_directory=True)
        (tmp_path / "a-hidden-alias").symlink_to(hidden, target_is_directory=True)
    except OSError:
        pytest.skip("Directory symlinks unavailable for this user/platform")
    assert [row["sample_id"] for row in specs(tmp_path)] == ["sample"]


@pytest.mark.parametrize("assembly", [None, "GCF_000000005.1"])
def test_distinct_roots_with_same_automatic_identity_fail_honestly(tmp_path, assembly):
    native(tmp_path / "run1" / "same", assembly)
    native(tmp_path / "run2" / "same", assembly)
    with pytest.raises(ValueError, match="duplicate suite sample_id.*explicit manifest"):
        discover_suite(tmp_path)


def test_metadata_overrides_filename_notes_and_binds_hashes_to_import(tmp_path):
    sample = native(tmp_path / "GCF_999999999.1_wrong_taxon_manual_notes", "GCF_000000006.2",
                    "Example bacterium isolate-A", "isolate-A; alternate")
    source = sample / "source.gbff"
    data = suite.load_suite(tmp_path)
    row = data["samples"][0]
    assert row["sample_id"] == "GCF_000000006.2"
    assert row["label"] == "Example bacterium isolate-A"
    assert row["genus"] == "Example" and row["species"] == "bacterium"
    assert row["display"] == {"taxon_label": "Example bacterium", "strain_label": "isolate-A"}
    assert source in data["inputs"]
    assert data["input_sha256"][str(source)] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert data["suite"]["samples"] == [{"sample_id": row["sample_id"], "directory": str(sample)}]
    assert any(role["kind"] == "sample_metadata" for role in data["input_roles"][str(source)])
    assert all(gene["sample_id"] == row["sample_id"] for gene in data["genes"])


def test_wrapped_source_name_and_multirecord_genbank_keep_one_sample(tmp_path):
    sample = native(tmp_path / "arbitrary", "GCF_000000007.1", "Example species subsp. subgroup strain X/Y", "X/Y")
    source = sample / "source.gbff"
    record = source.read_text(encoding="utf-8").replace(
        '/organism="Example species subsp. subgroup strain X/Y"',
        '/organism="Example species subsp.\n                     subgroup strain X/Y"')
    source.write_text(record + record, encoding="utf-8")
    row = specs(tmp_path)[0]
    assert row["label"] == "Example species subsp. subgroup X/Y"
    assert row["species"] == "species subsp. subgroup"
    assert row["display"]["strain_label"] == "X/Y"


@pytest.mark.parametrize("declaration", ["Alpha; Beta", "Alpha = Beta"])
@pytest.mark.parametrize("separator", [" = ", "; "])
def test_complete_declared_alias_chain_is_removed_without_losing_provenance(tmp_path, declaration, separator):
    organism = "Example species subsp. subgroup strain Alpha" + separator + "Beta"
    native(tmp_path / "arbitrary", "GCF_000000020.1", organism, declaration)
    item = discover_suite(tmp_path)["samples"][0]
    assert item["spec"]["display"] == {
        "taxon_label": "Example species subsp. subgroup", "strain_label": "Alpha"}
    qualifiers = item["metadata_sources"][0]["records"][0]["source"]
    assert qualifiers["organism"] == [organism]
    assert qualifiers["strain"] == [declaration]


@pytest.mark.parametrize("alias", ["ABC 123", "XYZ456"])
def test_explicit_organism_equivalence_is_anchored_by_declared_terminal_strain(tmp_path, alias):
    organism = "Example species subsp. subgroup " + alias + " = QRS 789"
    sample = native(tmp_path / "misleading folder", "GCF_000000021.1", organism, "QRS 789")
    chromosome = (sample / "source.gbff").read_text()
    # Source records may use either declared equivalent name or omit strain.
    plasmid = genbank("GCF_000000021.1", "Example species subsp. subgroup QRS 789", "QRS 789")
    other = genbank("GCF_000000021.1", "Example species subsp. subgroup " + alias, "").replace(
        '                     /strain=""\n', "")
    taxon_only = genbank("GCF_000000021.1", "Example species subsp. subgroup", "").replace(
        '                     /strain=""\n', "")
    (sample / "source.gbff").write_text(chromosome + plasmid + other + taxon_only)
    item = discover_suite(tmp_path)["samples"][0]
    assert item["spec"]["display"] == {
        "taxon_label": "Example species subsp. subgroup", "strain_label": alias}
    assert item["spec"]["species"] == "species subsp. subgroup"
    records = item["metadata_sources"][0]["records"]
    assert records[0]["source"]["organism"] == [organism]
    assert records[0]["source"]["strain"] == ["QRS 789"]
    assert "strain" not in records[2]["source"]


@pytest.mark.parametrize("organism", [
    "Example species unknown = QRS 789",
    "Example species ABC 123 = unknown = QRS 789",
])
def test_unknown_alias_chain_is_not_partially_stripped(tmp_path, organism):
    from defenseviz.io.suite_discovery import _organism_strain_suffix
    assert _organism_strain_suffix(organism, ["QRS 789"]) == (organism, [])
    native(tmp_path / "sample", "GCF_000000022.1", organism, "QRS 789")
    with pytest.raises(ValueError, match="ambiguous.*alias chain"):
        discover_suite(tmp_path)


def test_collection_shaped_species_text_is_not_removed_without_equivalence(tmp_path):
    native(tmp_path / "sample", "GCF_000000023.1", "Example species group ABC 123", "QRS 789")
    row = specs(tmp_path)[0]
    assert row["display"] == {"taxon_label": "Example species group ABC 123", "strain_label": "QRS 789"}


def test_organism_equivalence_requires_an_exact_declared_terminal_anchor():
    from defenseviz.io.suite_discovery import _organism_strain_suffix
    organism = "Example species ABC 123 = QRS 789"
    assert _organism_strain_suffix(organism, ["Other-1"]) == (organism, [])
    assert _organism_strain_suffix(organism, ["QRS 78"]) == (organism, [])


def test_explicit_alias_chain_does_not_reconcile_conflicting_strain(tmp_path):
    sample = native(tmp_path / "sample", "GCF_000000024.1", "Example species ABC 123 = QRS 789", "QRS 789")
    with (sample / "source.gbff").open("a") as stream:
        stream.write(genbank("GCF_000000024.1", "Example species", "Different-1"))
    with pytest.raises(ValueError, match="conflicting GenBank"):
        discover_suite(tmp_path)


def test_assembly_mentions_in_comments_are_not_sample_identities(tmp_path):
    sample = native(tmp_path / "arbitrary")
    text = genbank("", "Example bacterium", "S-1").replace(
        "FEATURES", "COMMENT     Comparison to GCF_000000009.1\nFEATURES")
    (sample / "source.gbff").write_text(text, encoding="utf-8")
    assert specs(tmp_path)[0]["sample_id"] == "arbitrary"


def test_missing_strain_on_a_contig_uses_only_exact_declared_suffix(tmp_path):
    sample = native(tmp_path / "sample", "GCF_000000014.1", "Example bacterium S-1", "S-1")
    source = sample / "source.gbff"
    record = source.read_text(encoding="utf-8")
    source.write_text(record + record.replace('                     /strain="S-1"\n', ""), encoding="utf-8")
    assert specs(tmp_path)[0]["label"] == "Example bacterium S-1"


def test_taxon_only_plasmid_inherits_only_same_taxon_unambiguous_strain(tmp_path):
    sample = native(tmp_path / "sample", "GCF_000000015.1", "Example bacterium S-1", "S-1")
    source = sample / "source.gbff"
    chromosome = source.read_text(encoding="utf-8")
    plasmid = genbank("GCF_000000015.1", "Example bacterium", "").replace('                     /strain=""\n', "")
    source.write_text(chromosome + plasmid, encoding="utf-8")
    item = discover_suite(tmp_path)["samples"][0]
    assert item["spec"]["label"] == "Example bacterium S-1"
    assert not item["metadata_sources"][0]["records"][1]["source"].get("strain")
    source.write_text(chromosome + plasmid.replace("Example bacterium", "Different taxon"), encoding="utf-8")
    with pytest.raises(ValueError, match="conflicting GenBank"):
        discover_suite(tmp_path)


@pytest.mark.parametrize("change", ["assembly", "taxon", "strain"])
def test_conflicting_multirecord_metadata_requires_explicit_selection(tmp_path, change):
    sample = native(tmp_path / "arbitrary", "GCF_000000008.1")
    other = genbank("GCF_000000010.1" if change == "assembly" else "GCF_000000008.1",
                    "Other species" if change == "taxon" else "Example bacterium",
                    "S-2" if change == "strain" else "S-1")
    with (sample / "source.gbff").open("a", encoding="utf-8") as stream:
        stream.write(other)
    with pytest.raises(ValueError, match="conflicting GenBank"):
        specs(tmp_path)


def test_explicit_manifest_preserves_custom_ids_labels_and_order(tmp_path):
    a = native(tmp_path / "first", "GCF_000000011.1")
    b = native(tmp_path / "second", "GCF_000000011.1")
    manifest = tmp_path / "selected.json"
    manifest.write_text(json.dumps({"schema_version": 1, "samples": [
        {"sample_id": "custom-B", "path": str(b), "label": "Manual B"},
        {"sample_id": "custom-A", "path": str(a), "label": "Manual A"}]}), encoding="utf-8")
    rows = suite.load_suite(tmp_path, manifest=manifest)["samples"]
    assert [(row["sample_id"], row["label"]) for row in rows] == [("custom-B", "Manual B"), ("custom-A", "Manual A")]


def test_source_mutation_between_discovery_and_import_is_rejected(tmp_path, monkeypatch):
    sample = native(tmp_path / "sample", "GCF_000000012.1")
    def mutate(root):
        inventory = discover_suite(root)
        with (sample / "source.gbff").open("a", encoding="utf-8") as stream:
            stream.write("\n")
        return inventory
    monkeypatch.setattr(suite, "discover_suite", mutate)
    with pytest.raises(ValueError, match="GenBank metadata changed since discovery"):
        suite.load_suite(tmp_path)


def test_old_accession_identity_is_stable_after_arbitrary_folder_rename(tmp_path):
    sample = native(tmp_path / "old decorated name", "GCA_000000013.3")
    before = suite.load_suite(sample)
    renamed = sample.rename(tmp_path / "new arbitrary name")
    after = suite.load_suite(renamed)
    assert before["samples"] == after["samples"]
    assert before["genes"] == after["genes"]


def test_empty_parent_does_not_silently_publish_zero_samples(tmp_path):
    with pytest.raises(ValueError, match="no native suite folders"):
        discover_suite(tmp_path)


def test_broken_legacy_annotation_fragments_never_become_genes_or_reference_hits(tmp_path):
    import openpyxl
    from defenseviz.io.suite_cache import attach_cached_evidence

    sample = native(tmp_path / "nested" / "native", "GCF_000000016.1")
    module = sample / "dnmb_module_rebasefinder"
    (module / "rebasefinder_input.tsv").write_text(
        "locus_tag\ttranslation\tproduct\ng1\tMAAA\tannotation\nfragment,with,commas\n", encoding="utf-8")
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "1.GenBank_table"
    sheet.append(["locus_tag", "translation", "product"])
    sheet.append(["g1", "MAAA", "annotation\nfragment,with,commas"])
    workbook.save(sample / "native_total.xlsx")
    (module / "DNMB_REBASEfinder_augmented_hits.tsv").write_text(
        "query\tfamily_id\thit_label\tblast_identity\n"
        "g1\tType II\tR.reference\t0.9\n"
        "annotation_fragment_1\tType II\tR.reference\t0.9\n", encoding="utf-8")
    (module / "query_sequences.fasta").write_text(
        ">g1\nMAAA\n>annotation_fragment_1\nCONSENSUS,ANNOTATION,TEXT\n", encoding="utf-8")
    (module / "blast_results.txt").write_text("".join(
        f"{query}\tR.reference\t90\t3\t0\t0\t1\t3\t1\t3\t1e-8\t40\n"
        for query in ("g1", "annotation_fragment_1")), encoding="utf-8")
    data = suite.load_suite(tmp_path)
    assert {g["gene_id"] for g in data["genes"]} == {"g1"}
    assert {c["gene_id"] for c in data["calls"]} == {"g1"}
    assert {d["gene_id"] for d in data["saved_decisions"]} == {"g1"}
    assert data["suite"]["unmapped_records"][0]["native_fields"]["query"] == "annotation_fragment_1"
    enriched = attach_cached_evidence(data, tmp_path)
    assert {r["gene_id"] for r in enriched["suite"]["cached_blast"]} == {"g1"}
    omissions = enriched["suite"]["native_cache_evidence"]["samples"][0]["omitted"]
    assert any(row["gene_id"] == "annotation_fragment_1" and row["reason"] == "saved_query_protein_not_bound" for row in omissions)

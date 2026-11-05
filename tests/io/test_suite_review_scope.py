"""Pinned raw evidence scope and binding, without scientific role inference."""
import csv
import hashlib
import json
from copy import deepcopy
import pytest
from defenseviz.io.suite import SuiteReader, attach_system_review_annotations, load_suite
from tests.io.test_suite import native, write_table


def domain(gene="g:a", sequence="MAAA", *, md5=None, length=4, start=1, end=4):
    return [gene, md5 if md5 is not None else hashlib.md5(sequence.encode()).hexdigest(), str(length),
            "Pfam", "PFsynthetic", "Synthetic domain", str(start), str(end), "1e-12", "T", "date", "-", "-", "-", "-"]


def interpro(root, rows, *, sites=False):
    path = root / "dnmb_module_interproscan" / ("interproscan_results.tsv.sites" if sites else "interproscan_results.tsv")
    path.parent.mkdir(exist_ok=True)
    path.write_text("".join("\t".join(row) + "\n" for row in rows))
    return path


def rewrite(path, mutate):
    sep = "," if path.suffix == ".csv" else "\t"
    with path.open() as f:
        reader = csv.DictReader(f, delimiter=sep)
        headers, rows = list(reader.fieldnames), list(reader)
    headers, rows = mutate(headers, rows)
    return write_table(path, headers, [[r.get(h, "") for h in headers] for r in rows])


def observed_fixture(tmp_path):
    root = native(tmp_path)
    # g_c remains only in native HMM evidence, absent from final/candidate sets.
    rewrite(root / "dnmb_module_padloc/padloc_query_proteins_padloc.csv",
            lambda h,rs: (h,[r for r in rs if r["target.name"] != "proteinC"]))
    rewrite(root / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
            lambda h,rs: (h,[r for r in rs if r["query"] != "g_c"]))
    rewrite(root / "dnmb_module_rebasefinder/rebasefinder_input.tsv", lambda h,rs:(h,rs+[
        dict(locus_tag="review_only",contig="chr B",start="1",end="12",direction="+",translation="MDDD",product="unknown")]))
    return root


def test_all_call_scope_includes_hmm_only_and_source_row_accounting(tmp_path):
    root = observed_fixture(tmp_path)
    interpro(root, [domain(), domain("g_c", "MCCC"), domain("review_only", "MDDD"), domain("unknown")])
    data = load_suite(root)
    profiles = [p for p in data["suite"]["profiles"] if p["namespace"] == "InterProScan"]
    assert {p["gene_id"] for p in profiles} == {"g:a", "g_c"}
    assert all(p["sequence_binding"] == "verified" and p["qualification"] == "unassessed" for p in profiles)
    assert data["effector_evidence"] == []
    c = next(c for c in data["suite"]["capabilities"] if "raw_interpro_parsed_rows" in c)
    assert (c["raw_interpro_parsed_rows"],c["raw_interpro_retained_rows"],c["raw_interpro_excluded_scope_rows"],c["raw_interpro_quarantined_rows"]) == (4,2,1,1)
    assert any(r.get("reason") == "unbound_interpro_query" for r in data["suite"]["unmapped_records"])
    assert c["raw_interpro_status"] == "partial"


def test_post_supplement_attachment_preserves_calls_and_deduplicates(tmp_path):
    root = observed_fixture(tmp_path)
    interpro(root, [domain(), domain("review_only", "MDDD")])
    data = load_suite(root)
    original = deepcopy(data["calls"])
    data["annotations"].append(dict(sample_id=root.name, gene_id="review_only", status="unassessed"))
    assert attach_system_review_annotations(data, root) is data
    first = deepcopy(data["suite"])
    attach_system_review_annotations(data, root)
    assert data["suite"] == first
    assert data["calls"] == original and data["effector_evidence"] == []
    rows = [p for p in data["suite"]["profiles"] if p["namespace"] == "InterProScan"]
    assert {p["gene_id"] for p in rows} == {"g:a", "review_only"}
    c = next(c for c in data["suite"]["capabilities"] if "raw_interpro_parsed_rows" in c)
    assert c["raw_interpro_retained_rows"] == 2 and c["raw_interpro_excluded_scope_rows"] == 0


@pytest.mark.parametrize("row,reason", [
    (domain(md5="0"*32),"interpro_sequence_md5_mismatch"),
    (domain(md5="not-md5"),"invalid_interpro_sequence_md5"),
    (domain(length=5),"interpro_sequence_length_mismatch"),
    (domain(start=0),"invalid_interpro_coordinates"),
    (domain(end=5),"invalid_interpro_coordinates"),
])
def test_stale_or_invalid_domains_are_quarantined_not_qualified(tmp_path,row,reason):
    root = native(tmp_path); interpro(root,[row])
    data = load_suite(root)
    assert not [r for r in data["suite"]["profiles"] if r["namespace"] == "InterProScan"]
    assert any(r.get("reason") == reason and r["native_fields"]["signature_accession"] == "PFsynthetic" for r in data["suite"]["unmapped_records"])
    assert data["effector_evidence"] == []
    with pytest.raises(ValueError,match=reason): load_suite(root,strict_inputs=True)


def test_missing_sequence_proof_remains_explicitly_unassessed(tmp_path):
    root = native(tmp_path); interpro(root,[domain(md5="-")])
    data = load_suite(root)
    r = next(r for r in data["suite"]["profiles"] if r["namespace"] == "InterProScan")
    assert r["sequence_binding"] == "unassessed" and r["qualification"] == "unassessed"


def test_post_attachment_honors_strict_input_mode(tmp_path):
    root = native(tmp_path)
    interpro(root, [domain(md5="0" * 32)])
    data = load_suite(root)
    with pytest.raises(ValueError, match="interpro_sequence_md5_mismatch"):
        attach_system_review_annotations(data, root, strict_inputs=True)


def test_site_records_validate_site_coordinates_and_never_assert_motif_function(tmp_path):
    root = native(tmp_path)
    row=["g:a",hashlib.md5(b"MAAA").hexdigest(),"4","SFLD","synthetic","1","4","1","M","1","1","site annotation"]
    bad=list(row);bad[9]="5";bad[10]="5"
    interpro(root,[row,bad],sites=True)
    data=load_suite(root)
    r=next(r for r in data["suite"]["profiles"] if r["namespace"]=="InterProScan")
    assert r["evidence_type"]=="site_annotation" and r["qualification"]=="unassessed"
    assert len([r for r in data["suite"]["unmapped_records"] if r.get("reason")=="invalid_interpro_coordinates"])==1
    assert data["effector_evidence"]==[]


def test_post_attachment_rejects_changed_source_and_relative_paths_bind(tmp_path):
    root=native(tmp_path);path=interpro(root,[domain()]);data=load_suite(root)
    source=next(s for s in data["suite"]["sources"] if s["module"]=="interproscan")
    source["resolved_path"]=str(path.relative_to(root))
    attach_system_review_annotations(data,root)
    path.write_text(path.read_text()+"\n")
    with pytest.raises(ValueError,match="pinned hash"):attach_system_review_annotations(data,root)


def test_exact_native_row_version_lift_and_conflicting_saved_version(tmp_path,monkeypatch):
    root=native(tmp_path)
    p=root/"dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv"
    rewrite(p,lambda h,rs:(h+["model_version"],[dict(r,model_version="synthetic-1") for r in rs]))
    data=load_suite(root)
    c=next(c for c in data["calls"] if c["caller"]=="DefenseFinder" and c["evidence_kind"]=="final_member")
    assert c["model_version"]=="synthetic-1" and c["model_version_source"]=="native_row.model_version"
    assert c["source_ref"]["sha256"]==hashlib.sha256(p.read_bytes()).hexdigest()
    monkeypatch.setattr(SuiteReader,"metadata",lambda *args:{"cache":{"signature":{"db_state":{"DefenseFinder":{"models_version":"different"}}}}})
    with pytest.raises(ValueError,match="conflicts with saved module metadata"):load_suite(root)


def test_manifest_sample_assessment_declaration_provenance(tmp_path):
    root=native(tmp_path)
    manifest=tmp_path/"manifest.json"
    manifest.write_text(json.dumps(dict(schema_version=1,samples=[dict(sample_id="bound",path=root.name)],
        options=dict(assessment_coverage={"motif":"partial"},assessment_coverage_by_sample={"bound":{"effector":"complete"}}))))
    data=load_suite(tmp_path,manifest=manifest)
    declarations=[r for r in data["source_metadata"] if r["source_kind"]=="assessment_coverage_declaration"]
    assert {r["declaration_scope"] for r in declarations}=={"global","sample"}
    assert all(r["sample_id"]=="bound" and r["source_ref"]["sha256"]==hashlib.sha256(manifest.read_bytes()).hexdigest() for r in declarations)


def test_domain_missing_source_and_explicit_empty_are_distinct(tmp_path):
    root=native(tmp_path)
    assert any(c.get("raw_interpro_status")=="missing" for c in load_suite(root)["suite"]["capabilities"])
    interpro(root,[])
    c=next(c for c in load_suite(root)["suite"]["capabilities"] if "raw_interpro_parsed_rows" in c)
    assert c["raw_interpro_status"]=="complete" and c["raw_interpro_parsed_rows"]==0


def test_duplicate_native_version_contradiction_rejected(tmp_path):
    root=native(tmp_path)
    p=root/"dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv"
    def versions(h,rows):
        out=[dict(r,model_version="v1") for r in rows]
        return h+["model_version"],out+[dict(out[0],model_version="v2")]
    rewrite(p,versions)
    with pytest.raises(ValueError,match="contradictory duplicate native model_version"):load_suite(root)


def test_database_version_fallback_remains_bound_to_native_row(tmp_path):
    root=native(tmp_path)
    p=root/"dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv"
    rewrite(p,lambda h,rs:(h+["database_version"],[dict(r,database_version="db1") for r in rs]))
    data=load_suite(root)
    c=next(c for c in data["calls"] if c["caller"]=="DefenseFinder" and c["evidence_kind"]=="final_member")
    assert c["model_version"]==c["database_version"]=="db1"
    assert c["model_version_source"]=="native_row.database_version"


def test_old_serialized_bundle_without_md5_cannot_claim_sequence_verification(tmp_path):
    root=native(tmp_path);interpro(root,[domain()]);data=load_suite(root)
    for g in data["genes"]:g.pop("protein_md5",None)
    attach_system_review_annotations(data,root)
    r=next(r for r in data["suite"]["profiles"] if r["namespace"]=="InterProScan")
    assert r["sequence_binding"]=="unassessed"
    assert r["qualification"]=="unassessed" and data["effector_evidence"]==[]

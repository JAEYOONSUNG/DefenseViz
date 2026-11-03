import hashlib
from copy import deepcopy
import pytest
from defenseviz.analysis.registry import assess_candidate, assess_registry, digest
from defenseviz.analysis.admission import apply_reference_anchors, entry_sha256
from defenseviz.analysis import evaluate, parse_evidence, parse_policy
from .test_integration import evidence


@pytest.fixture
def data(tmp_path):
    sequence = "ADCEFGHIKLMN"
    source = tmp_path / "synthetic_evidence.txt"
    source.write_text("Synthetic verified source fixture; no real biological inference.\n")
    sha = hashlib.sha256(source.read_bytes()).hexdigest()
    model = dict(accession="TEST.1",source="synthetic",model_sha256="1"*64,
        evidence_kind="nuclease_domain",role="nuclease",family_specific=True,min_model_coverage=.7)
    ref = dict(accession="REF",sequence=sequence,sequence_sha256=digest(sequence),
        numbering_origin="1-based full reference sequence",independent_mapping_methods=["method_A","method_B"],
        catalytic_sites=[dict(position=2,residue="D",role="nuclease",evidence_source="synthetic"),
                         dict(position=4,residue="E",role="nuclease",evidence_source="synthetic")])
    entry = dict(family_id="Synthetic nuclease",component="nuclease",substrate="DNA",domain_models=[model],references=[ref])
    registry = dict(schema_version="1.0",families=[entry])
    hit = dict(accession="TEST.1",source="synthetic",model_sha256="1"*64,
        query_sequence_sha256=digest(sequence),qstart=1,qend=len(sequence),hmm_start=1,hmm_end=12,hmm_length=12,
        target_length=12,independent_evalue=1e-30,
        qualification=dict(passed=True,rule="synthetic_qualified",source=source.name,source_sha256=sha))
    candidate = dict(candidate_id="any id",sample_id="arbitrary",query_id="unseen_gene",sequence=sequence,
        sequence_sha256=digest(sequence),domain_hits=[hit],site_maps=[dict(reference_accession="REF",reference_sequence_sha256=digest(sequence),
        query_sequence_sha256=digest(sequence),method=m,positions={"2":2,"4":4},source=source.name,source_sha256=sha) for m in ["method_A","method_B"]])
    rule=dict(family_id=entry["family_id"],RM_families=["Type II"],roles=["R"],entry_sha256=entry_sha256(entry),
        model_accession="TEST.1",model_sha256="1"*64,reference_accession="REF",required_anchor_subset=[2,4],source=source.name,source_sha256=sha)
    return registry,candidate,rule,tmp_path


def test_qualified_registry_anchor_is_not_activity_or_geometry(data):
    reg,c,rule,root=data
    out=assess_candidate(c,reg,root)
    f=out["families"][0]
    assert f["family_evidence"]["state"] == f["domain_evidence"]["state"] == f["catalytic_site_evidence"]["state"] == "supported"
    assert f["structure_evidence"]["state"] == "unassessed"
    assert out["activity"] == f["biological_activity"] == "not_assessed"
    assert out["source_files_verified"] is True


def test_actual_source_tampering_rejected(data):
    reg,c,rule,root=data
    (root/"synthetic_evidence.txt").write_text("Changed")
    with pytest.raises(ValueError,match="hash mismatch"):
        assess_candidate(c,reg,root)


def test_source_query_or_model_mismatch_cannot_support(data):
    reg,c,rule,root=data
    c["domain_hits"][0]["query_sequence_sha256"]="2"*64
    f=assess_candidate(c,reg,root)["families"][0]
    assert f["domain_evidence"]["state"] == "unassessed"
    c["domain_hits"][0]["model_sha256"]="3"*64
    assert assess_candidate(c,reg,root)["families"] == []


def test_independent_mapping_disagreement_is_conflict(data):
    reg,c,rule,root=data
    c["site_maps"][1]["positions"]["4"]=5
    f=assess_candidate(c,reg,root)["families"][0]
    assert f["catalytic_site_evidence"]["state"] == "conflict"


def test_missing_mapping_is_unassessed_not_negative(data):
    reg,c,rule,root=data
    c["site_maps"]=c["site_maps"][:1]
    f=assess_candidate(c,reg,root)["families"][0]
    assert f["catalytic_site_evidence"]["state"] == "unassessed"


def test_domains_cannot_stitch_sites_across_instances(data):
    reg,c,rule,root=data
    reg["families"][0]["domain_models"][0]["min_model_coverage"]=0
    h=c["domain_hits"][0];h["qend"]=3
    h2=deepcopy(h);h2["qstart"]=4;h2["qend"]=12
    c["domain_hits"].append(h2)
    f=assess_candidate(c,reg,root)["families"][0]
    assert f["catalytic_site_evidence"]["state"] == "conflict"
    assert "required_sites_stitched_across_distinct_domain_instances" in f["catalytic_site_evidence"]["references"][0]["failures"]


def test_unknown_registry_family_is_unassessed(data):
    reg,c,rule,root=data
    c["domain_hits"]=[]
    assert assess_candidate(c,reg,root)["disposition"].startswith("unassessed")


def test_rna_reference_not_dna_activity(data):
    reg,c,rule,root=data
    reg["families"][0]["substrate"]="RNA"
    f=assess_candidate(c,reg,root)["families"][0]
    assert f["DNA_catalytic_site_state"] == "not_applicable"


def test_registry_batch_duplicates_rejected(data):
    reg,c,rule,root=data
    with pytest.raises(ValueError,match="Duplicate"):
        assess_registry(reg,[c,c],registry_dir=root)


def records(c):
    e=evidence(old_context=True,curated_identity=.9,curated_query_coverage=.9,curated_reference_coverage=.9,
        curated_evalue=1e-20,curated_role_compatible=True,curated_family_compatible=True)
    return [dict(sample_id=c["sample_id"],gene_id=c["query_id"],evidence=e,provenance={"sequence_sha256":c["sequence_sha256"]})]


def test_generic_source_bound_admission_changes_only_eligible_evidence(data):
    reg,c,rule,root=data
    r=records(c)
    out,audit=apply_reference_anchors(r,[c],reg,[rule],registry_dir=root)
    assert "reference_nuclease_anchor_support" not in r[0]["evidence"]  # input preserved
    assert out[0]["evidence"]["reference_nuclease_anchor_profile_indices"] == [0]
    d=evaluate(parse_evidence(out[0]["evidence"]),parse_policy())
    assert d["reference_nuclease_anchor_motif_credit"]
    assert d["score_motif"] == 20
    assert d["score_structure"] == 0
    assert audit["records"][0]["status"] == "admitted"
    assert audit["thresholds_changed"] is False


def test_admission_requires_exact_sample_gene_and_sequence(data):
    reg,c,rule,root=data
    r=records(c);r[0]["provenance"]["sequence_sha256"]="f"*64
    with pytest.raises(ValueError,match="sequence mismatch"):
        apply_reference_anchors(r,[c],reg,[rule],registry_dir=root)
    r=records(c);r[0]["sample_id"]="different"
    out,audit=apply_reference_anchors(r,[c],reg,[rule],registry_dir=root)
    assert audit["records"][0]["status"] == "outside_evidence_pool"
    assert "reference_nuclease_anchor_support" not in out[0]["evidence"]


def test_locus_selector_and_unpinned_entry_forbidden(data):
    reg,c,rule,root=data
    rule["gene_id"]="unseen_gene"
    with pytest.raises(ValueError,match="whitelists"):
        apply_reference_anchors(records(c),[c],reg,[rule],registry_dir=root)
    del rule["gene_id"];rule["entry_sha256"]="0"*64
    with pytest.raises(ValueError,match="exact registered"):
        apply_reference_anchors(records(c),[c],reg,[rule],registry_dir=root)


def test_member_duplication_cannot_create_complete_module(data):
    reg,c,rule,root=data
    c["modules"]=[dict(family_id="Synthetic nuclease",members=[dict(query_id="same",protein_sha256="1"*64,component="A"),dict(query_id="same",protein_sha256="1"*64,component="B")])]
    with pytest.raises(ValueError,match="multiple components"):
        assess_candidate(c,reg,root)


def test_empty_or_unmatched_batch_still_validates_all_registry_and_rules(data):
    reg,c,rule,root=data
    bad=deepcopy(reg);bad["schema_version"]="unsupported"
    with pytest.raises(ValueError,match="schema"):
        apply_reference_anchors([],[],bad,[],registry_dir=root)
    with pytest.raises(ValueError,match="whitelists"):
        apply_reference_anchors([],[],reg,[dict(rule,gene_id="hidden")],registry_dir=root)
    bad_rule=dict(rule,reference_accession="not_registered")
    with pytest.raises(ValueError,match="reference"):
        apply_reference_anchors([],[],reg,[bad_rule],registry_dir=root)
    malformed=dict(sample_id="outside",gene_id="pool",candidate_id="broken")
    with pytest.raises(ValueError,match="candidate"):
        apply_reference_anchors([], [malformed], reg, [rule], registry_dir=root)


def test_external_reference_file_is_bound_even_empty_batch(data):
    reg,c,rule,root=data
    ref=reg["families"][0]["references"][0]
    (root/"reference.faa").write_text(">synthetic\n"+ref.pop("sequence")+"\n")
    ref["sequence_file"]="reference.faa"
    rule["entry_sha256"]=entry_sha256(reg["families"][0])
    expected=hashlib.sha256((root/"reference.faa").read_bytes()).hexdigest()
    assessed=assess_candidate(c,reg,root)
    assert assessed["source_sha256"][str((root/"reference.faa").resolve())]==expected
    _,audit=apply_reference_anchors([],[],reg,[rule],registry_dir=root)
    assert audit["source_sha256"][str((root/"reference.faa").resolve())]==expected


def test_repeated_admission_is_idempotent_and_linked_profile_tampering_rejected(data):
    reg,c,rule,root=data
    first,audit1=apply_reference_anchors(records(c),[c],reg,[rule],registry_dir=root)
    second,audit2=apply_reference_anchors(first,[c],reg,[rule],registry_dir=root)
    assert second==first
    assert audit2==audit1
    assert len(second[0]["evidence"]["profiles"])==1
    bad=deepcopy(first);bad[0]["evidence"]["profiles"][0]["evalue"]=1e-50
    with pytest.raises(ValueError,match="actual linked profile"):
        apply_reference_anchors(bad,[c],reg,[rule],registry_dir=root)

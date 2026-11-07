"""Generic component-architecture review and counterexamples, using no cohort IDs."""
from copy import deepcopy
from hashlib import md5, sha256

import pytest

from defenseviz.automatic_typeiv import apply_typeiv_component_review
from defenseviz.analysis.saved import analyze_saved


REF = dict(path="native.tsv", sha256="a" * 64, row=2)
IPR = dict(path="interproscan_results.tsv", sha256="b" * 64, row=1)
CALL = dict(path="native_members.tsv", sha256="c" * 64, row=2)
MAP = dict(path="identity.tsv", sha256="d" * 64, row=2)
SEQ = "M" + "A" * 799
SHA = sha256(SEQ.encode()).hexdigest()
MD5 = md5(SEQ.encode()).hexdigest()


def fixture(domains=("PF04471",), *, reference=True, native=None):
    gene = dict(sample_id="sample", gene_id="candidate", protein_sha256=SHA,
                protein_md5=MD5, protein_length=len(SEQ), source_refs=[REF])
    evidence = dict(family="Type IV", role="R", valid_sequence=True,
                    alternative_annotation=False, conflict=False, partial=False, hmm_role_ambiguous=False,
                    curated_role_compatible=reference, curated_family_compatible=reference)
    current = dict(**gene, state="not_retained", assignment_status="evaluated", evidence=evidence,
                   decision=dict(numeric_reference_pass=reference, score_total=35, keep=False))
    features = dict(sample_id="sample", gene_id="candidate", protein_sha256=SHA, source_ref=REF,
        native_fields=dict(query="candidate", curated_blast_family="Type IV", curated_blast_role="R",
            curated_blast_match="synthetic_reference", curated_blast_role_compatible=reference,
            curated_blast_family_compatible=reference))
    data = dict(samples=[dict(sample_id="sample", label="Synthetic isolate")], genes=[gene], calls=[], annotations=[],
        saved_decisions=[dict(**gene, state="not_retained", family="Type IV", role="R")],
        effector_evidence=[], groups=[], complexes=[], coverage=[],
        suite=dict(profiles=[], features=[features], identity_maps=[], sources=[dict(sample_id="sample", **r) for r in (REF, IPR, CALL, MAP)]),
        rm_evaluation=dict(records=[current]))
    set_domains(data, domains)
    if native:
        add_native(data, native)
    return data


def set_domains(data, domains):
    rows = []
    for i, accession in enumerate(domains):
        analysis = "PANTHER" if accession.startswith("PTHR") else "Pfam"
        rows.append(dict(sample_id="sample", gene_id="candidate", protein_sha256=SHA, namespace="InterProScan",
            source_ref={**IPR, "row": i + 1}, evidence_type="domain_hit", sequence_binding="verified", coordinate_binding="verified",
            native_fields=dict(protein_id="candidate", sequence_md5=MD5, protein_length=str(len(SEQ)),
                analysis=analysis, signature_accession=accession, start=str(10 + i * 120), end=str(99 + i * 120), status="T")))
    data["suite"]["profiles"] = rows


def add_native(data, caller="DefenseFinder", *, span=(1, 800)):
    native_id = "renamed_protein"
    fields = (dict(hit_id=native_id, subtype="RM_Type_IV", hit_seq_len="800", hit_begin_match=str(span[0]), hit_end_match=str(span[1]))
              if caller == "DefenseFinder" else {"target.name": native_id, "system": "RM_type_IV", "start": "1", "end": "800"})
    data["calls"].append(dict(sample_id="sample", gene_id="candidate", protein_sha256=SHA, caller=caller,
        family="RM_Type_IV", rm_type="Type IV", evidence_kind="final_member", selected=True,
        source_ref=CALL, native_fields=fields, system_id="native_system"))
    data["suite"]["identity_maps"].append(dict(sample_id="sample", gene_id="candidate", namespace=caller.lower(),
        source_id=native_id, source_ref=MAP))


def review(data):
    return apply_typeiv_component_review(data)["automatic_typeiv_review"]["records"][0]


@pytest.mark.parametrize("domains,rule,role", [
    (("PF12102", "PF07728"), "McrB_architecture", "Aux"),
    (("PF10117",), "McrC_domain", "R"),
    (("PF04471",), "Mrr_domain", "R"),
    (("PF03235", "PF07510"), "GmrSD_fused_architecture", "R"),
    (("PF09823", "PF04411"), "McrC_like_DUF2357_PDDEXK", "R"),
    (("PF07728", "PTHR37291:SF17"), "McrB_family_motor_variable_recognition_domain", "Aux"),
    (("PF13091", "PF04851", "PF00271"), "PLD_SF2_TypeIV_associated_architecture", "R"),
    (("PF01844",), "HNH_native_TypeIV_associated_component", "R"),
])
def test_all_eight_architectures(domains, rule, role):
    row = review(fixture(domains, native="DefenseFinder"))
    assert row["annotation_supported"] and row["rule_id"] == rule and row["role"] == role
    assert row["complete_complex"] == "not_inferred"
    assert row["catalytic_residue_validation"] == "not_performed"


def test_annotation_is_idempotent_and_preserves_strict_policy_and_all_other_evidence():
    data = fixture()
    before = deepcopy(data)
    after = apply_typeiv_component_review(data)
    assert data == before
    assert apply_typeiv_component_review(after) == after
    for key in data:
        if key != "annotations":
            assert after[key] == data[key]
    annotation = after["annotations"][0]
    assert annotation["role"] == "R" and annotation["include_in_union"]
    assert annotation["protocol_sha256"] == after["automatic_typeiv_review"]["protocol_sha256"]
    assert after["saved_decisions"][0]["state"] == "not_retained"
    assert after["rm_evaluation"]["records"][0]["decision"]["score_total"] == 35
    result = analyze_saved(after)
    assert result["global_summary"]["retained_CDS"] == 0
    assert result["global_summary"]["reviewed_CDS"] == 1
    assert result["global_summary"]["effector_CDS"] == 0
    assert result["global_summary"]["priority_support_qualified_CDS"] == 1


@pytest.mark.parametrize("field,value", [
    ("sequence_md5", "0" * 32), ("protein_length", "801"), ("protein_id", "another_gene"),
    ("start", "0"), ("end", "801"), ("start", "100.5"), ("end", "1"),
])
def test_domain_binding_and_coordinates_quarantine(field, value):
    data = fixture()
    data["suite"]["profiles"][0]["native_fields"][field] = value
    assert review(data)["status"] == "quarantined"


@pytest.mark.parametrize("field,value", [
    ("protein_sha256", "e" * 64), ("sequence_binding", "unassessed"), ("coordinate_binding", "reported_sequence_only")
])
def test_stale_protein_and_unbound_domain_quarantine(field, value):
    data = fixture()
    data["suite"]["profiles"][0][field] = value
    assert review(data)["status"] == "quarantined"


def test_hash_must_match_importer_pinned_inventory():
    data = fixture()
    data["suite"]["profiles"][0]["source_ref"]["sha256"] = "f" * 64
    assert review(data)["status"] == "quarantined"


@pytest.mark.parametrize("domains", [(), ("PF04411",), ("PF07728",), ("PF07510",), ("PF13091",), ("PF07510", "PF03235")])
def test_generic_or_incomplete_or_reversed_architectures_do_not_assign(domains):
    assert not review(fixture(domains))["annotation_supported"]


def test_overlapping_domains_and_failed_hits_are_not_complete_fusions():
    data = fixture(("PF03235", "PF07510"))
    data["suite"]["profiles"][1]["native_fields"]["start"] = "99"
    assert not review(data)["annotation_supported"]
    data = fixture()
    data["suite"]["profiles"][0]["native_fields"]["status"] = "F"
    assert not review(data)["annotation_supported"]


def test_specific_original_architecture_precedes_broader_same_family_rule():
    row = review(fixture(("PF12102", "PF07728", "PTHR37291:SF1")))
    assert row["annotation_supported"] and row["rule_id"] == "McrB_architecture"
    assert len(row["observed_architectures"]) == 1
    assert row["role"] == "Aux"


def test_distinct_component_families_require_explicit_review():
    row = review(fixture(("PF10117", "PF04471")))
    assert not row["annotation_supported"]
    assert "multiple_role_architectures_require_review" in row["reasons"]


def test_mrr_domain_requires_typeiv_corroboration_not_product_text():
    data = fixture(reference=False)
    data["genes"][0]["product"] = "Mrr Type IV nuclease PF04471"
    assert not review(data)["annotation_supported"]
    add_native(data)
    assert review(data)["annotation_supported"]
    data["rm_evaluation"]["records"][0]["evidence"]["family"] = "Type II"
    assert apply_typeiv_component_review(data)["automatic_typeiv_review"]["records"] == []


@pytest.mark.parametrize("domains", [("PF07728", "PTHR37291"), ("PF13091", "PF04851", "PF00271")])
def test_extended_reference_required_routes_cannot_use_only_native_membership(domains):
    assert not review(fixture(domains, reference=False, native="DefenseFinder"))["annotation_supported"]


def test_hnh_needs_exact_native_protein_domain_overlap_not_reference_or_genomic_coordinates():
    assert not review(fixture(("PF01844",)))["annotation_supported"]
    assert not review(fixture(("PF01844",), native="PADLOC"))["annotation_supported"]
    data = fixture(("PF01844",), reference=False)
    add_native(data, span=(100, 500))
    assert not review(data)["annotation_supported"]
    data["calls"][0]["native_fields"].update(hit_begin_match="10", hit_end_match="99")
    assert review(data)["annotation_supported"]


@pytest.mark.parametrize("model,fqn", [
    ("RM_Type_IV_1", "defense-finder-models/RM/RM/RM_Type_IV_1"),
    ("McrBC", "defense-finder-models/RM/McrBC/McrBC"),
    ("GmrSD_RM_Type_IV", "defense-finder-models/RM/GmrSD_RM_Type_IV/GmrSD_RM_Type_IV"),
    ("EcoKMcrA", "defense-finder-models/DefenseFinder/EcoKMcrA/EcoKMcrA"),
])
def test_exact_native_model_taxonomy_aliases_are_preserved(model, fqn):
    data = fixture(("PF01844",), reference=False, native="DefenseFinder")
    data["calls"][0]["native_fields"].update(subtype=model, model_fqn=fqn)
    assert review(data)["annotation_supported"]


@pytest.mark.parametrize("mutation", ["map", "native_id", "protein", "length", "span", "source", "type"])
def test_native_binding_failure_never_rescues_hnh(mutation):
    data = fixture(("PF01844",), native="DefenseFinder")
    call = data["calls"][0]
    if mutation == "map": data["suite"]["identity_maps"][0]["gene_id"] = "other"
    if mutation == "native_id": call["native_fields"]["hit_id"] = "other"
    if mutation == "protein": call["protein_sha256"] = "e" * 64
    if mutation == "length": call["native_fields"]["hit_seq_len"] = "801"
    if mutation == "span": call["native_fields"]["hit_end_match"] = "801"
    if mutation == "source": call["source_ref"] = {**CALL, "sha256": "e" * 64}
    if mutation == "type": call["native_fields"]["subtype"] = "RM_Type_I"
    assert review(data)["status"] == "quarantined"


@pytest.mark.parametrize("hold", ["alternative_annotation", "conflict", "partial", "hmm_role_ambiguous"])
def test_positive_or_unknown_source_holds_are_preserved(hold):
    data = fixture()
    data["rm_evaluation"]["records"][0]["evidence"][hold] = True
    assert not review(data)["annotation_supported"]
    data["rm_evaluation"]["records"][0]["evidence"][hold] = None
    assert not review(data)["annotation_supported"]


@pytest.mark.parametrize("retained_role", ["R", "M"])
def test_mcrb_aux_conflicting_retained_assignments_quarantine(retained_role):
    data = fixture(("PF12102", "PF07728"))
    data["saved_decisions"][0].update(state="retained", role=retained_role)
    assert review(data)["status"] == "quarantined"


def test_mcrb_aux_existing_noncatalytic_review_does_not_create_r_or_complete_complex():
    data = fixture(("PF12102", "PF07728"))
    data["annotations"].append(dict(**data["genes"][0], status="supported", assignment_kind="noncatalytic_component", rule_id="existing", role="Aux"))
    result = apply_typeiv_component_review(data)
    assert review(data)["role"] == "Aux"
    assert result["groups"] == [] and result["complexes"] == [] and result["effector_evidence"] == []
    assert result["annotations"][-1]["complete_complex"] == "not_inferred"


def test_existing_review_contradiction_is_never_silently_overridden():
    data = fixture()
    data["annotations"].append(dict(**data["genes"][0], status="supported", assignment_kind="reviewed_RM_component",
                                    family="Type II", role="M", rule_id="explicit"))
    after = apply_typeiv_component_review(data)
    assert review(data)["status"] == "quarantined"
    assert after["annotations"] == data["annotations"]


def test_source_quarantine_does_not_become_a_positive_from_other_rows():
    data = fixture()
    data["suite"]["unmapped_records"] = [dict(sample_id="sample", gene_id="candidate", namespace="InterProScan", reason="sequence_mismatch")]
    assert review(data)["status"] == "quarantined"


def test_no_partner_is_inferred_from_supported_single_component():
    data = fixture(("PF10117",))
    result = apply_typeiv_component_review(data)
    assert review(data)["annotation_supported"]
    assert result["groups"] == [] and result["complexes"] == []
    assert result["automatic_typeiv_review"]["records"][0]["observed_architectures"][0]["required_partner"] == "McrB"


def test_verified_terminal_stop_mapping_keeps_source_domains_and_native_span_bound():
    data = fixture(("PF01844",), native="DefenseFinder")
    raw = SEQ + "*"
    source_sha = sha256(raw.encode()).hexdigest()
    for row in [data["genes"][0], data["saved_decisions"][0], data["rm_evaluation"]["records"][0], data["calls"][0],
                data["suite"]["features"][0], data["suite"]["profiles"][0]]:
        row["protein_sha256"] = source_sha
    data["genes"][0].update(protein_md5=md5(raw.encode()).hexdigest(), protein_length=len(raw))
    data["suite"]["profiles"][0]["native_fields"].update(sequence_md5=md5(raw.encode()).hexdigest(), protein_length=str(len(raw)))
    data["rm_evaluation"]["records"][0]["protein_binding"] = dict(source_sha256=source_sha, canonical_sha256=source_sha,
        hmm_sha256=SHA, normalization="terminal_stop_removal_only", source_receipt="exact_canonical_translation")
    data["suite"]["features"][0]["native_fields"].update(source_protein_sha256=source_sha, canonical_protein_sha256=source_sha,
        hmm_protein_sha256=SHA, sequence_length_aa=800)
    assert review(data)["annotation_supported"]
    data["rm_evaluation"]["records"][0]["protein_binding"]["source_receipt"] = "unverified_normalization"
    assert not review(data)["annotation_supported"]


@pytest.mark.parametrize("decision", [None, {}])
def test_missing_policy_measurement_is_unassessed_not_a_zero_or_reference_pass(decision):
    data = fixture()
    data["rm_evaluation"]["records"][0]["decision"] = decision
    row = review(data)
    assert row["status"] == "unassessed" and row["reference_gate"] is False
    assert row["strict_policy_score"] is None
    assert not apply_typeiv_component_review(data)["annotations"]


@pytest.mark.parametrize("missing", ["evidence", "suite", "rm_evaluation"])
def test_missing_legacy_evidence_sections_do_not_create_positive_assignments(missing):
    data = fixture()
    if missing == "evidence":
        data["rm_evaluation"]["records"][0].update(family="Type IV", evidence=None, source_refs=None)
    else:
        data[missing] = None
    result = apply_typeiv_component_review(data)
    assert not result["annotations"]
    assert not any(r["annotation_supported"] for r in result["automatic_typeiv_review"]["records"])

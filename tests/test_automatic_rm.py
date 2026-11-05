"""Native scientific predicates and source binding, with synthetic proteins."""
from copy import deepcopy
from hashlib import sha256

import pytest

from defenseviz.automatic_rm import compile_rm_evidence, _motifs, _profiles, _ordered, _component

REF = {"path": "native.tsv", "sha256": "a" * 64, "row": 2}


def bundle(sequence="FAGAGDPPY" + "A" * 150, family="Type II", role="M", **feature):
    digest = sha256(sequence.encode()).hexdigest()
    gene = dict(sample_id="sample", gene_id="gene", protein_sha256=digest,
                protein_length=len(sequence), source_refs=[REF])
    fields = dict(query="gene", family_id=family, enzyme_role=role,
                  partial_status="full_length", **feature)
    annotation = dict(sample_id="sample", gene_id="gene", source_ref=REF,
                      native_fields={"translation": sequence, "product": "DNA methyltransferase",
                                     "Signature description_Pfam": "DNA methylase"})
    return dict(samples=[{"sample_id": "sample"}], genes=[gene], saved_decisions=[dict(sample_id="sample", gene_id="gene", source_ref=REF)],
                calls=[], suite={"features": [dict(sample_id="sample", gene_id="gene", source_ref=REF, native_fields=fields)],
                "annotations": [annotation], "profiles": [], "motifs": [], "contexts": [], "sources": []})


def record(data):
    return compile_rm_evidence(data)["records"][0]


def test_fresh_bound_sequence_features_are_not_detector_rerun_or_mutation():
    data = bundle()
    before = deepcopy(data)
    result = compile_rm_evidence(data)
    e = result["records"][0]["evidence"]
    assert e["legacy_mtase_motif"] is True and e["old_motif"] is True
    assert e["legacy_nuclease"] is False and e["legacy_motor"] is False
    assert result["motif_scans"] == 1 and result["engine_runs"] == 0
    assert data == before


@pytest.mark.parametrize("seq,text,expected", [
    ("HANAAAAAH", "hypothetical protein", False),
    ("HANAAAAAH", "HNH PF01844", True),
    ("HAKAAAD", "PLDc_2 PF13091", False),
    ("HAKAAADAAAAAHAKAAAD", "PLDc_2 PF13091", True),
    ("HAKAAADAAAAAHAKAAAD", "hypothetical protein", False),
    ("PDAAAAADAK", "DNA repair nuclease RecB", False),
    ("PDAAAAADAK", "restriction endonuclease", True),
])
def test_original_domain_gated_nuclease_rules(seq, text, expected):
    values, _ = _motifs([], seq, text, True)
    assert values["legacy_nuclease"] is expected


def test_motif_export_position_windows_do_not_change_original_feature_rule():
    seq = "PDAAAAADAK" + "A" * 180
    values, _ = _motifs([], seq, "restriction endonuclease", True)
    assert values["legacy_nuclease"] is True
    # Old Type III plotting windows cannot veto the actual native sequence rule.
    assert "PD-ExK" in values["raw_nuclease_classes"]


def test_methyltransferase_requires_pair_and_plausible_length():
    assert _motifs([], "FAGAGDPPY", "", True)[0]["legacy_mtase_motif"] is False
    assert _motifs([], "FAGAGDPPY" + "A" * 150, "", True)[0]["legacy_mtase_motif"] is True
    assert _motifs([], "FAGAG" + "A" * 150, "", True)[0]["legacy_mtase_motif"] is False


def test_motor_architecture_order_and_inclusive_span():
    assert _ordered({"a": {1}, "b": {300}, "c": {751}}, ["a", "b", "c"], 750)
    assert not _ordered({"a": {1}, "b": {300}, "c": {752}}, ["a", "b", "c"], 750)
    assert not _ordered({"a": {30}, "b": {20}, "c": {10}}, ["a", "b", "c"], 750)


def test_absent_source_axes_are_unknown_not_false():
    r = record(bundle())
    assert r["evidence"]["strict_structure"] is None
    assert r["evidence"]["old_strong"] is None
    assert r["evidence"]["old_context"] is None
    assert "strict_structure" in r["unknown_axes"] and "profiles" in r["unknown_axes"]


def test_structure_missing_export_is_not_geometry_negative():
    r = record(bundle(structure_supported="FALSE", structure_status="structure_missing"))
    assert r["evidence"]["strict_structure"] is None


def test_specificity_motif_is_not_applicable():
    r = record(bundle(family="Type I", role="S"))
    assert r["evidence"]["old_motif"] is False
    assert r["evidence"]["divergent_context"] is False


def test_type_iv_never_requires_methylation_partner():
    data = bundle(family="Type IV", role="R")
    data["suite"]["annotations"][0]["native_fields"]["product"] = "modification-dependent restriction endonuclease"
    r = record(data)
    assert r["family"] == "Type IV" and r["role"] == "R"
    assert r["evidence"]["old_context"] is None  # missing export, no artificial partner requirement
    assert r["evidence"]["divergent_context"] is False


def test_partial_can_remain_unknown_when_reference_rescue_missing():
    data = bundle()
    data["suite"]["features"][0]["native_fields"]["partial_status"] = "partial_or_short"
    r = record(data)
    assert r["evidence"]["partial"] is None
    data["suite"]["features"][0]["native_fields"]["partial_reason"] = "annotation suggests pseudogene"
    assert record(data)["evidence"]["partial"] is True


def test_explicit_stored_zero_factors_remain_false():
    r = record(bundle(factor_strong="FALSE", factor_context="FALSE", numeric_structure_admissible_structural_evidence="FALSE"))
    for key in ["old_strong", "old_context", "strict_structure"]:
        assert r["evidence"][key] is False


def test_conflicting_protein_annotation_blocks_assignment_and_is_not_pooled():
    data = bundle()
    data["suite"]["annotations"].append(dict(sample_id="sample", gene_id="gene", source_ref=REF,
        native_fields={"translation": "AAAA", "product": "RNA methyltransferase"}))
    r = record(data)
    assert r["assignment_status"] == "unresolved"
    assert r["protein_binding_status"] == "unassessed"
    assert "annotation_protein_hash_mismatch" in r["reason_codes"]
    assert r["evidence"]["alternative_annotation"] is False


def test_annotation_append_without_translation_cannot_erase_exact_binding():
    data = bundle()
    data["suite"]["annotations"].append(dict(sample_id="sample", gene_id="gene", source_ref=REF,
        native_fields={"Signature accession_Pfam": "PF01555"}))
    r = record(data)
    assert r["evidence"]["valid_sequence"] is True
    assert r["evidence"]["legacy_mtase_motif"] is True


def test_reference_alignment_does_not_invent_subject_length_or_coverage():
    r = record(bundle(blast_identity="1", blast_length="159", blast_evalue="0", hit_label="M.Example"))
    assert r["evidence"]["identity"] == 1
    assert r["evidence"]["reference_coverage"] is None
    assert r["evidence"]["old_strong"] is None


def test_blank_candidate_hit_is_not_proof_of_a_completed_search():
    r = record(bundle(blast_identity="NA", hit_label="NA"))
    assert r["evidence"]["old_strong"] is None
    assert r["evidence"]["old_moderate"] is None


def test_native_rm_calls_absent_from_saved_candidates_stay_in_scope():
    data = bundle()
    data["saved_decisions"] = []
    data["suite"]["features"] = []
    data["calls"] = [dict(sample_id="sample", gene_id="gene", caller="PADLOC", evidence_kind="final_member", rm_type="Type II", native_model="MTase_II_00001", source_ref=REF)]
    result = compile_rm_evidence(data)
    r = result["records"][0]
    assert r["source_scope"] == ["original_tool_rm_members"]
    assert (r["family"], r["role"]) == ("Type II", "M")


def test_exact_component_grammar_rejects_generic_labels_and_unknown_profiles():
    assert _component("PADLOC", "REase_II_00001") == ("Type II", "R")
    assert _component("PADLOC", "looks_like_REase_II_00001") is None
    assert _component("DefenseFinder", "RM__Type_I_REases_FAM_999.einsi_trimmed") is None
    assert _component("DefenseFinder", "RM_Type_IV__Type_IV_03") is None


def test_comparator_profile_diagnostics_do_not_receive_native_hmm_credit():
    data = bundle()
    data["suite"]["profiles"] = [dict(sample_id="sample", gene_id="gene", namespace="DefenseFinder", source_ref=REF,
        native_fields={"gene_name": "RM_Type_II__Type_II_MTases_FAM_0", "hit_sequence_length": 159, "i_eval": "1e-40", "hit_seq_cov": ".95", "hit_profile_cov": ".95"})]
    r = record(data)
    assert r["evidence"]["profiles"] == []
    assert r["axis_sources"]["profiles"]["comparator_profile_details"]


def test_same_model_domains_union_without_cross_model_or_cross_protein_borrowing():
    gene = {"protein_length": 200}
    def row(model, a, b, c, d, qlen=200):
        return dict(namespace="PADLOC", source_ref=REF, native_fields=dict(query_name=model, target_name="protein", query_length="100", target_length=str(qlen), ali_from=a, ali_to=b, hmm_from=c, hmm_to=d, evalue="1e-20"))
    ps = _profiles([row("REase_II_00001", 1, 100, 1, 50), row("REase_II_00001", 80, 180, 40, 90), row("REase_II_00002", 1, 20, 1, 99), row("REase_II_00003", 1, 300, 1, 100, 300)], gene)
    assert len(ps) == 2
    a = next(p for p in ps if p["model"] == "REase_II_00001")
    assert a["query_coverage"] == .9 and a["model_coverage"] == .9
    b = next(p for p in ps if p["model"] == "REase_II_00002")
    assert b["query_coverage"] == .1 and b["model_coverage"] == .99


def test_reference_alignment_requires_operation_protein_and_database_binding():
    data = bundle(blast_identity=".9", blast_length="159", blast_evalue="1e-50", hit_label="M.Example")
    row = dict(sample_id="sample", gene_id="gene", source_ref=REF,
               operation="fresh_reference_alignment", command_exit_code=0, protein_sha256=data["genes"][0]["protein_sha256"],
               reference_sequence_sha256="b" * 64, reference_metadata_source_refs=[{"path": "references.faa", "sha256": "c" * 64}],
               native_fields=dict(subject_id="M.Example", identity=.9, evalue=1e-50, qstart=1,qend=159,sstart=1,send=159,
                                  reference_length=159,reference_family="Type II",reference_role="M"))
    data["suite"]["reference_alignments"] = [row]
    r = record(data)
    assert r["evidence"]["old_strong"] is True
    assert r["evidence"]["reference_coverage"] == 1
    assert any(ref["path"] == "references.faa" for ref in r["axis_sources"]["old_strong"]["source_refs"])
    for changed in ({"operation": "unspecified"}, {"protein_sha256": "d" * 64}, {"reference_sequence_sha256": None}, {"reference_metadata_source_refs": []}):
        altered = deepcopy(data)
        altered["suite"]["reference_alignments"][0].update(changed)
        assert record(altered)["evidence"]["reference_coverage"] is None


def test_measured_poor_alignment_prevents_strong_credit_without_invented_coverage():
    r = record(bundle(blast_identity=".99", blast_length="50", blast_evalue="1e-50", hit_label="M.Example"))
    assert r["evidence"]["reference_coverage"] is None
    assert r["evidence"]["query_coverage"] == pytest.approx(50/159)
    assert r["evidence"]["old_strong"] is False
    assert r["evidence"]["old_moderate"] is False


def test_completed_exact_fresh_reference_nohit_is_observation_not_activity():
    data = bundle(blast_identity="1", blast_length="159", blast_evalue="0", hit_label="M.Example")
    row = dict(sample_id="sample", gene_id="gene", source_ref=REF, operation="fresh_reference_alignment",
               alignment_status="no_alignment", command_exit_code=0,
               protein_sha256=data["genes"][0]["protein_sha256"], reference_sequence_sha256="b" * 64,
               reference_metadata_source_refs=[{"path": "references.faa", "sha256": "c" * 64}],
               native_fields={"subject_id": "M.Example", "reference_length": 159, "reference_family": "Type II", "reference_role": "M"})
    data["suite"]["reference_alignments"] = [row]
    r = record(data)
    assert r["evidence"]["old_strong"] is None and r["evidence"]["old_moderate"] is None
    assert r["reference_measurement_diagnostic"]["alignment_status"] == "no_alignment"
    assert r["evidence"]["identity"] is None
    assert r["biological_activity"] == "not_inferred"
    row["command_exit_code"] = 1
    assert record(data)["evidence"]["old_strong"] is None


def test_fresh_reference_cannot_borrow_old_quality_or_spans():
    data = bundle(blast_identity="1", blast_length="159", blast_evalue="0", hit_label="M.Example", blast_alignment_quality="strong", factor_strong="TRUE", blast_query_coverage="1", blast_reference_coverage="1", blast_role_compatible="TRUE", blast_family_compatible="TRUE")
    data["suite"]["reference_alignments"] = [dict(sample_id="sample", gene_id="gene", source_ref=REF,
        operation="fresh_reference_alignment", command_exit_code=0,
        protein_sha256=data["genes"][0]["protein_sha256"], reference_sequence_sha256="b"*64,
        reference_metadata_source_refs=[{"path":"references.faa", "sha256":"c"*64}],
        native_fields={"subject_id":"M.Example", "identity":.9, "evalue":1e-20, "reference_length":159, "reference_family":"Type II", "reference_role":"M"})]
    e = record(data)["evidence"]
    assert e["old_strong"] is None and e["query_coverage"] is None and e["reference_coverage"] is None


def test_fresh_reference_unknown_role_does_not_borrow_candidate_role():
    data = bundle(blast_identity="1", blast_length="159", blast_evalue="0", hit_label="M.Example")
    data["suite"]["reference_alignments"] = [dict(sample_id="sample", gene_id="gene", source_ref=REF,
        operation="fresh_reference_alignment", command_exit_code=0,
        protein_sha256=data["genes"][0]["protein_sha256"], reference_sequence_sha256="b"*64,
        reference_metadata_source_refs=[{"path":"references.faa", "sha256":"c"*64}],
        native_fields={"subject_id":"M.Example", "identity":1, "evalue":0, "qstart":1,"qend":159,"sstart":1,"send":159,"reference_length":159})]
    e = record(data)["evidence"]
    assert e["role_compatible"] is None and e["family_compatible"] is None and e["old_strong"] is None


def test_selected_native_score_motif_and_annotation_diagnostic_are_distinct():
    data = bundle(factor_strong="FALSE", factor_moderate="FALSE", factor_context="FALSE", numeric_structure_admissible_structural_evidence="FALSE")
    data["suite"]["annotations"][0]["native_fields"].update(product="RNA methyltransferase", **{"Signature description_Pfam":"RNA methyltransferase"})
    r = record(data)
    e = r["evidence"]
    assert e["legacy_mtase_motif"] is True
    # Selected source 47c058... score_motif_support explicitly permits this
    # sequence evidence, but its diagnostic predicate and function do not.
    assert e["old_motif"] is True
    assert e["annotation"] is False and e["old_profile"] is False
    assert e["alternative_annotation"] is True
    assert r["axis_sources"]["diagnostic_motif_support"]["value"] is False
    assert r["biological_activity"] == "not_inferred"


def test_fresh_measured_weak_reference_overrides_legacy_strong_factor():
    data = bundle(blast_identity="1", blast_length="159", blast_evalue="0", hit_label="M.Example", blast_alignment_quality="strong", factor_strong="TRUE", factor_moderate="TRUE", blast_query_coverage="1", blast_reference_coverage="1", blast_role_compatible="TRUE", blast_family_compatible="TRUE")
    data["suite"]["reference_alignments"] = [dict(sample_id="sample", gene_id="gene", source_ref=REF,
        operation="fresh_reference_alignment", command_exit_code=0,
        protein_sha256=data["genes"][0]["protein_sha256"], reference_sequence_sha256="b"*64,
        reference_metadata_source_refs=[{"path":"references.faa", "sha256":"c"*64}],
        native_fields={"subject_id":"M.Example", "identity":.2, "evalue":1e-4, "qstart":1,"qend":50,"sstart":1,"send":50,"reference_length":159,"reference_family":"Type II","reference_role":"M"})]
    r = record(data); e = r["evidence"]
    assert e["identity"] is None and e["query_coverage"] is None
    assert e["old_strong"] is None and e["old_moderate"] is None
    assert r["reference_measurement_diagnostic"]["native_fields"]["identity"] == .2
    assert e["legacy_mtase_motif"] is True  # independent native axis preserved


def test_candidate_role_is_never_its_own_reference_compatibility():
    data = bundle(hit_label="M.Example", blast_identity=.9, blast_query_coverage=.95,
                  blast_reference_coverage=.95, blast_evalue=1e-80)
    e = record(data)["evidence"]
    assert e["role_compatible"] is None and e["family_compatible"] is None
    assert e["old_strong"] is None and e["old_moderate"] is None
    data["suite"]["features"][0]["native_fields"]["blast_alignment_quality"] = "strong"
    assert record(data)["evidence"]["old_strong"] is None


def test_independent_reference_metadata_can_establish_compatibility():
    data = bundle(hit_label="M.Example", blast_identity=.9, blast_query_coverage=.95,
                  blast_reference_coverage=.95, blast_evalue=1e-80,
                  reference_family="Type II", reference_role="M")
    assert record(data)["evidence"]["old_strong"] is True
    f = data["suite"]["features"][0]["native_fields"]
    del f["reference_family"], f["reference_role"]
    f.update(curated_blast_subject_id="M.Other", curated_blast_family="Type II", curated_blast_role="M")
    assert record(data)["evidence"]["role_compatible"] is None
    f["curated_blast_subject_id"] = "M.Example"
    assert record(data)["evidence"]["old_strong"] is True


def test_comparator_disagreement_cannot_veto_native_assignment():
    data = bundle()
    data["calls"] = [dict(sample_id="sample", gene_id="gene", caller="PADLOC", evidence_kind="final_member",
                          rm_type="Type I", native_model="MTase_I_00001", source_ref=REF)]
    r = record(data)
    assert (r["family"], r["role"], r["assignment_status"]) == ("Type II", "M", "supported")
    assert "comparator_assignment_disagreement_diagnostic" in r["reason_codes"]
    assert r["assignment_support"] == [{"family":"Type II", "role":"M", "rule":"saved_candidate_assignment"}]


@pytest.mark.parametrize("role", ["Aux", "A", "C", "F"])
def test_explicit_auxiliary_role_cannot_be_reinferred_as_restriction(role):
    data = bundle(family="Type IV", role=role)
    data["suite"]["annotations"][0]["native_fields"] = {"translation":"FAGAGDPPY"+"A"*150, "product":"restriction enzyme"}
    data["calls"] = [dict(sample_id="sample", gene_id="gene", caller="PADLOC", evidence_kind="final_member",
                          rm_type="Type IV", native_model="REase_IV_00001", source_ref=REF)]
    r = record(data)
    assert r["role"] == role and r["assignment_status"] == "unresolved"


@pytest.mark.parametrize("product,expected", [
    ("restriction endonuclease subunit S", "S"),
    ("DNA methyltransferase specificity subunit", "S"),
    ("restriction enzyme/methyltransferase specificity subunit", "RM"),
    ("DNA methyltransferase", "M"),
    ("restriction endonuclease", "R"),
])
def test_missing_role_precedence_matches_selected_native_source(product, expected):
    data = bundle(family="Type I", role="NA")
    data["suite"]["annotations"][0]["native_fields"] = {"translation":"FAGAGDPPY"+"A"*150, "product":product}
    r = record(data)
    assert r["role"] == expected and r["assignment_status"] == "supported"
    assert r["evidence"]["partial"] is None  # old threshold was for a missing role


@pytest.mark.parametrize("product", ["transposase", "DNA mismatch repair protein MutS"])
def test_alternative_keywords_do_not_create_an_assessed_scope_negative(product):
    from defenseviz.automatic import apply_automatic_rm
    data = bundle(family="NA", role="NA")
    data["saved_decisions"][0]["state"] = "unassessed"
    data["suite"]["annotations"][0]["native_fields"] = {"translation":"FAGAGDPPY"+"A"*150, "product":product}
    r = record(data)
    assert r["assignment_status"] == "unresolved"
    assert "alternative_annotation_diagnostic_not_scope_exclusion" in r["reason_codes"]
    result = apply_automatic_rm(data)
    assert result["saved_decisions"][0]["state"] == "unassessed"


def test_native_context_is_any_family_or_with_missing_exports_unknown():
    data = bundle(family="Type I", typei_operon_supported="FALSE", typeii_operon_supported="TRUE")
    assert record(data)["evidence"]["old_context"] is True
    data["suite"]["features"][0]["native_fields"].update(typeii_operon_supported="FALSE")
    assert record(data)["evidence"]["old_context"] is None
    data["suite"]["features"][0]["native_fields"].update(typeiii_operon_supported="FALSE")
    assert record(data)["evidence"]["old_context"] is False
    data["suite"]["contexts"] = [dict(sample_id="sample", gene_id="gene", family="Type III", source_ref=REF, native_fields={"operon_supported":"TRUE"})]
    assert record(data)["evidence"]["old_context"] is True


def test_motif_gating_uses_broader_native_annotation_selector():
    seq = "HANAAAAAH" + "A" * 150
    data = bundle(sequence=seq, role="R")
    data["suite"]["annotations"][0]["native_fields"] = {"translation":seq, "product":"hypothetical protein", "custom_domain_description":"HNH PF01844"}
    r = record(data)
    assert r["evidence"]["legacy_nuclease"] is True
    assert r["evidence"]["annotation"] is False  # curation selector stays distinct


def test_terminal_stop_uses_native_normalization_after_exact_raw_binding():
    seq = "FAGAGDPPY" + "A" * 150 + "*"
    r = record(bundle(sequence=seq))
    assert r["protein_sha256"] == sha256(seq.encode()).hexdigest()
    assert r["evidence"]["valid_sequence"] is True and r["evidence"]["legacy_mtase_motif"] is True
    assert r["motif_sequence_normalization"]["changed"] is True
    assert r["motif_sequence_normalization"]["normalized_sequence_sha256"] == sha256(seq[:-1].encode()).hexdigest()


@pytest.mark.parametrize("native,axis", [
    ("rease_nuclease_verified", "legacy_nuclease"), ("rease_motor_verified", "legacy_motor"),
    ("mtase_motif_verified", "legacy_mtase_motif"), ("factor_annotation", "annotation"),
    ("factor_motif", "old_motif"), ("factor_profile", "old_profile"), ("factor_context", "old_context"),
])
def test_explicit_native_feature_values_preserve_exported_parity(native, axis):
    for value in (False, True):
        assert record(bundle(**{native:value}))["evidence"][axis] is value


def test_runtime_pattern_hashes_are_computed_and_pinned():
    from defenseviz.automatic_rm import MOTIF_PATTERNS, ANNOTATION_PATTERNS, RULE_PROVENANCE, _pattern_digest, _EXPECTED_PATTERN_HASHES
    assert RULE_PROVENANCE["installed_motif_patterns_sha256"] == _pattern_digest(MOTIF_PATTERNS)
    assert RULE_PROVENANCE["installed_annotation_patterns_sha256"] == _pattern_digest(ANNOTATION_PATTERNS)
    assert all(RULE_PROVENANCE[k] == v for k,v in _EXPECTED_PATTERN_HASHES.items())


def test_exported_legacy_motif_flags_govern_score_after_rescan_override():
    assert record(bundle(mtase_motif_verified=False))["evidence"]["old_motif"] is False
    assert record(bundle(sequence="A"*159, mtase_motif_verified=True))["evidence"]["old_motif"] is True
    for nuc,motor,expected in [(True,True,True),(False,True,False),(True,False,False)]:
        r=record(bundle(role="R",rease_nuclease_verified=nuc,rease_motor_verified=motor))
        assert r['evidence']['old_motif'] is expected
    assert record(bundle(mtase_motif_verified=False,factor_motif=True))["evidence"]["old_motif"] is True


def test_native_profile_factor_reads_all_signature_prefix_columns():
    data=bundle()
    data['suite']['annotations'][0]['native_fields']={
        'translation':'FAGAGDPPY'+'A'*150,'product':'hypothetical protein',
        'Signature alternative_label_NewLibrary':'DNA methyltransferase'}
    r=record(data)
    assert r['evidence']['old_profile'] is True
    assert r['evidence']['annotation'] is False


def test_curated_subject_precedes_raw_hit_current_alignment():
    data=bundle(hit_label='M.Raw',curated_blast_resolved_subject_id='M.Curated')
    def alignment(subject,identity):
        return dict(sample_id='sample',gene_id='gene',source_ref=REF,operation='fresh_reference_alignment',command_exit_code=0,
            protein_sha256=data['genes'][0]['protein_sha256'],reference_sequence_sha256='b'*64,
            reference_metadata_source_refs=[{'path':'reference.faa','sha256':'c'*64}],
            native_fields=dict(subject_id=subject,identity=identity,evalue=1e-40,query_coverage=.95,reference_coverage=.95,
                               reference_family='Type II',reference_role='M',reference_length=159))
    data['suite']['reference_alignments']=[alignment('M.Raw',.2),alignment('M.Curated',.9)]
    r=record(data)
    assert r['evidence']['old_strong'] is True
    assert r['reference_measurement_diagnostic']['subject_id']=='M.Curated'
    data['suite']['reference_alignments']=[alignment('M.Raw',.9)]
    assert record(data)['evidence']['old_strong'] is None


def test_moderate_nominated_subject_does_not_rule_out_global_strong_reference():
    data=bundle(hit_label='M.Example')
    data['suite']['reference_alignments']=[dict(sample_id='sample',gene_id='gene',source_ref=REF,
        operation='fresh_reference_alignment',command_exit_code=0,protein_sha256=data['genes'][0]['protein_sha256'],
        reference_sequence_sha256='b'*64,reference_metadata_source_refs=[{'path':'reference.faa','sha256':'c'*64}],
        native_fields=dict(subject_id='M.Example',identity=.4,evalue=1e-40,query_coverage=.95,reference_coverage=.95,
                           reference_family='Type II',reference_role='M',reference_length=159))]
    r=record(data)
    assert r['evidence']['old_moderate'] is True and r['evidence']['old_strong'] is None


def test_merged_hsp_union_coverage_overrides_enclosing_coordinates():
    data=bundle(hit_label='M.Example')
    data['suite']['reference_alignments']=[dict(sample_id='sample',gene_id='gene',source_ref=REF,
        operation='fresh_reference_alignment',command_exit_code=0,protein_sha256=data['genes'][0]['protein_sha256'],
        reference_sequence_sha256='b'*64,reference_metadata_source_refs=[{'path':'reference.faa','sha256':'c'*64}],
        native_fields=dict(subject_id='M.Example',identity=.9,evalue=1e-40,qstart=1,qend=159,sstart=1,send=159,
            query_coverage=.75,reference_coverage=.75,dnmb_hsp_count=2,
            reference_family='Type II',reference_role='M',reference_length=159))]
    r=record(data)
    assert r['evidence']['query_coverage']==.75 and r['evidence']['reference_coverage']==.75
    assert r['evidence']['old_moderate'] is True and r['evidence']['old_strong'] is None

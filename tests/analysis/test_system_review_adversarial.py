"""Independent controls for generic saved-system review, using invented inputs.

These exercise scientific boundaries rather than duplicate the implementation.
The known Wadjet component names occur only in synthetic observations; no real
sequences or cohort results are included.
"""
from copy import deepcopy
import hashlib
import json

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.analysis.system_review import (
    review_all_systems,
    validate_system_review_rules,
)


def sha(value):
    return hashlib.sha256(value).hexdigest()


def bundle(tmp_path, records):
    """Make source-bound canonical observations with actual local artifacts."""
    records = deepcopy(records)
    calls, genes, samples = [], {}, {}
    for index, row in enumerate(records):
        sid = row.get("sample_id", "synthetic")
        gid = row.get("gene_id", f"protein-{index}")
        samples[sid] = {"sample_id": sid, "label": sid}
        seq = "M" + "ACDEFGHIKLMNPQRSTVWY" * (index + 1)
        genes.setdefault((sid, gid), dict(sample_id=sid, gene_id=gid,
            protein_sha256=sha(seq.encode()), contig="contig-1", start=1,
            end=len(seq)*3, strand="+"))
        call = dict(sample_id=sid, gene_id=gid, caller="DefenseFinder",
            family="FutureDefense", native_family="FutureDefense",
            native_model="FutureDefense__Effector", category="Unresolved",
            evidence_kind="final_member", selected=True, system_id="system-1",
            run_id="run-1", replicon_id="contig-1", model_version="synthetic-1",
            native_fields={"model_fqn": "defense-finder-models/DefenseFinder/Synthetic/FutureDefense"},
            protein_sha256=genes[(sid, gid)]["protein_sha256"])
        call.update({k: v for k, v in row.items() if k not in {"sample_id", "gene_id"}})
        calls.append(call)
    source = tmp_path / "native-observations.jsonl"
    source.write_text("".join(json.dumps(c, sort_keys=True) + "\n" for c in calls))
    ref = {"path": str(source), "sha256": sha(source.read_bytes())}
    for number, call in enumerate(calls, 1):
        call["source_ref"] = dict(ref, row=number)
    return dict(samples=list(samples.values()), genes=list(genes.values()), calls=calls,
        coverage=[dict(sample_id=sid, caller=provider, status="complete")
            for sid, provider in sorted({(c["sample_id"], c["caller"]) for c in calls})],
        saved_decisions=[], annotations=[], effector_evidence=[], groups=[], complexes=[],
        input_sha256={str(source): ref["sha256"]}, inputs=[str(source)])


def component_rule(rule_id="future-effector", *, model="FutureDefense__Effector",
                   component="FutureEffector", role="nuclease", **conclusion):
    return dict(rule_id=rule_id, version="1", route="native_component",
        endpoint="component_role", primary_sources=["synthetic:test-fixture-only"],
        match=dict(provider="DefenseFinder", exact_native_models=[model],
            evidence_kinds=["final_member"], selected=True,
            model_versions=["synthetic-1"]),
        requires=dict(protein_binding=True, source_binding=True),
        conclusion=dict(component=component, role=role, substrate="DNA", **conclusion))


def registry(*rules):
    return dict(schema_version=1, registry_id="synthetic-controls", version="1", rules=list(rules))


def reviewed(data, rules=None):
    accepted = analyze_saved(data)
    before_data, before_result = deepcopy(data), deepcopy(accepted)
    review = review_all_systems(data, accepted, rules)
    assert data == before_data, "Review must not rewrite canonical evidence"
    assert accepted == before_result, "Review must not mutate accepted analysis or figures"
    return accepted, review


def gene_review(review, gid="protein-0", sid="synthetic"):
    return next(x for x in review["gene_reviews"] if (x["sample_id"], x["gene_id"]) == (sid, gid))


def axis(review, name="component_role", gid="protein-0", sid="synthetic"):
    return gene_review(review, gid, sid)["axes"][name]


def architecture_rules(tmp_path, *, names=("A", "B"), family="FutureDefense"):
    components = [component_rule(f"component-{n}", model=f"FutureDefense__{n}",
                                component=n, role="module_component") for n in names]
    source = tmp_path / "synthetic-architecture-contract.txt"
    source.write_text("Synthetic all-component contract. Each distinct component requires one distinct protein.\n")
    complete = dict(rule_id="complete-future-core", version="1", route="native_architecture",
        endpoint="complex_core", primary_sources=["synthetic:test-fixture-only"],
        architecture_contract=dict(contract_id="synthetic-core", version="1",
            source_ref=dict(path=str(source), sha256=sha(source.read_bytes()))),
        match=dict(provider="DefenseFinder", exact_native_families=[family],
            evidence_kinds=["final_member"], selected=True, model_versions=["synthetic-1"]),
        requires=dict(protein_binding=True, source_binding=True, same_native_unit=True,
            all_components=[dict(component=n, component_rule_ids=[f"component-{n}"],
                                 minimum_members=1, maximum_members=1) for n in names]),
        conclusion=dict(role="nuclease_complex_core", substrate="DNA", individual_catalysis="unassigned"))
    return registry(*components, complete)


def test_unknown_family_is_accounted_without_invented_function(tmp_path):
    data = bundle(tmp_path, [{}])
    accepted, review = reviewed(data)
    assert review["review_execution"]["status"] == "completed"
    assert len(review["detection_accounting"]) == len(data["calls"])
    assert review["units"][0]["axes"]["native_detection"]["status"] == "supported"
    for endpoint in ("component_role", "nuclease_domain", "motif_site", "activity"):
        assert axis(review, endpoint)["status"] == "unassessed"
    assert review["pending_queue"]
    assert accepted["global_summary"]["native_union_CDS"] == 1
    assert accepted["global_summary"]["effector_CDS"] == 0


def test_exact_future_rule_supports_review_without_promoting_accepted_counts(tmp_path):
    data = bundle(tmp_path, [{}])
    _, unknown = reviewed(data)
    accepted, review = reviewed(data, registry(component_rule()))
    assert axis(unknown)["status"] == "unassessed"
    assert axis(review)["status"] == "supported"
    assert axis(review)["evaluated_rule_ids"]
    assert axis(review, "motif_site")["status"] == "unassessed"
    assert axis(review, "activity")["status"] == "unassessed"
    assert all(accepted["global_summary"][name] == 0
        for name in ("effector_CDS", "associated_CDS", "priority_CDS", "motif_CDS"))


@pytest.mark.parametrize("change", [
    {"native_model": "FutureDefense__Effector_extra"},
    {"native_model": "futuredefense__effector"},
    {"caller": "PADLOC"},
    {"model_version": "synthetic-2"},
    {"model_version": None, "software_version": "synthetic-1"},
    {"selected": False},
    {"evidence_kind": "hmm_hit"},
])
def test_near_match_missing_version_or_nonfinal_evidence_never_satisfies_rule(tmp_path, change):
    data = bundle(tmp_path, [change])
    _, review = reviewed(data, registry(component_rule()))
    assert axis(review)["status"] == "unassessed"


def test_model_digest_is_an_exact_alternative_not_a_substring(tmp_path):
    data = bundle(tmp_path, [{"model_version": None, "model_sha256": "1" * 64}])
    rule = component_rule()
    rule["match"].pop("model_versions")
    rule["match"]["model_sha256"] = "1" * 64
    _, good = reviewed(data, registry(rule))
    assert axis(good)["status"] == "supported"
    data["calls"][0]["model_sha256"] = "2" * 64
    _, wrong = reviewed(data, registry(rule))
    assert axis(wrong)["status"] == "unassessed"


def test_version_and_digest_constraints_are_both_required_when_supplied(tmp_path):
    data = bundle(tmp_path, [{"model_sha256": "1" * 64, "model_version": "wrong-version"}])
    rule = component_rule()
    rule["match"]["model_sha256"] = ["1" * 64]
    _, review = reviewed(data, registry(rule))
    assert axis(review)["status"] == "unassessed"


def test_stale_native_artifact_hash_cannot_support_new_review(tmp_path):
    data = bundle(tmp_path, [{}])
    accepted = analyze_saved(data)
    (tmp_path / "native-observations.jsonl").write_text("tampered after import\n")
    with pytest.raises(ValueError, match="hash|digest|changed|SHA"):
        review_all_systems(data, accepted, registry(component_rule()))


def test_call_protein_hash_cannot_be_borrowed_from_current_gene(tmp_path):
    data = bundle(tmp_path, [{}])
    accepted = analyze_saved(data)
    data["calls"][0]["protein_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="hash|protein|sequence"):
        review_all_systems(data, accepted, registry(component_rule()))


def test_missing_binding_is_reasoned_unassessed_not_positive(tmp_path):
    data = bundle(tmp_path, [{}])
    data["calls"][0].pop("source_ref")
    _, review = reviewed(data, registry(component_rule()))
    assert axis(review)["status"] == "unassessed"
    assert axis(review)["reason_codes"]


def test_native_manifest_relative_source_is_resolved_by_exact_pinned_source_record(tmp_path):
    data = bundle(tmp_path, [{}])
    ref = data["calls"][0]["source_ref"]
    resolved = ref["path"]
    ref["path"] = "runs/synthetic/native-observations.jsonl"
    data["suite"] = {"sources": [dict(path=ref["path"], sha256=ref["sha256"],
        resolved_path=resolved, sample_id="synthetic", module="defensefinder")]}
    _, review = reviewed(data, registry(component_rule()))
    assert axis(review)["status"] == "supported"


def test_same_relative_label_with_conflicting_bound_paths_is_not_guessed(tmp_path):
    data = bundle(tmp_path, [{}])
    ref = data["calls"][0]["source_ref"]
    actual = ref["path"]
    duplicate = tmp_path / "same-content-other-source.jsonl"
    duplicate.write_bytes((tmp_path / "native-observations.jsonl").read_bytes())
    ref["path"] = "ambiguous/native.jsonl"
    data["input_sha256"][str(duplicate)] = ref["sha256"]
    data["suite"] = {"sources": [dict(path=ref["path"], sha256=ref["sha256"],
        resolved_path=p, sample_id="synthetic", module="defensefinder") for p in (actual, str(duplicate))]}
    accepted = analyze_saved(data)
    with pytest.raises(ValueError, match="ambig|conflict|source|path"):
        review_all_systems(data, accepted, registry(component_rule()))


def test_every_call_is_accounted_even_with_nonexhaustive_context_groups(tmp_path):
    data = bundle(tmp_path, [
        {"gene_id": "grouped"},
        {"gene_id": "second-system", "system_id": "other-system"},
        {"gene_id": "ungrouped", "system_id": None},
        {"gene_id": "hmm-only", "system_id": None, "evidence_kind": "hmm_hit", "selected": False},
        {"gene_id": "unselected", "system_id": "unselected-system", "selected": False},
    ])
    g = data["genes"][0]
    data["groups"] = [dict(group_id="partial-context", sample_id="synthetic",
        provider="DefenseFinder", run_id="run-1", native_system_id="system-1",
        membership_kind="native_system", members=[dict(gene_id=g["gene_id"], protein_sha256=g["protein_sha256"])],
        source_refs=[data["calls"][0]["source_ref"]])]
    _, review = reviewed(data)
    assert len(review["detection_accounting"]) == 5
    assert len({x["detection_id"] for x in review["detection_accounting"]}) == 5
    assert {x["gene_id"] for x in review["gene_reviews"]} == {g["gene_id"] for g in data["genes"]}
    assert {u["unit_kind"] for u in review["units"]} >= {"native_system", "ungrouped_detection", "auxiliary_detection"}


def test_anti_defense_nuclease_name_does_not_enter_host_effector_set(tmp_path):
    data = bundle(tmp_path, [{"native_family": "Anti_CBASS", "family": "Anti_CBASS",
        "native_model": "acb1", "namespace": "ADF",
        "native_fields": {"model_fqn": "defense-finder-models/ADF/Anti_CBASS/acb1"}}])
    accepted, review = reviewed(data)
    assert accepted["global_summary"]["effector_CDS"] == 0
    assert axis(review)["status"] == "unassessed"
    assert len(review["detection_accounting"]) == 1


def test_other_mechanism_annotation_is_preserved_without_nuclease_promotion(tmp_path):
    data = bundle(tmp_path, [{"native_model": "FutureDefense__NucleotideModifier"}])
    rule = component_rule(model="FutureDefense__NucleotideModifier", component="Modifier", role="nucleotide_modification")
    accepted, review = reviewed(data, registry(rule))
    assert axis(review)["status"] == "supported"
    assert axis(review, "nuclease_domain")["status"] == "unassessed"
    assert accepted["global_summary"]["native_union_CDS"] == 1
    assert accepted["global_summary"]["effector_CDS"] == 0


def test_eggnog_transfer_and_broad_interpro_text_are_not_qualified_nuclease_evidence(tmp_path):
    data = bundle(tmp_path, [{}])
    ref = data["calls"][0]["source_ref"]
    data["suite"] = {"annotations": [dict(sample_id="synthetic", gene_id="protein-0", namespace="EggNOG", source_ref=ref,
        native_fields={"seed_ortholog": "synthetic.1", "Preferred_name": "nuclease", "PFAMs": "HNH", "Description": "active nuclease effector"})],
        "profiles": [dict(sample_id="synthetic", gene_id="protein-0", namespace="InterProScan", source_ref=ref,
        native_fields={"signature_accession": "SYNTHETIC_FOLD", "signature_description": "nuclease fold", "interpro_description": "nuclease", "status": "T"})]}
    accepted, review = reviewed(data)
    assert axis(review, "nuclease_domain")["status"] == "unassessed"
    assert axis(review, "motif_site")["status"] == "unassessed"
    assert axis(review)["status"] == "unassessed"
    assert accepted["global_summary"]["effector_CDS"] == 0


def test_missing_interpro_row_is_unassessed_not_not_applicable(tmp_path):
    data = bundle(tmp_path, [{}])
    data["suite"] = {"profiles": [], "capabilities": [dict(sample_id="synthetic", raw_domain_scope="observed_native_members_and_DNMB_candidates")]}
    _, review = reviewed(data)
    assert axis(review, "nuclease_domain")["status"] == "unassessed"
    assert axis(review, "motif_site")["status"] == "unassessed"


def qualified_domain_pack(tmp_path, data, configure=None):
    """Use the public evidence adapter; never manufacture assessed trust flags."""
    from defenseviz.system_review_inputs import apply_system_review_inputs
    sequence = "M" + "ACDEFGHIKLMNPQRSTVWY"
    source = tmp_path / "synthetic-domain-evidence.txt"
    source.write_text("Invented exact-model aligned-domain evidence for software tests only.\n")
    model = dict(accession="SYNTHETIC_NUCLEASE.1", source="synthetic-hmm",
        model_sha256="b" * 64, evidence_kind="nuclease_domain", role="nuclease",
        family_specific=True, min_model_coverage=.7)
    registry_data = dict(schema_version="1.0", families=[dict(family_id="SyntheticNovelDomain",
        component="SyntheticDomain", substrate="DNA", domain_models=[model], references=[])])
    candidate = dict(candidate_id="synthetic-domain-candidate", sample_id="synthetic", query_id="protein-0",
        sequence=sequence, sequence_sha256=sha(sequence.encode()),
        domain_hits=[dict(accession=model["accession"], source=model["source"], model_sha256=model["model_sha256"],
            query_sequence_sha256=sha(sequence.encode()), qstart=1, qend=len(sequence),
            hmm_start=1, hmm_end=len(sequence), hmm_length=len(sequence), target_length=len(sequence),
            independent_evalue=1e-30, qualification=dict(passed=True, rule="synthetic-qualified-domain",
                source=source.name, source_sha256=sha(source.read_bytes())))], site_maps=[])
    if configure is not None:
        configure(registry_data, candidate)
    reg_path, candidates_path, pack_path = (tmp_path / n for n in ("domain-registry.json", "domain-candidates.json", "domain-pack.json"))
    reg_path.write_text(json.dumps(registry_data) + "\n")
    candidates_path.write_text(json.dumps([candidate]) + "\n")
    pack_path.write_text(json.dumps(dict(schema_version=1, registry=reg_path.name, candidates=candidates_path.name)) + "\n")
    apply_system_review_inputs(data, evidence_path=pack_path)
    return source


def test_bound_qualified_domain_extends_unknown_review_but_not_catalysis_or_effector_counts(tmp_path):
    data = bundle(tmp_path, [{}])
    qualified_domain_pack(tmp_path, data)
    accepted, review = reviewed(data)
    assert axis(review, "nuclease_domain")["status"] == "supported"
    assert axis(review, "motif_site")["status"] == "unassessed"
    assert axis(review, "activity")["status"] == "unassessed"
    assert axis(review)["status"] == "unassessed"
    assert accepted["global_summary"]["effector_CDS"] == 0


def test_postevaluation_assessment_forgery_is_rejected(tmp_path):
    data = bundle(tmp_path, [{}])
    qualified_domain_pack(tmp_path, data)
    accepted = analyze_saved(data)
    family = data["system_review_evidence"]["assessments"][0]["families"][0]
    family["catalytic_site_evidence"]["state"] = "supported"
    with pytest.raises(ValueError, match="assessments|evaluation|evidence"):
        review_all_systems(data, accepted)


def test_qualified_evidence_source_is_rechecked_after_adapter_loading(tmp_path):
    data = bundle(tmp_path, [{}])
    source = qualified_domain_pack(tmp_path, data)
    accepted = analyze_saved(data)
    source.write_text("changed after qualification\n")
    with pytest.raises(ValueError, match="hash|changed|digest"):
        review_all_systems(data, accepted)


def test_claiming_verified_boolean_does_not_create_qualified_evidence(tmp_path):
    data = bundle(tmp_path, [{}])
    data["system_review_evidence"] = dict(schema_version=1, source_files_verified=True,
        registry_sha256="1" * 64, assessments=[dict(state="supported")])
    with pytest.raises(ValueError):
        review_all_systems(data, analyze_saved(data))


def test_wadjet_native_detection_does_not_impute_missing_jetd(tmp_path):
    data = bundle(tmp_path, [dict(gene_id=n, family="Wadjet II", native_family="Wadjet_II", native_model=f"Wadjet__Jet{n}_II") for n in "ABC"])
    accepted, review = reviewed(data)
    assert review["units"][0]["axes"]["native_detection"]["status"] == "supported"
    assert accepted["global_summary"]["effector_CDS"] == 0
    assert all(x["axes"]["component_role"]["status"] != "supported" for x in review["gene_reviews"])
    assert {x["gene_id"] for x in review["gene_reviews"]} == {"A", "B", "C"}


def test_complete_components_support_only_reviewed_complex_not_individual_catalysis(tmp_path):
    data = bundle(tmp_path, [{"gene_id": n, "native_model": f"FutureDefense__{n}"} for n in "AB"])
    accepted, review = reviewed(data, architecture_rules(tmp_path))
    assert review["units"][0]["axes"]["complex_core"]["status"] == "supported"
    assert review["units"][0]["axes"]["activity"]["status"] == "unassessed"
    assert accepted["global_summary"]["complex_core_CDS"] == 0
    assert accepted["global_summary"]["effector_CDS"] == 0


@pytest.mark.parametrize("difference", [
    {"system_id": "other-system"}, {"run_id": "other-run"},
    {"replicon_id": "other-contig"}, {"sample_id": "other-sample"},
])
def test_complex_cannot_be_completed_by_cross_unit_components(tmp_path, difference):
    data = bundle(tmp_path, [{"gene_id": "A", "native_model": "FutureDefense__A"},
        {"gene_id": "B", "native_model": "FutureDefense__B", **difference}])
    _, review = reviewed(data, architecture_rules(tmp_path))
    assert len(review["units"]) == 2
    assert all(u["axes"]["complex_core"]["status"] != "supported" for u in review["units"])


def test_duplicate_core_copies_do_not_satisfy_distinct_missing_component(tmp_path):
    data = bundle(tmp_path, [{"gene_id": n, "native_model": "FutureDefense__A"} for n in ("a1", "a2")])
    _, review = reviewed(data, architecture_rules(tmp_path))
    assert review["units"][0]["axes"]["complex_core"]["status"] != "supported"


def test_one_protein_cannot_implicitly_fill_two_distinct_core_slots(tmp_path):
    data = bundle(tmp_path, [{"gene_id": "fusion-claim", "native_model": f"FutureDefense__{n}"} for n in "AB"])
    _, review = reviewed(data, architecture_rules(tmp_path))
    assert review["units"][0]["axes"]["complex_core"]["status"] != "supported"


def test_compatible_multiple_roles_are_not_forced_into_conflict(tmp_path):
    data = bundle(tmp_path, [{}])
    rules = registry(component_rule("catalytic"), component_rule("regulatory", role="regulation"))
    _, review = reviewed(data, rules)
    assert axis(review)["status"] == "supported"


def test_explicit_exclusion_retains_both_claims_as_conflict(tmp_path):
    data = bundle(tmp_path, [{}])
    rules = registry(component_rule("catalytic", excludes_roles=["regulation"]), component_rule("regulatory", role="regulation"))
    _, review = reviewed(data, rules)
    assert axis(review)["status"] == "conflicting"
    assert axis(review)["conflicting_evidence_ids"]


def test_accepted_noncatalytic_review_is_not_silently_overridden_by_new_rule(tmp_path):
    data = bundle(tmp_path, [{}])
    gene = data["genes"][0]
    data["annotations"] = [dict(sample_id=gene["sample_id"], gene_id=gene["gene_id"],
        protein_sha256=gene["protein_sha256"], assignment_kind="noncatalytic_component",
        status="supported", rule_id="synthetic-explicit-noncatalytic-review", role="A",
        source_refs=[data["calls"][0]["source_ref"]])]
    accepted, review = reviewed(data, registry(component_rule()))
    assert axis(review)["status"] == "conflicting"
    assert gene_review(review)["preserved_assignment_ids"]
    assert accepted["global_summary"]["effector_CDS"] == 0


def test_complex_requires_bound_native_run_identity(tmp_path):
    data = bundle(tmp_path, [{"gene_id": n, "native_model": f"FutureDefense__{n}", "run_id": None} for n in "AB"])
    _, review = reviewed(data, architecture_rules(tmp_path))
    assert review["units"][0]["axes"]["complex_core"]["status"] == "unassessed"


def test_architecture_rule_cannot_claim_individual_catalysis(tmp_path):
    rules = architecture_rules(tmp_path)
    rules["rules"][-1]["conclusion"]["individual_catalysis"] = "confirmed"
    with pytest.raises(ValueError, match="catalysis"):
        validate_system_review_rules(rules)


def test_custom_rule_cannot_replace_bundled_identity_or_conclude_activity():
    reserved = component_rule("defensefinder_jetd")
    with pytest.raises(ValueError, match="reserved|Duplicate"):
        validate_system_review_rules(registry(reserved))
    claim = component_rule()
    claim["conclusion"]["activity"] = "confirmed"
    with pytest.raises(ValueError, match="unknown|fields"):
        validate_system_review_rules(registry(claim))


def test_shared_gene_does_not_transfer_complete_core_to_another_run(tmp_path):
    """Claude H3: a gene can have complex evidence without every unit having it."""
    data = bundle(tmp_path, [
        {"gene_id": "shared", "native_model": "FutureDefense__A"},
        {"gene_id": "partner", "native_model": "FutureDefense__B"},
        {"gene_id": "shared", "native_model": "FutureDefense__A", "run_id": "run-2",
         "family": "OtherDefense", "native_family": "OtherDefense", "system_id": "other-system"},
    ])
    _, review = reviewed(data, architecture_rules(tmp_path))
    complete = next(u for u in review["units"] if u["run_id"] == "run-1")
    other = next(u for u in review["units"] if u["run_id"] == "run-2")
    assert complete["axes"]["complex_core"]["status"] == "supported"
    assert gene_review(review, "shared")["axes"]["complex_core"]["status"] == "supported"
    assert other["axes"]["complex_core"]["status"] == "unassessed"
    assert not other["axes"]["complex_core"]["supporting_evidence_ids"]


def test_shared_gene_role_from_one_provider_is_not_native_unit_support_for_another(tmp_path):
    data = bundle(tmp_path, [{"gene_id": "shared"},
        {"gene_id": "shared", "caller": "PADLOC", "run_id": "pad-run", "system_id": "pad-system"}])
    _, review = reviewed(data, registry(component_rule()))
    assert gene_review(review, "shared")["axes"]["component_role"]["status"] == "supported"
    native = next(u for u in review["units"] if u["provider"] == "DefenseFinder")
    other = next(u for u in review["units"] if u["provider"] == "PADLOC")
    assert native["axes"]["component_role"]["status"] == "supported"
    assert other["axes"]["component_role"]["status"] == "unassessed"
    assert other["member_gene_evidence"], "Intrinsic evidence may be reported separately, not copied into unit claims"


def test_new_explicit_noncatalytic_role_conflicts_without_manual_excludes_list(tmp_path):
    data = bundle(tmp_path, [{}])
    _, review = reviewed(data, registry(component_rule("nuclease-claim"),
        component_rule("noncatalytic-claim", role="noncatalytic")))
    assert axis(review)["status"] == "conflicting"
    assert axis(review)["conflicting_evidence_ids"]


def qualified_rule(rule_id, endpoint, signature, family="SyntheticNovelDomain"):
    return dict(rule_id=rule_id, version="1", route="qualified_domain", endpoint=endpoint,
        primary_sources=["synthetic:test-fixture-only"],
        match=dict(family_ids=[family], namespace="synthetic-hmm", signature_accessions=[signature]),
        requires=dict(protein_binding=True, source_binding=True),
        conclusion=dict(family_id=family) if endpoint == "family_annotation" else dict(role="nuclease", substrate="DNA"))


def test_qualified_family_conclusion_cannot_rename_a_different_matched_family():
    rule = qualified_rule("family-mismatch", "family_annotation", "SYNTHETIC_NUCLEASE.1", family="RegisteredA")
    rule["conclusion"]["family_id"] = "DifferentB"
    with pytest.raises(ValueError, match="family"):
        validate_system_review_rules(registry(rule))


def test_multi_family_selector_cannot_transfer_a_conclusion_from_the_other_entry(tmp_path):
    data = bundle(tmp_path, [{}])
    qualified_domain_pack(tmp_path, data)
    rule = qualified_rule("wrong-entry-conclusion", "family_annotation", "SYNTHETIC_NUCLEASE.1")
    rule["match"]["family_ids"].append("DifferentB")
    rule["conclusion"]["family_id"] = "DifferentB"
    _, review = reviewed(data, registry(rule))
    evaluation = next(e for e in review["rule_evaluations"] if e["rule_id"] == rule["rule_id"])
    assert evaluation["status"] == "unassessed"
    assert "qualified_conclusion_family_mismatch" in evaluation["reason_codes"]
    assert axis(review, "family_annotation")["status"] == "supported", "Valid own-entry evidence is retained"


@pytest.mark.parametrize("endpoint,x_kind", [("family_annotation", "nuclease_domain"), ("nuclease_domain", "other_domain")])
def test_selected_signature_cannot_borrow_excluded_model_support(tmp_path, endpoint, x_kind):
    """Overall gene support from Y is valid; a rule scoped only to X cannot borrow it."""
    data = bundle(tmp_path, [{}])
    def configure(reg, candidate):
        x = reg["families"][0]["domain_models"][0]
        x.update(evidence_kind=x_kind, role="nuclease" if x_kind == "nuclease_domain" else "other",
                 family_specific=False)
        y = dict(x, accession="SYNTHETIC_SPECIFIC_Y.1", model_sha256="c" * 64,
                 evidence_kind="nuclease_domain", role="nuclease", family_specific=True)
        reg["families"][0]["domain_models"].append(y)
        candidate["domain_hits"].append(dict(candidate["domain_hits"][0],
            accession=y["accession"], model_sha256=y["model_sha256"]))
    qualified_domain_pack(tmp_path, data, configure)
    rules = registry(qualified_rule("only-X", endpoint, "SYNTHETIC_NUCLEASE.1"),
                     qualified_rule("only-Y", endpoint, "SYNTHETIC_SPECIFIC_Y.1"))
    _, review = reviewed(data, rules)
    evaluations = {e["rule_id"]: e for e in review["rule_evaluations"] if e["rule_id"] in {"only-X", "only-Y"}}
    assert evaluations["only-Y"]["status"] == "supported"
    assert evaluations["only-X"]["status"] == "unassessed"
    assert axis(review, endpoint)["status"] == "supported", "Valid independent Y evidence remains available"


def test_observed_but_unbound_component_is_not_reported_as_missing(tmp_path):
    data = bundle(tmp_path, [{"gene_id": n, "native_model": f"FutureDefense__{n}"} for n in "AB"])
    rules = architecture_rules(tmp_path)
    rules["rules"][1]["match"]["model_versions"] = ["different-model-version"]
    _, review = reviewed(data, rules)
    assessment = review["units"][0]["axes"]["complex_core"]
    assert assessment["status"] == "unassessed"
    assert "required_component_unbound" in assessment["reason_codes"]
    assert "required_component_unobserved" not in assessment["reason_codes"]


def test_matching_host_rule_does_not_reclassify_an_anti_defense_nuclease(tmp_path):
    data = bundle(tmp_path, [{"native_family": "Anti_CBASS", "family": "Anti_CBASS",
        "native_model": "acb1", "namespace": "ADF",
        "native_fields": {"model_fqn": "defense-finder-models/ADF/Anti_CBASS/acb1"}}])
    host = component_rule("host-acb1-claim", model="acb1")
    host["defense_scope"] = "host_defense"
    accepted, review = reviewed(data, registry(host))
    assert review["units"][0]["defense_scope"] == "anti_defense"
    assert axis(review)["status"] == "unassessed"
    assert accepted["global_summary"]["effector_CDS"] == 0
    explicit = deepcopy(host)
    explicit["rule_id"] = "explicit-anti-defense-annotation"
    explicit["defense_scope"] = "anti_defense"
    _, correctly_scoped = reviewed(data, registry(explicit))
    assert axis(correctly_scoped)["status"] == "supported"
    assert correctly_scoped["units"][0]["defense_scope"] == "anti_defense"


def test_retained_family_without_role_is_not_component_role_support(tmp_path):
    data = bundle(tmp_path, [{}])
    data["saved_decisions"] = [dict(sample_id="synthetic", gene_id="protein-0", state="retained",
        family="Type II", source_ref=data["calls"][0]["source_ref"])]
    accepted, review = reviewed(data)
    assert accepted["global_summary"]["retained_CDS"] == 1
    assert axis(review, "family_annotation")["status"] == "supported"
    assert axis(review)["status"] == "unassessed"


def test_review_coverage_reuses_resolved_analysis_and_preserves_annotation_capabilities(tmp_path):
    data = bundle(tmp_path, [{"native_model": "Wadjet__JetD_II"}])
    capability = dict(sample_id="synthetic", raw_interpro_status="partial",
        raw_interpro_excluded_rows=3, raw_interpro_quarantined_rows=1, raw_eggnog_status="complete")
    data["suite"] = {"capabilities": [capability]}
    accepted, review = reviewed(data)
    resolved = accepted["availability"]["assessment_coverage_by_sample"]["synthetic"]
    assert resolved["effector"] == "partial", "Positive observed evidence without a declaration is a lower bound"
    assert axis(review)["coverage"] == resolved["effector"]
    coverage = next(r for r in review["sample_coverage"] if r["sample_id"] == "synthetic")
    assert capability in coverage["annotation_capabilities"]


@pytest.mark.parametrize("declaration", [{"unknown-sample": {"effector": "complete"}}, {"synthetic": {"invented-axis": "complete"}}])
def test_standalone_review_rejects_invalid_scoped_coverage(declaration, tmp_path):
    data = bundle(tmp_path, [{}])
    accepted = analyze_saved(data)
    data["options"] = {"assessment_coverage_by_sample": declaration}
    with pytest.raises(ValueError, match="coverage|sample|assessment"):
        review_all_systems(data, accepted)


def test_native_rule_cannot_use_an_existing_but_unpinned_artifact(tmp_path):
    data = bundle(tmp_path, [{}])
    unrelated = tmp_path / "not-imported-native-input.txt"
    unrelated.write_text("Unpinned but readable; not an imported native result.\n")
    data["calls"][0]["source_ref"] = dict(path=str(unrelated), sha256=sha(unrelated.read_bytes()), row=1)
    _, review = reviewed(data, registry(component_rule()))
    assert axis(review)["status"] == "unassessed"
    assert "source_binding_unavailable" in axis(review)["reason_codes"]


def test_model_package_version_cannot_satisfy_database_release_predicate(tmp_path):
    data = bundle(tmp_path, [{}])
    rule = component_rule()
    rule["match"]["database_versions"] = rule["match"].pop("model_versions")
    _, review = reviewed(data, registry(rule))
    assert axis(review)["status"] == "unassessed"


def test_importer_attached_protein_hash_is_not_claimed_independent_producer_evidence(tmp_path):
    data = bundle(tmp_path, [{}])
    data["calls"][0].pop("protein_sha256")
    _, review = reviewed(data, registry(component_rule()))
    evaluation = next(e for e in review["rule_evaluations"] if e["rule_id"] == "future-effector")
    assert evaluation["protein_binding"]["independent_producer_binding"] is False
    assert evaluation["protein_binding"]["origin"] == "canonical_gene_identity"


@pytest.mark.parametrize("mutation", ["unknown_match_field", "gene_selector", "no_version_or_hash", "duplicate_rule", "unknown_requirement", "invalid_substrate"])
def test_unsupported_or_unscoped_rule_contracts_are_rejected(mutation):
    rule = component_rule()
    rules = registry(rule)
    if mutation == "unknown_match_field": rule["match"]["model_regex"] = ".*"
    elif mutation == "gene_selector": rule["match"]["gene_id"] = "protein-0"
    elif mutation == "no_version_or_hash": rule["match"].pop("model_versions")
    elif mutation == "duplicate_rule": rules["rules"].append(deepcopy(rule))
    elif mutation == "unknown_requirement": rule["requires"]["trust_annotation_text"] = True
    elif mutation == "invalid_substrate": rule["conclusion"]["substrate"] = "anything"
    with pytest.raises(ValueError):
        validate_system_review_rules(rules)

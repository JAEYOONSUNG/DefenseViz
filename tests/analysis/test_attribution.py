"""Synthetic source-bound scope controls, never manuscript outcome fixtures."""
from copy import deepcopy
import hashlib

import pytest

from defenseviz.analysis.attribution import (
    NATIVE_CONTEXT_CONTRACTS, apply_family_attribution, project_family_attribution,
)
from defenseviz.io.suite_curation import native_curation_record


def fixture(rule_index=0):
    rule = NATIVE_CONTEXT_CONTRACTS[rule_index]
    native_ref = dict(path="synthetic-native.tsv", sha256="a" * 64, row=2)
    gene_ref = dict(path="synthetic-genes.tsv", sha256="b" * 64, row=2)
    data = dict(genes=[], calls=[], saved_decisions=[], annotations=[], curation_records=[],
                input_sha256={native_ref["path"]: native_ref["sha256"], gene_ref["path"]: gene_ref["sha256"]})
    for index, component in enumerate(rule["required_components"], 1):
        gid = "synthetic_" + component
        digest = hashlib.sha256(component.encode()).hexdigest()
        gene = dict(sample_id="unseen_genome", gene_id=gid, protein_sha256=digest,
                    contig="chr", start=index*100, end=index*100+80, strand="+", source_refs=[gene_ref])
        data["genes"].append(gene)
        data["calls"].append(dict(sample_id=gene["sample_id"], gene_id=gid, protein_sha256=digest,
            caller=rule["provider"], native_family=rule["native_family"], family=rule["family"], category=rule["category"],
            evidence_kind="final_member", selected=True, run_id="synthetic-run", system_id="synthetic-unit",
            native_model="synthetic-model:"+component, source_ref=native_ref,
            native_fields={"system":rule["native_family"], "hmm.name":"synthetic-model:"+component,
                           "protein.name":component, "relative.position":str(index),
                           "seqid":"chr", "start":str(gene["start"]), "end":str(gene["end"]), "strand":"+"}))
        if component == rule["target_component"]:
            target = gene
    data["saved_decisions"].append(dict(sample_id=target["sample_id"], gene_id=target["gene_id"],
        protein_sha256=target["protein_sha256"], source_ref=native_ref, state="retained",
        family="Type II", rm_type="Type II", role="M", score_total=73, decision_reason="stored sequence evidence"))
    native = dict(query=target["gene_id"], curation_keep="TRUE", final_family="Type II", final_role="M",
                  rm_association_class="DNA_MTase_RM_association_unproven", curation_score="73")
    data["curation_records"].append(native_curation_record(target["sample_id"], target["gene_id"], native, native_ref, target))
    return data


def one(data):
    return apply_family_attribution(data)["family_attribution"][0]


@pytest.mark.parametrize("rule_index", range(3))
def test_supported_native_architecture_reassigns_scope_not_score(rule_index):
    data = fixture(rule_index); before = deepcopy(data)
    out = apply_family_attribution(data)
    assert data == before
    assert out["calls"] == before["calls"] and out["saved_decisions"] == before["saved_decisions"]
    row = out["family_attribution"][0]
    assert row["action"] == "assign_noncanonical_scope"
    assert row["target_family"] == NATIVE_CONTEXT_CONTRACTS[rule_index]["family"]
    assert row["activity"] == "unassessed" and row["preserved_component_role"] == "M"
    assert row["native_contexts"][0]["definition_source"]["sha256"]
    projected = project_family_attribution(out)
    result = projected["saved_decisions"][0]
    assert (result["state"], result["rm_type"], result["role"], result["score_total"]) == ("out_of_scope", None, None, 73)
    assert result["prior_saved_decision"] == before["saved_decisions"][0]
    assert projected["calls"] == before["calls"]
    assert project_family_attribution(out) == projected
    with pytest.raises(ValueError, match="unprojected"):
        project_family_attribution(projected)


@pytest.mark.parametrize("damage", ["missing", "split_run", "split_unit", "wrong_contig", "hash", "unselected", "raw_hit",
                                    "duplicate_component", "missing_ordinal", "far_ordinal", "reversed_coordinates", "multipart", "string_parts", "malformed_part", "native_coordinates"])
def test_partial_or_unbound_context_cannot_override(damage):
    data = fixture(); target = data["saved_decisions"][0]["gene_id"]
    other = next(c for c in data["calls"] if c["gene_id"] != target)
    gene = next(g for g in data["genes"] if g["gene_id"] == other["gene_id"])
    if damage == "missing": data["calls"].remove(other)
    elif damage == "split_run": other["run_id"] = "different"
    elif damage == "split_unit": other["system_id"] = "different"
    elif damage == "wrong_contig": gene["contig"] = "plasmid"
    elif damage == "hash": other["protein_sha256"] = "f"*64
    elif damage == "unselected": other["selected"] = False
    elif damage == "raw_hit": other["evidence_kind"] = "hmm_hit"
    elif damage == "duplicate_component": data["calls"].append(deepcopy(other))
    elif damage == "missing_ordinal": del other["native_fields"]["relative.position"]
    elif damage == "far_ordinal": other["native_fields"]["relative.position"] = "500"
    elif damage == "reversed_coordinates": gene["end"] = gene["start"]-1
    elif damage == "multipart": gene["location_parts"] = [[1,10,1],[100,120,1]]
    elif damage == "string_parts": gene["location_parts"] = "x"
    elif damage == "malformed_part": gene["location_parts"] = [[1,10,1]]
    elif damage == "native_coordinates": other["native_fields"]["start"] = "90000"
    assert one(data)["action"] == "retain_original_unassessed"
    assert project_family_attribution(data)["saved_decisions"] == data["saved_decisions"]


@pytest.mark.parametrize("damage", ["absent", "ambiguous", "unbound", "sha", "source_sha", "scope", "query", "projection"])
def test_no_transfer_without_exact_curation(damage):
    data = fixture(); row = data["curation_records"][0]
    if damage == "absent": data["curation_records"] = []
    elif damage == "ambiguous": data["curation_records"].append(deepcopy(row))
    elif damage == "unbound": row["usable_for_attribution"] = False
    elif damage == "sha": row["protein_sha256"] = "f"*64
    elif damage == "source_sha": row["source_ref"]["sha256"] = "f"*64
    elif damage == "scope": row["rm_association_class"] = row["native_fields"]["rm_association_class"] = None
    elif damage == "query": row["gene_id"] = "different_gene"
    elif damage == "projection": row["native_fields"]["rm_association_class"] = "classic_RM_supported"
    assert one(data)["action"] == "retain_original_unassessed"


@pytest.mark.parametrize("protection", ["R", "RM", "S", "coherent_context", "curation_conflict", "classic_RM", "independent_review", "native_RM"])
def test_real_RM_role_or_context_is_preserved(protection):
    data = fixture(); decision=data["saved_decisions"][0]; curation=data["curation_records"][0]
    if protection in {"R","RM","S"}: decision["role"] = protection
    elif protection == "coherent_context": decision["native_context"] = {"typeii_operon_supported":"TRUE"}
    elif protection == "curation_conflict":
        curation["evidence_conflict"] = True; curation["native_fields"]["evidence_conflict"] = "TRUE"
    elif protection == "classic_RM":
        curation["rm_association_class"] = curation["native_fields"]["rm_association_class"] = "classic_RM_supported"
    elif protection == "independent_review":
        data["annotations"] = [{**decision, "status":"supported", "assignment_kind":"reviewed_RM_component", "rule_id":"independent"}]
    elif protection == "native_RM":
        data["calls"].append({**decision, "evidence_kind":"final_member", "category":"R-M", "caller":"Independent"})
    assert one(data)["action"] == "retain_dual_or_conflicting"
    assert project_family_attribution(data)["saved_decisions"] == data["saved_decisions"]


def test_orphan_chemistry_review_and_explicit_memberships_project_with_provenance():
    data=fixture(); decision=data["saved_decisions"][0]
    data["annotations"]=[{**decision,"status":"supported","assignment_kind":"reviewed_RM_component","orphan":True}]
    data["options"]={"family_memberships":[{**decision,"family_id":"R-M","category":"R-M"},
                                            {**decision,"family_id":"SomethingElse","category":"Other"}]}
    before=deepcopy(data); out=project_family_attribution(data)
    assert data==before
    assert {r["family_id"] for r in out["options"]["family_memberships"]} == {"BREX","SomethingElse"}
    assert out["annotations"][0]["status"]=="out_of_scope"
    archive=out["family_attribution_projection"]["original_records"]
    assert archive["annotations"]==before["annotations"] and len(archive["family_memberships"])==1


def test_unreviewed_family_labels_and_injected_transition_do_not_override():
    data=fixture(); data["calls"]=[]
    data["family_attribution"]=[dict(action="assign_noncanonical_scope",target_family="BREX",**data["saved_decisions"][0])]
    projected=project_family_attribution(data)
    assert projected["saved_decisions"]==data["saved_decisions"]
    assert projected["family_attribution"]==[]


def test_general_fallback_or_partial_pbe_labels_are_not_scope_vetoes():
    for family in ("DMS_other","PDC-S02","HEC-06","PbeA"):
        data=fixture()
        for call in data["calls"]:
            call["family"]=family;call["native_family"]=family
        assert apply_family_attribution(data)["family_attribution"]==[]


def test_independent_saved_context_and_hidden_conflict_cannot_be_erased():
    data=fixture(); decision=data["saved_decisions"][0]
    data["groups"]=[dict(sample_id=decision["sample_id"],family="Type II",context_state="met",
                         members=[dict(gene_id=decision["gene_id"],protein_sha256=decision["protein_sha256"])],
                         source_ref=decision["source_ref"])]
    assert one(data)["action"]=="retain_dual_or_conflicting"
    data.pop("groups")
    data["curation_records"][0]["native_fields"]["evidence_conflict"]="TRUE"
    assert one(data)["action"]=="retain_original_unassessed"


def test_repeated_protein_hash_never_moves_curation_between_genes():
    data=fixture(); curation=data["curation_records"][0]
    curation["gene_id"]="different_locus_with_identical_protein"
    assert one(data)["action"]=="retain_original_unassessed"


def test_prohibited_component_blocks_disarm_contract():
    data=fixture(1); call=deepcopy(data["calls"][0]); gene=deepcopy(data["genes"][0])
    call["gene_id"]=gene["gene_id"]="prohibited"
    gene.update(start=500,end=580)
    call["native_fields"].update(start="500",end="580")
    call["native_fields"]["protein.name"]="DrmD"
    call["native_fields"]["relative.position"]="5"
    data["calls"].append(call); data["genes"].append(gene)
    row=one(data)
    assert row["action"]=="retain_original_unassessed"
    assert "prohibited_native_component" in row["native_contexts"][0]["reason_codes"]


def test_conflicting_complete_noncanonical_contexts_keep_original_scope():
    data=fixture(0); other=fixture(2)
    target=data["saved_decisions"][0]; old=other["saved_decisions"][0]["gene_id"]
    for row in other["genes"]+other["calls"]:
        if row["gene_id"]==old:
            row["gene_id"]=target["gene_id"]; row["protein_sha256"]=target["protein_sha256"]
        if "system_id" in row: row["system_id"]="other-context"
    # Preserve canonical coordinates for the shared target and order partners.
    target_gene=next(g for g in data["genes"] if g["gene_id"]==target["gene_id"])
    for gene in other["genes"]:
        if gene["gene_id"]==target["gene_id"]: continue
        offset=-100 if gene["gene_id"].endswith("HmaA") else 100
        gene.update(start=target_gene["start"]+offset,end=target_gene["end"]+offset)
        data["genes"].append(gene)
    by_gene={g["gene_id"]:g for g in data["genes"]}
    for call in other["calls"]:
        gene=by_gene[call["gene_id"]]
        call["native_fields"].update(start=str(gene["start"]),end=str(gene["end"]))
    data["calls"].extend(other["calls"])
    row=one(data)
    assert row["action"]=="retain_dual_or_conflicting"
    assert "conflicting_noncanonical_families" in row["reason_codes"]


def test_producer_separation_counts_intervening_proteins_not_coordinate_gap():
    data=fixture()
    for call in data["calls"][1:]:
        call["native_fields"]["relative.position"]=str(int(call["native_fields"]["relative.position"])+3)
    assert one(data)["action"]=="assign_noncanonical_scope"  # difference 4, three intervening proteins
    for call in data["calls"][1:]:
        call["native_fields"]["relative.position"]=str(int(call["native_fields"]["relative.position"])+1)
    assert one(data)["action"]=="retain_original_unassessed"

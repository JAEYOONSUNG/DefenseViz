from copy import deepcopy
from hashlib import sha256

import pytest

from defenseviz.automatic_neighborhood import cognate_subjects, divergent_geometry, enrich_neighborhood_evidence

GENE_REF = {"path":"genes.tsv","sha256":"a"*64}
FEATURE_REF = {"path":"curated.tsv","sha256":"b"*64}


def pair_data(gap=800, intervening=0):
    gids = ["r",*[f"spacer{i}" for i in range(intervening)],"m"]
    genes = []
    for i,gid in enumerate(gids):
        start = 1 if gid=="r" else 101+gap if gid=="m" else 110+10*i
        genes.append(dict(sample_id="s",gene_id=gid,contig="c",start=start,end=start+99 if gid in {"r","m"} else start+5,
                          strand="-" if gid=="r" else "+",protein_length=300,protein_sha256=sha256(gid.encode()).hexdigest(),
                          source_refs=[{**GENE_REF,"row":i+2}]))
    features=[];records=[]
    for gid,role,subject in (("r","R","R.Example"),("m","M","M.Example")):
        gene=next(g for g in genes if g["gene_id"]==gid)
        f=dict(family_id="Type II",enzyme_role=role,partial_status="full_length",hit_label=subject,
               curated_blast_resolved_subject_id=subject,curated_blast_family="Type II",curated_blast_role=role,
               curated_blast_role_compatible="TRUE",curated_blast_family_compatible="TRUE",
               curated_blast_identity=.5,curated_blast_query_coverage=.8,curated_blast_reference_coverage=.8,curated_blast_evalue=1e-5)
        features.append(dict(sample_id="s",gene_id=gid,source_ref=FEATURE_REF,native_fields=f))
        records.append(dict(sample_id="s",gene_id=gid,family="Type II",role=role,protein_sha256=gene["protein_sha256"],protein_binding_status="verified",
                           source_refs=[FEATURE_REF],unknown_axes=["divergent_context"],axis_sources={},
                           evidence=dict(family="Type II",role=role,conflict=False,alternative_annotation=False,valid_sequence=True,divergent_context=None)))
    data=dict(samples=[{"sample_id":"s"}],genes=genes,suite=dict(features=features,reference_alignments=[],
        sources=[dict(sample_id="s",role="gene_annotations",selected=True,row_count=len(genes),**GENE_REF)]))
    compiled=dict(schema_version=1,records=records,candidate_scope=[dict(sample_id="s",status="complete",gene_ids=["r","m"])])
    return data,compiled


def values(result):
    return {r["gene_id"]:r["evidence"]["divergent_context"] for r in result["records"]}


def test_exact_frozen_boundary_and_one_intervening_gene():
    data,compiled=pair_data(intervening=1)
    before=deepcopy(compiled)
    result=enrich_neighborhood_evidence(data,compiled)
    assert values(result)=={"r":True,"m":True}
    assert compiled==before
    assert len(result["neighborhood_evaluation"]["pairs"])==1
    pair=result["neighborhood_evaluation"]["pairs"][0]
    assert pair["gap_bp"]==800 and pair["intervening_genes"]==1
    assert pair["activity_inferred"] is False
    assert all("divergent_context" not in r["unknown_axes"] for r in result["records"])


@pytest.mark.parametrize("gap,intervening",[(801,0),(800,2)])
def test_outside_prespecified_geometry_has_no_credit(gap,intervening):
    data,compiled=pair_data(gap=gap,intervening=intervening)
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":False,"m":False}


@pytest.mark.parametrize("strand",["+","?",None])
def test_head_to_head_required(strand):
    data,compiled=pair_data()
    data["genes"][0]["strand"]=strand
    v=values(enrich_neighborhood_evidence(data,compiled))
    assert v==({"r":False,"m":False} if strand=="+" else {"r":None,"m":None})


def test_different_contigs_and_overlaps_do_not_form_pair():
    data,compiled=pair_data()
    data["genes"][1]["contig"]="other"
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":False,"m":False}
    data,compiled=pair_data(gap=-1)
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":False,"m":False}


@pytest.mark.parametrize("r,m,expected",[("R.Example","M.Example",True),("Example","M.Example",True),("R.Example","M.Other",False),("R.Example_1","M.Example_2",False),("","M.Example",False)])
def test_cognate_record_relationship_is_exact(r,m,expected):
    assert cognate_subjects(r,m) is expected


def test_no_arbitrary_nearby_methyltransferase_partner():
    data,compiled=pair_data()
    data["suite"]["features"][1]["native_fields"]["curated_blast_resolved_subject_id"]="M.Unrelated"
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":False,"m":False}


@pytest.mark.parametrize("field,value",[("curated_blast_identity",.499), ("curated_blast_query_coverage",.799),("curated_blast_reference_coverage",.799),("curated_blast_evalue",1.001e-5),("curated_blast_family","Type I"),("curated_blast_role","RM"),("curated_blast_role_compatible","FALSE")])
def test_same_alignment_all_original_prerequisites_required(field,value):
    data,compiled=pair_data()
    data["suite"]["features"][0]["native_fields"][field]=value
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":False,"m":False}


def test_missing_alignment_axis_is_unknown_when_pair_could_qualify():
    data,compiled=pair_data()
    del data["suite"]["features"][0]["native_fields"]["curated_blast_reference_coverage"]
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}


def test_literal_full_length_required_even_when_partial_penalty_rescued():
    data,compiled=pair_data()
    data["suite"]["features"][0]["native_fields"]["partial_status"]="partial_or_short"
    compiled["records"][0]["evidence"]["partial"]=False
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":False,"m":False}


def test_missing_partial_or_conflict_assessment_stays_unknown():
    for missing in ("partial_status","conflict","alternative_annotation"):
        data,compiled=pair_data()
        if missing=="partial_status":del data["suite"]["features"][0]["native_fields"][missing]
        else:compiled["records"][0]["evidence"][missing]=None
        assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}


@pytest.mark.parametrize("mutation",["row_count","row_ref","coordinate","hash"])
def test_incomplete_gene_inventory_cannot_assert_presence_or_absence(mutation):
    data,compiled=pair_data()
    if mutation=="row_count":data["suite"]["sources"][0]["row_count"]=3
    elif mutation=="row_ref":data["genes"][1]["source_refs"][0]["row"]=2
    elif mutation=="coordinate":data["genes"][1]["start"]=None
    else:data["genes"][0]["source_refs"][0]["sha256"]="d"*64
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}


def test_missing_gene_inventory_receipt_does_not_claim_no_neighbor():
    data,compiled=pair_data(gap=1000)
    data["suite"]["sources"]=[]
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}


def test_incomplete_candidate_scope_prevents_no_pair_claim_but_not_observed_pair():
    data,compiled=pair_data()
    compiled["candidate_scope"][0]["status"]="partial"
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":True,"m":True}
    data["genes"][1]["start"]+=100
    data["genes"][1]["end"]+=100
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}


def test_other_types_and_fused_or_specificity_roles_are_inapplicable():
    for family,role in (("Type I","R"),("Type III","M"),("Type IV","R"),("Type II","RM"),("Type II","S")):
        data,compiled=pair_data()
        compiled["records"][0].update(family=family,role=role)
        assert values(enrich_neighborhood_evidence(data,compiled))["r"] is False


def test_protein_hash_mismatch_cannot_provide_neighborhood_evidence():
    data,compiled=pair_data()
    compiled["records"][0]["protein_sha256"]="d"*64
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}


def test_fresh_alignment_does_not_borrow_old_curated_strength():
    data,compiled=pair_data()
    for gene,role,subject in (("r","R","R.Example"),("m","M","M.Example")):
        sha=next(g["protein_sha256"] for g in data["genes"] if g["gene_id"]==gene)
        data["suite"]["reference_alignments"].append(dict(sample_id="s",gene_id=gene,operation="fresh_reference_alignment",command_exit_code=0,
            protein_sha256=sha,reference_sequence_sha256="e"*64,source_ref={"path":"fresh.tsv","sha256":"f"*64},
            reference_metadata_source_refs=[{"path":"db.faa","sha256":"1"*64}],
            native_fields=dict(subject_id=subject,identity=.9,evalue=1e-20,qstart=1,qend=300,sstart=1,send=300,reference_length=300,reference_family="Type II",reference_role=role)))
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":True,"m":True}
    data["suite"]["reference_alignments"][0]["native_fields"]["identity"]=.1
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}
    data["suite"]["reference_alignments"][0]["native_fields"].pop("reference_role")
    data["suite"]["reference_alignments"][0]["native_fields"]["identity"]=.9
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}


def add_fresh_alignment(data,gid="r",subject="R.Example",**native):
    gene=next(g for g in data["genes"] if g["gene_id"]==gid)
    row=dict(sample_id="s",gene_id=gid,operation="fresh_reference_alignment",command_exit_code=0,
        search_scope="saved_nominated_subjects",search_is_exhaustive=False,
        protein_sha256=gene["protein_sha256"],reference_sequence_sha256="e"*64,
        source_ref={"path":"fresh.tsv","sha256":"f"*64,"rows":[1,2]},
        reference_metadata_source_refs=[{"path":"db.faa","sha256":"1"*64}],
        native_fields=dict(subject_id=subject,identity=.9,evalue=1e-20,qstart=1,qend=300,
            sstart=1,send=300,query_length=300,reference_length=300,query_coverage=1.,
            reference_coverage=1.,dnmb_hsp_count=2,reference_family="Type II",reference_role="R",**native))
    data["suite"]["reference_alignments"].append(row)
    return row


def test_frozen_curated_subject_is_not_replaced_by_different_raw_hit():
    data,compiled=pair_data()
    data["suite"]["features"][0]["native_fields"]["hit_label"]="R.Other"
    add_fresh_alignment(data,subject="R.Other")
    result=enrich_neighborhood_evidence(data,compiled)
    assert values(result)=={"r":True,"m":True}
    assert result["neighborhood_evaluation"]["pairs"][0]["r_subject"]=="R.Example"
    own=result["neighborhood_evaluation"]["candidates"][0]["eligibility"]
    assert own["source_kind"]=="stored_curated_alignment"


def test_fresh_curated_subject_is_selected_instead_of_raw_subject_alignment():
    data,compiled=pair_data()
    data["suite"]["features"][0]["native_fields"]["hit_label"]="R.Other"
    add_fresh_alignment(data,subject="R.Other")
    add_fresh_alignment(data)
    result=enrich_neighborhood_evidence(data,compiled)
    assert values(result)=={"r":True,"m":True}
    own=result["neighborhood_evaluation"]["candidates"][0]["eligibility"]
    assert own["source_kind"]=="fresh_reference_alignment"
    assert own["subject"]=="R.Example"


@pytest.mark.parametrize("failure",["no_alignment","low_identity","wrong_role","wrong_family"])
def test_subject_limited_failure_cannot_exclude_frozen_neighborhood(failure):
    data,compiled=pair_data()
    row=add_fresh_alignment(data)
    if failure=="no_alignment":row["alignment_status"]="no_alignment"
    elif failure=="low_identity":row["native_fields"]["identity"]=.2
    elif failure=="wrong_role":row["native_fields"]["reference_role"]="M"
    else:row["native_fields"]["reference_family"]="Type IV"
    result=enrich_neighborhood_evidence(data,compiled)
    assert values(result)=={"r":None,"m":None}
    own=result["neighborhood_evaluation"]["candidates"][0]["eligibility"]
    assert own["subject"] is None
    assert own["measured_subject"]=="R.Example"


def test_failed_raw_subject_cannot_use_noncognate_name_to_exclude_pair():
    data,compiled=pair_data()
    f=data["suite"]["features"][0]["native_fields"]
    for name in list(f):
        if name.startswith("curated_blast_"):f.pop(name)
    f["hit_label"]="R.Other"
    row=add_fresh_alignment(data,subject="R.Other")
    row["native_fields"]["identity"]=.1
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}


def test_positive_unselected_nominee_cannot_exclude_other_cognate_references():
    data,compiled=pair_data()
    f=data["suite"]["features"][0]["native_fields"]
    for name in list(f):
        if name.startswith("curated_blast_"):f.pop(name)
    f["hit_label"]="R.Other"
    add_fresh_alignment(data,subject="R.Other")
    result=enrich_neighborhood_evidence(data,compiled)
    assert values(result)=={"r":None,"m":None}
    assert result["neighborhood_evaluation"]["candidates"][0]["eligibility"]["reference_selection"]=="nominated_subject"


@pytest.mark.parametrize("coverage,expected",[(.8,True),(.79,None),(None,None)])
def test_multiple_hsp_coverage_uses_interval_union_not_enclosing_span(coverage,expected):
    data,compiled=pair_data()
    row=add_fresh_alignment(data)
    if coverage is None:
        row["native_fields"].pop("query_coverage")
        row["native_fields"].pop("reference_coverage")
    else:
        row["native_fields"].update(query_coverage=coverage,reference_coverage=coverage)
    result=enrich_neighborhood_evidence(data,compiled)
    assert values(result)=={"r":expected,"m":expected}


def test_threshold_override_is_rejected():
    with pytest.raises(ValueError,match="frozen"):
        divergent_geometry({}, {}, 0, 1500)


def add_exact_translations(data,compiled,sequence="A"*300,product="restriction enzyme"):
    data["suite"]["annotations"]=[]
    for gene in data["genes"]:
        gene["protein_sha256"]=sha256(sequence.encode()).hexdigest()
        gene["protein_length"]=len(sequence)
        rec=next((r for r in compiled["records"] if r["gene_id"]==gene["gene_id"]),None)
        if rec:
            rec["protein_sha256"]=gene["protein_sha256"]
            rec["evidence"].update(partial=None,old_strong=True,query_coverage=1.,reference_coverage=1.)
            rec["unknown_axes"].append("partial")
        data["suite"]["annotations"].append(dict(sample_id="s",gene_id=gene["gene_id"],source_ref=GENE_REF,native_fields=dict(translation=sequence,product=product,note="")))


def test_missing_legacy_partial_is_freshly_evaluated_from_exact_sequence():
    data,compiled=pair_data()
    add_exact_translations(data,compiled)
    for row in data["suite"]["features"]:row["native_fields"].pop("partial_status")
    before=deepcopy(data)
    result=enrich_neighborhood_evidence(data,compiled)
    assert values(result)=={"r":True,"m":True}
    assert all(r["evidence"]["partial"] is False for r in result["records"])
    assert all(r["sequence_partial_assessment"]["status"]=="full_length" for r in result["records"])
    assert data==before


def test_short_sequence_rescue_changes_penalty_not_literal_neighborhood_status():
    data,compiled=pair_data()
    add_exact_translations(data,compiled,sequence="A"*149)
    for row in data["suite"]["features"]:row["native_fields"].pop("partial_status")
    result=enrich_neighborhood_evidence(data,compiled)
    assert values(result)=={"r":False,"m":False}
    assert all(r["evidence"]["partial"] is False for r in result["records"])
    assert all(r["sequence_partial_assessment"]["status"]=="partial_or_short" for r in result["records"])


def test_explicit_partial_annotation_is_not_rescued_by_full_reference():
    data,compiled=pair_data()
    add_exact_translations(data,compiled,product="partial restriction enzyme")
    for row in data["suite"]["features"]:row["native_fields"].pop("partial_status")
    result=enrich_neighborhood_evidence(data,compiled)
    assert all(r["evidence"]["partial"] is True for r in result["records"])
    assert values(result)=={"r":False,"m":False}


def test_short_sequence_unknown_reference_rescue_keeps_penalty_unknown():
    data,compiled=pair_data()
    add_exact_translations(data,compiled,sequence="A"*149)
    for row in data["suite"]["features"]:row["native_fields"].pop("partial_status")
    for row in compiled["records"]:row["evidence"]["old_strong"]=None
    result=enrich_neighborhood_evidence(data,compiled)
    assert all(r["evidence"]["partial"] is None for r in result["records"])


def test_partial_reconstruction_rejects_wrong_protein_source():
    data,compiled=pair_data()
    add_exact_translations(data,compiled)
    for row in data["suite"]["features"]:row["native_fields"].pop("partial_status")
    data["suite"]["annotations"][0]["native_fields"]["translation"]="WRONG"
    result=enrich_neighborhood_evidence(data,compiled)
    assert result["records"][0]["evidence"]["partial"] is None


def test_native_partial_length_strips_terminal_stop_only_after_raw_hash_binding():
    data,compiled=pair_data()
    add_exact_translations(data,compiled,sequence="A"*199+"*")
    for row in data["suite"]["features"]:row["native_fields"].pop("partial_status")
    result=enrich_neighborhood_evidence(data,compiled)
    methylase=next(r for r in result["records"] if r["role"]=="M")
    assessment=methylase["sequence_partial_assessment"]
    assert assessment["length_aa"]==199
    assert assessment["status"]=="partial_or_short"
    assert assessment["normalized_sequence_sha256"]==sha256(("A"*199).encode()).hexdigest()
    assert methylase["protein_sha256"]==sha256(("A"*199+"*").encode()).hexdigest()
    data["suite"]["annotations"][1]["native_fields"]["translation"]="A"*199
    rejected=enrich_neighborhood_evidence(data,compiled)
    assert rejected["records"][1]["evidence"]["partial"] is None


def test_native_full_length_under_unknown_old_role_is_not_borrowed():
    data,compiled=pair_data()
    add_exact_translations(data,compiled,sequence="A"*149)
    for row in data["suite"]["features"]:
        row["native_fields"].update(family_id="NA",enzyme_role="NA")
    result=enrich_neighborhood_evidence(data,compiled)
    assert all(r["sequence_partial_assessment"]["status"]=="partial_or_short" for r in result["records"])
    assert all(r["sequence_partial_assessment"]["native_hypothesis_compatible"] is False for r in result["records"])
    assert values(result)=={"r":False,"m":False}


def test_invalid_reference_identifier_is_not_a_scientific_negative():
    data,compiled=pair_data()
    data["suite"]["features"][0]["native_fields"]["curated_blast_resolved_subject_id"]=123
    assert values(enrich_neighborhood_evidence(data,compiled))=={"r":None,"m":None}
    assert cognate_subjects(123,"M.123") is False

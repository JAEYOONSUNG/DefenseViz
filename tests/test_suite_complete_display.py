"""Observed annotation counts are not hypothetical biological upper bounds."""
from copy import deepcopy
import json

import pytest

from defenseviz.analysis.saved import analyze_saved
from defenseviz.plot.publication import render_publication
from defenseviz.plot.publication.atlas import _type_display
from defenseviz.plot.publication.drawing import number
from defenseviz.plot.publication.manuscript_model import enrich_manuscript_model
from tests.plot.test_manuscript_model import fixture as bundle_fixture, REF
from tests.plot.test_manuscript_atlas import synthetic_result, render_pages, texts
from tests.plot.test_publication import fixture as publication_fixture


def computed_scope():
    return [dict(sample_id="sample",family=typ,status="complete",
                 scope="current_retained_groups_from_completed_context_scan",
                 source_refs=[REF]) for typ in ("Type I","Type III")]


@pytest.mark.parametrize("value,state,expected", [
    (0,"partial","0"), (12,"partial","12"), (0,"complete","0"),
    (12,"observed_lower_bound","12"), (0,"unassessed","—"),
    (None,"complete","—"), (0,"not_applicable","—"),
])
def test_numbers_show_observed_counts_without_repeated_qualifiers(value,state,expected):
    assert number(value,state)==expected


def test_computed_group_inventory_does_not_require_activity_or_manual_layout_coverage():
    bundle=bundle_fixture()
    bundle["options"]={}
    # An undecided candidate outside this exact observed context does not make
    # its three retained members or the count of computed groups unknowable.
    bundle["saved_decisions"][-1]["state"]="unassessed"
    result=analyze_saved(bundle)
    assert result["sample_summary"][0]["retained_counts_status"]=="partial"
    bundle["rm_evaluation"]={"current_context_review":{"observed_inventory_coverage":computed_scope()}}
    original=deepcopy(result)
    model=enrich_manuscript_model(bundle,result)["manuscript_atlas"]
    assert result==original
    comps=model["rows"][0]["compositions"]
    assert {key:(r["count"],r["status"]) for key,r in comps.items()}=={
        "I_MS":(0,"complete"),"I_RMS":(1,"complete"),
        "III_M":(0,"complete"),"III_MR":(0,"complete")}
    assert comps["I_RMS"]["R_support_all"] is None
    assert comps["I_RMS"]["inventory_coverage"]["source_refs"]==[REF]


@pytest.mark.parametrize("defect", ["duplicate","wrong_scope","wrong_sample","wrong_family","missing_provenance"])
def test_computed_inventory_requires_valid_analysis_scope(defect):
    bundle=bundle_fixture(); result=analyze_saved(bundle); records=computed_scope()
    if defect=="duplicate": records.append(deepcopy(records[0]))
    elif defect=="wrong_scope": records[0]["scope"]="whole_genome_activity"
    elif defect=="wrong_sample": records[0]["sample_id"]="not_in_inputs"
    elif defect=="wrong_family": records[0]["family"]="Type IV"
    else: records[0]["source_refs"]=[]
    bundle["rm_evaluation"]={"current_context_review":{"observed_inventory_coverage":records}}
    with pytest.raises(ValueError):enrich_manuscript_model(bundle,result)


def test_missing_context_scan_is_not_an_observed_zero():
    bundle=bundle_fixture();bundle["options"]={};result=analyze_saved(bundle)
    records=computed_scope();records[1]["status"]="missing"
    bundle["rm_evaluation"]={"current_context_review":{"observed_inventory_coverage":records}}
    comp=enrich_manuscript_model(bundle,result)["manuscript_atlas"]["rows"][0]["compositions"]["III_MR"]
    assert comp["count"]==0 and number(comp["count"],comp["status"])=="—"


def test_atlas_retains_counts_and_support_geometry_without_question_or_plus_marks():
    result=synthetic_result(samples=1,categories=[("Nuclease",2)])
    row=result["manuscript_atlas"]["rows"][0]
    for record in row["rm_roles"].values(): record["status"]="partial"
    for record in row["compositions"].values():record.update(count=0,status="complete")
    cells=list(row["cells"].values())
    cells[0]["metric_status"].update(associated="partial",effector="partial")
    cells[1]["metric_status"].update(associated="unassessed",effector="unassessed")
    pages=render_pages(result)
    for page in pages:
        labels=texts(page)
        assert not {"NA","?","+","0+","1+"}.intersection(labels)
        assert not any("unassessed" in t for t in labels)
        assert "M(S)-only groups" in labels and "Orphan M(S)" not in labels
        assert labels.count("0")>=4
        glyphs=page.ledger["glyphs"]
        assert glyphs[0]["inner_area_fraction"]==pytest.approx(2/8)
        assert glyphs[1]["inner_radius_mm"]==0
        assert glyphs[1]["outer_radius_mm"]>0
        page.fig.clear()


def test_orphan_parentheses_only_describe_documented_positive_subsets():
    row=publication_fixture(samples=1)["type_tables"][0]
    row.update(orphan_M=None,orphan_state={"M":"unassessed"})
    assert _type_display(row,"M")[0]=="1"
    row.update(orphan_M=1,orphan_state={"M":"complete"})
    assert _type_display(row,"M")[0]=="1 (1)"


def test_publication_distinguishes_missing_from_zero_with_one_scope_note(tmp_path):
    result=publication_fixture(samples=1,families=1)
    for row in result["inventory_rows"]["families"]+result["inventory_rows"]["categories"]:
        row["metric_status"]["effector"]="unassessed"
        row["effector_CDS"]=0
        row["source_count_status"]["DefenseFinder"]="partial"
    before=deepcopy(result)
    render_publication(result,tmp_path,sections=["comparison","type_benchmarks"])
    assert result==before
    layout=json.loads((tmp_path/"publication_layout.json").read_text())
    labels=[t["text"] for p in layout["pages"] for t in p["texts"]]
    assert not any("unassessed" in t or "(NA)" in t for t in labels)
    assert not {"NA","?","+","0+"}.intersection(labels)
    effectors=[t for p in layout["pages"] for t in p["texts"] if (t.get("source_ref") or {}).get("field")=="effector"]
    assert effectors and all(t["text"]=="—" for t in effectors)
    assert any("Counts: recorded annotations" in t for t in labels)
    caption=(tmp_path/"publication_caption.md").read_text()
    assert "unassessed" not in caption
    assert "M(S)-only composition does not establish orphanhood" in caption


def test_missing_overlap_is_explicit_and_does_not_render_zero_counts(tmp_path):
    result=publication_fixture(samples=1,families=1)
    for field in ("overlap_regions","priority_overlap_regions"):
        result[field].update(coverage_status="unassessed",sample_ids=[],excluded_sample_ids=["sample/0"])
    render_publication(result,tmp_path,sections=["comparison"])
    layout=json.loads((tmp_path/"publication_layout.json").read_text())
    labels=[t for p in layout["pages"] for t in p["texts"]]
    assert [t["text"] for t in labels if t["role"]=="missing_overlap"]==["Source overlap: unavailable"]
    assert not any(t["role"]=="venn_count" for t in labels)


def policy_display_bundle(supported=True):
    bundle=bundle_fixture()
    gene=next(g for g in bundle['genes'] if g['gene_id']=='r')
    saved=next(r for r in bundle['saved_decisions'] if r['gene_id']=='r')
    saved.update(protein_sha256=gene['protein_sha256'],origin='full_native_evidence_evaluation',
                 evaluation_id='e'*64,policy_sha256='b'*64,integrated_R_component_support=supported)
    result=analyze_saved(bundle)
    bundle['derived_artifact_sha256']={REF['path']:REF['sha256']}
    bundle['rm_evaluation']=dict(operation='fresh_native_evidence_selected_policy_evaluation',
        evaluation_id='e'*64,policy_sha256='b'*64,records=[dict(sample_id='sample',gene_id='r',
        protein_sha256=gene['protein_sha256'],state='retained',assignment_status='evaluated',
        evidence=dict(family='Type I',role='R'),decision=dict(keep=True,integrated_R_component_support=supported),
        source_refs=[REF])])
    return bundle,result


@pytest.mark.parametrize('supported',[True,False])
def test_exact_current_policy_predicate_supplies_R_badge_without_activity_or_reference_claim(supported):
    bundle,result=policy_display_bundle(supported);before=deepcopy(result)
    model=enrich_manuscript_model(bundle,result)['manuscript_atlas'];row=model['rows'][0]
    assert result==before
    assert row['compositions']['I_RMS']['R_support_all'] is supported
    assert row['compositions']['I_RMS']['R_support_status']=='complete'
    assert row['compositions']['I_RMS']['reference_any'] is None
    assert row['compositions']['I_MS']['R_support_all'] is None
    assert row['rm_roles']['I_R']['support_count']==int(supported)
    assert row['rm_roles']['I_R']['support_status']=='complete'
    group=next(g for g in model['group_records'] if g['composition_class']=='I_RMS')
    assert group['evidence_refs']==[REF]
    assert group['R_support_scope']=='all_current_retained_R_members_meet_selected_policy_predicate'
    assert result['global_summary']['effector_CDS']==0


@pytest.mark.parametrize('defect',['protein','source_hash','evaluation','predicate','role'])
def test_current_policy_badge_requires_bound_identity_and_exact_predicate(defect):
    bundle,result=policy_display_bundle();record=bundle['rm_evaluation']['records'][0]
    if defect=='protein':record['protein_sha256']='0'*64
    elif defect=='source_hash':bundle['derived_artifact_sha256'][REF['path']]='0'*64
    elif defect=='evaluation':bundle['rm_evaluation']['evaluation_id']='0'*64
    elif defect=='predicate':record['decision']['integrated_R_component_support']=False
    else:record['evidence']['role']='M'
    with pytest.raises(ValueError,match='R-support display'):enrich_manuscript_model(bundle,result)


def test_comparison_observed_zero_needs_explicit_analysis_count_availability(tmp_path):
    result=publication_fixture(samples=1,families=1)
    for row in result['inventory_rows']['categories']+result['inventory_rows']['families']:
        row['effector_CDS']=0;row['metric_status']['effector']='unassessed'
    result['inventory_rows']['categories'][0]['observed_metric_status']={'effector':'complete'}
    before=deepcopy(result)
    render_publication(result,tmp_path,sections=['comparison'])
    assert result==before
    layout=json.loads((tmp_path/'publication_layout.json').read_text())
    values={t['source_ref']['table']:t['text'] for p in layout['pages'] for t in p['texts']
            if t['role']=='numeric' and t['source_ref']['field']=='effector'}
    assert values=={'categories':'0','families':'—'}
    assert result['inventory_rows']['categories'][0]['metric_status']['effector']=='unassessed'

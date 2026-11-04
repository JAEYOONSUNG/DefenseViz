"""Synthetic checks for the original five-card vector workflow contract."""
import math
from copy import deepcopy
import xml.etree.ElementTree as ET

import matplotlib
import pytest

from defenseviz.plot.publication.drawing import Canvas, Fonts
from defenseviz.plot.publication.manuscript_workflow import draw_manuscript_workflow


def canvas(width=447, height=55.53):
    return Canvas(width, height, Fonts(), "workflow", "workflow-test")


def full_config():
    return {"manuscript_workflow": {"parameters": {"structure": 25, "max_divergent_gap_bp": 800,
            "retain_cutoff": 50}, "pt_pairs": [["SspBCD", "SspE / SspFGH"]]}}


def test_original_five_card_artwork_and_live_label_contract(tmp_path):
    c = canvas()
    meta = draw_manuscript_workflow(c, {}, full_config(), [0, 0, 447, 55.53])
    assert [r["number"] for r in meta["cards"]] == [1, 2, 3, 4, 5]
    assert len(meta["labels"]) == 48 and len(meta["glyphs"]) == 16
    assert len(meta["arrows"]) == 18 and all(a["compound_path"] for a in meta["arrows"])
    assert {g["kind"] for g in meta["glyphs"]} >= {
        "database_layered_cylinder", "CDS_directional_arrow", "domain_profile_and_motif_slots",
        "protein_ribbon_schematic", "native_residue_pair_distance", "final_score_threshold_tick",
        "DNA_duplex", "CRISPR_array_context", "cell_context"}
    titles = [r["text"] for r in meta["labels"] if r["text"][:1] in "12345"]
    assert titles == ["1  R-M evidence review", "2  Proposed validation", "3  Modification-associated review",
                      "4  Nuclease evidence", "5  Complementary defense priorities"]
    c.validate()
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        c.fig.savefig(tmp_path / "workflow.svg")
    tree = ET.parse(tmp_path / "workflow.svg")
    ns = {"s": "http://www.w3.org/2000/svg"}
    assert not tree.findall(".//s:image", ns)
    assert len(tree.findall(".//s:text", ns)) == 48
    assert len(c.ax.patches) > 40


def test_original_transforms_compaction_and_uniform_physical_scaling():
    c = canvas(500, 100)
    m = draw_manuscript_workflow(c, {}, full_config(), [4, 3, 279.4, 55.53*279.4/447])
    scale = 279.4/447
    assert m["scale"] == pytest.approx(scale)
    first = m["cards"][0]["bounds_mm"]
    assert first == pytest.approx([4+8*447/417.5*scale, 3+16*53/63*scale,
                                   126*447/417.5*scale, 45*53/63*scale])
    heading = next(r for r in m["labels"] if r["text"].startswith("1  "))
    assert heading["font_pt"] == pytest.approx(9.5*scale)
    assert heading["anchor_mm"][0] == pytest.approx(4+(8+3*126/156)*447/417.5*scale)
    structure_route = next(r for r in m["arrows"] if r["id"] == "rm_context_to_structure")
    assert len(structure_route["route_mm"]) == 4
    assert structure_route["route_mm"][-1][0] < structure_route["route_mm"][-2][0]
    assert m["arrows"][-4]["dashed"] and m["arrows"][-3]["dashed"]
    assert not m["arrows"][-2]["dashed"] and not m["arrows"][-1]["dashed"]


def test_missing_policy_does_not_invent_manuscript_criteria_or_pair_presence():
    c = canvas()
    m = draw_manuscript_workflow(c, {"samples": [{"sample_id": "arbitrary"}]}, {}, [0, 0, 447, 55.53])
    labels = {r["text"] for r in m["labels"]}
    assert {"Score (saved)", "Neighborhood review", "Structure evidence (saved)"} <= labels
    assert not {"Score >=50", "Neighborhood <=800 bp", "Structure evidence +25/0", "SspBCD"} & labels
    assert not m["activity_assessed"] and not m["assay_executed"] and not m["new_engine_runs"]
    assert not m["sample_counts_displayed"]
    assert len(m["cards"]) == 5
    c.validate()


def test_parameters_come_from_explicit_saved_values_and_configuration():
    c = canvas()
    result = {"selected_policy": {"parameters": {"structure": 7, "max_divergent_gap_bp": 123}},
              "policy": {"retain_cutoff": 9}}
    m = draw_manuscript_workflow(c, result, {"manuscript_workflow": {"parameters": {"structure": 7}}}, [0, 0, 447, 55.53])
    labels = {r["text"] for r in m["labels"]}
    assert {"Score >=9", "Neighborhood <=123 bp", "Structure evidence +7/0"} <= labels
    assert m["parameter_origins"]["structure"] == "selected_policy.parameters"


def test_config_can_fill_absent_criteria_without_replacing_saved_zero():
    result = {"selected_policy": {"parameters": {"structure": 0, "retain_cutoff": None}},
              "policy": {"parameters": {"structure": 0, "max_divergent_gap_bp": 123}}}
    m = draw_manuscript_workflow(canvas(), result,
        {"manuscript_workflow": {"parameters": {"structure": 0, "retain_cutoff": 9}}}, [0, 0, 447, 55.53])
    assert m["parameters"] == {"structure": 0, "max_divergent_gap_bp": 123, "retain_cutoff": 9}
    assert m["parameter_origins"]["retain_cutoff"] == "configuration"


@pytest.mark.parametrize("result,parameters,key", [
    ({"selected_policy": {"parameters": {"structure": 7}}}, {"structure": 8}, "structure"),
    ({"policy": {"max_divergent_gap_bp": 123}}, {"max_divergent_gap_bp": 456}, "max_divergent_gap_bp"),
    ({"selected_policy": {"parameters": {"structure": 7}}, "policy": {"structure": 8}}, {}, "structure"),
    ({"policy": {"structure": 7, "parameters": {"structure": 8}}}, {}, "structure"),
    ({"policy": {"retain_cutoff": 7}, "retain_cutoff": [8]}, {}, "retain_cutoff"),
    ({"retain_cutoff": 7}, {"retain_cutoff": 8}, "retain_cutoff"),
])
def test_conflicting_saved_or_configured_criteria_are_rejected(result, parameters, key):
    with pytest.raises(ValueError, match=f"Conflicting workflow {key}.*saved criteria"):
        draw_manuscript_workflow(canvas(), result,
            {"manuscript_workflow": {"parameters": parameters}}, [0, 0, 447, 55.53])


def test_matching_saved_values_from_all_supported_forms_are_accepted():
    result = {"selected_policy": {"parameters": {"retain_cutoff": 7}},
              "policy": {"retain_cutoff": 7.0}, "retain_cutoff": [7]}
    meta = draw_manuscript_workflow(canvas(), result,
        {"manuscript_workflow": {"parameters": {"retain_cutoff": 7}}}, [0, 0, 447, 55.53])
    assert meta["parameters"]["retain_cutoff"] == 7


def test_actual_original_provider_set_preserves_the_exact_reference_label():
    result = {"sources": [{"id": "REBASE", "kind": "reference_annotation"},
              {"id": "DefenseFinder"}, {"id": "DefensePredictor"}, {"id": "PADLOC"}]}
    meta = draw_manuscript_workflow(canvas(), result, full_config(), [0, 0, 447, 55.53])
    assert len(meta["labels"]) == 48
    label = next(x for x in meta["labels"] if x["text"] == "PAD | DF | DP | REBASE")
    assert label["native_anchor"] == [26, 44] and label["font_pt"] == 6.7
    assert [s["id"] for s in meta["providers"]] == ["PADLOC", "DefenseFinder", "DefensePredictor", "REBASE"]


def test_enriched_provider_selection_precedes_other_declared_sources():
    result = {"sources": [{"id": "Unselected"}, {"id": "HMMscan", "kind": "hmm_evidence"}],
              "manuscript_atlas": {"source_order": ["CustomNative"], "reference_sources": ["CustomRef"]}}
    c = canvas()
    meta = draw_manuscript_workflow(c, result, {}, [0, 0, 447, 55.53])
    assert meta["provider_label_lines"] == ["CustomNative | CustomRef"]
    assert meta["provider_origin"] == "manuscript_atlas"
    assert [s["id"] for s in meta["providers"]] == ["CustomNative", "CustomRef"]
    c.validate()


def test_hmm_only_sources_are_excluded_from_the_declared_provider_label():
    result = {"sources": [{"id": "NativeA", "kind": "native_final_member"},
              {"id": "ProfileOnly", "kind": "hmm_evidence"}, {"id": "ReferenceB", "kind": "reference_annotation"}]}
    meta = draw_manuscript_workflow(canvas(), result, {}, [0, 0, 447, 55.53])
    assert [s["id"] for s in meta["providers"]] == ["NativeA", "ReferenceB"]
    assert not any("ProfileOnly" in x["text"] for x in meta["labels"])


def test_missing_source_metadata_does_not_invent_the_four_reference_providers():
    for result in ({"samples": []}, {"sources": []},
                   {"sources": [{"id": "HMMscan", "kind": "hmm_evidence"}]}):
        meta = draw_manuscript_workflow(canvas(), result, {}, [0, 0, 447, 55.53])
        assert meta["providers"] == [] and meta["provider_label_lines"] == ["Sources not declared"]
    reference = draw_manuscript_workflow(canvas(), {}, {}, [0, 0, 447, 55.53])
    assert reference["provider_origin"] == "schematic_reference"
    assert reference["provider_label_lines"] == ["PAD | DF | DP | REBASE"]


def test_long_source_identity_has_a_bound_alias_and_keeps_the_full_identity():
    source = "A very long explicitly declared native defense provider identifier"
    c = canvas()
    meta = draw_manuscript_workflow(c, {"callers": [source]}, {}, [0, 0, 447, 55.53])
    assert meta["providers"][0]["id"] == source
    assert len(meta["providers"][0]["label"]) == 5
    c.validate()


def test_explicit_aliases_wrap_within_the_original_provider_label_area():
    c = canvas()
    result = {"sources": [{"id": "A"}, {"id": "B"}]}
    config = {"aliases": {"A": "First synthetic source", "B": "Second synthetic source"}}
    meta = draw_manuscript_workflow(c, result, config, [0, 0, 447, 55.53])
    assert len(meta["provider_label_lines"]) == 2
    assert [s["label"] for s in meta["providers"]] == list(config["aliases"].values())
    c.validate()


def test_unrenderable_provider_aliases_require_shorter_explicit_labels():
    with pytest.raises(ValueError, match="provide shorter source aliases"):
        draw_manuscript_workflow(canvas(), {"sources": ["A"]},
            {"aliases": {"A": "Very long provider annotation " * 15}}, [0, 0, 447, 55.53])


@pytest.mark.parametrize("parameters", [{"structure": math.nan}, {"retain_cutoff": True}, {"max_divergent_gap_bp": -1}])
def test_invalid_explicit_criteria_are_rejected(parameters):
    with pytest.raises(ValueError, match="explicit nonnegative finite"):
        draw_manuscript_workflow(canvas(), {}, {"manuscript_workflow": {"parameters": parameters}}, [0, 0, 447, 55.53])


@pytest.mark.parametrize("bounds", [[0, 0, -1, 5], [0, 0, 447, math.nan], [0, 0, 1]])
def test_invalid_bounds_are_rejected(bounds):
    with pytest.raises(ValueError, match="bounds"):
        draw_manuscript_workflow(canvas(), {}, {}, bounds)


def pair_rows(group="Explicit pair", left="Modifier (PT)", right="Effector (PT)"):
    return [{"id": group + "-left", "label": left, "N": 2, "visible": True,
             "category": "Modification", "display_group": "Modification", "partner_group": group, "partner_order": 0},
            {"id": group + "-right", "label": right, "N": 1, "visible": True,
             "category": "Nuclease", "display_group": "Modification", "partner_group": group, "partner_order": 1}]


def result_with_families(rows):
    return {"sources": [], "inventory_rows": {"families": rows}}


def test_observed_canonical_pt_display_metadata_restores_original_mini_boxes():
    rows = pair_rows("PT Ssp", "SspBCD (PT)", "SspE / SspFGH (PT)")
    rows[0]["id"], rows[1]["id"] = "SspABCD (PT)", "SspFGH (PT)"
    absent = pair_rows("PT Dnd", "DndABCDE", "DndFGH (PT)")
    absent[1].update(N=0, visible=False)
    result = result_with_families([rows[1], *absent, rows[0]])
    before = deepcopy(result)
    c = canvas()
    meta = draw_manuscript_workflow(c, result, {}, [0, 0, 447, 55.53])
    assert result == before
    assert meta["pt_pairs"] == [["SspBCD", "SspE / SspFGH"]]
    selection = meta["pt_pair_selection"]
    assert selection["origin"] == "observed_family_display_metadata"
    assert selection["selected"][0]["family_ids"] == ["SspABCD (PT)", "SspFGH (PT)"]
    assert selection["selected"][0]["labels"] == ["SspBCD (PT)", "SspE / SspFGH (PT)"]
    assert selection["omitted"][0]["partner_group"] == "PT Dnd"
    assert not any("Dnd" in label["text"] for label in meta["labels"])
    assert [box[2] for box in meta["pt_pair_label_boxes_mm"]] == pytest.approx([21*447/417.5, 25*447/417.5])
    for label in selection["selected"]:
        assert all(w <= available for w, available in zip(label["widths_mm"], label["available_widths_mm"]))
    assert not meta["activity_assessed"] and not meta["sample_counts_displayed"]
    c.validate()


def test_default_pairing_requires_authored_metadata_not_familiar_family_names():
    rows = [{"id": "SspBCD (PT)", "label": "SspBCD (PT)", "N": 3, "category": "Modification"},
            {"id": "SspFGH (PT)", "label": "SspE / SspFGH (PT)", "N": 2, "category": "Nuclease"}]
    meta = draw_manuscript_workflow(canvas(), result_with_families(rows), {}, [0, 0, 447, 55.53])
    assert meta["pt_pairs"] == []


@pytest.mark.parametrize("value", [0, None, math.nan, -1, True])
def test_default_pairing_never_invents_zero_or_unobserved_partner(value):
    rows = pair_rows()
    rows[1]["N"] = value
    meta = draw_manuscript_workflow(canvas(), result_with_families(rows), {}, [0, 0, 447, 55.53])
    assert meta["pt_pairs"] == []
    assert meta["pt_pair_selection"]["omitted"][0]["reason"] == "requires_two_observed_visible_families"


@pytest.mark.parametrize("edit,reason", [
    ({"visible": False}, "requires_two_observed_visible_families"),
    ({"partner_order": None}, "missing_or_ambiguous_partner_order"),
    ({"partner_order": 0}, "missing_or_ambiguous_partner_order"),
    ({"display_group": "Other mechanisms"}, "not_declared_modification_display_pair"),
    ({"label": None}, "missing_display_label"),
])
def test_incomplete_or_ambiguous_display_pairs_are_omitted_with_reason(edit, reason):
    rows = pair_rows(); rows[1].update(edit)
    meta = draw_manuscript_workflow(canvas(), result_with_families(rows), {}, [0, 0, 447, 55.53])
    assert meta["pt_pairs"] == []
    assert meta["pt_pair_selection"]["omitted"][0]["reason"] == reason


def test_two_arbitrary_authored_pairs_fit_original_template_without_function_inference():
    result = result_with_families(pair_rows("A", "Markase (PT)", "Cutase (PT)") +
                                  pair_rows("B", "Modifier B", "Partner B"))
    c = canvas()
    meta = draw_manuscript_workflow(c, result, {}, [0, 0, 279.4, 55.53*279.4/447])
    assert meta["pt_pairs"] == [["Markase", "Cutase"], ["Modifier B", "Partner B"]]
    assert "not system membership or activity" in meta["pt_pair_selection"]["scope"]
    c.validate()


def test_pair_group_order_follows_explicit_display_order_before_lexical_order():
    first = pair_rows("Z", "First mark", "First partner")
    second = pair_rows("A", "Next mark", "Next partner")
    for row in first: row["display_order"] = 1
    for row in second: row["display_order"] = 2
    meta = draw_manuscript_workflow(canvas(), result_with_families(second + first), {}, [0, 0, 447, 55.53])
    assert meta["pt_pairs"] == [["First mark", "First partner"], ["Next mark", "Next partner"]]


def test_duplicate_family_identity_cannot_make_a_pair_with_itself():
    rows = pair_rows(); rows[1]["id"] = rows[0]["id"]
    meta = draw_manuscript_workflow(canvas(), result_with_families(rows), {}, [0, 0, 447, 55.53])
    assert meta["pt_pairs"] == []
    assert meta["pt_pair_selection"]["omitted"][0]["reason"] == "duplicate_family_identity"


def test_more_than_two_authored_groups_does_not_select_a_silent_subset():
    rows = sum((pair_rows(group) for group in ("A", "B", "C")), [])
    result = result_with_families(rows); before = deepcopy(result)
    meta = draw_manuscript_workflow(canvas(), result, {}, [0, 0, 447, 55.53])
    assert result == before and meta["pt_pairs"] == []
    assert len(meta["pt_pair_selection"]["omitted"]) == 3
    assert {r["reason"] for r in meta["pt_pair_selection"]["omitted"]} == {"more_than_two_display_pairs_no_subset_selected"}
    assert "pt_pair_label_boxes_mm" not in meta


@pytest.mark.parametrize("config,result_pairs,expected,origin", [
    ({"pt_pairs": []}, [["Result left", "Result right"]], [], "configuration"),
    ({}, [], [], "result.workflow_pt_pairs"),
    ({"pt_pairs": [["Configured", "Partner"]]}, [["Result left", "Result right"]], [["Configured", "Partner"]], "configuration"),
    ({}, [["Result left", "Result right"]], [["Result left", "Result right"]], "result.workflow_pt_pairs"),
])
def test_explicit_pairs_and_empty_disable_precede_default_metadata(config, result_pairs, expected, origin):
    result = result_with_families(pair_rows()); result["workflow_pt_pairs"] = result_pairs
    meta = draw_manuscript_workflow(canvas(), result, {"manuscript_workflow": config}, [0, 0, 447, 55.53])
    assert meta["pt_pairs"] == expected and meta["pt_pair_selection"]["origin"] == origin


def test_derived_long_labels_are_not_truncated_or_drawn_outside_boxes():
    rows = pair_rows(left="A deliberately long unshortened modification-family label (PT)")
    meta = draw_manuscript_workflow(canvas(), result_with_families(rows), {}, [0, 0, 447, 55.53])
    assert meta["pt_pairs"] == []
    omitted = meta["pt_pair_selection"]["omitted"][0]
    assert omitted["reason"] == "labels_exceed_mini_boxes"
    assert omitted["labels"][0] == rows[0]["label"]
    assert omitted["display_labels"][0] == rows[0]["label"].removesuffix(" (PT)")


def test_explicit_labels_that_cannot_fit_require_shorter_display_configuration():
    with pytest.raises(ValueError, match="pt_pairs labels exceed.*shorter explicit display labels"):
        draw_manuscript_workflow(canvas(), {}, {"manuscript_workflow": {
            "pt_pairs": [["An explicitly configured label much too long for the miniature box", "Partner"]]}}, [0, 0, 447, 55.53])

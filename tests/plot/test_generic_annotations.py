"""Small-cohort presentation must not invent study-specific annotations."""
from copy import deepcopy
import math

from matplotlib.colors import to_rgba
from matplotlib.patches import Rectangle
import pytest

from defenseviz.plot.publication import DEFAULTS, _config
from defenseviz.plot.publication.comparison import comparison_pages, ordered_families
from defenseviz.plot.publication.drawing import Canvas, Fonts, display_label
from tests.plot.test_manuscript_atlas import render_pages, synthetic_result, texts
from tests.plot.test_publication import fixture as publication_fixture


def colored_rectangles(page, color):
    return [patch for patch in page.ax.patches if isinstance(patch, Rectangle)
            and patch.get_facecolor() == to_rgba(color)]


@pytest.mark.parametrize("sample_count", [2, 3])
def test_small_cohorts_keep_design_without_implied_target_or_phenotype(sample_count):
    result = synthetic_result(samples=sample_count)
    before = deepcopy(result)
    page = render_pages(result)[0]
    assert result == before
    assert DEFAULTS["highlight_label"] is None and DEFAULTS["sample_marker"] is None
    assert "Target strain" not in texts(page) and "Transformable" not in texts(page)
    # Legacy phenotype metadata is deliberately present in this fixture; it
    # must not cause an undeclared stripe or label in the generic renderer.
    assert not colored_rectangles(page, "#526A75")
    rows = [cell for cell in page.ledger["cells"] if cell["kind"] == "strain"]
    assert len(rows) == sample_count
    assert [row["highlight"] for row in rows] == [i in (0, 2) for i in range(sample_count)]
    assert [card["number"] for card in page.ledger["manuscript_workflow"]["cards"]] == [1, 2, 3, 4, 5]
    assert len(page.ledger["panel_bounds"]) == 14
    assert len(page.ledger["glyphs"]) == 48 * sample_count
    page.validate()


def test_explicit_highlight_label_only_describes_highlights_on_its_page():
    result = synthetic_result(samples=3)
    pages = render_pages(result, max_rows_per_page=1, highlight_label="Chosen isolates")
    assert ["Chosen isolates" in texts(page) for page in pages] == [True, False, True]
    assert all("Target strain" not in texts(page) for page in pages)
    for page in pages:
        page.validate()


def test_generic_marker_requires_literal_true_and_has_matching_current_page_key():
    result = synthetic_result(samples=5)
    values = [True, False, "yes", 1, None]
    rows = result["manuscript_atlas"]["rows"]
    for row, value in zip(rows, values):
        row["display"]["available_in_collection"] = value
    before = deepcopy(result)
    marker = {"field": "available_in_collection", "label": "Collection isolate", "color": "#31546E"}
    pages = render_pages(result, max_rows_per_page=1, sample_marker=marker)
    assert result == before
    assert ["Collection isolate" in texts(page) for page in pages] == [True, False, False, False, False]
    for index, page in enumerate(pages):
        row = next(cell for cell in page.ledger["cells"] if cell["kind"] == "strain")
        _, y, _, h = row["bounds_mm"]
        patches = colored_rectangles(page, marker["color"])
        body = [patch for patch in patches if math.isclose(patch.get_y(), y)
                and math.isclose(patch.get_height(), h)]
        keys = [patch for patch in patches if patch.get_y() > y + h]
        assert len(body) == len(keys) == (1 if index == 0 else 0)
        page.validate()


def test_legacy_phenotype_field_can_be_used_only_via_explicit_generic_marker():
    result = synthetic_result(samples=2)
    marker = {"field": "transformable", "label": "Supplied phenotype", "color": "#526A75"}
    page = render_pages(result, sample_marker=marker)[0]
    assert "Supplied phenotype" in texts(page)
    assert "Transformable" not in texts(page)
    assert len(colored_rectangles(page, marker["color"])) == 2  # one row and its key
    page.validate()


@pytest.mark.parametrize("options", [
    {"highlight_label": ""},
    {"highlight_label": "  "},
    {"highlight_label": False},
    {"sample_marker": {}},
    {"sample_marker": "phenotype"},
    {"sample_marker": {"field": "x", "label": "X", "color": "not a color"}},
    {"sample_marker": {"field": "", "label": "X", "color": "#31546E"}},
    {"sample_marker": {"field": "x", "label": "  ", "color": "#31546E"}},
    {"sample_marker": {"field": "x", "label": "X", "color": "#31546E", "unknown": True}},
])
def test_annotation_options_reject_ambiguous_or_invalid_configuration(options):
    with pytest.raises(ValueError, match="highlight_label|sample_marker"):
        _config("Synthetic result", options)


def test_missing_taxonomy_and_categories_do_not_create_empty_legend_sections():
    result = synthetic_result(samples=2, categories=[("Nuclease", 1)])
    for row in result["manuscript_atlas"]["rows"] + result["samples"]:
        row.pop("genus", None)
        row.pop("species", None)
    page = render_pages(result)[0]
    assert not {"Genus", "Species", "Other", "Mechanism not assigned yet"}.intersection(texts(page))
    assert "Family evidence" in texts(page)
    page.validate()
    # The category-key absence above must not be implemented by suppressing
    # those categories unconditionally.
    populated = render_pages(synthetic_result(samples=2))[0]
    assert {"Genus", "Species", "Other", "Mechanism not assigned yet"}.issubset(texts(populated))


def test_absence_key_requires_a_visible_complete_zero_on_the_current_page():
    result = synthetic_result(samples=3, categories=[("Nuclease", 1)])
    model = result["manuscript_atlas"]
    model["absence_label"] = "Example feature"
    family = model["families"][0]
    family["N"] = 8
    for row, value, state in zip(model["rows"], [0, 0, 8], ["complete", "partial", "complete"]):
        row["display"]["absence_family_ids"] = [family["id"]]
        cell = row["cells"][family["id"]]
        cell.update(N=value, associated_CDS=0, effector_CDS=0, counts_status=state)
        cell["metric_status"]["inventory"] = state
    pages = render_pages(result, max_rows_per_page=1)
    assert ["No Example feature" in texts(page) for page in pages] == [True, False, False]
    for page in pages:
        page.validate()


@pytest.mark.parametrize("labels", [
    ["Shared label", "Shared label"],
    ["Shared label", "Shared label", "Shared label"],
    ["Shared label", "Shared label", "pool-a"],
])
def test_colliding_family_labels_preserve_distinct_pools_and_resolve_secondary_collisions(labels):
    rows = [dict(id=f"pool-{chr(97+i)}", label=label, category="Nuclease",
                 N=i+2, source_counts={"Source": i+1}) for i, label in enumerate(labels)]
    rows.append(dict(id="pool-unaffected", label="Independent label", category="Nuclease", N=11))
    before = deepcopy(rows)
    rendered = ordered_families(rows)
    assert rows == before
    assert [row["id"] for row in rendered] == [row["id"] for row in before]
    assert [display_label(row) for row in rendered] == [row["id"] for row in before[:-1]] + ["Independent label"]
    assert [{key: value for key, value in row.items() if key != "_display_label"}
            for row in rendered] == before


def test_atlas_and_comparison_show_the_same_unique_family_identity_after_label_collision():
    result = publication_fixture(samples=2, families=3, sources=3)
    families = result["inventory_rows"]["families"]
    for row, label in zip(families, ["Shared label", "Shared label", "family/0"]):
        row.update(label=label, category="Nuclease")
    for cell in result["cells"]:
        family = next(row for row in families if row["id"] == cell["family_id"])
        cell.update(label=family["label"], category=family["category"])
    before = deepcopy(result)
    atlas = render_pages(result)[0]
    fonts = Fonts()
    pages = comparison_pages(result, fonts, _config("Synthetic comparison", {}),
                             lambda width, height, section: Canvas(width, height, fonts, section, 1))
    assert result == before
    expected = {row["id"]: row["id"] for row in families}
    atlas_labels = {meta["text"] for _, meta in atlas.text_artists if meta["rotation"] == 90}
    assert set(expected.values()) <= atlas_labels and "Shared label" not in atlas_labels
    shown = {meta["source_ref"]["row_id"]: meta["text"] for page in pages for _, meta in page.text_artists
             if (meta.get("source_ref") or {}).get("field") == "label"
             and meta["source_ref"]["row_id"] in expected}
    assert shown == expected
    glyphs = {(glyph["sample_id"], glyph["family_id"]): glyph["total"] for glyph in atlas.ledger["glyphs"]}
    assert glyphs == {(cell["sample_id"], cell["family_id"]): cell["N"] for cell in result["cells"]}
    comparison = {row["row_id"]: row for page in pages for row in page.ledger["cells"]
                  if row.get("table") == "families"}
    for family in families:
        assert comparison[family["id"]]["count"] == family["N"]
        assert comparison[family["id"]]["source_counts"] == family["source_counts"]
    atlas.validate()
    for page in pages:
        page.validate()

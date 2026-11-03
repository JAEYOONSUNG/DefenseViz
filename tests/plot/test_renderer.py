import copy
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from defenseviz.plot import render


def fixture(n=4, f=6, callers=3):
    names = [f"caller-{j}" for j in range(callers)]
    samples = [{"sample_id": f"sample/{i}", "label": f"Organism strain {i}", "highlight": i == 1} for i in range(n)]
    families = [{"family": f"Family {j}", "category": "R-M" if j == 0 else "Unresolved"} for j in range(f)]
    summaries = [{"sample_id": s["sample_id"], "union_cds": 10 + i, "rm_pool_cds": 5,
                  "retained_cds": 2, "unassessed_cds": 3,
                  "caller_counts": {name: i + j for j, name in enumerate(names)}} for i, s in enumerate(samples)]
    cells = [{"sample_id": s["sample_id"], **fam, "cds_count": (i + j) % 5, "system_count": 0}
             for i, s in enumerate(samples) for j, fam in enumerate(families)]
    components = [{"sample_id": s["sample_id"], "rm_type": typ, "role": role, "count": 2, "supported_count": 1}
                  for s in samples for typ in ("Type I", "Type II") for role in ("R", "M", "S")]
    return {"samples": samples, "callers": names, "families": families, "sample_summary": summaries,
            "cells": cells, "rm_components": components,
            "coverage": [{"sample_id": s["sample_id"], "caller": c, "status": "complete"} for s in samples for c in names],
            "decisions": [], "policy": {"retain_cutoff": 50, "min_capped_evidence_axes": 1}}


def read_layout(out):
    return json.loads((out / "figure2_layout.json").read_text(encoding="utf-8"))


def test_single_page_live_text_counts_geometry(tmp_path):
    result = fixture()
    before = copy.deepcopy(result)
    paths = render(result, tmp_path)
    assert len(paths) == 3
    assert result == before
    layout = read_layout(tmp_path)
    assert layout["page_count"] == 1
    bars = layout["pages"][0]["panel_b"]["bars"]
    assert {b["bar_height_mm"] for b in bars} == {1.65}
    for bar in bars:
        width = bar["track_bounds_mm"][2] - bar["track_bounds_mm"][0]
        assert bar["filled_width_mm"] == pytest.approx(width * bar["value"] / bar["scale_max"])
    xml = ET.parse(tmp_path / "figure2_ab.svg").getroot()
    texts = xml.findall(".//{http://www.w3.org/2000/svg}text")
    assert len(texts) > 50
    assert not xml.findall(".//{http://www.w3.org/2000/svg}image")
    assert any(t.text == "Organism strain 1" for t in texts)
    assert "#f3bcd2" in (tmp_path / "figure2_ab.svg").read_text(encoding="utf-8")


def test_rectangular_pagination_never_drops_samples_or_families(tmp_path):
    result = fixture(n=9, f=11)
    render(result, tmp_path, max_rows_per_page=4, max_families_per_page=5)
    layout = read_layout(tmp_path)
    assert layout["page_count"] == 9
    observed = {(c["sample_id"], c["family"]) for p in layout["pages"] for c in p["panel_b"]["family_cells"]}
    assert observed == {(s["sample_id"], f["family"]) for s in result["samples"] for f in result["families"]}
    assert len(list(tmp_path.glob("figure2_ab_page*.svg"))) == 9
    assert not (tmp_path / "figure2_ab.svg").exists()


def test_missing_caller_is_na_not_zero(tmp_path):
    data = fixture(n=1, f=2, callers=1)
    data["coverage"][0]["status"] = "missing"
    data["sample_summary"][0]["caller_counts"]["caller-0"] = None
    render(data, tmp_path)
    panel = read_layout(tmp_path)["pages"][0]["panel_b"]
    assert panel["bars"][0]["value"] is None
    assert all(c["count_status"] == "unassessed" for c in panel["family_cells"])


def test_explicit_zero_saved_native_calls(tmp_path):
    data = fixture(n=1, f=1, callers=1)
    render(data, tmp_path)
    cell = read_layout(tmp_path)["pages"][0]["panel_b"]["family_cells"][0]
    assert cell["cds_count"] == 0
    assert cell["count_status"] == "observed"


def test_empty_dataset(tmp_path):
    render(fixture(n=0, f=0, callers=0), tmp_path)
    assert read_layout(tmp_path)["page_count"] == 1


def test_long_labels_unicode_and_no_clipping(tmp_path):
    data = fixture(n=1, f=2)
    data["samples"][0]["label"] = "Organism " + "abcdefghijklmno " * 12 + "αβγ"
    data["families"][1]["family"] = "Family with a lengthy label for test " * 3
    data["cells"][1]["family"] = data["families"][1]["family"]
    render(data, tmp_path)
    page = read_layout(tmp_path)["pages"][0]
    w, h = page["size_mm"]
    for t in page["text"]:
        x1, y1, x2, y2 = t["bounds_mm"]
        assert x1 >= -.2 and y1 >= -.2 and x2 <= w + .2 and y2 <= h + .2


def test_many_callers_and_components_all_covered(tmp_path):
    data = fixture(n=1, f=2, callers=11)
    data["rm_components"] = [{"sample_id": "sample/0", "rm_type": f"Type {i}", "role": "R", "count": 1, "supported_count": 0} for i in range(19)]
    render(data, tmp_path)
    pages = read_layout(tmp_path)["pages"]
    assert set(c for p in pages for c in p["panel_b"]["callers"]) == set(data["callers"])
    assert len({tuple(c) for p in pages for c in p["panel_b"]["rm_components"]}) == 19


@pytest.mark.parametrize("kwargs", [{"max_rows_per_page": 0}, {"max_families_per_page": True}, {"title": ""}])
def test_rejects_invalid_settings(tmp_path, kwargs):
    with pytest.raises(ValueError):
        render(fixture(), tmp_path, **kwargs)
    assert not list(tmp_path.iterdir())


def test_rejects_duplicate_identity_and_negative_counts(tmp_path):
    data = fixture()
    data["cells"].append(data["cells"][0])
    with pytest.raises(ValueError, match="Duplicate"):
        render(data, tmp_path)
    data = fixture()
    data["cells"][0]["cds_count"] = -1
    with pytest.raises(ValueError, match="nonnegative"):
        render(data, tmp_path)


def test_no_overwrite(tmp_path):
    render(fixture(n=1, f=1), tmp_path)
    prior = (tmp_path / "figure2_ab.pdf").read_bytes()
    with pytest.raises(FileExistsError):
        render(fixture(n=1, f=1), tmp_path)
    assert prior == (tmp_path / "figure2_ab.pdf").read_bytes()


def test_wide_glyph_sample_label_stays_in_sample_column(tmp_path):
    data = fixture(n=1, f=2)
    data['samples'][0]['label'] = 'W' * 150
    render(data, tmp_path)
    page = read_layout(tmp_path)['pages'][0]
    labels = [t for t in page['text'] if t['role'] == 'sample_label']
    assert len(labels) == 1
    assert labels[0]['bounds_mm'][2] < 71
    body = page['panel_b']['rows'][0]['bounds_mm']
    assert labels[0]['bounds_mm'][1] >= body[1]
    assert labels[0]['bounds_mm'][3] <= body[3]


def test_bottom_labels_do_not_overlap_last_data_row(tmp_path):
    render(fixture(n=2, f=6), tmp_path)
    page = read_layout(tmp_path)['pages'][0]
    body_end = page['panel_b']['rows'][-1]['bounds_mm'][3]
    for t in page['text']:
        if t['role'] in {'family_label', 'rm_type_label'}:
            assert t['bounds_mm'][1] > body_end


def test_large_integer_counts_remain_legible_in_their_cells(tmp_path):
    data = fixture(n=1, f=1, callers=1)
    data['sample_summary'][0]['caller_counts']['caller-0'] = 100000
    data['sample_summary'][0]['union_cds'] = 100000
    data['cells'][0]['cds_count'] = 100000
    render(data, tmp_path)
    page = read_layout(tmp_path)['pages'][0]
    family_cell = page['panel_b']['family_cells'][0]['bounds_mm']
    number = next(t for t in page['text'] if t['role'] == 'family_count')['bounds_mm']
    assert number[0] >= family_cell[0] and number[2] <= family_cell[2]


def test_native_caller_headers_do_not_overlap(tmp_path):
    data = fixture(n=1, f=1, callers=3)
    names = ['DefenseFinder', 'PADLOC', 'DefensePredictor']
    data['callers'] = names
    data['sample_summary'][0]['caller_counts'] = {c: 1 for c in names}
    data['coverage'] = [{'sample_id': 'sample/0', 'caller': c, 'status': 'complete'} for c in names]
    render(data, tmp_path)
    page = read_layout(tmp_path)['pages'][0]
    bounds = [t['bounds_mm'] for t in page['text'] if t['role'] == 'caller_header']
    for a, b in zip(bounds, bounds[1:]):
        assert a[2] + .5 < b[0]
    assert page['panel_b']['rm_components'] == [['Type I', 'R'], ['Type I', 'M'], ['Type I', 'S'], ['Type II', 'R'], ['Type II', 'M'], ['Type II', 'S']]


def test_long_title_wraps_and_expands_header_without_collisions(tmp_path):
    title = 'A comparative defense landscape with a long publication title ' * 6
    render(fixture(n=1, f=1), tmp_path, title=title)
    page = read_layout(tmp_path)['pages'][0]
    text = {r['role']: r for r in page['text'] if r['role'] in {'title', 'subtitle', 'page_number'}}
    title_bounds = text['title']['bounds_mm']
    assert '\n' in text['title']['text']
    assert title_bounds[2] < text['page_number']['bounds_mm'][0] - 2
    assert title_bounds[3] < text['subtitle']['bounds_mm'][1]
    assert text['subtitle']['bounds_mm'][3] < page['panel_a']['bounds_mm'][1]
    assert page['panel_b']['bounds_mm'][1] > 60


def test_custom_policy_visible_label(tmp_path):
    data = fixture(n=1, f=1)
    data['policy_metadata'] = {'is_selected_preset': False}
    render(data, tmp_path)
    text = [r['text'] for r in read_layout(tmp_path)['pages'][0]['text']]
    assert '3  DNMB policy (custom)' in text


def test_labels_with_dollars_and_latex_are_literal_svg_text(tmp_path):
    literal = r'$bad_{$ \alpha'
    title = 'Literal title: ' + literal
    data = fixture(n=1, f=1)
    data['samples'][0]['label'] = literal
    data['families'][0]['family'] = literal
    data['cells'][0]['family'] = literal
    render(data, tmp_path, title=title)
    xml = ET.parse(tmp_path / 'figure2_ab.svg').getroot()
    texts = [''.join(t.itertext()) for t in xml.findall('.//{http://www.w3.org/2000/svg}text')]
    assert texts.count(literal) == 2
    assert title in texts
    layout = read_layout(tmp_path)
    labels = [t['text'] for t in layout['pages'][0]['text'] if t['role'] in {'sample_label', 'family_label', 'title'}]
    assert labels == [title, literal, literal]


def test_layout_write_explicit_utf8_under_legacy_default_encoding(tmp_path, monkeypatch):
    original = Path.write_text

    def legacy_write(path, data, *args, **kwargs):
        # Emulate Windows' legacy locale default without changing CI globally.
        kwargs.setdefault('encoding', 'cp1252')
        return original(path, data, *args, **kwargs)

    monkeypatch.setattr(Path, 'write_text', legacy_write)
    data = fixture(n=1, f=1)
    data['samples'][0]['label'] = 'Isolate αβγ'
    render(data, tmp_path)
    raw = (tmp_path / 'figure2_layout.json').read_bytes()
    assert 'Isolate αβγ'.encode('utf-8') in raw
    assert json.loads(raw.decode('utf-8'))['page_count'] == 1

"""Saved DP classification and explicit threshold provenance, with synthetic CDSs."""
import csv
import hashlib
import json

import pytest

from defenseviz.io.suite import load_suite


def write_table(path, headers, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(rows)


def fixture(tmp_path, rows, *, flags=False, overrides=None, metadata=None):
    folder = tmp_path / "isolate"
    write_table(folder / "dnmb_module_rebasefinder/rebasefinder_input.tsv",
                ["locus_tag", "translation"], [["a", "MAAA"], ["b", "MBBB"], ["c", "MCCC"]])
    headers = ["locus_tag", "mean_log_odds"] + (["predicted_general_defense"] if flags else [])
    scores = tmp_path / "scores.tsv"
    write_table(scores, headers, rows)
    module = {"calls": "scores.tsv", **(overrides or {})}
    if metadata is not None:
        (tmp_path / "status.json").write_text(json.dumps(metadata), encoding="utf-8")
        module["decision_metadata"] = "status.json"
    manifest = tmp_path / "suite.json"
    manifest.write_text(json.dumps({"schema_version": 1, "samples": [{
        "sample_id": "isolate", "path": "isolate", "modules": {"defensepredictor": module},
    }]}), encoding="utf-8")
    return manifest, scores


def selected(data):
    return {c["gene_id"] for c in data["calls"] if c["caller"] == "DefensePredictor"}


def dp_metadata(data):
    return next(x for x in data["source_metadata"] if x["provider"] == "DefensePredictor")


@pytest.mark.parametrize("operator, expected", [(">", {"c"}), (">=", {"b", "c"})])
def test_explicit_operator_preserves_threshold_equality(tmp_path, operator, expected):
    manifest, _ = fixture(tmp_path, [["a", 3.9], ["b", 4], ["c", 4.1]],
                          overrides={"score_threshold": 4, "score_operator": operator})
    data = load_suite(tmp_path, manifest=manifest)
    assert selected(data) == expected
    assert all(r["score_rule_explicit"] and r["selection_source"] == "explicit_score_rule"
               for r in data["suite"]["scores"])
    assert dp_metadata(data)["decision_rule"]["score_operator_source"] == "manifest"


def test_saved_flags_are_authoritative_even_without_operator_metadata(tmp_path):
    manifest, scores = fixture(tmp_path, [["a", -3, "True"], ["b", 8, "False"], ["c", 4, "False"]], flags=True)
    original = scores.read_bytes()
    data = load_suite(tmp_path, manifest=manifest)
    assert selected(data) == {"a"}
    assert [r["selected"] for r in data["suite"]["scores"]] == [True, False, False]
    assert all(r["selection_source"] == "native_predicted_general_defense" for r in data["suite"]["scores"])
    assert data["suite"]["conflicts"] == []  # An assumed rule is not verified producer truth.
    assert scores.read_bytes() == original
    assert data["input_sha256"][str(scores)] == hashlib.sha256(original).hexdigest()
    assert dp_metadata(data)["decision_rule"]["selection_precedence"] == "saved_flag_then_score_rule"


@pytest.mark.parametrize("flag", ["", "NA", "null", "-", "maybe", "2"])
def test_missing_or_invalid_saved_flag_is_not_silently_negative(tmp_path, flag):
    manifest, _ = fixture(tmp_path, [["a", 4, flag]], flags=True,
                          overrides={"score_threshold": 4, "score_operator": ">"})
    with pytest.raises(ValueError, match="predicted_general_defense"):
        load_suite(tmp_path, manifest=manifest)


def test_status_operator_and_hash_are_preserved(tmp_path):
    manifest, _ = fixture(tmp_path, [["a", 4], ["b", 5]], metadata={"threshold": 4, "operator": ">", "exit_code": 0})
    data = load_suite(tmp_path, manifest=manifest)
    assert selected(data) == {"b"}
    rule = dp_metadata(data)["decision_rule"]
    assert rule["threshold_source"] == rule["score_operator_source"] == "decision_metadata"
    ref = rule["decision_metadata_source_ref"]
    status = tmp_path / "status.json"
    assert ref["sha256"] == data["input_sha256"][str(status)] == hashlib.sha256(status.read_bytes()).hexdigest()
    assert any(s["role"] == "decision_metadata" and not s["native_result"] for s in data["suite"]["sources"])


@pytest.mark.parametrize("overrides", [{"score_threshold": 5}, {"score_operator": ">="}])
def test_source_rule_cannot_be_overwritten_by_disagreeing_manifest(tmp_path, overrides):
    manifest, _ = fixture(tmp_path, [["a", 4, "False"]], flags=True, overrides=overrides,
                          metadata={"threshold": 4, "operator": ">"})
    with pytest.raises(ValueError, match="contradicts.*decision_metadata"):
        load_suite(tmp_path, manifest=manifest)


@pytest.mark.parametrize("strict", [False, True])
def test_explicit_score_rule_disagreement_never_replaces_saved_flag(tmp_path, strict):
    manifest, _ = fixture(tmp_path, [["a", 5, "False"], ["b", 1, "True"]], flags=True,
                          metadata={"threshold": 4, "operator": ">"})
    if strict:
        with pytest.raises(ValueError, match="saved_dp_flag_score_rule_mismatch.*strict"):
            load_suite(tmp_path, manifest=manifest, strict_inputs=True)
        return
    data = load_suite(tmp_path, manifest=manifest)
    assert selected(data) == {"b"}
    conflicts = data["suite"]["conflicts"]
    assert len(conflicts) == 2
    assert all(c["code"] == "saved_dp_flag_score_rule_mismatch" and c["action"] == "native_saved_flag_preserved" for c in conflicts)
    assert [r["score_rule_selected"] for r in data["suite"]["scores"]] == [True, False]
    assert [r["selected"] for r in data["suite"]["scores"]] == [False, True]


def test_legacy_fallback_is_explicitly_an_assumption(tmp_path):
    manifest, _ = fixture(tmp_path, [["a", 4]])
    data = load_suite(tmp_path, manifest=manifest)
    assert selected(data) == {"a"}  # Preserve older DNMB wrapper behavior.
    score = data["suite"]["scores"][0]
    assert score["selection_source"] == "legacy_score_rule_assumption"
    assert score["threshold_source"] == "legacy_default_assumption_4.0"
    assert score["score_operator_source"] == "legacy_inclusive_assumption"
    assert score["score_rule_explicit"] is False
    assert any("legacy threshold/operator assumptions" in w for w in data["suite"]["warnings"])


@pytest.mark.parametrize("metadata", [{}, {"threshold": 4}, {"operator": ">"}, [],
                                     {"threshold": True, "operator": ">"},
                                     {"threshold": "4", "operator": ">"},
                                     {"threshold": 4, "operator": "=="}])
def test_malformed_decision_metadata_cannot_trigger_fallback(tmp_path, metadata):
    manifest, _ = fixture(tmp_path, [["a", 4]], metadata=metadata)
    with pytest.raises(ValueError, match="decision_metadata|score_threshold|score_operator"):
        load_suite(tmp_path, manifest=manifest)


@pytest.mark.parametrize("overrides", [{"score_operator": op} for op in ("<", "=", "GT", None, True, 1)] +
                         [{"score_threshold": value} for value in (True, None, "4")])
def test_invalid_explicit_rule_is_rejected_even_with_saved_flags(tmp_path, overrides):
    manifest, _ = fixture(tmp_path, [["a", 4, "False"]], flags=True, overrides=overrides)
    with pytest.raises(ValueError, match="score_threshold|score_operator"):
        load_suite(tmp_path, manifest=manifest)


def test_saved_flag_and_consistent_explicit_metadata_do_not_change_source_membership(tmp_path):
    manifest, _ = fixture(tmp_path, [["a", 4, "False"], ["b", 5, "True"], ["c", -2, "False"]], flags=True,
                          overrides={"score_threshold": 4.0, "score_operator": ">"},
                          metadata={"threshold": 4, "operator": ">"})
    data = load_suite(tmp_path, manifest=manifest, strict_inputs=True)
    assert selected(data) == {"b"}
    assert data["suite"]["conflicts"] == []
    assert data["suite"]["scores"][0]["selected"] is False
    assert dp_metadata(data)["decision_rule"]["saved_flag_present"] is True

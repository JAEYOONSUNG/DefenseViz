"""Only absent optional columns in valid 13-field domain TSVs are padded."""
from copy import deepcopy
import hashlib
import json

import pytest

from defenseviz.io.suite import attach_system_review_annotations, load_suite
from tests.io.test_suite import native
from tests.io.test_suite_review_scope import domain, interpro


def test_mixed_legacy13_and_current15_rows_preserve_all_core_observations(tmp_path):
    root = native(tmp_path)
    old = domain()[:13]
    current = domain("g_b", "MBBB")
    current[-2:] = ["GO:synthetic", "pathway:synthetic"]
    path = interpro(root, [old, current])
    data = load_suite(root)
    rows = [r for r in data["suite"]["profiles"] if r.get("namespace") == "InterProScan"]
    assert len(rows) == 2
    assert list(rows[0]["native_fields"].values())[:13] == old
    assert list(rows[0]["native_fields"].values())[13:] == ["", ""]
    assert rows[0]["absent_optional_fields"] == ["go_annotations", "pathway_annotations"]
    assert rows[0]["native_column_count"] == 13
    assert list(rows[1]["native_fields"].values()) == current
    assert "absent_optional_fields" not in rows[1]
    assert all(r["sequence_binding"] == "verified" and r["qualification"] == "unassessed" for r in rows)
    assert all(r["source_ref"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest() for r in rows)
    capability = next(c for c in data["suite"]["capabilities"] if "raw_interpro_parsed_rows" in c)
    assert capability["raw_interpro_input_column_counts"] == {"13": 1, "15": 1}
    assert capability["raw_interpro_legacy13_rows"] == 1
    assert capability["raw_interpro_status"] == "complete"
    transform, = data["suite"]["table_transforms"]
    assert transform["affected_rows"] == 1 and transform["original_core_fields_unchanged"]


@pytest.mark.parametrize("count", [11, 12, 14, 16])
def test_ambiguous_or_malformed_column_counts_are_not_repaired(tmp_path, count):
    root = native(tmp_path)
    row = (domain() + ["extra"])[:count]
    interpro(root, [row])
    with pytest.raises(ValueError, match=f"unsupported InterPro column count {count}"):
        load_suite(root)


@pytest.mark.parametrize("row,reason", [
    (domain(md5="0" * 32)[:13], "interpro_sequence_md5_mismatch"),
    (domain(length=5)[:13], "interpro_sequence_length_mismatch"),
    (domain(start=0)[:13], "invalid_interpro_coordinates"),
    (domain(end=5)[:13], "invalid_interpro_coordinates"),
    (domain(gene="unknown")[:13], "unbound_interpro_query"),
])
def test_legacy_format_does_not_weaken_identity_or_coordinate_checks(tmp_path, row, reason):
    root = native(tmp_path)
    interpro(root, [row])
    data = load_suite(root)
    assert not [r for r in data["suite"]["profiles"] if r.get("namespace") == "InterProScan"]
    record = next(r for r in data["suite"]["unmapped_records"] if r.get("namespace") == "InterProScan")
    assert record["reason"] == reason
    assert record["native_column_count"] == 13
    assert list(record["native_fields"].values())[:13] == row
    with pytest.raises(ValueError, match=reason):
        load_suite(root, strict_inputs=True)


def test_sites_still_require_their_exact_schema(tmp_path):
    root = native(tmp_path)
    site = ["g:a", hashlib.md5(b"MAAA").hexdigest(), "4", "SFLD", "synthetic", "1", "4", "1", "M", "1", "1", "site annotation"]
    interpro(root, [site + ["extra"]], sites=True)
    with pytest.raises(ValueError, match="unsupported InterPro column count 13"):
        load_suite(root)


def test_rereading_for_expanded_scope_is_idempotent_and_keeps_hash_gate(tmp_path):
    root = native(tmp_path)
    path = interpro(root, [domain()[:13]])
    data = load_suite(root)
    expected = deepcopy(data["suite"])
    attach_system_review_annotations(data, root)
    assert data["suite"]["profiles"] == expected["profiles"]
    assert data["suite"]["table_transforms"] == expected["table_transforms"]
    assert sorted(json.dumps(x, sort_keys=True) for x in data["suite"]["capabilities"]) == sorted(
        json.dumps(x, sort_keys=True) for x in expected["capabilities"])
    stable = deepcopy(data["suite"])
    attach_system_review_annotations(data, root)
    assert data["suite"] == stable
    path.write_text(path.read_text().replace("PFsynthetic", "PFchanged"))
    with pytest.raises(ValueError, match="changed from its pinned hash"):
        attach_system_review_annotations(data, root)

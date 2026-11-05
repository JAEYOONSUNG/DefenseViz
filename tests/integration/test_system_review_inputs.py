"""Real-file proof loading and failure isolation for optional system evidence."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from defenseviz.system_review_inputs import apply_system_review_inputs


def proof_fixture(tmp_path, evidence_newline=None):
    source = Path(__file__).resolve().parents[2] / "examples/reference_extension"
    dest = tmp_path / "proof"
    shutil.copytree(source, dest)
    evidence = dest / "synthetic_evidence.txt"
    if evidence_newline is not None:
        evidence.write_bytes(evidence.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", evidence_newline))
    candidates = json.loads((dest / "candidates.json").read_text())
    if evidence_newline is not None:
        # Only deliberately converted fixtures get new byte-specific hashes.
        # The unmodified public example must retain its checked-in bindings.
        for item in candidates:
            refs = [hit["qualification"] for hit in item["domain_hits"]] + item["site_maps"]
            for ref in refs:
                ref["source_sha256"] = hashlib.sha256((dest / ref["source"]).read_bytes()).hexdigest()
        (dest / "candidates.json").write_bytes((json.dumps(candidates) + "\n").encode("utf-8"))
    candidate = candidates[0]
    data = {"genes": [{"sample_id": candidate["sample_id"], "gene_id": candidate["query_id"],
                       "protein_sha256": candidate["sequence_sha256"]}],
            "inputs": [], "input_sha256": {}}
    pack = dest / "system-evidence.json"
    pack.write_text(json.dumps({"schema_version": 1, "registry": "registry.json", "candidates": "candidates.json"}))
    return data, pack


def test_unmodified_public_reference_example_keeps_checked_in_hashes(tmp_path):
    data, path = proof_fixture(tmp_path)
    original = Path(__file__).resolve().parents[2] / "examples/reference_extension"
    for name in ("candidates.json", "registry.json", "reference.faa", "synthetic_evidence.txt"):
        assert (path.parent / name).read_bytes() == (original / name).read_bytes()
    loaded = apply_system_review_inputs(data, evidence_path=path)
    assert loaded["system_review_evidence"]["source_files_verified"] is True


@pytest.mark.parametrize("evidence_newline", [b"\n", b"\r\n"], ids=["LF", "CRLF"])
def test_optional_domain_and_sites_are_verified_and_all_files_pinned(tmp_path, evidence_newline):
    data, path = proof_fixture(tmp_path, evidence_newline)
    loaded = apply_system_review_inputs(data, evidence_path=path)
    proof = loaded["system_review_evidence"]
    assert proof["source_files_verified"] is True
    family = proof["assessments"][0]["families"][0]
    assert family["domain_evidence"]["state"] == "supported"
    assert family["catalytic_site_evidence"]["state"] == "supported"
    assert "calls" not in loaded and "effector_evidence" not in loaded
    for name, sha in loaded["input_sha256"].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == sha
    assert {p.name for p in loaded["inputs"]} >= {
        "system-evidence.json", "registry.json", "candidates.json", "reference.faa", "synthetic_evidence.txt"}


@pytest.mark.parametrize("defect", ["wrong_gene", "wrong_protein", "changed_evidence", "changed_line_endings", "unknown_pack_field"])
def test_invalid_optional_proof_does_not_partly_attach(tmp_path, defect):
    data, path = proof_fixture(tmp_path)
    if defect == "wrong_gene":
        data["genes"][0]["gene_id"] = "different"
    elif defect == "wrong_protein":
        data["genes"][0]["protein_sha256"] = "0" * 64
    elif defect == "changed_evidence":
        (path.parent / "synthetic_evidence.txt").write_text("changed source")
    elif defect == "changed_line_endings":
        evidence = path.parent / "synthetic_evidence.txt"
        original = evidence.read_bytes()
        changed = original.replace(b"\r\n", b"\n") if b"\r\n" in original else original.replace(b"\n", b"\r\n")
        assert changed != original
        evidence.write_bytes(changed)
    else:
        pack = json.loads(path.read_text())
        pack["trust_unverified"] = True
        path.write_text(json.dumps(pack))
    before = deepcopy(data)
    with pytest.raises(ValueError):
        apply_system_review_inputs(data, evidence_path=path)
    assert data == before


def test_review_cli_exposes_extensions_without_changing_default_calls():
    from defenseviz.cli import parser
    args = parser().parse_args(["inspect-suite", "native", "--system-review-rules", "rules.json",
                               "--system-review-evidence", "proof.json"])
    assert args.system_review_rules == Path("rules.json")
    assert args.system_review_evidence == Path("proof.json")

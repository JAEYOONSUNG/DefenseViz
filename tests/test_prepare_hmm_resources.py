"""Resource preparation preserves explicit provider evidence and exact HMMs."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/prepare_hmm_resources.py"
spec = importlib.util.spec_from_file_location("prepare_hmm_resources", SCRIPT)
preparer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preparer)


def record(name, length=10, ga="GA 12.3 7.8;\n"):
    return f"HMMER3/f\nNAME {name}\nLENG {length}\n{ga}//\n"


def inputs(tmp_path, *, metadata=None, library=None, off="excluded\n"):
    hmm = tmp_path / "source.hmm"
    hmm.write_text(library or record("unhelpful_name") + record("fusion") + record("excluded"))
    meta = tmp_path / "restriction_metadata.csv"
    meta.write_text(metadata or "Name,System,Enzyme_type\nunhelpful_name,RM_Type_IIM,RE\nfusion,RM_Type_IIG,IIG\n")
    excluded = tmp_path / "off_target.txt"; excluded.write_text(off)
    return hmm, meta, excluded


def test_exact_library_explicit_metadata_and_output_are_bound(tmp_path):
    from defenseviz.resources import _profile_binding
    source = inputs(tmp_path)
    before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source}
    out = tmp_path / "prepared"
    receipt = preparer.prepare(*source, out, argv=["test", "--out", str(out)])
    rows = {r["profile_name"]: r for r in csv.DictReader((out/"profile_roles.tsv").open(), delimiter="\t")}
    assert rows["unhelpful_name"]["roles"] == "R" and rows["unhelpful_name"]["type_hints"] == "IIM"
    assert rows["fusion"]["roles"] == "M;R" and rows["fusion"]["type_hints"] == "IIG"
    assert rows["excluded"]["off_target"] == "true" and rows["excluded"]["roles"] == rows["excluded"]["type_hints"] == ""
    assert (out/"RM_HMMs.hmm").read_bytes() == source[0].read_bytes()
    assert not (out/"RM_HMMs.hmm").is_symlink() and not (out/"RM_HMMs.hmm").samefile(source[0])
    assert all(r["library_sha256"] == before[str(source[0])] and r["profile_length"] == "10" for r in rows.values())
    assert {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source} == before
    assert {r["path"]: r["sha256"] for r in receipt["inputs"].values()} == before
    assert receipt["inputs_unchanged"] and receipt["original_GA_cutoffs_preserved"] and receipt["detector_runs"] == 0
    assert receipt["argv"] == ["test", "--out", str(out)]
    for name, expected in receipt["outputs"].items():
        assert hashlib.sha256((out/name).read_bytes()).hexdigest() == expected
    assert _profile_binding(out/"RM_HMMs.hmm", out/"profile_roles.tsv")["profiles"] == 3
    assert json.loads((out/"preparation.json").read_text())["passed"] is True


@pytest.mark.parametrize("rows,match", [
    ("model,RM_Type_II,unknown\n", "unknown or ambiguous source Enzyme_type"),
    ("model,RM_Type_II,MT\nmodel,RM_Type_II,RE\n", "ambiguous source Enzyme_type"),
    ("model,RM_Type_II,RE\nmodel,RM_Type_III,RE\n", "ambiguous source System"),
    ("model,unrecognized_system,RE\n", "unknown or ambiguous source System"),
    ("other,RM_Type_II,RE\n", "Missing, unknown or ambiguous source Enzyme_type"),
    ("model,RM_Type_II,RE,unexpected\n", "Malformed source metadata row"),
])
def test_unmapped_or_ambiguous_nonexcluded_profiles_fail_without_inference(tmp_path, rows, match):
    source = inputs(tmp_path, metadata="Name,System,Enzyme_type\n" + rows, library=record("model"), off="")
    with pytest.raises(ValueError, match=match):
        preparer.prepare(*source, tmp_path/"out")
    assert not (tmp_path/"out").exists()


def test_duplicate_identical_metadata_and_unused_provider_profiles_are_audited(tmp_path):
    source = inputs(tmp_path, metadata="Name,System,Enzyme_type\nmodel,RM_Type_II,RE\nmodel,RM_Type_II,RE\nunused,RM_Type_I,MT\n",
                    library=record("model"), off="unrepresented_exclusion\n")
    result = preparer.prepare(*source, tmp_path/"out")
    rows = list(csv.DictReader((tmp_path/"out/profile_roles.tsv").open(), delimiter="\t"))
    assert len(rows) == 1 and rows[0]["metadata_lines"] == "2;3"
    assert result["metadata_profiles_not_in_library"] == ["unused"]
    assert result["off_target_names_not_in_library"] == ["unrepresented_exclusion"]


def test_explicit_type_prefix_alias_has_the_same_meaning(tmp_path):
    source = inputs(tmp_path, library=record("any_name"),
                    metadata="Name,System,Enzyme_type\nany_name,Type_IV,RE\n", off="")
    preparer.prepare(*source, tmp_path/"out")
    rows = list(csv.DictReader((tmp_path/"out/profile_roles.tsv").open(), delimiter="\t"))
    assert rows[0]["roles"] == "R" and rows[0]["type_hints"] == "IV"


@pytest.mark.parametrize("library,match", [
    (record("model") + record("model"), "unique names"),
    (record("model", ga=""), "native GA"),
    (record("model", ga="GA nan 1;\n"), "native GA"),
    (record("model", length=0), "positive integer LENG"),
    (record("model").replace("//\n", ""), "terminator"),
])
def test_invalid_hmm_contract_is_rejected_without_repairing_cutoffs(tmp_path, library, match):
    source = inputs(tmp_path, library=library, metadata="Name,System,Enzyme_type\nmodel,RM_Type_II,RE\n", off="")
    with pytest.raises(ValueError, match=match):
        preparer.prepare(*source, tmp_path/"out")
    assert not (tmp_path/"out").exists()


def test_existing_output_is_preserved(tmp_path):
    source = inputs(tmp_path)
    out = tmp_path/"out"; out.mkdir(); sentinel = out/"keep.txt"; sentinel.write_text("original")
    with pytest.raises(FileExistsError):
        preparer.prepare(*source, out)
    assert sentinel.read_text() == "original"


def test_standalone_cli_needs_no_database_download_or_package_install(tmp_path):
    source = inputs(tmp_path)
    out = tmp_path/"out"
    run = subprocess.run([sys.executable, str(SCRIPT), "--hmm", str(source[0]), "--metadata", str(source[1]),
                          "--off-target", str(source[2]), "--out", str(out)], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert json.loads(run.stdout) == dict(output=str(out), profiles=3, off_target_profiles=1, passed=True)

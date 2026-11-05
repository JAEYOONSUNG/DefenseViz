"""Make invented DNMBsuite exports for an unknown defense review demonstration.

No original biological sequences, model classifications, or activity claims are
provided. The rule and proteins intentionally exist only for software testing.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil


def table(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_example(output):
    output = Path(output).resolve()
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Example output already exists: {output}")
    output.mkdir(parents=True)
    samples = [
        ("example-supported", "SyntheticFutureDefense__Effector", "synthetic-models-1"),
        ("example-unversioned", "SyntheticFutureDefense__Effector", ""),
        ("example-unknown", "SyntheticUnreviewedDefense__Protein", "synthetic-models-1"),
    ]
    manifest = dict(schema_version=1, samples=[])
    for number, (sid, model, version) in enumerate(samples):
        folder = output / "runs" / sid
        protein = "M" + "ACDEFGHIKLMNPQRSTVWY" * (number + 2)
        coords = dict(locus_tag="synthetic-protein", contig="synthetic-contig",
                      start=1, end=len(protein)*3, direction="+")
        table(folder / "dnmb_module_rebasefinder/rebasefinder_input.tsv",
            ["locus_tag", "contig", "start", "end", "direction", "translation", "product"],
            [dict(coords, translation=protein, product="Invented demonstration protein")])
        table(folder / "dnmb_module_defensefinder/defensefinder_id_map.tsv",
            ["defensefinder_id", "locus_tag", "contig", "start", "end", "direction"],
            [dict(coords, defensefinder_id="native-synthetic-protein")])
        table(folder / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv",
            ["hit_id", "sys_id", "type", "subtype", "gene_name", "model_version", "model_fqn"],
            [dict(hit_id="native-synthetic-protein", sys_id="synthetic-system",
                type="SyntheticFutureDefense", subtype="SyntheticFutureDefense",
                gene_name=model, model_version=version,
                model_fqn="defense-finder-models/DefenseFinder/Synthetic/SyntheticFutureDefense")])
        manifest["samples"].append(dict(sample_id=sid, path=f"runs/{sid}",
                                       label=f"Invented isolate {number + 1}"))
    (output / "cohort.json").write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.copyfile(Path(__file__).with_name("rules.json"), output / "review-rules.json")
    return dict(root=output / "runs", manifest=output / "cohort.json", rules=output / "review-rules.json")


def verify_example(inputs, output):
    """Exercise native import and additive review with and without the rule."""
    from copy import deepcopy
    from defenseviz.analysis.saved import analyze_saved
    from defenseviz.analysis.system_review import review_all_systems
    from defenseviz.io.suite import load_suite
    from defenseviz.system_review_inputs import apply_system_review_inputs

    data = load_suite(inputs["root"], manifest=inputs["manifest"])
    baseline = analyze_saved(data)
    unknown = review_all_systems(data, baseline)
    apply_system_review_inputs(data, rules_path=inputs["rules"])
    before = deepcopy(baseline)
    assessed = review_all_systems(data, baseline, data["system_review_rules"])
    assert baseline == before, "Additive review mutated accepted scientific memberships"
    key = lambda r: (r["sample_id"], r["gene_id"])
    old = {key(r): r for r in unknown["gene_reviews"]}
    new = {key(r): r for r in assessed["gene_reviews"]}
    assert len(new) == len(old) == 3
    assert all(r["axes"]["component_role"]["status"] == "unassessed" for r in old.values())
    assert new[("example-supported", "synthetic-protein")]["axes"]["component_role"]["status"] == "supported"
    missing = new[("example-unversioned", "synthetic-protein")]["axes"]["component_role"]
    assert missing["status"] == "unassessed"
    assert "producer_version_unavailable" in missing["reason_codes"]
    assert new[("example-unknown", "synthetic-protein")]["axes"]["component_role"]["status"] == "unassessed"
    assert baseline["global_summary"]["native_union_CDS"] == 3
    count_keys = ("effector_CDS", "associated_CDS", "priority_CDS", "motif_CDS")
    assert all(baseline["global_summary"][k] == 0 for k in count_keys)
    assert len(assessed["detection_accounting"]) == len(data["calls"]) == 3
    output = Path(output)
    (output / "review-before.json").write_text(json.dumps(unknown, indent=2) + "\n")
    (output / "review-after.json").write_text(json.dumps(assessed, indent=2) + "\n")
    receipt = dict(result="SYSTEM_REVIEW_EXAMPLE_PASS", synthetic_only=True,
        analysis_kind="evaluation_of_saved_system_evidence", detector_runs=0,
        canonical_detections=len(data["calls"]), reviewed_genes=len(new),
        accepted_counts_unchanged={k: baseline["global_summary"][k] for k in count_keys},
        component_role_by_sample={sid: row["axes"]["component_role"]["status"] for (sid, _), row in new.items()},
        input_sha256=data["input_sha256"],
        rule_file_sha256=hashlib.sha256(inputs["rules"].read_bytes()).hexdigest())
    (output / "verification.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Fresh output directory")
    parser.add_argument("--verify", action="store_true", help="Run source-bound review and assert the three expected outcomes")
    args = parser.parse_args()
    inputs = make_example(args.output)
    if args.verify:
        receipt = verify_example(inputs, args.output)
        print(receipt["result"])
    else:
        print(json.dumps({k: str(v) for k, v in inputs.items()}, indent=2))

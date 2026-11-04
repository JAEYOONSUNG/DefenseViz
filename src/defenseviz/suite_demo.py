"""Entirely synthetic native-format DNMBsuite example; no manuscript data."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


def _table(path, fields, rows, delimiter="\t"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_suite_project(output_dir):
    """Write three invented suite folders, including a completed empty export."""
    out = Path(output_dir).absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Demo input already exists: {out}")
    out.mkdir(parents=True)
    manifest = {"schema_version": 1, "samples": []}
    protein_hashes = {}
    # Synthetic component IDs demonstrate native membership vs role evidence.
    specs = {
        "isolate-a": [
            ("r", "rm", "RM_Type_I", "Synthetic_R", "Type I", "R"),
            ("m", "rm", "RM_Type_I", "Synthetic_M", "Type I", "M"),
            ("s", "rm", "RM_Type_I", "Synthetic_S", "Type I", "S"),
            ("jet-a", "wadjet", "Wadjet_II", "Wadjet__JetA_II", None, None),
            ("jet-b", "wadjet", "Wadjet_II", "Wadjet__JetB_II", None, None),
            ("jet-c", "wadjet", "Wadjet_II", "Wadjet__JetC_II", None, None),
            ("jet-d", "wadjet", "Wadjet_II", "Wadjet__JetD_II", None, None),
            ("orphan-m", None, None, None, "Type III", "M"),
        ],
        "isolate-b": [
            ("type-iv", "rm-iv", "RM_Type_IV", "Synthetic_IV_R", "Type IV", "R"),
            ("ssp-f", "pt", "SspBCDE", "SspBCDE__SspF", None, None),
            ("ssp-g", "pt", "SspBCDE", "SspBCDE__SspG", None, None),
            ("ssp-h", "pt", "SspBCDE", "SspBCDE__SspH", None, None),
            ("reference-only", None, None, None, None, None),
        ],
        "isolate-c": [],
    }
    for number, (sid, records) in enumerate(specs.items()):
        folder = out / "runs" / sid
        manifest["samples"].append({"sample_id": sid, "path": f"runs/{sid}",
                                    "label": f"Example isolate {chr(65 + number)}",
                                    "highlight": number == 0, "cohort": "synthetic"})
        genes, mapping, df, pad, dp, rebase = [], [], [], [], [], []
        for i, (gid, system, family, model, rm_type, role) in enumerate(records):
            protein = "M" + "ACDEFGHIKLMNPQRSTVWY" * (i + 2) + "G" * number
            protein_hashes[(sid, gid)] = hashlib.sha256(protein.encode("ascii")).hexdigest()
            coordinates = {"locus_tag": gid, "contig": "synthetic-contig", "start": 100 + i * 1000,
                           "end": 100 + i * 1000 + 3 * len(protein) - 1, "direction": "+"}
            genes.append({**coordinates, "translation": protein, "product": "Synthetic demonstration protein"})
            mapping.append(coordinates)
            dp.append({"locus_tag": gid, "mean_log_odds": 6 if i % 2 == 0 and system else 1})
            if system:
                df.append({"hit_id": f"df-{gid}", "sys_id": system, "type": family,
                           "subtype": family, "gene_name": model})
                if gid not in {"jet-c", "ssp-h"}:
                    pad_family = {"Wadjet_II": "wadjet_type_II", "SspBCDE": "PT_SspABCD"}.get(family, family)
                    pad.append({"system.number": system, "seqid": "synthetic-contig", "system": pad_family,
                                "target.name": f"pad-{gid}", "hmm.name": f"Synthetic_{gid}", "protein.name": gid})
            if rm_type or gid == "reference-only":
                rebase.append({"query": gid, "curation_keep": "TRUE" if rm_type else "FALSE",
                               "final_family": rm_type or "Unresolved", "final_role": role or "",
                               "blast_match": f"Synthetic_reference_{gid}", "blast_identity": "80",
                               "curation_tier": "synthetic", "curation_score": "75" if rm_type else "10"})
        fields = ["locus_tag", "contig", "start", "end", "direction"]
        _table(folder / "dnmb_module_rebasefinder/rebasefinder_input.tsv", fields + ["translation", "product"], genes)
        _table(folder / "dnmb_module_defensefinder/defensefinder_id_map.tsv", ["defensefinder_id", *fields],
               [{"defensefinder_id": f"df-{r['locus_tag']}", **r} for r in mapping])
        _table(folder / "dnmb_module_padloc/padloc_query_proteins_id_map.tsv", ["query_id", *fields],
               [{"query_id": f"pad-{r['locus_tag']}", **r} for r in mapping])
        _table(folder / "dnmb_module_defensepredictor/defense_predictor_query_id_map.tsv", ["query_id", *fields],
               [{"query_id": f"dp-{r['locus_tag']}", **r} for r in mapping])
        _table(folder / "dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv",
               ["hit_id", "sys_id", "type", "subtype", "gene_name"], df)
        _table(folder / "dnmb_module_padloc/padloc_query_proteins_padloc.csv",
               ["system.number", "seqid", "system", "target.name", "hmm.name", "protein.name"], pad, ",")
        _table(folder / "dnmb_module_defensepredictor/defense_predictor_output.csv",
               ["locus_tag", "mean_log_odds"], dp, ",")
        _table(folder / "dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv",
               ["query", "curation_keep", "final_family", "final_role", "blast_match", "blast_identity", "curation_tier", "curation_score"], rebase)
    (out / "cohort.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    review = out / "synthetic-review.txt"
    review.write_text("Synthetic example only. The isolate-a r component has an invented explicit effector/site review. No experimental or biological observation is represented.\n", encoding="utf-8")
    native = out / "runs/isolate-a/dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv"
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    supplement = {"schema_version": 1,
                  "input_sha256": {native.relative_to(out).as_posix(): sha(native)},
                  "effector_evidence": [{"sample_id": "isolate-a", "gene_id": "r",
                     "protein_sha256": protein_hashes[("isolate-a", "r")], "status": "supported",
                     "assignment_kind": "individual_nuclease_effector", "rule_id": "synthetic_review",
                     "motif_state": "supported", "motif_grade": "synthetic_mapped_reference_site",
                     "source_refs": [{"path": review.name, "sha256": sha(review)}]}]}
    (out / "review.json").write_text(json.dumps(supplement, indent=2) + "\n", encoding="utf-8")
    layout = {"preview_png": True, "type_samples": ["isolate-a", "isolate-b"]}
    (out / "layout.json").write_text(json.dumps(layout, indent=2) + "\n", encoding="utf-8")
    (out / "README.txt").write_text("Entirely synthetic DNMBsuite-format data. No genomes, sequences or quantitative results from a manuscript are included. Reproduce with: defenseviz from-suite runs --manifest cohort.json --supplement review.json --layout layout.json --out ../rerun\n", encoding="utf-8")
    return {"root": out / "runs", "manifest": out / "cohort.json", "supplement": out / "review.json", "layout": out / "layout.json"}


def demo_suite(output_dir):
    from .suite_pipeline import run_suite
    out = Path(output_dir).absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Demo output already exists: {out}")
    inputs = write_suite_project(out / "inputs")
    return run_suite(inputs.pop("root"), out / "results", **inputs,
                     title="Synthetic defense annotation comparison", argv=["defenseviz", "demo", "--suite", "--out", str(out)])

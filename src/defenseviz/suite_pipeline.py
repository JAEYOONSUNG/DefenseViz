"""Transactional publication workflow for saved DNMBsuite directories."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import json
import shutil
import sys
import tempfile
import time

from . import __version__
from .pipeline import digest, write_json, _write_table, software_state, provenance_path


def _serializable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        if any(not isinstance(k, str) for k in value):
            raise ValueError("JSON object keys must be strings; encode gene identities as arrays")
        return {k: _serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serializable(v) for v in value]
    return value


def _portable(value, paths):
    """Relocate actual imported paths, preserving unrelated upstream metadata."""
    if isinstance(value, dict):
        return {paths.get(k, k): _portable(v, paths) for k, v in value.items()}
    if isinstance(value, list):
        return [_portable(v, paths) for v in value]
    return paths.get(value, value) if isinstance(value, str) else value


RESULT_TABLES = {
    "gene_ledger": ("sample_id", "gene_id", "priority", "effector", "associated"),
    "priority_ledger": ("sample_id", "gene_id", "selection_route"),
    "cells": ("sample_id", "family_id", "N", "associated_CDS", "effector_CDS", "counts_status"),
    "rm_components": ("sample_id", "rm_type", "role", "count"),
    "type_tables": ("sample_id", "rm_type"),
    "groups_selected": ("sample_id", "group_id"),
}


def _table_rows(value, name):
    if not isinstance(value, list) or any(not isinstance(row, dict) for row in value):
        raise ValueError(f"{name} must be a list of record objects")
    return value


def _review_summary(review):
    """Flat review index for spreadsheet inspection; detailed proof stays in JSON."""
    genes = {(r["sample_id"], r["gene_id"]): r for r in review["gene_reviews"]}
    rows = []
    for unit in review["units"]:
        members = [genes[tuple(k)] for k in unit["member_keys"]]
        row = {"sample_id": unit["sample_id"], "provider": unit["provider"],
               "native_families": "; ".join(unit["native_families"]),
               "native_defense_scope": unit["defense_scope"],
               "unit_kind": unit["unit_kind"], "coding_sequences": len(members),
               "selected_coding_sequences": len(unit["selected_member_keys"]),
               "accepted_effector_CDSs": sum(g["accepted_membership"]["effector"] for g in members),
               "members_with_role_evidence": sum(g["axes"]["component_role"]["status"] == "supported" for g in members),
               "members_with_nuclease_domain_evidence": sum(g["axes"]["nuclease_domain"]["status"] == "supported" for g in members),
               "members_with_motif_site_evidence": sum(g["axes"]["motif_site"]["status"] == "supported" for g in members),
               "review_execution": unit["review_execution"]["status"], "unit_id": unit["unit_id"]}
        for name, axis in unit["axes"].items():
            row[name + "_status"] = axis["status"]
            row[name + "_reasons"] = "; ".join(axis["reason_codes"])
        rows.append(row)
    return rows


def suite_warnings(data, result):
    """Expose import conflicts/coverage warnings as well as analysis warnings."""
    return list(dict.fromkeys([*data.get("warnings", []),
                               *data.get("suite", {}).get("warnings", []),
                               *result.get("warnings", [])]))


def prepare_suite(root, *, manifest=None, supplement=None, rm_evidence=None, rscript=None, strict_inputs=False,
                  system_review_rules=None, system_review_evidence=None):
    from .io.suite import load_suite
    data = load_suite(Path(root), manifest=manifest, rscript=rscript, strict_inputs=strict_inputs)
    if rm_evidence is not None:
        from .rm_evaluation import apply_rm_evidence
        data = apply_rm_evidence(data, rm_evidence)
    if supplement is not None:
        from .supplement import apply_supplement
        data = apply_supplement(data, supplement)
    if supplement is not None or rm_evidence is not None:
        from .io.suite import attach_system_review_annotations
        data = attach_system_review_annotations(data, root, strict_inputs=strict_inputs)
    if system_review_rules is not None or system_review_evidence is not None:
        from .system_review_inputs import apply_system_review_inputs
        data = apply_system_review_inputs(data, rules_path=system_review_rules, evidence_path=system_review_evidence)
    return data


def run_suite(root, output_dir, *, manifest=None, supplement=None, rm_evidence=None, rscript=None,
              strict_inputs=False, title="DefenseViz", layout=None, argv=None,
              system_review_rules=None, system_review_evidence=None):
    """Import native exports, aggregate saved evidence, publish a fresh directory."""
    from .analysis.saved import analyze_saved
    from .plot.publication import render_publication
    from .io.project import Reader

    started = time.monotonic()
    source_root = Path(root).resolve()
    out = Path(output_dir).absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output already exists; choose a fresh directory: {out}")
    if not isinstance(title, str) or not title.strip():
        raise ValueError("Title must be nonempty text")
    software_before = software_state()
    data = prepare_suite(source_root, manifest=manifest, supplement=supplement, rm_evidence=rm_evidence,
                         rscript=rscript, strict_inputs=strict_inputs,
                         system_review_rules=system_review_rules, system_review_evidence=system_review_evidence)
    import_seconds = time.monotonic() - started
    layout_options = {}
    if layout is not None:
        layout_path = Path(layout).resolve()
        reader = Reader(layout_path)
        layout_options = reader.json(layout_path)
        if not isinstance(layout_options, dict):
            raise ValueError("Layout configuration must be an object")
        if "title" in layout_options:
            raise ValueError("Set the figure title with --title, not the layout file")
        data.setdefault("inputs", []).extend(reader.inputs)
        data.setdefault("input_sha256", {}).update(reader.hashes)
    paths = sorted({Path(p).resolve() for p in data.get("inputs", [])})
    pinned = {str(p): digest(p) for p in paths}
    for p, sha in data.get("input_sha256", {}).items():
        if pinned.get(str(Path(p).resolve())) != sha:
            raise ValueError(f"Input changed during suite import: {p}")
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{out.name}.staging-", dir=out.parent))
    try:
        tick = time.monotonic()
        result = analyze_saved(data)
        if data.get("rm_evaluation"):
            result.setdefault("provenance", {})["rm_evaluation"] = data["rm_evaluation"]
            result["provenance"]["scoring_runs"] = data["rm_evaluation"]["scoring_runs"]
        from .plot.publication.manuscript_model import enrich_manuscript_model
        result = enrich_manuscript_model(data, result)
        from .knowledge import describe_systems
        result['system_information'] = describe_systems(data, result)
        analysis_seconds = time.monotonic() - tick
        path_map = {str(p): provenance_path(p, source_root) for p in paths}
        path_map[str(source_root)] = "."
        native_export = _portable(_serializable(data), path_map)
        analysis_export = _portable(_serializable(result), path_map)
        required = {*RESULT_TABLES, "inventory_rows", "availability", "global_summary", "warnings"}
        if not required <= set(result):
            raise ValueError(f"Analysis is missing required outputs: {sorted(required - set(result))}")
        if not isinstance(result["availability"], dict):
            raise ValueError("Analysis availability must be an object")
        write_json(stage / "native_bundle.json", native_export)
        write_json(stage / "analysis.json", analysis_export)
        review = analysis_export["system_review"]
        write_json(stage / "system_review.json", review)
        if "rm_evaluation" in native_export:
            write_json(stage / "rm_evaluation.json", native_export["rm_evaluation"])
        table_dir = stage / "tables"
        table_dir.mkdir()
        from .knowledge_output import write_system_information
        write_system_information(stage, analysis_export['system_information'])
        for name, headers in RESULT_TABLES.items():
            _write_table(table_dir / f"{name}.tsv", _table_rows(analysis_export[name], name), headers)
        review_tables = {
            "units": ("unit_id", "unit_kind", "sample_id"),
            "gene_reviews": ("sample_id", "gene_id"),
            "rule_evaluations": ("rule_id",),
            "pending_queue": ("queue_id", "sample_id", "target_kind", "target", "axis", "queue_class"),
            "detection_accounting": ("detection_id", "unit_id"),
            "sample_coverage": ("sample_id",),
            "quarantine": (),
        }
        for name, headers in review_tables.items():
            _write_table(table_dir / f"system_review_{name}.tsv",
                         _table_rows(review[name], "system_review." + name), headers)
        _write_table(table_dir / "system_review_summary.tsv", _review_summary(review),
                     ("sample_id", "provider", "native_families", "unit_kind", "coding_sequences"))
        if "context_assessments" in analysis_export:
            _write_table(table_dir / "context_assessments.tsv",
                         _table_rows(analysis_export["context_assessments"], "context_assessments"),
                         ("sample_id", "group_ids", "original_context_state", "context_state"))
        inventory = analysis_export["inventory_rows"]
        if not isinstance(inventory, dict) or not {"categories", "families"} <= set(inventory):
            raise ValueError("Analysis inventory_rows requires categories and families")
        for name in ("categories", "families"):
            _write_table(table_dir / f"{name}.tsv", _table_rows(inventory[name], name),
                         ("id", "label", "N", "associated_CDS", "effector_CDS", "counts_status"))
        _write_table(table_dir / "availability.tsv",
                     [{"endpoint": key, "state": value} for key, value in analysis_export["availability"].items()],
                     ("endpoint", "state"))
        write_json(table_dir / "availability.json", analysis_export["availability"])
        for name in ("calls", "coverage", "saved_decisions", "annotations", "effector_evidence", "groups", "complexes"):
            _write_table(table_dir / f"imported_{name}.tsv", _table_rows(native_export.get(name, []), name),
                         ("sample_id", "caller", "status") if name == "coverage" else ("sample_id", "gene_id"))
        tick = time.monotonic()
        figure_dir = stage / "figures"
        figure_dir.mkdir()
        # Reference context is exported alongside the analysis; it is not a
        # rendering input and must not change the publication model identity.
        figure_model = {key: value for key, value in analysis_export.items()
                        if key != 'system_information'}
        figures = render_publication(figure_model, figure_dir, title=title, **layout_options)
        render_seconds = time.monotonic() - tick
        resolved_figures = [Path(p).resolve() for p in figures]
        if not resolved_figures or len(set(resolved_figures)) != len(resolved_figures):
            raise RuntimeError("Publication renderer returned no figures or duplicate paths")
        for p in resolved_figures:
            if not p.is_file() or not p.is_relative_to(figure_dir.resolve()):
                raise RuntimeError(f"Renderer output must exist inside its figures directory: {p}")
        if not {".pdf", ".svg"} <= {p.suffix.lower() for p in resolved_figures}:
            raise RuntimeError("Publication renderer must produce both PDF and editable SVG")
        for p in paths:
            if digest(p) != pinned[str(p)]:
                raise ValueError(f"Input changed during suite publication: {p}")
        software_after = software_state()
        if software_before["source_sha256"] != software_after["source_sha256"]:
            raise ValueError("DefenseViz implementation changed during execution; run again from a stable checkout")
        receipt = {
            "schema_version": 1, "tool": "DefenseViz", "version": __version__,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "operation": ("saved_suite_evidence_reevaluation_aggregation_and_render" if rm_evidence is not None
                          else "saved_suite_result_import_aggregation_and_render"),
            "scoring_runs": data.get("rm_evaluation", {}).get("scoring_runs", 0),
            "detector_runs": 0, "structure_prediction_runs": 0, "threshold_fitting": False,
            "argv": argv if argv is not None else sys.argv, "exit_code": 0,
            "elapsed_seconds": time.monotonic() - started,
            "stages": [{"name": "native_import_and_rm_reevaluation" if rm_evidence is not None else "native_import", "elapsed_seconds": import_seconds},
                       {"name": "saved_evidence_aggregation", "elapsed_seconds": analysis_seconds},
                       {"name": "publication_render", "elapsed_seconds": render_seconds}],
            "input_path_base": "supplied suite root",
            "inputs": [{"path": provenance_path(p, source_root), "sha256": pinned[str(p)]} for p in paths],
            "software": software_before, "software_after": software_after,
            "source_metadata": native_export.get("source_metadata", []),
            "options": {"manifest": provenance_path(Path(manifest).resolve(), source_root) if manifest else None,
                        "supplement": provenance_path(Path(supplement).resolve(), source_root) if supplement else None,
                        "rm_evidence": provenance_path(Path(rm_evidence).resolve(), source_root) if rm_evidence else None,
                        "system_review_rules": provenance_path(Path(system_review_rules).resolve(), source_root) if system_review_rules else None,
                        "system_review_evidence": provenance_path(Path(system_review_evidence).resolve(), source_root) if system_review_evidence else None,
                        "rscript": rscript, "strict_inputs": strict_inputs, "title": title, "layout": layout_options},
            "conflicts": native_export.get("suite", {}).get("conflicts", []),
            "external_processes": _serializable(data.get("suite", {}).get("commands", [])),
            "command_path_scope": "Exact executed argv; machine-specific paths are retained for execution provenance.",
            "endpoint_definitions": result.get("endpoint_definitions", {}),
            "availability": result.get("availability", {}),
            "system_review": {"accounting": review["accounting"], "provenance": review["provenance"]},
            "system_information": {"accounting": analysis_export['system_information']['accounting'],
                                   "provenance": analysis_export['system_information']['provenance']},
            "warnings": suite_warnings(data, result),
            "scientific_scope": "Saved annotation/evidence only. Missing is unassessed; neither activity nor independent performance is inferred.",
            "outputs": {p.relative_to(stage).as_posix(): digest(p) for p in sorted(stage.rglob("*")) if p.is_file()},
        }
        write_json(stage / "receipt.json", receipt)
        if out.exists() or out.is_symlink():
            raise FileExistsError(f"Output appeared during execution: {out}")
        relative_figures = [p.relative_to(stage.resolve()) for p in resolved_figures]
        response = _serializable({"output": str(out), "samples": len(data["samples"]),
                                  "figures": [str(out / p) for p in relative_figures],
                                  "summary": result["global_summary"], "warnings": receipt["warnings"]})
        json.dumps(response, allow_nan=False)
        stage.rename(out)
        return response
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise

"""Transactional publication workflow for saved DNMBsuite directories."""
from __future__ import annotations

from datetime import datetime, timezone
from copy import deepcopy
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
                  system_review_rules=None, system_review_evidence=None, evaluation="auto",
                  reference_cache=None, evidence_workspace=None, engine_config=None):
    from .io.suite import load_suite
    if evaluation not in {"auto", "stored", "full"}:
        raise ValueError("evaluation must be auto, stored or full")
    from .native_pipeline import engine_config_path, read_engine_config
    configured = engine_config_path(engine_config) if evaluation != "stored" else None
    full_run = evaluation == "full" or (evaluation == "auto" and configured is not None and rm_evidence is None)
    if full_run:
        if configured is None:
            raise ValueError("Full evaluation requires --engine-config or DEFENSEVIZ_ENGINE_CONFIG with the DNMB runtime and reference/HMM resources")
        if rm_evidence is not None or reference_cache is not None:
            raise ValueError("Full execution uses its configured references; do not combine with rm_evidence or reference_cache")
        if evidence_workspace is None: raise ValueError("Full evaluation requires a private evidence workspace")
        resolved_engine = read_engine_config(configured)
        if rscript is None: rscript = resolved_engine["rscript"]
    data = load_suite(Path(root), manifest=manifest, rscript=rscript, strict_inputs=strict_inputs)
    original_import = deepcopy(data)
    if evaluation != "stored" and rm_evidence is None:
        from .io.suite_cache import attach_cached_evidence
        data = attach_cached_evidence(data, root, rscript=rscript)
    if rm_evidence is not None:
        from .rm_evaluation import apply_rm_evidence
        data = apply_rm_evidence(data, rm_evidence)
    if supplement is not None:
        from .supplement import apply_supplement
        data = apply_supplement(data, supplement)
    if supplement is not None or rm_evidence is not None or evaluation in {"auto", "full"}:
        from .io.suite import attach_system_review_annotations
        data = attach_system_review_annotations(data, root, strict_inputs=strict_inputs)
    if system_review_rules is not None or system_review_evidence is not None:
        from .system_review_inputs import apply_system_review_inputs
        data = apply_system_review_inputs(data, rules_path=system_review_rules, evidence_path=system_review_evidence)
    if full_run:
        from .native_pipeline import execute_full_pipeline
        from .full_evaluation import apply_full_evaluation
        engine_started = time.monotonic()
        full_stages = execute_full_pipeline(data, root, configured, Path(evidence_workspace) / "full_pipeline")
        full_execution = json.loads((Path(evidence_workspace) / "full_pipeline" / "execution.json").read_text())
        full_execution["elapsed_seconds"] = time.monotonic() - engine_started
        policy_started = time.monotonic()
        data = apply_full_evaluation(data, full_stages)
        data["rm_evaluation"]["elapsed_seconds"] = time.monotonic() - policy_started
        data["suite"]["full_execution"] = full_execution
        data["automatic_evaluation"] = True
    if reference_cache is not None:
        if evaluation != "auto" or rm_evidence is not None:
            raise ValueError("reference_cache requires automatic evaluation without explicit rm_evidence")
        if evidence_workspace is None:
            raise ValueError("Fresh reference evaluation requires an evidence workspace")
        from .automatic_reference import enrich_reference_evidence
        data = enrich_reference_evidence(data, reference_cache, evidence_workspace, rscript=rscript)
        data["derived_artifact_sha256"] = {
            str(Path(p).resolve()): digest(Path(p))
            for p in data.get("suite", {}).get("automatic_reference", {}).get("artifacts", [])}
    if evaluation == "auto" and rm_evidence is None and not full_run:
        from .automatic import apply_automatic_rm
        data = apply_automatic_rm(data)
        from .automatic_context import apply_current_context_review
        data = apply_current_context_review(data)
        from .automatic import refresh_decision_audit
        data = refresh_decision_audit(data)
        data["automatic_evaluation"] = True
    data["evaluation_execution"] = {
        "requested_mode": evaluation,
        "rm_mode": "full_native_execution" if full_run else "explicit_complete_evidence" if rm_evidence is not None else "automatic" if evaluation == "auto" else "stored_replay",
        "automatic_system_review": "enabled" if data.get("automatic_evaluation") else "skipped_explicit_review_path" if rm_evidence is not None else "skipped_stored_replay",
        "explicit_reviewed_outputs_preserved": rm_evidence is not None,
        "fresh_full_execution": full_run,
    }
    data["_native_import_bundle"] = original_import
    return data


def analyze_suite(data):
    """Run system review and scoped aggregation after optional R-M evaluation."""
    from .analysis.saved import analyze_saved
    result = analyze_saved(data)
    if data.get("automatic_evaluation"):
        from .automatic_systems_output import apply_system_review
        data, result = apply_system_review(data, result)
    return data, result


def run_suite(root, output_dir, *, manifest=None, supplement=None, rm_evidence=None, rscript=None,
              strict_inputs=False, title="DefenseViz", layout=None, argv=None,
              system_review_rules=None, system_review_evidence=None, evaluation="auto", reference_cache=None, engine_config=None):
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
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{out.name}.staging-", dir=out.parent))
    try:
        data = prepare_suite(source_root, manifest=manifest, supplement=supplement, rm_evidence=rm_evidence,
                             rscript=rscript, strict_inputs=strict_inputs,
                             system_review_rules=system_review_rules, system_review_evidence=system_review_evidence,
                             evaluation=evaluation, reference_cache=reference_cache, engine_config=engine_config,
                             evidence_workspace=stage / "evidence" / "reference")
        original_import = data.pop("_native_import_bundle", None) or dict(data)
        import_seconds = time.monotonic() - started
        layout_options = {"preview_png": True}
        if layout is not None:
            layout_path = Path(layout).resolve()
            reader = Reader(layout_path)
            supplied_layout = reader.json(layout_path)
            if not isinstance(supplied_layout, dict):
                raise ValueError("Layout configuration must be an object")
            if "title" in supplied_layout:
                raise ValueError("Set the figure title with --title, not the layout file")
            layout_options.update(supplied_layout)
            data.setdefault("inputs", []).extend(reader.inputs)
            data.setdefault("input_sha256", {}).update(reader.hashes)
        paths = sorted({Path(p).resolve() for p in data.get("inputs", [])})
        pinned = {str(p): digest(p) for p in paths}
        for p, sha in data.get("input_sha256", {}).items():
            if pinned.get(str(Path(p).resolve())) != sha:
                raise ValueError(f"Input changed during suite import: {p}")
        tick = time.monotonic()
        data, result = analyze_suite(data)
        if data.get("rm_evaluation"):
            # Execution timings belong to the evaluation receipt, not the
            # deterministic scientific model consumed by the renderer.
            result.setdefault("provenance", {})["rm_evaluation"] = {
                k: v for k, v in data["rm_evaluation"].items() if k != "elapsed_seconds"}
            result["provenance"]["scoring_runs"] = data["rm_evaluation"]["scoring_runs"]
        from .plot.publication.manuscript_model import enrich_manuscript_model
        result = enrich_manuscript_model(data, result)
        from .knowledge import describe_systems
        result['system_information'] = describe_systems(data, result)
        analysis_seconds = time.monotonic() - tick
        path_map = {str(p): provenance_path(p, source_root) for p in paths}
        path_map[str(source_root)] = "."
        for source in data.get("suite", {}).get("sources", []):
            label, resolved = source.get("path"), source.get("resolved_path")
            if label and resolved and str(Path(resolved).resolve()) in pinned:
                if pinned[str(Path(resolved).resolve())] != source.get("sha256"):
                    raise ValueError("Source registry hash conflicts with imported input")
                portable = provenance_path(Path(resolved).resolve(), source_root)
                if label in path_map and path_map[label] != portable:
                    raise ValueError("Ambiguous source label in publication provenance")
                path_map[label] = portable
        for artifact in stage.rglob("*"):
            if artifact.is_file():
                path_map[str(artifact.resolve())] = artifact.relative_to(stage).as_posix()
        native_export = _portable(_serializable(data), path_map)
        imported_export = _portable(_serializable(original_import), path_map)
        analysis_export = _portable(_serializable(result), path_map)
        required = {*RESULT_TABLES, "inventory_rows", "availability", "global_summary", "warnings"}
        if not required <= set(result):
            raise ValueError(f"Analysis is missing required outputs: {sorted(required - set(result))}")
        if not isinstance(result["availability"], dict):
            raise ValueError("Analysis availability must be an object")
        write_json(stage / "native_bundle.json", imported_export)
        write_json(stage / "analysis_bundle.json", native_export)
        write_json(stage / "analysis.json", analysis_export)
        review = analysis_export["system_review"]
        write_json(stage / "system_review.json", review)
        if "rm_evaluation" in native_export:
            write_json(stage / "rm_evaluation.json", native_export["rm_evaluation"])
        table_dir = stage / "tables"
        table_dir.mkdir()
        if "automatic_system_review" in analysis_export:
            automatic_review = analysis_export["automatic_system_review"]
            write_json(stage / "automatic_system_review.json", automatic_review)
            _write_table(table_dir / "automatic_gene_assessments.tsv", automatic_review["gene_assessments"],
                         ("sample_id", "gene_id", "review_execution"))
        if native_export.get("automatic_evaluation"):
            evaluation_rows = []
            for row in native_export["rm_evaluation"]["records"]:
                decision = row.get("decision") or {}
                facts = row.get("evidence") or {}
                evaluation_rows.append({
                    "sample_id": row["sample_id"], "gene_id": row["gene_id"],
                    "state": row["state"], "reason": row["reason"],
                    "rm_type": facts.get("family"), "role": facts.get("role"),
                    "score_min": decision.get("score_min"), "score_max": decision.get("score_max"),
                    "score_total": decision.get("score_total"),
                    "R_component_support": decision.get("integrated_R_component_support"),
                    "unknown_axes": row.get("unknown_axes", []),
                    "protein_sha256": row.get("protein_sha256"),
                    "source_refs": row.get("source_refs", []),
                    "assessment": decision.get("assessment"),
                })
            _write_table(table_dir / "rm_candidate_evaluation.tsv", evaluation_rows,
                         ("sample_id", "gene_id", "state", "reason"))
            for state in ("retained", "not_retained", "unassessed", "out_of_scope"):
                _write_table(table_dir / ("rm_" + state + ".tsv"),
                             [r for r in evaluation_rows if r["state"] == state],
                             ("sample_id", "gene_id", "state", "reason"))
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
            _write_table(table_dir / f"imported_{name}.tsv", _table_rows(imported_export.get(name, []), name),
                         ("sample_id", "caller", "status") if name == "coverage" else ("sample_id", "gene_id"))
        for name in ("saved_decisions", "annotations", "effector_evidence", "groups", "complexes", "coverage"):
            export_rows = _table_rows(native_export.get(name, []), name)
            if name == "saved_decisions" and native_export.get("automatic_evaluation"):
                scopes = {(r["sample_id"],r["gene_id"]):r.get("source_scope",[])
                          for r in native_export["rm_evaluation"]["records"]}
                export_rows = [dict(r,source_scope=scopes.get((r["sample_id"],r["gene_id"]),[])) for r in export_rows]
            _write_table(table_dir / f"analysis_{name}.tsv", export_rows,
                         ("sample_id", "gene_id"))
        from .suite_exports import export_suite_workbook
        annotation_exports = export_suite_workbook(native_export, analysis_export, stage)
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
        for path, expected in data.get("derived_artifact_sha256", {}).items():
            if digest(Path(path)) != expected:
                raise ValueError(f"Derived evidence changed during publication: {path}")
        software_after = software_state()
        if software_before["source_sha256"] != software_after["source_sha256"]:
            raise ValueError("DefenseViz implementation changed during execution; run again from a stable checkout")
        rm_meta = data.get("rm_evaluation", {})
        ref_meta = data.get("suite", {}).get("automatic_reference", {})
        full_meta = data.get("suite", {}).get("full_execution", {})
        rm_seconds = rm_meta.get("elapsed_seconds", 0)
        ref_seconds = ref_meta.get("elapsed_seconds", 0)
        full_seconds = full_meta.get("elapsed_seconds", 0)
        policy_decisions = rm_meta.get("policy_decisions",rm_meta.get("policy_evaluations",rm_meta.get("scoring_runs",0)))
        stages = [{"name": "native_import_and_context_preparation", "elapsed_seconds": max(0, import_seconds-rm_seconds-ref_seconds-full_seconds)}]
        if full_meta:
            stages.append({"name": "fresh_DNMB_HMM_BLAST_context_and_feature_execution", "elapsed_seconds": full_seconds,
                           "hmm_search_runs": full_meta["hmm_search_runs"], "homology_search_runs": full_meta["homology_search_runs"]})
        if ref_meta:
            stages.append({"name": "fresh_reference_alignment", "elapsed_seconds": ref_seconds,
                           "homology_search_runs": ref_meta.get("homology_search_runs",0)})
        if rm_meta:
            stages.append({"name": "current_native_evidence_policy_evaluation" if full_meta else "stored_feature_collection_and_policy_evaluation" if policy_decisions else "stored_feature_collection_and_unresolved_policy_bounds" if rm_meta.get("policy_evaluations") else "stored_feature_collection_and_native_decision_replay",
                           "elapsed_seconds": rm_seconds, "policy_evaluations": rm_meta.get("policy_evaluations",0),
                           "policy_decisions": policy_decisions,
                           "preserved_decisions": rm_meta.get("preserved_decisions",0)})
        stages += [{"name": "system_review_and_aggregation", "elapsed_seconds": analysis_seconds},
                   {"name": "publication_render", "elapsed_seconds": render_seconds}]
        receipt = {
            "schema_version": 1, "tool": "DefenseViz", "version": __version__,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "operation": ("full_DNMB_RM_execution_and_publication" if data.get("evaluation_execution",{}).get("fresh_full_execution") else
                          "saved_suite_with_fresh_reference_evaluation_and_render" if data.get("suite", {}).get("automatic_reference", {}).get("homology_search_runs") else
                          "saved_suite_evidence_reevaluation_aggregation_and_render" if policy_decisions
                          else "saved_suite_unresolved_evidence_bounds_and_replay_render" if rm_meta.get("policy_evaluations")
                          else "saved_suite_result_import_aggregation_and_render"),
            "scoring_runs": data.get("rm_evaluation", {}).get("scoring_runs", 0),
            "detector_runs": sum(full_meta.get("detector_runs",{}).values()), "structure_prediction_runs": 0, "threshold_fitting": False,
            "original_defense_tool_runs": 0,
            "hmm_search_runs": full_meta.get("hmm_search_runs",0),
            "homology_search_runs": full_meta.get("homology_search_runs",0) + ref_meta.get("homology_search_runs",0),
            "full_execution": _portable(_serializable(full_meta),path_map),
            "sequence_feature_evaluations": data.get("rm_evaluation", {}).get("compiler", {}).get("motif_scans", 0),
            "argv": argv if argv is not None else sys.argv, "exit_code": 0,
            "elapsed_seconds": time.monotonic() - started,
            "stages": stages,
            "evaluation_execution": data.get("evaluation_execution", {}),
            "policy_evaluations": rm_meta.get("policy_evaluations",0),
            "policy_decisions": policy_decisions,
            "stored_vs_auto_summary": native_export.get("rm_evaluation",{}).get("stored_vs_auto_summary"),
            "preserved_native_decisions": rm_meta.get("preserved_decisions",0),
            "output_bundles": {"native_bundle.json":"original imported channels", "analysis_bundle.json":"current source-bound derived analysis input"},
            "input_path_base": "supplied suite root",
            "inputs": [{"path": provenance_path(p, source_root), "sha256": pinned[str(p)]} for p in paths],
            "software": software_before, "software_after": software_after,
            "source_metadata": native_export.get("source_metadata", []),
            "options": {"evaluation": evaluation, "engine_config": str(engine_config) if engine_config else None, "reference_cache": str(Path(reference_cache).resolve()) if reference_cache else None, "manifest": provenance_path(Path(manifest).resolve(), source_root) if manifest else None,
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
            "scientific_scope": ("Fresh DNMB R-M HMM search of all translated CDSs, full configured reference search of native annotation/HMM nominees, context and fixed-policy evaluation. Original DF/PAD/DP calls preserved; biological activity and independent performance are not inferred." if full_meta else "Saved calls plus automatic stored/sequence evidence evaluation; optional explicitly configured fresh homology search. Missing measurements remain explicit; neither activity nor independent performance is inferred."),
            "outputs": {p.relative_to(stage).as_posix(): digest(p) for p in sorted(stage.rglob("*")) if p.is_file()},
        }
        write_json(stage / "receipt.json", receipt)
        if out.exists() or out.is_symlink():
            raise FileExistsError(f"Output appeared during execution: {out}")
        relative_figures = [p.relative_to(stage.resolve()) for p in resolved_figures]
        response = _serializable({"output": str(out), "samples": len(data["samples"]),
                                  "annotation_exports": {name: str(out / Path(path).relative_to(stage))
                                                         for name, path in annotation_exports.items()},
                                  "figures": [str(out / p) for p in relative_figures],
                                  "summary": result["global_summary"], "warnings": receipt["warnings"]})
        json.dumps(response, allow_nan=False)
        stage.rename(out)
        return response
    except BaseException as exc:
        full_work = stage / "evidence" / "reference" / "full_pipeline"
        if full_work.exists():
            write_json(stage / "failure.json", {"error": str(exc), "published": False,
                       "reason": "Full execution logs preserved; no successful publication claimed"})
            failed = stage.with_name(stage.name.replace(".staging-", ".failed-", 1))
            stage.rename(failed)
            if isinstance(exc, Exception):
                raise RuntimeError(str(exc) + "; full execution evidence preserved at " + str(failed)) from exc
        else:
            shutil.rmtree(stage, ignore_errors=True)
        raise

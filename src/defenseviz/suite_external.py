"""Analyze independent type-table examples without extending the main cohort."""
from __future__ import annotations

from pathlib import Path
import time

from .suite_project import ANALYSIS_OPTIONS, bind_project_inputs, resolve_suite_project


def prepare_external_type_cases(entries, primary_ids, workspace):
    from .suite_pipeline import prepare_suite, analyze_suite, _portable, _serializable
    from .plot.publication.display_scope import model_digest
    if not isinstance(entries, list):
        raise ValueError('External type projects must be a list')
    cases, execution, inputs, derived = [], [], {}, {}
    selected_ids = set(primary_ids)
    for index, entry in enumerate(entries):
        started = time.monotonic()
        if not isinstance(entry, dict) or set(entry) - {'project', 'sample_ids'} or 'project' not in entry:
            raise ValueError('External type project requires project and optional sample_ids')
        root, settings, project_refs = resolve_suite_project(entry['project'])
        if settings.get('external_type_projects'):
            raise ValueError('External type projects cannot recursively include more projects')
        kwargs = {k: v for k, v in settings.items() if k in ANALYSIS_OPTIONS}
        data = prepare_suite(root, **kwargs, evidence_workspace=Path(workspace) / f'case{index + 1:03d}')
        data.pop('_native_import_bundle', None)
        bind_project_inputs(data, project_refs)
        data, result = analyze_suite(data)
        available = {s['sample_id'] for s in result['samples']}
        selected = entry.get('sample_ids')
        selected = list(selected) if selected is not None else [s['sample_id'] for s in result['samples']]
        if not selected or len(selected) != len(set(selected)) or not set(selected) <= available:
            raise ValueError('External type sample selection must be distinct and present in its result')
        if set(selected) & selected_ids:
            raise ValueError('External type samples must not duplicate primary or other displayed samples')
        selected_ids.update(selected)
        refs = []
        for raw_path, expected in data.get('input_sha256', {}).items():
            path = str(Path(raw_path).resolve())
            if path in inputs and inputs[path] != expected:
                raise ValueError('External type projects have conflicting source hashes')
            inputs[path] = expected
            refs.append({'path': path, 'sha256': expected})
        if not refs:
            raise ValueError('External type result requires actual hashed inputs')
        # Resolve each project's relative source labels before the primary
        # project's path map is applied. Identical labels can name different files.
        aliases = {}
        for source in data.get('suite', {}).get('sources', []):
            label, resolved = source.get('path'), source.get('resolved_path')
            if not label or not resolved:
                continue
            path = str(Path(resolved).resolve())
            if inputs.get(path) != source.get('sha256'):
                raise ValueError('External source registry does not match its input hash')
            if label in aliases and aliases[label] != path:
                raise ValueError('Ambiguous external source label')
            aliases[label] = path
        result = _portable(_serializable(result), aliases)
        content_hash = model_digest(result)
        cases.append({'result': result, 'sample_ids': selected, 'source_refs': refs,
                      'result_sha256': content_hash})
        rm = data.get('rm_evaluation', {})
        full = data.get('suite', {}).get('full_execution', {})
        reference = data.get('suite', {}).get('automatic_reference', {})
        execution.append({'sample_ids': selected, 'analyzed_sample_ids': sorted(available),
                          'scope': 'independent type-table examples; excluded from primary aggregates',
                          'result_sha256': content_hash, 'source_refs': refs,
                          'evaluation_execution': data.get('evaluation_execution', {}),
                          'elapsed_seconds': time.monotonic() - started,
                          'policy_evaluations': rm.get('policy_evaluations', 0),
                          'policy_decisions': rm.get('policy_decisions', rm.get('policy_evaluations', rm.get('scoring_runs', 0))),
                          'preserved_native_decisions': rm.get('preserved_decisions', 0),
                          'sequence_feature_evaluations': rm.get('compiler', {}).get('motif_scans', 0),
                          'full_execution': full,
                          'external_processes': data.get('suite', {}).get('commands', []),
                          'scoring_runs': rm.get('scoring_runs', 0),
                          'hmm_search_runs': full.get('hmm_search_runs', 0),
                          'homology_search_runs': full.get('homology_search_runs', 0) + reference.get('homology_search_runs', 0),
                          'detector_runs': sum(full.get('detector_runs', {}).values()),
                          'original_defense_tool_runs': 0})
        for path, expected in data.get('derived_artifact_sha256', {}).items():
            if path in derived and derived[path] != expected:
                raise ValueError('Conflicting derived external evidence hash')
            derived[path] = expected
    return cases, execution, inputs, derived

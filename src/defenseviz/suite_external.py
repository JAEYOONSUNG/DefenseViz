"""Analyze independent type-table examples without extending the main cohort."""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import re
import time

from .suite_project import ANALYSIS_OPTIONS, bind_project_inputs, resolve_suite_project



_ASSEMBLY = re.compile(r'(?<![A-Z0-9])GC[AF]_\d{9,}(?:\.\d+)?(?!\d)', re.I)
_DIGEST = re.compile(r'[0-9a-f]{64}')


def cohort_identities(data):
    """Bound specimen identities; no display-label or individual-protein matching.

    A protein digest covers the complete imported multiset, preserving duplicate
    sequences. It does not claim that the importer covered an entire genome.
    Empty or incompletely sequenced imports cannot establish this identity axis.
    """
    samples = {s['sample_id']: s for s in data.get('samples', [])}
    identities = {sid: {'sample_id': sid, 'assembly_accessions': set(),
                        'sequence_artifact_sha256': set(), 'specimen_paths': set(),
                        'protein_set_sha256': None, 'protein_count': 0} for sid in samples}
    def assembly(sid, value):
        if sid in identities and isinstance(value, str):
            identities[sid]['assembly_accessions'].update(m.group().upper().split('.')[0]
                                                          for m in _ASSEMBLY.finditer(value))
    def metadata(sid, value):
        if isinstance(value, dict):
            for k, v in value.items():
                if k in {'sample_id', 'assembly_accession', 'assembly_id', 'assembly', 'accession', 'source_directory'}:
                    assembly(sid, v)
                elif isinstance(v, dict):
                    metadata(sid, v)
    for sid, sample in samples.items():
        metadata(sid, sample)
    for row in data.get('source_metadata', []):
        metadata(row.get('sample_id'), row)
    for row in data.get('suite', {}).get('samples', []):
        assembly(row.get('sample_id'), Path(row.get('directory', '')).name)
    genes = {sid: [] for sid in samples}
    for gene in data.get('genes', []):
        if gene.get('sample_id') in genes:
            genes[gene['sample_id']].append(gene.get('protein_sha256'))
    hashes = {str(Path(p).resolve()): h for p, h in data.get('input_sha256', {}).items()}
    def source(sid, path, role, module):
        if sid not in identities or not path:
            return
        path = str(Path(path).resolve())
        if path not in hashes:
            raise ValueError('Specimen identity source is missing its input hash')
        row = identities[sid]
        row['specimen_paths'].add(path)
        has_protein = any(isinstance(h, str) and _DIGEST.fullmatch(h) for h in genes[sid])
        sequence_source = ((role == 'protein_sequences' and has_protein) or module == 'genbank'
                           or Path(path).suffix.lower() in {'.gb', '.gbk', '.gbff'})
        if sequence_source and hashes[path] != hashlib.sha256(b'').hexdigest():
            row['sequence_artifact_sha256'].add(hashes[path])
            assembly(sid, Path(path).name)
    for path, roles in data.get('input_roles', {}).items():
        for role in roles:
            source(role.get('sample_id'), path, role.get('kind'), role.get('module'))
    for row in data.get('suite', {}).get('sources', []):
        source(row.get('sample_id'), row.get('resolved_path'), row.get('role'), row.get('module'))
    for sid, proteins in genes.items():
        if proteins and all(isinstance(h, str) and _DIGEST.fullmatch(h) for h in proteins):
            identities[sid]['protein_set_sha256'] = hashlib.sha256(
                json.dumps(sorted(proteins), separators=(',', ':')).encode('utf-8')).hexdigest()
            identities[sid]['protein_count'] = len(proteins)
        row = identities[sid]
        row['identity_status'] = ('assembly_or_sequence_bound' if row['assembly_accessions'] or row['sequence_artifact_sha256'] or row['protein_set_sha256'] else 'unassessed_no_assembly_or_sequence_binding')
    return list(identities.values())


def reject_cohort_overlap(incoming, existing):
    for row in incoming:
        for prior in existing:
            axes = [name for name in ('assembly_accessions', 'sequence_artifact_sha256', 'specimen_paths')
                    if row[name] & prior[name]]
            if row['sample_id'] == prior['sample_id']:
                axes.append('sample_id')
            if row['protein_set_sha256'] and row['protein_set_sha256'] == prior['protein_set_sha256']:
                axes.append('imported_protein_multiset')
            if axes:
                raise ValueError('External type samples must not duplicate primary or other external specimens: '
                                 + row['sample_id'] + ' / ' + prior['sample_id'] + ' (' + ', '.join(axes) + ')')


def _execution_overrides(overrides):
    overrides = dict(overrides or {})
    if set(overrides) - {'rscript', 'strict_inputs'}:
        raise ValueError('External execution overrides support only rscript and strict_inputs')
    if 'strict_inputs' in overrides and type(overrides['strict_inputs']) is not bool:
        raise ValueError('External strict_inputs override must be boolean')
    if 'rscript' in overrides and (not isinstance(overrides['rscript'], str) or not overrides['rscript'].strip()):
        raise ValueError('External rscript override must be nonempty text')
    return overrides


def prepare_external_type_cases(entries, primary_data, workspace, *, execution_overrides=None):
    from .suite_pipeline import prepare_suite, analyze_suite, _portable, _serializable
    from .plot.publication.display_scope import model_digest
    if not isinstance(entries, list):
        raise ValueError('External type projects must be a list')
    cases, execution, inputs, derived = [], [], {}, {}
    occupied = cohort_identities(primary_data)
    execution_overrides = _execution_overrides(execution_overrides)
    for index, entry in enumerate(entries):
        started = time.monotonic()
        if not isinstance(entry, dict) or set(entry) - {'project', 'sample_ids'} or 'project' not in entry:
            raise ValueError('External type project requires project and optional sample_ids')
        root, settings, project_refs = resolve_suite_project(entry['project'], execution_overrides)
        if settings.get('external_type_projects'):
            raise ValueError('External type projects cannot recursively include more projects')
        kwargs = {k: v for k, v in settings.items() if k in ANALYSIS_OPTIONS}
        data = prepare_suite(root, **kwargs, evidence_workspace=Path(workspace) / f'case{index + 1:03d}',
                             _input_guard=lambda native: reject_cohort_overlap(cohort_identities(native), occupied))
        identities = cohort_identities(data)
        reject_cohort_overlap(identities, occupied)
        data.pop('_native_import_bundle', None)
        bind_project_inputs(data, project_refs)
        data, result = analyze_suite(data)
        available = {s['sample_id'] for s in result['samples']}
        selected = entry.get('sample_ids')
        selected = list(selected) if selected is not None else [s['sample_id'] for s in result['samples']]
        if not selected or len(selected) != len(set(selected)) or not set(selected) <= available:
            raise ValueError('External type sample selection must be distinct and present in its result')
        occupied.extend(identities)
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
                          'pre_portable_result_sha256': content_hash, 'source_refs': refs,
                          'execution_settings': {'rscript': kwargs.get('rscript'), 'strict_inputs': kwargs.get('strict_inputs', False)},
                          'explicit_execution_overrides': execution_overrides,
                          'identity_check': [{k: sorted(v) if isinstance(v, set) else v
                                              for k, v in row.items()} for row in identities],
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

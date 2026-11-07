"""Portable settings for a repeatable DNMBsuite analysis command.

Projects contain input paths and display/evaluation settings, never executable
code or expected result counts. Explicit CLI options take precedence. Relative
paths are resolved against the project file, independently of the working dir.
"""
from __future__ import annotations

from pathlib import Path

from .io.project import Reader


PROJECT_NAME = 'defenseviz_project.json'
PATH_OPTIONS = frozenset({'manifest', 'supplement', 'rm_evidence', 'layout',
                         'engine_config', 'reference_cache', 'system_review_rules',
                         'system_review_evidence'})
ANALYSIS_OPTIONS = frozenset({'manifest', 'supplement', 'rm_evidence', 'engine_config',
                             'reference_cache', 'system_review_rules',
                             'system_review_evidence', 'evaluation', 'strict_inputs', 'rscript'})
SETTINGS = PATH_OPTIONS | {'evaluation', 'strict_inputs', 'rscript', 'title', 'external_type_projects'}


def resolve_suite_project(root, overrides=None):
    """Resolve one optional project and return settings plus its hashed source.

    A directory without the exact project filename keeps normal recursive
    suite discovery. A JSON file passed explicitly must satisfy this contract.
    This does not import observations, run an engine, or write anything.
    """
    location = Path(root).resolve()
    project = location if location.is_file() else location / PROJECT_NAME
    settings = {}
    refs = []
    source_root = location
    if project.is_file():
        reader = Reader(project)
        doc = reader.json(project)
        if not isinstance(doc, dict) or type(doc.get('schema_version')) is not int or doc['schema_version'] != 1:
            raise ValueError('Suite project requires schema_version: 1')
        unknown = set(doc) - SETTINGS - {'schema_version', 'root', 'description'}
        if unknown:
            raise ValueError(f'Unknown suite project fields: {sorted(unknown)}')
        source_root = reader.path(doc.get('root', '.'))
        if not source_root.is_dir():
            raise ValueError('Suite project root must be an existing input directory')
        for name in SETTINGS:
            if name not in doc:
                continue
            value = doc[name]
            if name in PATH_OPTIONS:
                if not isinstance(value, str) or not value.strip():
                    raise ValueError(f'Suite project {name} must be a nonempty path')
                settings[name] = reader.path(value)
            elif name == 'external_type_projects':
                if not isinstance(value, list):
                    raise ValueError('external_type_projects must be a list')
                cases = []
                for item in value:
                    entry = {'project': item} if isinstance(item, str) else item
                    if not isinstance(entry, dict) or set(entry) - {'project', 'sample_ids'}:
                        raise ValueError('External type project supports project and optional sample_ids')
                    if not isinstance(entry.get('project'), str) or not entry['project'].strip():
                        raise ValueError('External type project requires a path')
                    selected = entry.get('sample_ids')
                    if selected is not None and (not isinstance(selected, list) or not selected
                            or any(not isinstance(s, str) or not s.strip() for s in selected)
                            or len(set(selected)) != len(selected)):
                        raise ValueError('External sample_ids must contain distinct nonempty IDs')
                    cases.append({'project': reader.path(entry['project']), 'sample_ids': selected})
                settings[name] = cases
            elif name == 'evaluation':
                if value not in {'auto', 'stored', 'full'}:
                    raise ValueError('Invalid suite project evaluation mode')
                settings[name] = value
            elif name == 'strict_inputs':
                if type(value) is not bool:
                    raise ValueError('Suite project strict_inputs must be boolean')
                settings[name] = value
            else:
                if not isinstance(value, str) or not value.strip():
                    raise ValueError(f'Suite project {name} must be nonempty text')
                if name == 'rscript' and ('/' in value or '\\' in value):
                    value = str(reader.path(value))
                settings[name] = value
        refs = [{'path': str(project), 'sha256': reader.hashes[str(project)],
                 'role': 'suite_project_settings'}]
    elif not location.is_dir():
        raise ValueError('Suite input must be an existing directory or project JSON')
    for name, value in (overrides or {}).items():
        if name not in SETTINGS:
            raise ValueError(f'Unknown suite setting override: {name}')
        if value is not None:
            settings[name] = value
    settings.setdefault('evaluation', 'auto')
    settings.setdefault('strict_inputs', False)
    return source_root, settings, refs


def bind_project_inputs(data, refs):
    """Include settings in the same input-integrity audit as native artifacts."""
    from .pipeline import digest
    for ref in refs:
        path = Path(ref['path']).resolve()
        if digest(path) != ref['sha256']:
            raise ValueError('Suite project settings changed during import')
        previous = data.setdefault('input_sha256', {}).get(str(path))
        if previous is not None and previous != ref['sha256']:
            raise ValueError('Conflicting suite project input hash')
        data['input_sha256'][str(path)] = ref['sha256']
        if str(path) not in {str(p) for p in data.setdefault('inputs', [])}:
            data['inputs'].append(path)
    if refs:
        data.setdefault('suite', {})['project_sources'] = list(refs)

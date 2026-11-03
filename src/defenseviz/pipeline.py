"""Transactional offline integration, rendering, tables and reproducibility receipts."""
from __future__ import annotations
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from . import __version__


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def _write_table(path: Path, rows: list[dict], required: tuple[str, ...] = ()) -> None:
    fields = list(dict.fromkeys([*required, *(key for row in rows for key in row)]))
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False, allow_nan=False) if isinstance(value, (dict, list, tuple)) else '' if value is None else value for key, value in row.items()})


def software_state() -> dict:
    base = Path(__file__).resolve().parent
    state = {'defenseviz': __version__, 'python': platform.python_version(), 'platform': platform.platform(),
             'source_sha256': {str(p.relative_to(base)): digest(p) for p in sorted(base.rglob('*')) if p.is_file() and p.suffix in {'.py', '.json'}},
             'dependencies': {name: importlib.metadata.version(name) for name in ('matplotlib', 'numpy')}}
    # A wheel install has no repository; do not attribute it to an unrelated parent git tree.
    repo = base.parent.parent
    if (repo / '.git').exists():
        git = {}
        for key, args in [('commit', ['rev-parse', 'HEAD']), ('working_tree', ['status', '--porcelain=v1'])]:
            try:
                p = subprocess.run(['git', '-C', str(repo), *args], capture_output=True, text=True, timeout=10)
                git[key] = {'exit_code': p.returncode, 'value': p.stdout.strip()}
            except (OSError, subprocess.TimeoutExpired) as exc:
                git[key] = {'unavailable': str(exc)}
        state['git'] = git
    else:
        state['git'] = {'status': 'not_a_source_checkout', 'identity': 'package version and source hashes above'}
    return state


def provenance_path(path, base):
    """Relative when possible; Windows cross-drive sources stay explicit."""
    try:
        return os.path.relpath(path, base)
    except ValueError:
        return Path(path).absolute().as_posix()


def portable_provenance(value, base):
    """Use config-relative presentation pointers; receipts retain exact input hashes."""
    if isinstance(value, list):
        return [portable_provenance(v, base) for v in value]
    if not isinstance(value, dict):
        return value
    result = {}
    for key, child in value.items():
        if key == 'source_sha256' and isinstance(child, dict):
            result[key] = {provenance_path(path, base) if Path(path).is_absolute() else path: sha for path, sha in child.items()}
        elif key == 'source' and isinstance(child, str) and Path(child).is_absolute():
            result[key] = provenance_path(child, base)
        else:
            result[key] = portable_provenance(child, base)
    return result


def validate_options(options):
    if not isinstance(options, dict):
        raise ValueError('options must be an object')
    unknown = set(options) - {'policy', 'title', 'max_rows_per_page', 'max_families_per_page'}
    if unknown:
        raise ValueError('Unknown options: ' + ', '.join(sorted(unknown)))
    if 'title' in options and (not isinstance(options['title'], str) or not options['title'].strip()):
        raise ValueError('options.title must be nonempty text')
    for key in ('max_rows_per_page', 'max_families_per_page'):
        if key in options and (type(options[key]) is not int or options[key] <= 0):
            raise ValueError(f'options.{key} must be a positive integer')


def evaluate_project(data):
    from .analysis import integrate
    options = data.get('options', {})
    validate_options(options)
    records = data['evidence_records']
    audit = None
    extension = data.get('evidence_extension')
    if extension:
        from .analysis.admission import apply_reference_anchors
        records, audit = apply_reference_anchors(records, extension['candidates'], extension['registry'], extension['admission_rules'],
            registry_dir=extension['registry_dir'], source_root=extension['source_root'], policy=options.get('policy'))
    result = integrate(data['samples'], data['calls'], records, coverage=data.get('coverage'), genes=data.get('genes'), policy=options.get('policy'))
    if audit is not None:
        result['reference_anchor_admission'] = audit
    return result


def run(config_path: str | Path, output_dir: str | Path, *, argv: list[str] | None = None) -> dict:
    from .io import load_project
    from .analysis import integrate
    from .plot import render
    started = time.monotonic()
    config = Path(config_path).resolve()
    out = Path(output_dir).absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Output already exists; choose a fresh directory: {out}')
    data = load_project(config)
    options = data.get('options', {})
    validate_options(options)
    source_paths = sorted({config, *(Path(p).resolve() for p in data['inputs'])})
    source_hashes = {str(p): digest(p) for p in source_paths}
    for key, value in data.get('input_sha256', {}).items():
        if source_hashes[str(Path(key).resolve())] != value:
            raise ValueError(f'Input changed while loading: {key}')
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f'.{out.name}.staging-', dir=out.parent))
    try:
        t = time.monotonic()
        result = evaluate_project(data)
        for path, sha in result.get('reference_anchor_admission', {}).get('source_sha256', {}).items():
            resolved = str(Path(path).resolve())
            if resolved in source_hashes and source_hashes[resolved] != sha:
                raise ValueError(f'Input changed during evidence admission: {path}')
            source_hashes[resolved] = sha
        source_paths = sorted(Path(p) for p in source_hashes)
        integration_seconds = time.monotonic() - t
        result = portable_provenance(result, config.parent)
        write_json(stage / 'analysis.json', result)
        tables = stage / 'tables'; tables.mkdir()
        for field in ('sample_summary', 'cells', 'decisions', 'overlap', 'rm_components', 'native_system_counts', 'coverage', 'evidence_calls'):
            _write_table(tables / f'{field}.tsv', result.get(field, []), ('sample_id',))
        _write_table(tables / 'normalized_calls.tsv', data['calls'], ('sample_id', 'gene_id', 'caller', 'family'))
        t = time.monotonic()
        figures = render(result, stage, title=options.get('title', 'DefenseViz'), max_rows_per_page=options.get('max_rows_per_page', 40), max_families_per_page=options.get('max_families_per_page', 24))
        render_seconds = time.monotonic() - t
        if not figures or not (stage / 'figure2_ab.pdf').is_file():
            raise RuntimeError('Renderer did not produce the required Figure 2 PDF')
        for p in source_paths:
            if digest(p) != source_hashes[str(p)]:
                raise ValueError(f'Input changed during execution: {p}')
        policy = result.get('policy', {})
        receipt = {'schema_version': 1, 'tool': 'DefenseViz', 'version': __version__, 'created_at': datetime.now(timezone.utc).isoformat(),
                   'operation': 'saved_call_evaluation_and_render', 'detector_runs': 0, 'structure_prediction_runs': 0,
                   'argv': argv if argv is not None else sys.argv, 'exit_code': 0, 'elapsed_seconds': time.monotonic()-started,
                   'stages': [{'name': 'saved_call_evaluation', 'elapsed_seconds': integration_seconds}, {'name': 'figure_render', 'elapsed_seconds': render_seconds}],
                   'input_path_base': 'directory containing supplied config', 'config_name': config.name,
                   'inputs': [{'path': provenance_path(p, config.parent), 'sha256': source_hashes[str(p)]} for p in source_paths],
                   'software': software_state(), 'policy_sha256': hashlib.sha256(json.dumps(policy, sort_keys=True, allow_nan=False).encode()).hexdigest(),
                   'outputs': {str(p.relative_to(stage)): digest(p) for p in sorted(stage.rglob('*')) if p.is_file()},
                   'scientific_scope': 'Annotation evidence only. Missing/unreviewed is unassessed. No biological activity or broad classifier-performance inference.',
                   'warnings': result.get('warnings', [])}
        write_json(stage / 'receipt.json', receipt)
        # Existing destinations are rejected; unrelated concurrent writers are unsupported.
        if out.exists():
            raise FileExistsError(f'Output appeared during execution: {out}')
        stage.rename(out)
        return {'output': str(out), 'figures': [str(out / Path(p).relative_to(stage)) for p in figures], 'samples': len(data['samples']), 'warnings': result.get('warnings', [])}
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise

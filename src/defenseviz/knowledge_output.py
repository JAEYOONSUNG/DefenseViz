"""Offline catalogue exports; reference context does not change analysis."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import shutil
import tempfile
import time

from . import __version__
from .pipeline import _write_table, digest, software_state, write_json


_HEADERS = {
    'units': ('unit_id', 'sample_id', 'reference_titles', 'reference_mechanism_summaries',
              'reference_source_mechanism_statement_states', 'reference_wiki_urls'),
    'detections': ('detection_id', 'sample_id', 'gene_id', 'identity_status',
                   'version_relation', 'reason_codes'),
    'wiki_systems': ('id', 'title', 'mechanism_summary', 'source_mechanism_statement_state',
                     'sensor', 'activator', 'effector', 'pfam', 'citation_dois', 'wiki_url'),
    'model_definitions': ('definition_id', 'model_name', 'source_family', 'defense_scope',
                          'package_version', 'source_commit', 'source_url'),
    'model_components': ('definition_id', 'native_component_id', 'relationship',
                         'declared_presence', 'parent_component_id'),
    'links': ('link_id', 'wiki_id', 'definition_id'),
    'sources': ('source_id', 'repository', 'commit'),
    'package_configurations': ('source_id',),
}


def _table(path, name, rows):
    if name not in _HEADERS:
        raise ValueError('Invalid system information table name')
    _write_table(path, rows, _HEADERS[name])


def write_system_information(directory, information):
    from .knowledge import information_tables

    directory = Path(directory)
    write_json(directory / 'system_information.json', information)
    tables = directory / 'tables'
    tables.mkdir(exist_ok=True)
    rows = information_tables(information)
    for name, records in rows.items():
        _table(tables / f'system_information_{name}.tsv', name, records)
    _table(tables / 'system_information.tsv', 'units', rows['units'])


def export_catalog(output_dir, *, query=None, argv=None):
    """Publish a fresh searchable-reference dataset from bundled snapshots."""
    from .knowledge import catalog_tables, search_catalog

    started = time.monotonic()
    out = Path(output_dir).absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Output already exists; choose a fresh directory: {out}')
    before = software_state()
    catalog = search_catalog(query)
    tables = catalog_tables(catalog)
    out.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f'.{out.name}.staging-', dir=out.parent))
    try:
        write_json(stage / 'catalog.json', catalog)
        for name, rows in tables.items():
            _table(stage / f'{name}.tsv', name, rows)
        (stage / 'README.md').write_text(
            '# Defense system reference information\n\n'
            'These offline reference records describe the pinned official DefenseFinder wiki, '
            'DefenseFinder model definitions and CasFinder definitions. Wiki pages and model '
            'definitions have different scopes; their row counts are not numbers of unique '
            'biological systems or detected CDSs.\n\n'
            'Use `wiki_systems.tsv` for system-level information and citations, '
            '`model_definitions.tsv` for source model identities, and `model_components.tsv` '
            'for native component metadata. `catalog.json` preserves the structured detail. '
            'Source-unknown, missing and ambiguous information remain explicit. '
            'For sensor, activator, effector and pfam, read each value together '
            'with its `_state` column: `unknown` is stated by the source, while '
            '`not_provided` means the source omits that field. '
            '`source_mechanism_statement_state` describes the source statement, '
            'not confidence in a query protein or a biological experiment. '
            'Search results are candidates for browsing, not automatic identity assignments.\n\n'
            'Display-family labels in figures can pool native systems; they are '
            'not catalogue identities and must not be used as definition or wiki IDs. '
            'Join detected-system information through its original detection/unit IDs.\n\n'
            'Reference metadata does not establish activity, query-protein function, or '
            'completeness of a detected system. Model versions and source hashes must be '
            'distinguished from the versions that produced a saved detection.\n\n'
            '[DefenseViz documentation](https://github.com/JAEYOONSUNG/DefenseViz)\n',
            encoding='utf-8')
        after = software_state()
        if before['source_sha256'] != after['source_sha256']:
            raise ValueError('Catalogue implementation changed during export')
        receipt = dict(schema_version=1, tool='DefenseViz', version=__version__,
                       operation='offline_reference_catalogue_export',
                       created_at=datetime.now(timezone.utc).isoformat(), query=query,
                       argv=list(argv or []), exit_code=0,
                       elapsed_seconds=time.monotonic()-started,
                       detector_runs=0, scoring_runs=0, software=before,
                       row_counts={name: len(rows) for name, rows in tables.items()},
                       outputs={p.relative_to(stage).as_posix(): digest(p)
                                for p in sorted(stage.rglob('*')) if p.is_file()})
        write_json(stage / 'receipt.json', receipt)
        if out.exists() or out.is_symlink():
            raise FileExistsError(f'Output appeared during export: {out}')
        stage.rename(out)
        return dict(output=str(out), query=query, row_counts=receipt['row_counts'])
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise

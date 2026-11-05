"""Offline reference exports preserve provenance and publish atomically."""
import csv
import hashlib
import json

import pytest

from defenseviz.knowledge_output import export_catalog


def test_catalogue_export_is_source_bound_and_readable(tmp_path):
    out = tmp_path / 'catalogue'
    result = export_catalog(out, query='Wadjet', argv=['catalog', 'Wadjet'])
    assert result['row_counts']['wiki_systems'] == 1
    with (out / 'wiki_systems.tsv').open(newline='') as stream:
        records = list(csv.DictReader(stream, delimiter='\t'))
    assert records[0]['title'] == 'Wadjet'
    assert records[0]['mechanism_summary'] and records[0]['source_sha256']
    receipt = json.loads((out / 'receipt.json').read_text())
    assert receipt['detector_runs'] == receipt['scoring_runs'] == 0
    assert receipt['argv'] == ['catalog', 'Wadjet']
    for name, digest in receipt['outputs'].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest
    with pytest.raises(FileExistsError):
        export_catalog(out)


def test_empty_search_has_readable_empty_tables(tmp_path):
    out = tmp_path / 'none'
    result = export_catalog(out, query='no_such_synthetic_system_9e33')
    assert result['row_counts']['wiki_systems'] == 0
    for name in ('wiki_systems', 'model_definitions', 'model_components', 'links'):
        with (out / (name + '.tsv')).open(newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            assert reader.fieldnames and list(reader) == []


def test_catalogue_mutation_during_export_leaves_no_partial_result(tmp_path, monkeypatch):
    import defenseviz.knowledge_output as output
    states = iter([{'source_sha256': {'file': 'before'}},
                   {'source_sha256': {'file': 'after'}}])
    monkeypatch.setattr(output, 'software_state', lambda: next(states))
    with pytest.raises(ValueError, match='changed during export'):
        export_catalog(tmp_path / 'result')
    assert not list(tmp_path.iterdir())

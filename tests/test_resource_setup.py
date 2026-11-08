"""Preparation preserves paired inputs and never claims readiness on failure."""
import hashlib
import json
from pathlib import Path
import sys

import pytest

from defenseviz import resource_setup as setup


def fixture(tmp_path, monkeypatch):
    source = tmp_path / 'source'; source.mkdir()
    cache = source / 'cache'; cache.mkdir()
    for name in setup._REFERENCE_REQUIRED:
        (cache / name).write_text(name)
    (cache / 'rebase_reference_state.rds').write_text('provider state')
    (cache / 'snapshot.json').write_text('{"version":"synthetic"}')
    hmm = source / 'source.hmm'
    hmm.write_text('HMMER3/f\nNAME  test\nLENG  60\nGA  20 20;\n//\n')
    roles = source / 'roles.tsv'
    roles.write_text('profile_name\toff_target\troles\ttype_hints\tprofile_length\tlibrary_sha256\n'
                     'test\tfalse\tR\tII\t60\t' + setup.native.digest(hmm) + '\n')
    libs = source / 'r-library'; libs.mkdir()
    cfg = source / 'engine.json'
    cfg.write_text(json.dumps(dict(schema_version=1, reference_cache='cache', hmm_database='source.hmm',
        profile_metadata='roles.tsv', r_libraries=['r-library'], rscript=sys.executable,
        blastp=sys.executable, blastdbcmd=sys.executable, hmmsearch=sys.executable,
        cath_resolve_hits=sys.executable)))
    calls = []
    def metadata(config, work, env):
        calls.append(('metadata', config))
        return dict(passed=True, metadata_matches_FASTA=True)
    def validate(config, out):
        calls.append(('resources', json.loads(config.read_text())))
        out.mkdir()
        (out / 'resource_check.json').write_text('{"ready":true}')
        return dict(ready=True)
    monkeypatch.setattr(setup, '_validate_metadata', metadata)
    monkeypatch.setattr(setup, '_hmmstat_path', lambda cfg: sys.executable)
    monkeypatch.setattr(setup, '_validate_hmm', lambda *args: dict(passed=True, complete_models_parsed=True))
    monkeypatch.setattr(setup.resources, 'verify_resources', validate)
    return source, cfg, cache, hmm, roles, calls


def test_copy_is_independent_preserves_provenance_and_resolves_template_paths(tmp_path, monkeypatch):
    source, cfg, cache, hmm, roles, calls = fixture(tmp_path, monkeypatch)
    before = setup._hashes(p for p in source.rglob('*') if p.is_file())
    result = setup.prepare_resources(cfg, tmp_path / 'out', argv=['example'], source_label='synthetic snapshot v1')
    assert result['ready'] and result['inputs_unchanged'] and result['exit_code'] == 0
    assert result['resource_downloads'] == result['detector_runs'] == 0
    assert result['provider_snapshot'] == 'synthetic snapshot v1'
    assert result['argv'] == ['example']
    assert setup._hashes(p for p in source.rglob('*') if p.is_file()) == before
    out = tmp_path / 'out'
    for source_file in cache.iterdir():
        copied = out / 'rebase' / source_file.name
        assert copied.read_bytes() == source_file.read_bytes()
        assert not copied.samefile(source_file) and not copied.is_symlink()
    new = json.loads((out / 'engine.json').read_text())
    assert new['r_libraries'] == [str(source / 'r-library')]
    assert new['reference_cache'] == str(out / 'rebase')
    assert [call[0] for call in calls] == ['metadata', 'resources']
    assert all(setup.native.digest(out / name) == value for name, value in result['outputs'].items())


def test_unknown_provider_version_is_explicit(tmp_path, monkeypatch):
    _, cfg, *_ = fixture(tmp_path, monkeypatch)
    result = setup.prepare_resources(cfg, tmp_path / 'out')
    assert result['provider_snapshot'].startswith('unknown provider date/version')


@pytest.mark.parametrize('fault', ['hash', 'missing_index', 'existing_output', 'inside_cache', 'inside_library'])
def test_unsafe_or_mismatched_inputs_fail_before_writing(tmp_path, monkeypatch, fault):
    source, cfg, cache, hmm, roles, _ = fixture(tmp_path, monkeypatch)
    out = tmp_path / 'out'
    if fault == 'hash':
        roles.write_text(roles.read_text().replace(setup.native.digest(hmm), '0' * 64))
    elif fault == 'missing_index':
        (cache / 'rebase_db.fasta.psq').unlink()
    elif fault == 'existing_output':
        out.mkdir(); (out / 'keep').write_text('keep')
    elif fault == 'inside_cache':
        out = cache / 'new'
    elif fault == 'inside_library':
        out = source / 'r-library/new'
    with pytest.raises((ValueError, FileExistsError, FileNotFoundError)):
        setup.prepare_resources(cfg, out)
    assert not out.exists() or fault == 'existing_output'
    if fault == 'existing_output':
        assert (out / 'keep').read_text() == 'keep'


def test_dangling_output_symlink_is_not_followed(tmp_path, monkeypatch):
    _, cfg, *_ = fixture(tmp_path, monkeypatch)
    out = tmp_path / 'out'; out.symlink_to(tmp_path / 'missing')
    with pytest.raises(FileExistsError):
        setup.prepare_resources(cfg, out)
    assert out.is_symlink() and not (tmp_path / 'missing').exists()


@pytest.mark.parametrize('kwargs', [
    {'reference_cache': 'one', 'rebase_source': 'two'},
    {'hmm_metadata': 'one'},
    {'hmm': 'one', 'hmm_metadata': 'two'},
    {'hmm': 'one', 'hmm_metadata': 'two', 'off_target': 'three', 'profile_metadata': 'four'},
])
def test_conflicting_modes_fail_without_output(tmp_path, monkeypatch, kwargs):
    _, cfg, *_ = fixture(tmp_path, monkeypatch)
    with pytest.raises(ValueError):
        setup.prepare_resources(cfg, tmp_path / 'out', **kwargs)
    assert not (tmp_path / 'out').exists()


def test_metadata_mismatch_retains_failure_not_readiness(tmp_path, monkeypatch):
    _, cfg, *_ = fixture(tmp_path, monkeypatch)
    def fail(*args):
        raise ValueError('Reference FASTA proteins do not match reference rows')
    monkeypatch.setattr(setup, '_validate_metadata', fail)
    out = tmp_path / 'out'
    with pytest.raises(RuntimeError, match='proteins do not match'):
        setup.prepare_resources(cfg, out)
    result = json.loads((out / 'preparation.json').read_text())
    assert result['ready'] is False and result['exit_code'] == 1
    assert not (out / 'validation').exists()


def test_source_mutation_cannot_be_reported_unchanged(tmp_path, monkeypatch):
    _, cfg, cache, *_ = fixture(tmp_path, monkeypatch)
    def changed(*args):
        (cache / 'rebase_data.rds').write_text('mutated source by test')
        return dict(passed=True, metadata_matches_FASTA=True)
    monkeypatch.setattr(setup, '_validate_metadata', changed)
    with pytest.raises(RuntimeError, match='Source changed during preparation'):
        setup.prepare_resources(cfg, tmp_path / 'out')
    result = json.loads((tmp_path / 'out/preparation.json').read_text())
    assert result['ready'] is False and 'inputs_unchanged' not in result


def test_provider_hmm_build_keeps_cutoffs_and_audit(tmp_path, monkeypatch):
    source, cfg, _, hmm, _, _ = fixture(tmp_path, monkeypatch)
    metadata = source / 'restriction_metadata.csv'
    metadata.write_text('Name,System,Enzyme_type\ntest,RM_Type_II,RE\n')
    excluded = source / 'off_target.txt'; excluded.write_text('')
    result = setup.prepare_resources(cfg, tmp_path / 'out', hmm=hmm, hmm_metadata=metadata, off_target=excluded)
    assert result['hmm_mode'] == 'explicit_provider_metadata'
    assert (tmp_path / 'out/rm-hmm/RM_HMMs.hmm').read_bytes() == hmm.read_bytes()
    assert json.loads((tmp_path / 'out/rm-hmm/preparation.json').read_text())['original_GA_cutoffs_preserved']


def test_hmmer_parse_failure_never_reaches_reference_check_or_readiness(tmp_path, monkeypatch):
    _, cfg, _, _, _, calls = fixture(tmp_path, monkeypatch)
    def fail(*args):
        raise RuntimeError('hmmstat_execution failed')
    monkeypatch.setattr(setup, '_validate_hmm', fail)
    with pytest.raises(RuntimeError, match='hmmstat_execution failed'):
        setup.prepare_resources(cfg, tmp_path / 'out')
    result = json.loads((tmp_path / 'out/preparation.json').read_text())
    assert result['ready'] is False and result['exit_code'] == 1
    assert not calls


@pytest.mark.parametrize('inventory', ['', '1 wrong - 2 2.0 60 1 1\n', '1 test - 2 2.0 59 1 1\n',
                                        '2 test - 2 2.0 60 1 1\n', 'unexpected text\n'])
def test_successful_hmmer_exit_requires_exact_model_inventory(tmp_path, monkeypatch, inventory):
    hmm = tmp_path / 'model.hmm'
    hmm.write_text('HMMER3/f\nNAME test\nLENG 60\nGA 20 20;\n//\n')
    monkeypatch.setattr(setup, '_run', lambda *a, **kw: dict(stdout=inventory, executable_sha256='hash'))
    with pytest.raises(ValueError):
        setup._validate_hmm({'hmm_database':str(hmm)}, tmp_path, {}, sys.executable)

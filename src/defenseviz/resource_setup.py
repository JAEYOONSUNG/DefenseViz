"""Prepare an isolated, verified DNMB resource set from explicit local inputs.

No provider download, detector search, or decision-policy change occurs here.
Raw REBASE inputs require the already installed pinned R helper namespace.
"""
from __future__ import annotations

from datetime import datetime, timezone
from importlib.metadata import version, PackageNotFoundError
from importlib.resources import files
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from . import hmm_preparation, native_pipeline as native, resources

_ALLOWED = {'schema_version', 'rscript', 'r_libraries', 'reference_cache', 'hmm_database',
            'profile_metadata', 'hmmsearch', 'cath_resolve_hits', 'blastp', 'blastdbcmd',
            'cpu', 'runtime_dependencies'}
_DB_KEYS = {'reference_cache', 'hmm_database', 'profile_metadata'}
_REFERENCE_REQUIRED = ('rebase_data.rds', 'rebase_db.fasta', 'rebase_db.fasta.phr',
                       'rebase_db.fasta.pin', 'rebase_db.fasta.psq')
_REFERENCE_OPTIONAL = ('REBASE_protein_seqs.txt', 'protein_seqs.txt', 'rebase_reference_state.rds')


def _path(value, parent):
    p = Path(value).expanduser()
    return (p if p.is_absolute() else parent / p).resolve()


def _executable(value, parent):
    if not isinstance(value, str) or not value:
        raise ValueError('Executable must be a nonempty string')
    located = shutil.which(value) if not any(x in value for x in ('/', '\\')) else None
    p = Path(located).resolve() if located else _path(value, parent)
    if not p.is_file() or not os.access(p, os.X_OK):
        raise FileNotFoundError('Required resource preparation executable: ' + str(p))
    return str(p)


def _template(path):
    path = Path(path).expanduser().resolve()
    cfg = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(cfg, dict) or cfg.get('schema_version') != 1 or set(cfg) - _ALLOWED:
        raise ValueError('Engine configuration needs schema_version 1 and supported resource fields')
    for key in ('reference_cache', 'hmm_database', 'profile_metadata', 'cath_resolve_hits'):
        if key in cfg:
            if not isinstance(cfg[key], str) or not cfg[key]:
                raise ValueError('Invalid ' + key)
            cfg[key] = str(_path(cfg[key], path.parent))
    for key, default in (('rscript', 'Rscript'), ('blastp', 'blastp'), ('hmmsearch', 'hmmsearch')):
        cfg[key] = _executable(cfg.get(key, default), path.parent)
    if 'cath_resolve_hits' not in cfg:
        raise ValueError('Engine configuration requires cath_resolve_hits')
    cfg['cath_resolve_hits'] = _executable(cfg['cath_resolve_hits'], path.parent)
    for key in ('r_libraries', 'runtime_dependencies'):
        values = cfg.get(key, [])
        if not isinstance(values, list) or any(not isinstance(v, str) or not v for v in values):
            raise ValueError(key + ' must be a list of nonempty paths')
        cfg[key] = [str(_path(v, path.parent)) for v in values]
        if any(not Path(v).exists() for v in cfg[key]):
            raise FileNotFoundError('Missing ' + key + ' resource')
    cfg['blastdbcmd'] = native._resolve_blastdbcmd(dict(cfg, _config_path=str(path)))
    cpu = cfg.get('cpu', 2)
    if type(cpu) is not int or not 1 <= cpu <= 256:
        raise ValueError('cpu must be an integer from 1 to 256')
    cfg['cpu'] = cpu
    return path, cfg


def _reference_files(cache):
    cache = Path(cache).expanduser().resolve()
    for name in _REFERENCE_REQUIRED:
        if not (cache / name).is_file():
            raise FileNotFoundError('Incomplete REBASE cache; missing ' + str(cache / name))
    names = set(_REFERENCE_REQUIRED)
    names.update(p.name for p in cache.glob('rebase_db.fasta.*') if p.is_file())
    names.update(name for name in _REFERENCE_OPTIONAL if (cache / name).is_file())
    # Provider/preparation receipts remain attached to the snapshot, when supplied.
    names.update(p.name for p in cache.glob('*.json') if p.is_file())
    return [cache / name for name in sorted(names)]


def _hashes(paths):
    return {str(Path(p).resolve()): native.digest(p) for p in paths}


def _unchanged(inputs):
    for path, expected in inputs.items():
        if not Path(path).is_file() or native.digest(path) != expected:
            raise ValueError('Source changed during preparation: ' + path)


def _copy(source, target, expected):
    shutil.copyfile(source, target)
    if native.digest(target) != expected:
        raise ValueError('Copied resource differs from source: ' + str(source))


def _run(argv, work, name, env, timeout=None):
    start = time.monotonic()
    receipt = dict(argv=[str(a) for a in argv], executable_sha256=native.digest(argv[0]))
    try:
        p = subprocess.run(argv, cwd=work, env=env, capture_output=True, text=True, timeout=timeout)
        receipt.update(exit_code=p.returncode, stdout=p.stdout, stderr=p.stderr)
    except (OSError, subprocess.TimeoutExpired) as error:
        receipt.update(exit_code=None, error=str(error))
        raise
    finally:
        receipt['elapsed_seconds'] = time.monotonic() - start
        native.write_json(work / (name + '.json'), receipt)
    if p.returncode:
        raise RuntimeError(name + ' failed; inspect ' + str(work / (name + '.json')))
    return receipt


def _hmmstat_path(cfg):
    return _executable(str(Path(cfg['hmmsearch']).with_name('hmmstat')), Path.cwd())


def _validate_hmm(cfg, work, env, executable):
    """Have HMMER parse every complete model without searching query sequences."""
    version = _run([executable, '-h'], work, 'hmmstat_version', env, timeout=30)
    parsed = _run([executable, cfg['hmm_database']], work, 'hmmstat_execution', env, timeout=120)
    observed = []
    for line in parsed['stdout'].splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        if len(fields) < 6 or not fields[0].isdigit() or not fields[5].isdigit():
            raise ValueError('Unexpected hmmstat model inventory format')
        if int(fields[0]) != len(observed) + 1:
            raise ValueError('HMMER model inventory indices are not consecutive')
        observed.append((fields[1], int(fields[5])))
    expected = [(row['profile_name'], int(row['profile_length']))
                for row in hmm_preparation.hmm_profiles(cfg['hmm_database'])]
    if not observed or observed != expected:
        raise ValueError('HMMER model inventory differs from HMM header names/lengths')
    result = dict(passed=True, profiles=len(observed), complete_models_parsed=True,
                  detector_searches=0, database_sha256=native.digest(cfg['hmm_database']),
                  executable_sha256=parsed['executable_sha256'],
                  version_record='hmmstat_version.json', execution_record='hmmstat_execution.json')
    native.write_json(work / 'hmm_model_validation.json', result)
    return result


def _validate_metadata(cfg, work, env):
    """Bind the RDS to exact encoded FASTA identifiers and protein sequences."""
    script = work / 'validate_reference_metadata.R'
    policy = Path(str(files('defenseviz').joinpath('data/native/native_metadata.R')))
    script.write_text('''a <- commandArgs(TRUE)
.libPaths(c(a[-seq_len(3L)], .libPaths()))
source(a[1], local=TRUE)
ref <- readRDS(file.path(a[2], 'rebase_data.rds'))
required <- c('enzyme_name','sequence','enz_type','rm_type','subunit','seq_length','rec_seq','is_gold_standard')
stopifnot(is.data.frame(ref), nrow(ref)>0L, all(required %in% names(ref)))
idx <- .dnmb_reference_index(ref, file.path(a[2], 'rebase_db.fasta'))
stopifnot(length(idx)>0L)
jsonlite::write_json(list(passed=TRUE, reference_rows=nrow(ref), encoded_rows=length(idx),
 metadata_matches_FASTA=TRUE, reference_metadata_sha256=digest::digest(file=file.path(a[2],'rebase_data.rds'),algo='sha256'),
 fasta_sha256=digest::digest(file=file.path(a[2],'rebase_db.fasta'),algo='sha256'),
 validator_sha256=digest::digest(file=a[1],algo='sha256')), a[3],auto_unbox=TRUE,pretty=TRUE)
''', encoding='utf-8')
    output = work / 'reference_metadata_validation.json'
    _run([cfg['rscript'], '--vanilla', str(script), str(policy), cfg['reference_cache'],
          str(output), *cfg['r_libraries']], work, 'reference_metadata_execution', env)
    result = json.loads(output.read_text())
    if result.get('passed') is not True or result.get('metadata_matches_FASTA') is not True:
        raise ValueError('Reference metadata validation did not confirm exact pairing')
    return result


def _raw_reference(raw, cfg, target, work, env, makeblastdb):
    target.mkdir()
    _copy(raw, target / 'REBASE_protein_seqs.txt', native.digest(raw))
    script = work / 'prepare_rebase.R'
    script.write_text('''a <- commandArgs(TRUE)
.libPaths(c(a[-seq_len(3L)], .libPaths()))
suppressPackageStartupMessages(library(dplyr))
stopifnot(requireNamespace('DefenseViz',quietly=TRUE))
ref <- DefenseViz::parse_rebase_sequences(a[1])
required <- c('enzyme_name','sequence','enz_type','rm_type','subunit','seq_length','rec_seq','is_gold_standard')
stopifnot(is.data.frame(ref), nrow(ref)>0L, all(required %in% names(ref)))
saveRDS(ref, file.path(a[2],'rebase_data.rds'))
DefenseViz::write_fasta_for_blast(ref,file.path(a[2],'rebase_db.fasta'),id_col='enzyme_name',seq_col='sequence')
ns <- asNamespace('DefenseViz')
functions <- c('parse_rebase_sequences','write_fasta_for_blast')
hashes <- setNames(lapply(functions,function(n)digest::digest(list(formals(get(n,ns)),body(get(n,ns))),algo='sha256')),functions)
jsonlite::write_json(list(provider_version='unknown unless supplied in source_label',reference_rows=nrow(ref),
 helper_namespace='DefenseViz',helper_version=as.character(utils::packageVersion('DefenseViz')),
 helper_path=getNamespaceInfo(ns,'path'),helper_function_sha256=hashes),a[3],auto_unbox=TRUE,pretty=TRUE)
''', encoding='utf-8')
    _run([cfg['rscript'], '--vanilla', str(script), str(raw), str(target),
          str(work / 'rebase_helper.json'), *cfg['r_libraries']], work, 'rebase_parse_execution', env)
    _run([makeblastdb, '-version'], work, 'makeblastdb_version', env)
    fasta = str(target / 'rebase_db.fasta')
    _run([makeblastdb, '-in', fasta, '-dbtype', 'prot', '-out', fasta, '-title', 'REBASE'],
         work, 'makeblastdb_execution', env)


def _source_state():
    result = {'module_sha256': {p.name: native.digest(p) for p in
              (Path(__file__), Path(hmm_preparation.__file__), Path(resources.__file__), Path(native.__file__))}}
    try:
        result['package_version'] = version('defenseviz')
    except PackageNotFoundError:
        result['package_version'] = 'source checkout; not installed'
    root = Path(__file__).resolve().parents[2]
    if (root / '.git').exists():
        commands = []
        for label, args in [('commit', ['rev-parse', 'HEAD']), ('working_tree', ['status', '--porcelain=v1'])]:
            start = time.monotonic()
            run = subprocess.run(['git', '-C', str(root), *args], capture_output=True, text=True)
            commands.append(dict(argv=['git', '-C', str(root), *args], exit_code=run.returncode,
                                 elapsed_seconds=time.monotonic()-start))
            result[label] = run.stdout.strip() if not run.returncode else 'unavailable'
        result['commands'] = commands
    else:
        result['commit'] = 'not available in installed distribution; use module hashes'
        result['working_tree'] = 'not a source checkout'
    return result


def prepare_resources(engine_config, out, *, reference_cache=None, rebase_source=None,
                      hmm=None, profile_metadata=None, hmm_metadata=None, off_target=None,
                      makeblastdb=None, source_label=None, argv=None):
    """Copy or build paired local databases, then validate the complete runtime.

    The output must be new. Failed preparations retain a failure receipt and
    logs, but are never reported ready. Paths in the template are resolved before
    writing the new configuration. Original inputs are never modified.
    """
    start = time.monotonic()
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError('Output must be a new directory: ' + str(out))
    out = out.resolve()
    config_path, cfg = _template(engine_config)
    hmmstat = _hmmstat_path(cfg)
    if reference_cache is not None and rebase_source is not None:
        raise ValueError('Choose reference_cache or rebase_source, not both')
    raw_hmm = hmm_metadata is not None or off_target is not None
    if raw_hmm and (hmm is None or hmm_metadata is None or off_target is None):
        raise ValueError('Raw HMM preparation requires hmm, hmm_metadata and off_target together')
    if raw_hmm and profile_metadata is not None:
        raise ValueError('Raw HMM metadata cannot be combined with prepared profile_metadata')
    hmm_path = Path(hmm or cfg.get('hmm_database', '')).expanduser().resolve()
    raw = Path(rebase_source).expanduser().resolve() if rebase_source is not None else None
    cache = None if raw else Path(reference_cache or cfg.get('reference_cache', '')).expanduser().resolve()
    if raw is None and not reference_cache and 'reference_cache' not in cfg:
        raise ValueError('Supply reference_cache or rebase_source')
    if hmm is None and 'hmm_database' not in cfg:
        raise ValueError('Supply an HMM database')
    sources = [config_path, hmm_path, Path(hmmstat)]
    if raw:
        sources.append(raw)
        blast_builder = _executable(makeblastdb or str(Path(cfg['blastp']).with_name('makeblastdb')), config_path.parent)
        sources.append(Path(blast_builder))
        reference_sources = []
    else:
        reference_sources = _reference_files(cache)
        sources.extend(reference_sources)
        blast_builder = None
    if raw_hmm:
        source_metadata = Path(hmm_metadata).expanduser().resolve()
        exclusions = Path(off_target).expanduser().resolve()
        sources.extend([source_metadata, exclusions])
    else:
        if profile_metadata is None and 'profile_metadata' not in cfg:
            raise ValueError('Supply prepared profile_metadata or raw HMM metadata')
        roles = Path(profile_metadata or cfg['profile_metadata']).expanduser().resolve()
        sources.append(roles)
        resources._profile_binding(hmm_path, roles)
    hmm_receipts = []
    if not raw_hmm and hmm_path.parent == roles.parent:
        for name in ('preparation.json', 'source_preparation.json'):
            candidate = hmm_path.parent / name
            if candidate.is_file():
                value = json.loads(candidate.read_text(encoding='utf-8'))
                outputs = value.get('outputs', {}) if isinstance(value, dict) else {}
                if (outputs.get(hmm_path.name) == native.digest(hmm_path)
                        and outputs.get(roles.name) == native.digest(roles)):
                    hmm_receipts.append(candidate)
        sources.extend(hmm_receipts)
    for protected in [cache, *map(Path, cfg['r_libraries'])]:
        if protected is not None and out.is_relative_to(protected):
            raise ValueError('Output cannot be inside a source database or R library: ' + str(protected))
    inputs = _hashes(sources)
    before = _source_state()
    out.mkdir(parents=True)
    work = out / 'preparation_logs'; work.mkdir()
    receipt = dict(schema_version=1, operation='local_resource_preparation', ready=False,
                   created_at=datetime.now(timezone.utc).isoformat(), argv=list(argv or sys.argv),
                   provider_snapshot=source_label or 'unknown provider date/version; exact input hashes identify this snapshot',
                   inputs=inputs, source=before, resource_downloads=0, detector_runs=0,
                   full_analysis_executed=False, original_GA_cutoffs_preserved=True,
                   reference_mode='raw_provider_file' if raw else 'existing_paired_cache',
                   hmm_mode='explicit_provider_metadata' if raw_hmm else 'existing_paired_library')
    try:
        env = native._runtime_environment(cfg, work / 'cache')
        rebase_out = out / 'rebase'; hmm_out = out / 'rm-hmm'
        if raw:
            _raw_reference(raw, cfg, rebase_out, work, env, blast_builder)
        else:
            rebase_out.mkdir()
            for p in reference_sources:
                _copy(p, rebase_out / p.name, inputs[str(p.resolve())])
        if raw_hmm:
            hmm_preparation.prepare(hmm_path, source_metadata, exclusions, hmm_out, argv=receipt['argv'])
        else:
            hmm_out.mkdir()
            for source, name in [(hmm_path, 'RM_HMMs.hmm'), (roles, 'profile_roles.tsv')]:
                _copy(source, hmm_out / name, inputs[str(source.resolve())])
            for source in hmm_receipts:
                _copy(source, hmm_out / source.name, inputs[str(source.resolve())])
        cfg.update(reference_cache=str(rebase_out), hmm_database=str(hmm_out / 'RM_HMMs.hmm'),
                   profile_metadata=str(hmm_out / 'profile_roles.tsv'))
        output_config = out / 'engine.json'
        native.write_json(output_config, cfg)
        checked_cfg = native.read_engine_config(output_config)
        receipt['hmm_models'] = _validate_hmm(checked_cfg, work, env, hmmstat)
        receipt['reference_metadata'] = _validate_metadata(checked_cfg, work, env)
        receipt['validation'] = resources.verify_resources(output_config, out / 'validation')
        _unchanged(inputs)
        receipt.update(ready=True, inputs_unchanged=True, exit_code=0,
                       engine_config=str(output_config), elapsed_seconds=time.monotonic()-start,
                       outputs={p.relative_to(out).as_posix(): native.digest(p)
                                for directory in (rebase_out, hmm_out) for p in directory.rglob('*') if p.is_file()})
        receipt['outputs']['engine.json'] = native.digest(output_config)
        native.write_json(out / 'preparation.json', receipt)
        return receipt
    except Exception as error:
        receipt.update(error=str(error), exit_code=1, elapsed_seconds=time.monotonic()-start)
        native.write_json(out / 'preparation.json', receipt)
        raise RuntimeError('Resource preparation failed: ' + str(error) + '; logs: ' + str(out)) from error


def add_arguments(parser):
    parser.add_argument('--engine-config', type=Path, required=True, help='Runtime configuration or template; relative paths resolve here')
    parser.add_argument('--out', type=Path, required=True, help='New directory for resources, engine.json and verification receipts')
    parser.add_argument('--reference-cache', type=Path, help='Existing paired REBASE metadata, FASTA and BLAST cache')
    parser.add_argument('--rebase-source', type=Path, help='Authorized provider protein_seqs.txt to parse and index locally')
    parser.add_argument('--hmm', type=Path, help='Exact source HMM library; defaults to configured hmm_database')
    parser.add_argument('--profile-metadata', type=Path, help='Prepared metadata bound to the HMM library')
    parser.add_argument('--hmm-metadata', type=Path, help='Provider restriction_metadata.csv for raw HMM preparation')
    parser.add_argument('--off-target', type=Path, help='Provider off_target.txt; required with --hmm-metadata')
    parser.add_argument('--makeblastdb', help='BLAST index builder; defaults to the sibling of configured blastp')
    parser.add_argument('--source-label', help='Provider release/date/source description; omitted means unknown')
    return parser


def from_args(args, argv=None):
    return prepare_resources(args.engine_config, args.out,
        **{key: getattr(args, key, None) for key in ('reference_cache', 'rebase_source', 'hmm',
          'profile_metadata', 'hmm_metadata', 'off_target', 'makeblastdb', 'source_label')}, argv=argv)

"""Read-only readiness checks for configured DNMB analysis resources."""
from __future__ import annotations

import csv
import json
from importlib.resources import files
from pathlib import Path
import re
import subprocess
import tempfile
import time

from . import native_pipeline as native


def _profile_binding(database, metadata):
    """Validate the exact named GA models and their declared role mapping."""
    database=Path(database);metadata=Path(metadata)
    database_hash=native.digest(database)
    records=[];lengths={};name=None;ga=False;length=None
    for line in database.read_text().splitlines():
        if line.startswith('NAME '):
            if name is not None:raise ValueError('HMM record lacks // terminator before '+line)
            name=line.split()[1];ga=False;length=None
        elif line.startswith('LENG '):
            fields=line.split()
            if name is None or length is not None or len(fields)!=2 or not fields[1].isdigit() or int(fields[1])<=0:
                raise ValueError('HMM profile needs one positive integer LENG: '+str(name))
            length=int(fields[1])
        elif line.startswith('GA '):
            try:
                values=[float(x.rstrip(';')) for x in line.split()[1:]]
                ga=len(values)==2 and all(float('-inf')<v<float('inf') for v in values)
            except ValueError:ga=False
        elif line.strip()=='//':
            if name is None:raise ValueError('HMM record lacks a NAME')
            if not ga:raise ValueError('HMM profile lacks valid sequence/domain GA cutoffs: '+name)
            if length is None:raise ValueError('HMM profile lacks LENG: '+name)
            lengths[name]=length
            records.append(name);name=None
    if name is not None:raise ValueError('Final HMM record lacks // terminator')
    if not records or len(records)!=len(set(records)):
        raise ValueError('HMM database must contain nonempty unique profile names')
    with metadata.open(newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        required={'profile_name','off_target','roles','type_hints','profile_length','library_sha256'}
        missing=required-set(reader.fieldnames or [])
        if missing:raise ValueError('Profile metadata lacks columns: '+', '.join(sorted(missing)))
        if len(reader.fieldnames)!=len(set(reader.fieldnames)):
            raise ValueError('Profile metadata has duplicate columns')
        rows=list(reader)
        if any(None in row or None in row.values() for row in rows):
            raise ValueError('Profile metadata row/header widths differ')
    names=[row['profile_name'] for row in rows]
    if len(names)!=len(set(names)) or set(names)!=set(records):
        raise ValueError('HMM profile names and role metadata do not match exactly')
    for row in rows:
        if row['off_target'].lower() not in {'true','false','1','0'}:
            raise ValueError('Invalid off_target value for '+row['profile_name'])
        if row['roles'] and not set(row['roles'].split(';')) <= {'R','M','S','unknown'}:
            raise ValueError('Invalid role hint for '+row['profile_name'])
        recorded=row['library_sha256']
        if recorded!=database_hash:
            raise ValueError('Profile metadata library_sha256 does not match configured HMM database')
        declared=row['profile_length']
        if not declared.isdigit() or int(declared)!=lengths[row['profile_name']]:
            raise ValueError('Profile metadata profile_length does not match HMM LENG: '+row['profile_name'])
    excluded=sum(row['off_target'].lower() in {'true','1'} for row in rows)
    return dict(profiles=len(records),off_target_profiles=excluded,reportable_profiles=len(records)-excluded,
                database_sha256=database_hash,metadata_sha256=native.digest(metadata),
                all_profiles_have_GA=True,exact_profile_name_binding=True)


def _origin(path,cfg,r_home):
    path=Path(path).resolve()
    if path.is_relative_to((Path(r_home)/'library').resolve()):return 'R_installation'
    if any(path.is_relative_to(Path(p).resolve()) for p in cfg['r_libraries']):return 'configured_library'
    return 'other_installed_library'


def _verify(engine_config,work,retained):
    started=time.monotonic()
    cfg=native.read_engine_config(engine_config)
    cfg['blastdbcmd']=native._resolve_blastdbcmd(cfg)
    env=native._runtime_environment(cfg,work/'cache')
    cache=Path(cfg['reference_cache']);prefix=cache/'rebase_db.fasta'
    required=[cache/'rebase_data.rds',prefix,*[Path(str(prefix)+ext) for ext in ('.phr','.pin','.psq')]]
    for path in required:
        if not path.is_file():raise FileNotFoundError('Incomplete REBASE cache; missing '+str(path))
    inputs={str(Path(cfg[k]).resolve()):native.digest(cfg[k]) for k in
            ('_config_path','hmm_database','profile_metadata','rscript','blastp','blastdbcmd','hmmsearch','cath_resolve_hits')}
    for p in [*required,*sorted(cache.glob('rebase_db.fasta.*')),*map(Path,cfg['runtime_dependencies'])]:
        inputs[str(p.resolve())]=native.digest(p)
    native_path=Path(str(files('defenseviz').joinpath('data/native')))
    inputs.update(native._verify_bundled_sources(native_path))
    inputs[str(Path(__file__).resolve())]=native.digest(__file__)
    profiles=_profile_binding(cfg['hmm_database'],cfg['profile_metadata'])
    versions={}
    for name,flags in [('rscript',['--version']),('blastp',['-version']),('blastdbcmd',['-version']),
                       ('hmmsearch',['-h']),('cath_resolve_hits',['--version'])]:
        argv=[cfg[name],*flags];start=time.monotonic()
        p=subprocess.run(argv,capture_output=True,text=True,env=env,timeout=30)
        versions[name]=dict(argv=argv,exit_code=p.returncode,elapsed_seconds=time.monotonic()-start,
            stdout=p.stdout,stderr=p.stderr,executable_sha256=native.digest(cfg[name]))
        native.write_json(work/'tool_versions.json',versions)
        if p.returncode:raise RuntimeError(name+' version check failed; inspect tool_versions.json')
    runtime,runtime_hashes=native._capture_r_runtime(cfg,work,env)
    inputs.update(runtime_hashes)
    reference=native._validate_reference_index(prefix,cfg['blastdbcmd'],work/'reference_index_validation.json',env)
    for path,expected in inputs.items():
        if native.digest(path)!=expected:raise ValueError('Read-only resource changed during validation: '+path)
    packages=[dict(name=p['name'],version=p['version'],path=p['path'],
                   origin=_origin(p['path'],cfg,runtime['R_home'])) for p in runtime['packages']]
    result=dict(schema_version=1,ready=True,execution='resource_validation_only',
        full_analysis_executed=False,resource_downloads=0,inputs_unchanged=True,
        engine_config=dict(path=cfg['_config_path'],sha256=cfg['_config_sha256']),
        reference=dict(protein_entries=reference['entries'],FASTA_sha256=native.digest(prefix),
                       metadata_sha256=native.digest(cache/'rebase_data.rds'),index_matches_FASTA=True),
        hmm=profiles,R=dict(version=runtime['R_version'],home=runtime['R_home'],packages=packages),
        tools=versions,elapsed_seconds=time.monotonic()-started,report_retained=retained)
    native.write_json(work/'resource_input_hashes.json',inputs)
    if retained:
        result['report_directory']=str(work)
        result['artifact_sha256']={p.relative_to(work).as_posix():native.digest(p)
                                   for p in work.rglob('*') if p.is_file()}
    native.write_json(work/'resource_check.json',result)
    return result


def verify_resources(engine_config,out=None):
    """Check resources without running HMM/BLAST analysis or changing inputs.

    ``out`` must be a new directory. Failures preserve logs there and raise an
    actionable RuntimeError. Without ``out``, only the returned readiness record
    remains; temporary logs are removed on return.
    """
    temporary=None
    if out is None:
        temporary=tempfile.TemporaryDirectory(prefix='defenseviz-resource-check-')
        work=Path(temporary.name)
    else:
        work=Path(out).expanduser().resolve();work.mkdir(parents=True,exist_ok=False)
    try:
        return _verify(engine_config,work,out is not None)
    except Exception as error:
        native.write_json(work/'resource_check.json',dict(schema_version=1,ready=False,
            full_analysis_executed=False,error=str(error),engine_config=str(engine_config)))
        detail='; logs: '+str(work) if out is not None else ''
        raise RuntimeError('DNMB resources are not ready: '+str(error)+detail) from error
    finally:
        if temporary is not None:temporary.cleanup()

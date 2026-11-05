"""Execute the complete DNMB R-M evidence stages in a private output workspace.

Native DF/PAD/DP calls remain input observations. R-M HMM nomination searches
all translated CDSs; native annotation/context nominees are then searched against
the complete configured REBASE resource. No reference update or activity claim.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from importlib.resources import files
from types import SimpleNamespace

from .pipeline import digest as _digest, write_json


def digest(path):
    return _digest(Path(path))


def engine_config_path(explicit=None):
    if explicit is not None:
        return Path(explicit).expanduser().resolve()
    if os.environ.get('DEFENSEVIZ_ENGINE_CONFIG'):
        return Path(os.environ['DEFENSEVIZ_ENGINE_CONFIG']).expanduser().resolve()
    path = Path(os.environ.get('XDG_CONFIG_HOME', str(Path.home()/'.config')))/'defenseviz'/'engine.json'
    return path if path.is_file() else None


def read_engine_config(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text(encoding='utf-8'))
    allowed = {'schema_version','rscript','r_libraries','reference_cache','hmm_database',
               'profile_metadata','hmmsearch','cath_resolve_hits','blastp','blastdbcmd','cpu','runtime_dependencies'}
    if not isinstance(config,dict) or config.get('schema_version') != 1 or set(config)-allowed:
        raise ValueError('Engine configuration needs schema_version 1 and supported resource fields')
    for key in ('reference_cache','hmm_database','profile_metadata','cath_resolve_hits'):
        if not isinstance(config.get(key),str) or not config[key]:
            raise ValueError('Engine configuration requires '+key)
    out = dict(config)
    for key in ('reference_cache','hmm_database','profile_metadata','cath_resolve_hits'):
        p = Path(config[key]).expanduser()
        if not p.is_absolute(): p = path.parent/p
        out[key] = str(p.resolve())
        if not p.exists(): raise FileNotFoundError(str(p))
    for key, default in [('rscript','Rscript'),('blastp','blastp'),('hmmsearch','hmmsearch')]:
        command = config.get(key,default)
        if not isinstance(command,str) or not command: raise ValueError('Invalid '+key)
        candidate = Path(command).expanduser()
        path_like = candidate.is_absolute() or '/' in command or '\\' in command
        located = None if path_like else shutil.which(command)
        if located is None:
            p = candidate if candidate.is_absolute() else path.parent/candidate
            located = str(p) if p.is_file() and os.access(p,os.X_OK) else None
        if located is None: raise FileNotFoundError('Required full evaluation executable: '+command)
        out[key] = str(Path(located).resolve())
    for key in ('r_libraries','runtime_dependencies'):
        entries = config.get(key,[])
        if not isinstance(entries,list) or any(not isinstance(x,str) for x in entries):
            raise ValueError(key+' must be a list of paths')
        out[key] = [str((Path(x).expanduser() if Path(x).expanduser().is_absolute() else path.parent/x).resolve()) for x in entries]
        if any(not Path(x).exists() for x in out[key]): raise FileNotFoundError('Missing '+key+' resource')
    cpu = config.get('cpu',2)
    if type(cpu) is not int or not 1 <= cpu <= 256: raise ValueError('cpu must be an integer from 1 to 256')
    out['cpu'] = cpu
    out['_config_path'] = str(path)
    out['_config_sha256'] = digest(path)
    return out


def _runtime_environment(cfg, cache):
    """Resolve native prechecks and bootstrap configured R libraries identically."""
    env = dict(os.environ)
    bins = [str(Path(cfg[k]).parent) for k in ('blastp','hmmsearch','cath_resolve_hits','rscript') if cfg.get(k)]
    env['PATH'] = os.pathsep.join(dict.fromkeys([*bins, *env.get('PATH','').split(os.pathsep)]))
    libraries = [*cfg.get('r_libraries', []), *env.get('R_LIBS','').split(os.pathsep)]
    env['R_LIBS'] = os.pathsep.join(dict.fromkeys(x for x in libraries if x))
    env.update(DNMB_BLASTP=cfg['blastp'], DNMB_CACHE_ROOT=str(cache),
               DEFENSEVIZ_CACHE_DIR=str(Path(cache)/'db_modules/rebasefinder/cache'),
               XDG_CACHE_HOME=str(cache), PYTHONDONTWRITEBYTECODE='1')
    return env


def _verify_bundled_sources(native):
    """Check shipped numerical/context sources against their frozen manifest."""
    native = Path(native)
    provenance = json.loads((native/'provenance.json').read_text())
    expected = provenance.get('source_files_sha256', {})
    actual = {p.relative_to(native/'source').as_posix() for p in (native/'source').rglob('*.R')}
    if not expected or actual != set(expected):
        raise ValueError('Bundled native source inventory differs from frozen provenance')
    for name, value in expected.items():
        if not re.fullmatch('[0-9a-f]{64}', str(value)) or digest(native/'source'/name) != value:
            raise ValueError('Bundled native source differs from frozen provenance: '+name)
    return {str(p.resolve()):digest(p) for p in native.rglob('*') if p.is_file()}


def _resolve_blastdbcmd(cfg):
    explicit = cfg.get('blastdbcmd')
    if explicit is not None:
        if not isinstance(explicit, str) or not explicit:
            raise ValueError('Invalid blastdbcmd executable')
        candidate = Path(explicit).expanduser()
        path_like = candidate.is_absolute() or '/' in explicit or '\\' in explicit
        located = None if path_like else shutil.which(explicit)
        if not candidate.is_absolute(): candidate = Path(cfg['_config_path']).parent/candidate
        if located is None and candidate.is_file() and os.access(candidate, os.X_OK): located = str(candidate)
    else:
        sibling = Path(cfg['blastp']).with_name('blastdbcmd')
        located = str(sibling) if sibling.is_file() and os.access(sibling, os.X_OK) else shutil.which('blastdbcmd')
    if located is None:
        raise FileNotFoundError('Full reference validation needs blastdbcmd beside blastp or in engine configuration')
    return str(Path(located).resolve())


def _reference_entries(path):
    name = None; parts = []
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(parts)
                name = line[1:].strip(); parts = []
                if not name: raise ValueError('Reference FASTA has an empty identifier')
            elif line.strip():
                if name is None: raise ValueError('Reference sequence precedes a FASTA identifier')
                parts.append(''.join(line.split()).upper())
    if name is not None: yield name, ''.join(parts)


def _validate_reference_index(prefix, executable, output, env):
    """Stream every indexed sequence and bind it to the corresponding FASTA entry."""
    prefix=Path(prefix); output=Path(output)
    command=[str(executable),'-db',str(prefix),'-entry','all','-outfmt','%t\t%s']
    receipt=dict(argv=command,executable_sha256=digest(executable),fasta_sha256=digest(prefix),
                 normalization='uppercase and whitespace removal; one terminal <> marker removed only')
    started=time.monotonic(); count=0; repaired=0; raw_length=0; residue_length=0
    observed=hashlib.sha256(); entries=iter(_reference_entries(prefix)); seen=set(); process=None
    try:
        with output.with_suffix('.stderr.log').open('w') as stderr:
            process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=stderr,text=True,env=env)
            for line in process.stdout:
                observed.update(line.encode())
                fields=line.rstrip('\r\n').split('\t')
                if len(fields)!=2: raise ValueError('Unexpected blastdbcmd sequence inventory format')
                try: name, raw=next(entries)
                except StopIteration: raise ValueError('BLAST index has entries absent from its FASTA') from None
                if name in seen: raise ValueError('Reference FASTA identifier is duplicated: '+name)
                seen.add(name)
                clean=raw[:-2] if raw.endswith('<>') else raw
                repaired+=int(clean!=raw)
                if not re.fullmatch('[A-Z*]+',clean):
                    raise ValueError('Unexplained non-protein reference sequence characters: '+name)
                if fields != [name,clean]:
                    raise ValueError('BLAST index and reference FASTA disagree at entry '+str(count+1)+': '+name)
                count+=1; raw_length+=len(raw); residue_length+=len(clean)
            code=process.wait(); receipt['exit_code']=code
            if code: raise RuntimeError('blastdbcmd failed; see '+str(output.with_suffix('.stderr.log')))
            try: next(entries)
            except StopIteration: pass
            else: raise ValueError('Reference FASTA has entries absent from its BLAST index')
            if not count: raise ValueError('Complete reference database is empty')
        receipt.update(passed=True,entries=count,terminal_markers_removed=repaired,
                       raw_character_count=raw_length,indexed_residue_count=residue_length,
                       indexed_stream_sha256=observed.hexdigest())
    except Exception as error:
        if process is not None and process.poll() is None: process.kill(); process.wait()
        receipt.update(passed=False,exit_code=process.returncode if process is not None else None,error=str(error))
        raise
    finally:
        if process is not None and process.stdout is not None: process.stdout.close()
        entries.close()
        receipt['elapsed_seconds']=time.monotonic()-started
        write_json(output,receipt)
    return receipt


def _capture_r_runtime(cfg, work, env):
    """Record loaded namespace locations, versions and actual executable assets."""
    work=Path(work); probe=work/'r_runtime_probe.R'; output=work/'r_runtime.json'
    probe.write_text('''a<-commandArgs(TRUE)
.libPaths(c(a[-1],.libPaths()))
for(p in c("DNMB","DefenseViz","dplyr","jsonlite","digest","Biostrings","Peptides","gtools","openxlsx","plyr","reshape2","stringr","tidyr","withr"))loadNamespace(p)
invisible(La.svd(matrix(1)))
packages<-lapply(sort(loadedNamespaces()),function(p){d<-normalizePath(if(p=="base")file.path(R.home("library"),"base")else getNamespaceInfo(asNamespace(p),"path"));list(name=p,path=d,version=as.character(utils::packageVersion(p,lib.loc=dirname(d))))})
dll<-lapply(getLoadedDLLs(),function(d)d[["path"]])
jsonlite::write_json(list(R_version=R.version.string,R_home=R.home(),library_paths=.libPaths(),packages=packages,loaded_DLLs=dll,blastp=unname(Sys.which("blastp"))),a[1],auto_unbox=TRUE,pretty=TRUE)
''')
    command=[cfg['rscript'],'--vanilla',str(probe),str(output),*cfg['r_libraries']]
    started=time.monotonic()
    process=subprocess.run(command,capture_output=True,text=True,env=env,cwd=work)
    receipt=dict(argv=command,exit_code=process.returncode,elapsed_seconds=time.monotonic()-started,
                 stdout=process.stdout,stderr=process.stderr,executable_sha256=digest(cfg['rscript']))
    write_json(work/'r_runtime_execution.json',receipt)
    if process.returncode or not output.is_file():
        raise RuntimeError('Configured R runtime did not load required packages; see '+str(work/'r_runtime_execution.json'))
    runtime=json.loads(output.read_text())
    if Path(runtime.get('blastp','')).resolve()!=Path(cfg['blastp']).resolve():
        raise ValueError('R preflight resolved a different blastp executable')
    package_names={row['name'] for row in runtime['packages']}
    if not {'DNMB','DefenseViz','jsonlite','digest'}<=package_names:
        raise ValueError('R preflight omitted required loaded namespaces')
    paths={probe,output,work/'r_runtime_execution.json'}
    for package in runtime['packages']:
        directory=Path(package['path'])
        for name in ('DESCRIPTION','NAMESPACE','Meta/package.rds'):
            if (directory/name).is_file(): paths.add(directory/name)
        for name in ('R','libs','data','extdata'):
            if (directory/name).is_dir(): paths.update(p for p in (directory/name).rglob('*') if p.is_file())
    paths.update(Path(p) for p in runtime.get('loaded_DLLs',{}).values() if p and Path(p).is_file())
    r_binary=Path(runtime['R_home'])/'bin/exec/R'
    if r_binary.is_file(): paths.add(r_binary)
    r_libraries=Path(runtime['R_home'])/'lib'
    if r_libraries.is_dir(): paths.update(p for p in r_libraries.rglob('*') if p.is_file())
    return runtime,{str(p.resolve()):digest(p) for p in paths}


def _verify_native_runtime(path, expected, resources):
    """Bind namespaces loaded by the actual engine to the measured preflight."""
    path = Path(path)
    actual = json.loads(path.read_text())
    if (actual.get('R_version') != expected.get('R_version')
            or Path(actual.get('R_home', '')).resolve() != Path(expected.get('R_home', '')).resolve()):
        raise ValueError('Native engine loaded a different R runtime')
    def packages(runtime):
        rows = runtime.get('packages')
        if not isinstance(rows, list) or not rows:
            raise ValueError('Native R runtime lacks loaded namespace records')
        result = {}
        for row in rows:
            if not isinstance(row, dict) or any(not row.get(k) for k in ('name', 'path', 'version')):
                raise ValueError('Native R namespace record is incomplete')
            if row['name'] in result:
                raise ValueError('Native R namespace record is duplicated')
            result[row['name']] = (str(Path(row['path']).resolve()), row['version'])
        return result
    before, after = packages(expected), packages(actual)
    if not {'DNMB', 'DefenseViz', 'jsonlite', 'digest'} <= after.keys():
        raise ValueError('Native engine omitted required loaded namespaces')
    for name, record in after.items():
        if before.get(name) != record:
            raise ValueError('Native engine loaded an unpinned or different R namespace: ' + name)
    actual_dlls = actual.get('loaded_DLLs')
    expected_dlls = expected.get('loaded_DLLs')
    if not isinstance(actual_dlls, dict) or not isinstance(expected_dlls, dict):
        raise ValueError('Native R runtime lacks loaded DLL records')
    for name, location in actual_dlls.items():
        previous = expected_dlls.get(name)
        if name == 'base' and location == previous == 'base':
            continue
        resolved = str(Path(location).resolve()) if isinstance(location, str) and location else ''
        if (not previous or resolved != str(Path(previous).resolve())
                or resolved not in resources or digest(resolved) != resources[resolved]):
            raise ValueError('Native engine loaded an unpinned or different R DLL: ' + name)
    return dict(passed=True, R_version=actual['R_version'],
                loaded_namespaces=len(after), loaded_DLLs=len(actual_dlls),
                runtime_path=str(path.resolve()), runtime_sha256=digest(path))


def _table(path,rows,columns):
    with Path(path).open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,columns,delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def _sequence(value):
    return ''.join(str(value).split()).upper()


def _nomination_union(data, sid, fresh_path, output):
    """Preserve qualified tool nominations with exact canonical query identity."""
    from .automatic_rm import _profiles, _qualifies
    with Path(fresh_path).open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        columns = list(reader.fieldnames or [])
        fresh = list(reader)
    profiles = {}
    for row in data.get('suite', {}).get('profiles', []):
        if row.get('sample_id') == sid:
            profiles.setdefault(row['gene_id'], []).append(row)
    translations = {}
    for row in data.get('suite', {}).get('annotations', []):
        fields = row.get('native_fields', {})
        if (row.get('sample_id') == sid and fields.get('locus_tag') == row.get('gene_id')
                and fields.get('translation') not in {None, '', 'NA'}):
            translations.setdefault(row['gene_id'], set()).add(_sequence(fields['translation']))
    nominated = []
    for gene in data['genes']:
        if gene['sample_id'] != sid: continue
        local_profiles = profiles.get(gene['gene_id'], [])
        if not local_profiles: continue
        profile_genes = [gene]
        nominee_hash = gene.get('protein_sha256')
        nominee_length = gene.get('protein_length')
        normalization = 'none'
        donor = translations.get(gene['gene_id'], set())
        if donor:
            if len(donor) != 1:
                raise ValueError('Conflicting exact protein binding for comparator nomination')
            sequence = next(iter(donor))
            if (hashlib.sha256(sequence.encode()).hexdigest() != gene.get('protein_sha256')
                    or len(sequence) != gene.get('protein_length')):
                raise ValueError('Comparator nomination source protein binding mismatch')
            clean = sequence.rstrip('*')
            if not clean:
                raise ValueError('Comparator nominee lacks amino acids after stop removal')
            nominee_hash = hashlib.sha256(clean.encode()).hexdigest()
            nominee_length = len(clean)
            if clean != sequence:
                normalization = 'terminal_stop_removal_only'
                # Profile coordinates retain their reported protein length. A
                # raw query with a terminal stop and its stop-free counterpart
                # are both exact representations of this canonical sequence.
                # Any other length remains rejected by _profiles.
                profile_genes.append(dict(gene, protein_sha256=nominee_hash,
                                          protein_length=nominee_length))
        seen = set()
        for profile_gene in profile_genes:
            for profile in _profiles(local_profiles, profile_gene):
                if not _qualifies(profile): continue
                identity = repr((profile['model'], profile['family'], profile['role'],
                                 profile['evalue'], profile['query_coverage'],
                                 profile['model_coverage'], profile['source_refs']))
                if identity in seen: continue
                seen.add(identity)
                nominated.append(dict(locus_tag=gene['gene_id'], profile_name=profile['model'],
                    role_hint=profile['role'] if profile['role'] in {'R','M','S'} else 'unknown',
                    protein_sha256=nominee_hash, sequence_length=nominee_length,
                    source='saved_DNMBsuite_HMM', evalue=profile['evalue'],
                    query_coverage=profile['query_coverage'], model_coverage=profile['model_coverage'],
                    source_protein_sha256=gene.get('protein_sha256'),
                    source_sequence_length=gene.get('protein_length'),
                    sequence_normalization=normalization))
    rows = fresh + nominated
    columns = list(dict.fromkeys([*columns, *(k for r in nominated for k in r)]))
    _table(output, rows, columns)
    return dict(fresh_HMM_nominees=len({r['locus_tag'] for r in fresh}),
                saved_qualified_nominees=len({r['locus_tag'] for r in nominated}),
                union_nominees=len({r['locus_tag'] for r in rows}),
                scoring_credit='native DNMB profiles only; saved comparator profiles nominate candidates')


def _sample_path(data,sid,root):
    # load_suite records the selected directory independently of labels.
    candidates = [s.get('directory') or s.get('path') for s in data.get('suite',{}).get('samples',[]) if s.get('sample_id')==sid]
    candidates = [Path(p) for p in candidates if p]
    if len(candidates)==1:
        return candidates[0] if candidates[0].is_absolute() else Path(root)/candidates[0]
    # Source registry contains resolved identities, so this fallback never
    # guesses genome identity from a folder label or a fuzzy accession.
    parents=set()
    for s in data.get('suite',{}).get('sources',[]):
        if s.get('sample_id')!=sid or not s.get('resolved_path'): continue
        p=Path(s['resolved_path'])
        for parent in [p.parent,*p.parents]:
            if parent==Path(root).parent: break
            if list(parent.glob('*.gbff')) or list(parent.glob('*.gbk')) or list(parent.glob('*.gb')):
                parents.add(parent);break
    if len(parents)!=1: raise ValueError('Full evaluation requires one exact GenBank source directory for '+sid)
    return parents.pop()


def _write_saved_rm_nominees(data, sid, genes, sequences, out):
    """Source-bound prior identities nominate searches, never roles or points."""
    registry = {}
    for source in data.get('suite', {}).get('sources', []):
        if source.get('sample_id') == sid and source.get('resolved_path'):
            registry[source.get('path')] = source
    hashes = data.get('input_sha256', {})
    prior = {}; excluded = []
    for record in data.get('saved_decisions', []):
        if record.get('sample_id') != sid: continue
        gid = record['gene_id']; gene = genes.get(gid)
        if gene is None: raise ValueError('Saved R-M nominee is absent from canonical genes: '+gid)
        if gid not in sequences or not gene.get('protein_sha256'):
            excluded.append(dict(gene_id=gid, reason='no_translated_protein'))
            continue
        if record.get('protein_sha256') != gene['protein_sha256']:
            raise ValueError('Saved R-M nominee protein binding mismatch: '+gid)
        refs = []
        for ref in [record.get('source_ref'), *record.get('source_refs', [])]:
            if not isinstance(ref, dict): continue
            source = registry.get(ref.get('path'), {})
            name = source.get('resolved_path', ref.get('path'))
            if not isinstance(name, str) or not Path(name).is_absolute():
                raise ValueError('Saved R-M nominee source path is not exactly bound: '+gid)
            expected = ref.get('sha256')
            if not expected or hashes.get(name) != expected or digest(name) != expected:
                raise ValueError('Saved R-M nominee source hash mismatch: '+gid)
            bound = dict(ref, path=name)
            if bound not in refs: refs.append(bound)
        if not refs: raise ValueError('Saved R-M nominee has no source evidence: '+gid)
        clean = sequences[gid].rstrip('*')
        row = prior.setdefault(gid, dict(sample_id=sid, locus_tag=gid,
            source_protein_sha256=gene['protein_sha256'],
            normalized_protein_sha256=hashlib.sha256(clean.encode()).hexdigest(),
            aa_length=len(clean), source_refs=[]))
        for ref in refs:
            if ref not in row['source_refs']: row['source_refs'].append(ref)
    rows = [dict(row, source_refs=json.dumps(row['source_refs'], sort_keys=True)) for _, row in sorted(prior.items())]
    _table(out/'saved_rm_nominees.tsv', rows, ['sample_id','locus_tag','source_protein_sha256',
        'normalized_protein_sha256','aa_length','source_refs'])
    write_json(out/'saved_rm_nomination.json', dict(sample_id=sid, source_bound_prior_nominees=len(rows),
        excluded=excluded, role_type_score_inherited=False, nomination_scope='saved_RM_candidate_identities_only'))
    return rows


def _prepare_sample(data,sid,directory,out):
    genes = {g['gene_id']:g for g in data['genes'] if g['sample_id']==sid}
    rows = [r['native_fields'] for r in data.get('suite',{}).get('annotations',[])
            if r.get('sample_id')==sid and r.get('native_fields',{}).get('locus_tag')==r.get('gene_id')
            and r.get('native_fields',{}).get('translation') not in {None,'','NA'}]
    if not rows: raise ValueError('Full evaluation needs exact protein-bound native annotations for '+sid)
    sequences={}
    for row in rows:
        name=row['locus_tag']; seq=_sequence(row['translation'])
        if hashlib.sha256(seq.encode()).hexdigest()!=genes[name].get('protein_sha256'):
            raise ValueError('Native annotation/protein binding mismatch: '+sid+'/'+name)
        if name in sequences and sequences[name]!=seq: raise ValueError('Conflicting exact protein binding')
        sequences[name]=seq
    expected={k for k,g in genes.items() if g.get('protein_sha256')}
    if set(sequences)!=expected:
        raise ValueError('Full HMM scope omits imported proteins: '+sid)
    _write_saved_rm_nominees(data,sid,genes,sequences,out)
    columns=list(dict.fromkeys(k for r in rows for k in r if not k.startswith('REBASEfinder_')
                              and k not in {'nt_seq','nt_seq_rev_comp','rearranged_nt_seq'}))
    _table(out/'annotations.tsv',rows,columns)
    fasta=out/'proteins.faa'
    fasta.write_text(''.join('>'+k+'\n'+v.rstrip('*')+'\n' for k,v in sorted(sequences.items())),encoding='utf-8')
    maps=[dict(sample_id=sid,query_id=k,locus_tag=k,aa_length=len(v.rstrip('*')),
               protein_sha256=hashlib.sha256(v.rstrip('*').encode()).hexdigest()) for k,v in sorted(sequences.items())]
    _table(out/'gene_map.tsv',maps,['sample_id','query_id','locus_tag','aa_length','protein_sha256'])
    genomes=sorted({*directory.glob('*.gbff'),*directory.glob('*.gbk'),*directory.glob('*.gb')})
    if len(genomes)!=1: raise ValueError('Full evaluation needs exactly one native GenBank file for '+sid)
    write_json(out/'protein_binding.json',dict(sample_id=sid,genbank=str(genomes[0]),genbank_sha256=digest(genomes[0]),
        protein_count=len(maps),normalization='uppercase whitespace-free; terminal stop removed for HMM only',
        source_protein_sha256={k:genes[k]['protein_sha256'] for k in sorted(sequences)}))
    return genomes[0]


def execute_full_pipeline(data,root,config_path,workspace):
    """Return completed per-sample stages, or fail before publishing any figure."""
    from .full_hmm import run as run_hmm
    cfg=read_engine_config(config_path)
    work=Path(workspace).resolve();work.mkdir(parents=True,exist_ok=False)
    cfg['blastdbcmd']=_resolve_blastdbcmd(cfg)
    env=_runtime_environment(cfg,work/'cache')
    resources={str(Path(cfg[k])):digest(cfg[k]) for k in ('_config_path','hmm_database','profile_metadata','rscript','blastp','blastdbcmd','hmmsearch','cath_resolve_hits')}
    resources[str(Path(__file__).resolve())] = digest(__file__)
    for p in cfg['runtime_dependencies']:resources[p]=digest(p)
    native=Path(str(files('defenseviz').joinpath('data/native')))
    resources.update(_verify_bundled_sources(native))
    actual_r_runtime,r_runtime_resources=_capture_r_runtime(cfg,work,env)
    resources.update(r_runtime_resources)
    versions={}
    for key,args in [('rscript',['--version']),('blastp',['-version']),('blastdbcmd',['-version']),('hmmsearch',['-h']),('cath_resolve_hits',['--version'])]:
        proc=subprocess.run([cfg[key],*args],capture_output=True,text=True,timeout=30,env=env)
        versions[key]=dict(argv=[cfg[key],*args],exit_code=proc.returncode,stdout=proc.stdout,stderr=proc.stderr)
    write_json(work/'tool_versions.json',versions)
    resources[str(work/'tool_versions.json')]=digest(work/'tool_versions.json')
    # Copy resource files once per pipeline. The native runner cannot update
    # original reference resources even if an upstream package has cache hooks.
    source=Path(cfg['reference_cache']);private=work/'cache'/'db_modules'/'rebasefinder'/'cache';private.mkdir(parents=True)
    reference=[source/'rebase_data.rds',source/'rebase_db.fasta',*sorted(source.glob('rebase_db.fasta.*'))]
    for name in ('rebase_db.fasta.phr','rebase_db.fasta.pin','rebase_db.fasta.psq'):
        if not (source/name).is_file():raise FileNotFoundError('Complete BLAST database missing '+name)
    private_hashes={}
    for p in reference:
        if not p.is_file():raise FileNotFoundError(p)
        resources[str(p.resolve())]=digest(p)
        shutil.copyfile(p,private/p.name)
        if digest(private/p.name)!=resources[str(p.resolve())]:raise ValueError('Reference copy hash mismatch')
        private_hashes[str((private/p.name).resolve())]=resources[str(p.resolve())]
    resources.update(private_hashes)
    database_validation=_validate_reference_index(private/'rebase_db.fasta',cfg['blastdbcmd'],work/'reference_index_validation.json',env)
    resources[str(work/'reference_index_validation.json')]=digest(work/'reference_index_validation.json')
    database_binding=dict(prefix=str(private/'rebase_db.fasta'),files=private_hashes,
        validation_receipt=dict(path=str(work/'reference_index_validation.json'),sha256=digest(work/'reference_index_validation.json')),
        entries=database_validation['entries'])
    runtime={k:v for k,v in cfg.items() if not k.startswith('_')}
    runtime['reference_cache']=str(private)
    write_json(work/'runtime.json',runtime)
    resources[str(work/'runtime.json')]=digest(work/'runtime.json')
    stages=[]
    for sample in data['samples']:
        sid=sample['sample_id'];key=hashlib.sha256(sid.encode()).hexdigest()[:16]
        sample_dir=work/key;sample_dir.mkdir();inputs=sample_dir/'inputs';inputs.mkdir()
        directory=_sample_path(data,sid,root)
        gb=_prepare_sample(data,sid,directory,inputs);resources[str(gb.resolve())]=digest(gb)
        # Genbank_organizer writes an adjacent FAA even with save_output=FALSE.
        # Keep that side effect in the private workspace, not the input folder.
        private_genbank_dir=inputs/'genbank';private_genbank_dir.mkdir()
        private_gb=private_genbank_dir/gb.name
        shutil.copyfile(gb,private_gb)
        gb_sha256=resources[str(gb.resolve())]
        if digest(gb)!=gb_sha256 or digest(private_gb)!=gb_sha256:
            raise ValueError('GenBank copy hash mismatch: '+str(gb))
        resources[str(private_gb.resolve())]=gb_sha256
        hmm=sample_dir/'hmm'
        run_hmm(SimpleNamespace(sample=sid,fasta=str(inputs/'proteins.faa'),gene_map=str(inputs/'gene_map.tsv'),
            hmm=cfg['hmm_database'],profiles=cfg['profile_metadata'],hmmsearch=cfg['hmmsearch'],
            cath=cfg['cath_resolve_hits'],cpu=cfg['cpu'],runtime_dependency=cfg['runtime_dependencies'],
            output=str(hmm),reuse=False))
        hmm_bindings={str(p.resolve()):digest(p) for p in hmm.rglob('*') if p.is_file()}
        if str(hmm/'execution.json') not in hmm_bindings:
            raise ValueError('HMM computation lacks its execution receipt')
        nominees=inputs/'hmm_nominees.tsv'
        nomination=_nomination_union(data,sid,hmm/'nominees.tsv',nominees)
        write_json(inputs/'nomination.json',nomination)
        output=sample_dir/'native';output.mkdir()
        cmd=[cfg['rscript'],'--vanilla',str(native/'full_rm.R'),str(work/'runtime.json'),str(native),
             str(private_gb),str(inputs/'annotations.tsv'),str(nominees),str(work/'cache'),str(output),str(cfg['cpu'])]
        input_bindings={**resources,**data.get('input_sha256',{}),**hmm_bindings,
            **{str(p.resolve()):digest(p) for p in inputs.rglob('*') if p.is_file()}}
        for path,expected in input_bindings.items():
            if digest(path)!=expected:raise ValueError('Full native input changed before execution: '+path)
        start=time.monotonic()
        with (sample_dir/'stdout.log').open('w') as stdout,(sample_dir/'stderr.log').open('w') as stderr:
            proc=subprocess.run(cmd,cwd=output,stdout=stdout,stderr=stderr,env=env)
        receipt=dict(sample_id=sid,argv=cmd,exit_code=proc.returncode,elapsed_seconds=time.monotonic()-start,
                     executable_sha256=digest(cfg['rscript']),native_source_sha256=dict(resources),
                     input_sha256=input_bindings,blast_database=database_binding,
                     r_runtime=dict(path=str(work/'r_runtime.json'),sha256=digest(work/'r_runtime.json')),
                     environment={k:env[k] for k in ('PATH','R_LIBS','DNMB_BLASTP')})
        write_json(sample_dir/'execution.json',receipt)
        if proc.returncode:raise RuntimeError('DNMB full R-M stage failed for '+sid+'; see '+str(sample_dir/'stderr.log'))
        completion=json.loads((output/'completion.json').read_text())
        if completion.get('passed') is not True or completion.get('full_reference_search') is not True:
            raise ValueError('Native R-M stage did not complete its declared search scope')
        runtime_check = _verify_native_runtime(output/'runtime_packages.json', actual_r_runtime, r_runtime_resources)
        runtime_check['preflight'] = dict(path=str(work/'r_runtime.json'), sha256=digest(work/'r_runtime.json'))
        write_json(output/'runtime_validation.json', runtime_check)
        for path,expected in input_bindings.items():
            if digest(path)!=expected:raise ValueError('Full native input changed during execution: '+path)
        input_bindings[str(sample_dir/'execution.json')]=digest(sample_dir/'execution.json')
        write_json(output/'stage_manifest.json',dict(schema_version=1,sample_id=sid,exit_code=0,
            artifacts={p.relative_to(output).as_posix():digest(p) for p in output.rglob('*') if p.is_file()},
            input_sha256=input_bindings,blast_database=database_binding,
            native_execution=dict(path=str(sample_dir/'execution.json'),sha256=digest(sample_dir/'execution.json')),
            hmm_execution=dict(path=str(hmm/'execution.json'),sha256=hmm_bindings[str(hmm/'execution.json')]),
            r_runtime=dict(path=str(work/'r_runtime.json'),sha256=digest(work/'r_runtime.json')),
            capture='immediately_after_successful_process_exit'))
        stages.append(dict(sample_id=sid,native_dir=str(output),hmm_dir=str(hmm)))
    for path,expected in resources.items():
        if digest(path)!=expected:raise ValueError('Full pipeline input changed: '+path)
    searches=[json.loads((Path(s['native_dir'])/'search_commands.json').read_text()) for s in stages]
    write_json(work/'execution.json',dict(schema_version=1,stages=stages,input_sha256=resources,
        completed=True,execution='fresh_RM_HMM_and_full_reference_search_then_native_context_and_policy_features',
        detector_runs=dict(RM_HMM=len(stages),REBASE_primary=sum(bool(s) for s in searches)),
        hmm_search_runs=len(stages),homology_search_runs=sum(len(s) for s in searches),
        reference_index_validation_runs=1,reference_entries=database_validation['entries'],
        actual_R_version=actual_r_runtime['R_version'],original_defense_tool_runs=0,source_inputs_unchanged=True))
    return stages

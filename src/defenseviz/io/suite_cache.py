"""Attach exact-protein-bound saved DNMB evidence without replacing live calls.

Old stage caches can contain richer alignments than the live TSV while belonging
to an earlier candidate scope. Their observations remain a separate auxiliary
channel; no cache count, type label or missing row is a new retention decision.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from importlib.resources import files
from pathlib import Path
import shutil
import subprocess

_HSP_COLUMNS = ('query_id','subject_id','pident','length','mismatch','gapopen',
                'qstart','qend','sstart','send','evalue','bitscore')


def _hash(path, algorithm='sha256'):
    before = path.stat()
    h = hashlib.new(algorithm)
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f'{path}: saved evidence changed while hashing')
    return h.hexdigest()


def _read_cache(path, rscript):
    helper = Path(str(files('defenseviz').joinpath('data/read_suite_evidence.R')))
    command = [str(rscript), '--vanilla', str(helper), str(path)]
    done = subprocess.run(command, text=True, encoding='utf-8', capture_output=True)
    if done.returncode:
        raise ValueError(f'{path}: cannot read saved DNMB evidence: {done.stderr.strip()}')
    try:
        rows = [json.loads(line) for line in done.stdout.splitlines() if line.strip()]
    except (ValueError, TypeError) as exc:
        raise ValueError(f'{path}: malformed saved evidence projection') from exc
    if not rows or rows[0].get('schema_version') != 1:
        raise ValueError(f'{path}: unsupported saved evidence projection')
    return rows, {'command': command, 'exit_code': done.returncode,
                  'reader_sha256': _hash(helper), 'stderr': done.stderr.strip()}


def _directory(data, sid, root):
    # Resolved source paths avoid guessing how a manifest-relative label was
    # originally based. Every module in a single suite sample shares one root.
    dirs = set()
    for row in data.get('suite', {}).get('sources', []):
        if row.get('sample_id') != sid or not row.get('resolved_path'):
            continue
        p = Path(row['resolved_path']).resolve()
        for parent in p.parents:
            if parent.name.startswith('dnmb_module_'):
                dirs.add(parent.parent)
                break
    if len(dirs) > 1:
        raise ValueError(f'{sid}: multiple native sample directories for cache evidence')
    if dirs:
        return next(iter(dirs))
    for row in data.get('source_metadata', []):
        if row.get('sample_id') == sid and row.get('source_directory'):
            p = Path(row['source_directory'])
            return (p if p.is_absolute() else root / p).resolve()
    return None


def _fasta(path):
    result, name, sequence = {}, None, []
    def flush():
        if name is not None:
            value = ''.join(sequence).upper()
            if name in result and result[name] != value:
                raise ValueError(f'{path}: conflicting saved query protein {name}')
            result[name] = value
    with path.open(encoding='utf-8-sig') as stream:
        for line in stream:
            if line.startswith('>'):
                flush()
                parts = line[1:].split()
                if not parts:
                    raise ValueError(f'{path}: empty saved query identity')
                name, sequence = parts[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(f'{path}: sequence before FASTA header')
                sequence.append(''.join(line.split()))
    flush()
    return result


def _sequence_hash(sequence):
    return hashlib.sha256(''.join(str(sequence).split()).upper().encode()).hexdigest()


def attach_cached_evidence(data, root, rscript=None):
    """Return an enriched copy; original calls, decisions and fields are intact.

    ``suite.cached_features`` and ``suite.cached_blast`` are auxiliary rows with
    exact current protein hashes and native source rows. Consumers must still
    match the measured subject before using alignment values and cannot treat
    cache absence as a completed search or native type as adjudicated truth.
    """
    root = Path(root).resolve()
    result = dict(data)
    suite = dict(data.get('suite', {}))
    result['suite'] = suite
    if suite.get('native_cache_evidence'):
        raise ValueError('saved cache evidence is already attached')
    result['inputs'] = list(data.get('inputs', []))
    result['input_sha256'] = dict(data.get('input_sha256', {}))
    result['input_roles'] = {p:list(v) for p,v in data.get('input_roles', {}).items()}
    suite['sources'] = list(suite.get('sources', []))
    features, hsps, audits = [], [], []
    genes = {(g['sample_id'],g['gene_id']):g for g in data.get('genes', [])}
    executable = shutil.which(str(rscript)) if rscript else shutil.which('Rscript')
    if rscript and executable is None:
        raise ValueError('Configured Rscript was not found for saved evidence')

    def pin(path, sid, kind):
        path = path.resolve()
        digest = _hash(path)
        previous = result['input_sha256'].get(str(path))
        if previous and previous != digest:
            raise ValueError(f'{path}: source changed since native import')
        result['input_sha256'][str(path)] = digest
        if str(path) not in {str(p) for p in result['inputs']}:
            result['inputs'].append(str(path))
        role = dict(sample_id=sid, module='rebasefinder', kind=kind, native_result=True)
        roles = result['input_roles'].setdefault(str(path), [])
        if role not in roles: roles.append(role)
        ref = dict(path=str(path), sha256=digest)
        if not any(s.get('resolved_path') == str(path) and s.get('role') == kind for s in suite['sources']):
            suite['sources'].append(dict(source_id='cache-source-'+hashlib.sha256((sid+str(path)+kind).encode()).hexdigest()[:20],
                sample_id=sid,module='rebasefinder',role=kind,selected=True,
                resolved_path=str(path),native_result=True,**ref))
        return ref

    for sample in data.get('samples', []):
        sid = sample['sample_id']
        directory = _directory(data, sid, root)
        audit = dict(sample_id=sid, status='missing', cache_candidate_scope='auxiliary_only',
                     search_scope='historical_reported_hsps', source_refs=[], tables=[],
                     omitted=[], commands=[], genome_binding='not_available')
        audits.append(audit)
        if directory is None:
            audit['reason'] = 'sample_source_directory_unavailable'
            continue
        module = directory / 'dnmb_module_rebasefinder'
        cache = directory / '.dnmb_module_stage_cache.rds'
        cache_rows, cache_ref = [], None
        if cache.is_file():
            cache_ref = pin(cache,sid,'saved_cache_auxiliary')
            audit['source_refs'].append(cache_ref)
            if executable:
                projection, command = _read_cache(cache,executable)
                if _hash(cache) != cache_ref['sha256']:
                    raise ValueError(f'{cache}: source changed during evidence extraction')
                audit['commands'].append(command)
                metadata = projection[0]
                audit.update(status=metadata['status'], native_module_status=metadata.get('module_status'),
                             written_at=metadata.get('written_at'))
                audit['tables'] = [r for r in projection if r.get('kind') == 'table']
                cache_rows = [r for r in projection if r.get('kind') == 'row']
                signature = metadata.get('genbank_files') or []
                if isinstance(signature,dict): signature=[signature]
                checked=[]
                for g in signature:
                    name=g.get('file')
                    if not isinstance(name,str) or Path(name).name!=name or not g.get('md5'):
                        continue
                    path=directory/name
                    checked.append(path.is_file() and _hash(path,'md5')==g['md5'])
                    if path.is_file():audit['source_refs'].append(pin(path,sid,'cache_genome_binding'))
                if checked: audit['genome_binding']='matched' if all(checked) else 'mismatched'
            else:
                audit.update(status='reader_unavailable',reason='Rscript required to inspect the saved stage cache')
        bindings={}
        fasta_bindings={}
        # The source query FASTA proves that raw alignment rows refer to the
        # current whole protein. No hash is inferred from an alignment length.
        for relative in ('query_sequences.fasta','supplemental_rebase_blast/supplemental_queries.faa'):
            path=module/relative
            if not path.is_file(): continue
            ref=pin(path,sid,'saved_query_proteins'); audit['source_refs'].append(ref)
            fasta_bindings[relative]={}
            for gid,sequence in _fasta(path).items():
                gene=genes.get((sid,gid)); digest=_sequence_hash(sequence)
                if gene and gene.get('protein_sha256')==digest:
                    bindings.setdefault(gid,[]).append(ref)
                    fasta_bindings[relative][gid]=[ref]
                else:
                    audit['omitted'].append(dict(gene_id=gid,reason='saved_query_protein_not_bound',source_ref=ref))
            if _hash(path)!=ref['sha256']:raise ValueError(f'{path}: saved query FASTA changed')
        pending=[]
        for row in cache_rows:
            native=row.get('fields',{}); table=row['table']; gid=native.get('locus_tag',native.get('query',native.get('query_id')))
            ref={**cache_ref,'rds_table':'module_results/REBASEfinder/module_result/'+table,'row':row['row']}
            sequence=native.get('translation')
            gene=genes.get((sid,gid))
            if sequence and gene and gene.get('protein_sha256')==_sequence_hash(sequence):
                bindings.setdefault(gid,[]).append(ref)
            pending.append((gid,native,ref,table))
        def accept(gid,native,ref,channel,required_bindings=None):
            gene=genes.get((sid,gid))
            effective_bindings=bindings if required_bindings is None else required_bindings
            if not gene or not gene.get('protein_sha256') or gid not in effective_bindings:
                audit['omitted'].append(dict(gene_id=gid,reason='cached_observation_not_protein_bound',source_ref=ref))
                return
            if native.get('translation') and _sequence_hash(native['translation'])!=gene['protein_sha256']:
                audit['omitted'].append(dict(gene_id=gid,reason='cached_row_protein_mismatch',source_ref=ref))
                return
            common=dict(sample_id=sid,gene_id=gid,protein_sha256=gene['protein_sha256'],
                        evidence_kind='saved_cache_auxiliary',source_refs=[ref],
                        protein_source_refs=list({json.dumps(r,sort_keys=True):r for r in effective_bindings[gid]}.values()))
            if all(k in native for k in _HSP_COLUMNS):
                try:
                    values={k:float(native[k]) for k in _HSP_COLUMNS[2:]}
                    valid=(all(math.isfinite(v) for v in values.values()) and bool(native.get('subject_id')) and
                           0<=values['pident']<=100 and values['length']>0 and
                           min(values[k] for k in ('mismatch','gapopen','evalue','bitscore'))>=0 and
                           min(values[k] for k in ('qstart','qend','sstart','send'))>=1)
                    if gene.get('protein_length'):
                        valid=valid and max(values['qstart'],values['qend'])<=gene['protein_length']
                    if not valid:raise ValueError('HSP bounds')
                except (ValueError,TypeError):
                    audit['omitted'].append(dict(gene_id=gid,reason='invalid_saved_alignment',source_ref=ref));return
                normalized=dict(query_id=gid,subject_id=str(native['subject_id']),**values)
                for name in ('length','mismatch','gapopen','qstart','qend','sstart','send'):
                    if normalized[name]!=int(normalized[name]):
                        raise ValueError(f'{ref}: noninteger saved alignment coordinate/count')
                    normalized[name]=int(normalized[name])
                hsps.append(dict(**common,source_channel=channel,native_fields=normalized,
                                 search_scope='historical_reported_hsps',search_is_exhaustive=False))
            else:
                fields={k:v for k,v in native.items() if k!='translation'}
                features.append(dict(**common,source_channel=channel,native_fields=fields,
                                     genome_binding=audit['genome_binding']))
        for gid,native,ref,table in pending: accept(gid,native,ref,table)
        for relative in ('blast_results.txt','supplemental_rebase_blast/supplemental_blast_results.txt',
                         'supplemental_rebase_blast/supplemental_best_hits.tsv'):
            path=module/relative
            if not path.is_file():continue
            ref=pin(path,sid,'saved_blast_auxiliary');audit['source_refs'].append(ref)
            with path.open(encoding='utf-8-sig',newline='') as stream:
                if path.suffix=='.tsv':
                    reader=csv.DictReader(stream,delimiter='\t');start=2
                else:
                    reader=csv.reader(stream,delimiter='\t');start=1
                for n,row in enumerate(reader,start):
                    if not row:continue
                    if isinstance(row,list):
                        if len(row)!=len(_HSP_COLUMNS):
                            raise ValueError(f'{path}:{n}: expected 12-column saved BLAST output')
                        row=dict(zip(_HSP_COLUMNS,row))
                    query_file=('supplemental_rebase_blast/supplemental_queries.faa' if relative.startswith('supplemental_') else 'query_sequences.fasta')
                    accept(row.get('query_id'),row,{**ref,'row':n},relative,fasta_bindings.get(query_file,{}))
            if _hash(path)!=ref['sha256']:raise ValueError(f'{path}: saved BLAST changed')
        live={r['gene_id'] for r in data.get('suite',{}).get('features',[]) if r.get('sample_id')==sid}
        cached={gid for gid,_,_,table in pending if table=='hits'}
        audit.update(current_candidate_genes=len(live),cached_hit_genes=len(cached),
                     cached_hit_scope_differs=bool(cached and cached!=live),
                     bound_cached_feature_observations=sum(r['sample_id']==sid for r in features),
                     bound_cached_hsp_observations=sum(r['sample_id']==sid for r in hsps))
    def dedup(rows):
        merged={}
        for row in rows:
            key=json.dumps({k:v for k,v in row.items() if k not in {'source_refs','protein_source_refs','source_channel'}},sort_keys=True)
            if key not in merged:
                merged[key]={**row,'source_channels':[row['source_channel']]}
                continue
            previous=merged[key]
            for field in ('source_refs','protein_source_refs'):
                previous[field]=list({json.dumps(r,sort_keys=True):r for r in previous[field]+row[field]}.values())
            previous['source_channels']=sorted(set(previous['source_channels']+[row['source_channel']]))
        return list(merged.values())
    suite['cached_features']=dedup(features)
    suite['cached_blast']=dedup(hsps)
    for audit in audits:
        sid=audit['sample_id']
        audit['bound_cached_features']=sum(r['sample_id']==sid for r in suite['cached_features'])
        audit['bound_cached_hsps']=sum(r['sample_id']==sid for r in suite['cached_blast'])
    suite['native_cache_evidence']=dict(schema_version=1,samples=audits,detector_runs=0,
                                      homology_search_runs=0,retention_decisions=0)
    return result

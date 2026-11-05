#!/usr/bin/env python3
"""DNMB all-translated-CDS HMM nomination; no annotation/threshold fitting."""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def table(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_table(path, rows, fields):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def sequences(path):
    result, name, chunks = {}, None, []
    def save():
        if name is not None:
            if name in result:
                raise ValueError('Duplicate FASTA ID: ' + name)
            seq = ''.join(chunks).upper().rstrip('*')
            if not seq:
                raise ValueError('Empty sequence: ' + name)
            result[name] = seq
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            save()
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            if name is None:
                raise ValueError('FASTA sequence before header')
            chunks.append(''.join(line.split()))
    save()
    return result


def validate_mapping(faa, mapping, sample):
    seqs = sequences(faa)
    rows = table(mapping)
    ids = [r['query_id'] for r in rows]
    if len(ids) != len(set(ids)) or set(ids) != set(seqs):
        raise ValueError('Gene map and FASTA must have identical unique IDs')
    loci = [r['locus_tag'] for r in rows]
    duplicate_loci = {v for v in loci if v and loci.count(v) > 1}
    mapped = {}
    for row in rows:
        q = row['query_id']
        if row.get('sample_id') != sample:
            raise ValueError('Wrong sample in gene map')
        digest = hashlib.sha256(seqs[q].encode()).hexdigest()
        if digest != row['protein_sha256'] or len(seqs[q]) != int(row['aa_length']):
            raise ValueError('Sequence identity mismatch: ' + q)
        row['mapping_status'] = ('missing_locus' if not row['locus_tag'] else
                                 'ambiguous_locus' if row['locus_tag'] in duplicate_loci else 'unique')
        mapped[q] = row
    return seqs, mapped


DOMAIN_FIELDS = ['sample_id', 'query_id', 'locus_tag', 'protein_sha256', 'mapping_status',
                 'profile_name', 'profile_accession', 'sequence_length', 'model_length',
                 'sequence_evalue', 'sequence_score', 'domain_number', 'domain_count',
                 'conditional_evalue', 'domain_evalue', 'domain_score', 'hmm_from', 'hmm_to',
                 'ali_from', 'ali_to', 'env_from', 'env_to', 'accuracy', 'off_target',
                 'role_hints', 'type_hints', 'resolved_selected', 'resolved_boundaries',
                 'database_sha256', 'source', 'command_id']
NOMINEE_FIELDS = DOMAIN_FIELDS + ['role_hint']


def parse_domains(path, sample, mapping, profiles, database_sha):
    rows = []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line.strip() or line.startswith('#'):
            continue
        a = line.split(maxsplit=22)
        if len(a) < 22:
            raise ValueError('Malformed domtblout line ' + str(number))
        q, p = a[0], a[3]
        if q not in mapping or p not in profiles:
            raise ValueError('Unmapped query/profile: ' + q + '/' + p)
        m, prof = mapping[q], profiles[p]
        if int(a[2]) != int(m['aa_length']):
            raise ValueError('HMM output sequence length mismatch')
        r = dict(sample_id=sample, query_id=q, locus_tag=m['locus_tag'],
                 protein_sha256=m['protein_sha256'], mapping_status=m['mapping_status'],
                 profile_name=p, profile_accession=a[4], sequence_length=int(a[2]),
                 model_length=int(a[5]), sequence_evalue=a[6], sequence_score=a[7],
                 domain_number=int(a[9]), domain_count=int(a[10]), conditional_evalue=a[11],
                 domain_evalue=a[12], domain_score=a[13], hmm_from=int(a[15]),
                 hmm_to=int(a[16]), ali_from=int(a[17]), ali_to=int(a[18]),
                 env_from=int(a[19]), env_to=int(a[20]), accuracy=a[21],
                 off_target=prof['off_target'], role_hints=prof['roles'],
                 type_hints=prof['type_hints'], resolved_selected=False, resolved_boundaries='',
                 database_sha256=database_sha, source='DNMB_hmm_nomination', command_id='hmmsearch')
        rows.append(r)
    return rows


def resolve_domains(path, domains):
    index = {}
    for r in domains:
        k = (r['query_id'], r['profile_name'], f"{r['env_from']}-{r['env_to']}")
        index.setdefault(k, []).append(r)
    selected = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        a = line.split()
        if len(a) != 7:
            raise ValueError('Unexpected CATH output schema: ' + line)
        matches = index.get((a[0], a[1], a[3]), [])
        matches = [r for r in matches if float(r['domain_score']) == float(a[2])
                   and float(r['domain_evalue']) == float(a[6])]
        if len(matches) != 1:
            raise ValueError('CATH domain does not map uniquely to raw evidence: ' + line)
        r = matches[0]
        if r['resolved_selected']:
            raise ValueError('Duplicate resolved domain: ' + line)
        r['resolved_selected'], r['resolved_boundaries'] = True, a[4]
        selected.append(r)
    return selected


def nominations(selected):
    result = []
    for r in selected:
        if str(r['off_target']).lower() in ('true', '1') or r['mapping_status'] != 'unique':
            continue
        roles = sorted(set(r['role_hints'].split(';'))) if r['role_hints'] else ['unknown']
        if not set(roles) <= {'R', 'M', 'S', 'unknown'}:
            raise ValueError('Invalid role hint')
        result.extend(dict(r, role_hint=role) for role in roles)
    return result


def execute(cmd, prefix, output, env):
    start = time.perf_counter()
    executable_sha256 = sha(cmd[0])
    with (output / (prefix + '.stdout.log')).open('w') as stdout, (output / (prefix + '.stderr.log')).open('w') as stderr:
        proc = subprocess.run(cmd, stdout=stdout, stderr=stderr, env=env)
    unchanged = sha(cmd[0]) == executable_sha256
    receipt = dict(command=cmd, executable_sha256=executable_sha256, executable_unchanged=unchanged,
                   exit_code=proc.returncode, elapsed_seconds=time.perf_counter()-start)
    (output / (prefix + '.json')).write_text(json.dumps(receipt, indent=2)+'\n')
    if proc.returncode:
        raise RuntimeError(prefix + ' failed; see saved logs')
    if not unchanged:
        raise RuntimeError(prefix + ' executable changed during the search')
    return receipt


def run(args):
    start = time.perf_counter()
    paths = {k: Path(getattr(args, k)).resolve() for k in ('fasta', 'gene_map', 'hmm', 'profiles', 'hmmsearch', 'cath')}
    inputs = {k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()}
    for dep in args.runtime_dependency:
        p = Path(dep).resolve()
        inputs['runtime_dependency:' + str(p)] = dict(path=str(p), sha256=sha(p))
    inputs['wrapper'] = dict(path=str(Path(__file__).resolve()), sha256=sha(__file__))
    signature = dict(inputs=inputs, sample_id=args.sample, cpu=args.cpu, threshold='model_GA',
                     overlap_policy='CATH_default', off_target_policy='resolved_non_offtarget_only')
    cache_key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
    output = Path(args.output).resolve()
    if output.exists():
        if not args.reuse:
            raise FileExistsError('Output exists; use explicit --reuse with hash verification')
        old = json.loads((output/'execution.json').read_text())
        if old.get('exit_code') != 0 or old.get('cache_key') != cache_key:
            raise ValueError('Cache signature mismatch or unsuccessful prior run')
        for name, digest in old['artifact_hashes'].items():
            if sha(output/name) != digest:
                raise ValueError('Cached artifact modified: '+name)
        reuse = dict(run_kind='validated_cache_reuse', cache_key=cache_key, exit_code=0,
                     original_receipt=str(output/'execution.json'), elapsed_seconds=time.perf_counter()-start)
        (output/'reuse.json').write_text(json.dumps(reuse, indent=2)+'\n')
        return reuse
    output.mkdir(parents=True)
    record = dict(run_kind='fresh_HMM_candidate_search', cache_key=cache_key, signature=signature,
                  started_utc=dt.datetime.now(dt.timezone.utc).isoformat(), host_argv=list(sys.argv),
                  invocation=dict(sample=args.sample,cpu=args.cpu,**{k:str(v) for k,v in paths.items()}),
                  source_hmm_collection='configured profile database; exact source bound by signature, not independent comparator evidence')
    try:
        seqs, mapping = validate_mapping(paths['fasta'], paths['gene_map'], args.sample)
        profile_rows = table(paths['profiles'])
        profiles = {r['profile_name']: r for r in profile_rows}
        names = [line.split()[1] for line in paths['hmm'].read_text().splitlines() if line.startswith('NAME ')]
        if len(names) != len(set(names)) or set(names) != set(profiles) or len(profiles) != len(profile_rows):
            raise ValueError('Profile mapping does not match exact HMM names')
        env = os.environ.copy()
        env['LC_ALL'] = 'C'
        record['tool_versions'] = {}
        for key, flags in [('hmmsearch',['-h']),('cath',['--version'])]:
            command=[str(paths[key]),*flags]
            version=subprocess.run(command,capture_output=True,text=True,env=env,timeout=30)
            record['tool_versions'][key]=dict(command=command,exit_code=version.returncode,
                stdout=version.stdout,stderr=version.stderr,executable_sha256=inputs[key]['sha256'])
        record['hmm'] = execute([str(paths['hmmsearch']), '--cut_ga', '--cpu', str(args.cpu),
                                '--domtblout', str(output/'domains.domtblout'), str(paths['hmm']),
                                str(paths['fasta'])], 'hmmsearch', output, env)
        dom = parse_domains(output/'domains.domtblout', args.sample, mapping, profiles, inputs['hmm']['sha256'])
        if dom:
            record['resolution'] = execute([str(paths['cath']), '--input-format', 'hmmer_domtblout',
                                           str(output/'domains.domtblout'), '--hits-text-to-file',
                                           str(output/'resolved.hits')], 'cath', output, env)
        else:
            (output/'resolved.hits').write_text('# Zero reported HMM domains; no CATH invocation needed\n')
            record['resolution'] = dict(run_kind='not_needed_zero_domains', exit_code=0, elapsed_seconds=0)
        selected = resolve_domains(output/'resolved.hits', dom)
        nominees = nominations(selected)
        write_table(output/'all_domains.tsv', dom, DOMAIN_FIELDS)
        write_table(output/'resolved_domains.tsv', selected, DOMAIN_FIELDS)
        write_table(output/'nominees.tsv', nominees, NOMINEE_FIELDS)
        record.update(proteins=len(seqs), profiles=len(profiles), raw_domains=len(dom),
                      resolved_domains=len(selected), nominee_rows=len(nominees),
                      nominated_CDS=len({r['query_id'] for r in nominees}),
                      unresolved_mapping_CDS=len({r['query_id'] for r in selected if r['mapping_status']!='unique'}),
                      exit_code=0)
        for v in inputs.values():
            if sha(v['path']) != v['sha256']:
                raise RuntimeError('Input mutated during search: '+v['path'])
        record['inputs_unchanged'] = True
        record['artifact_hashes'] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
    except Exception as e:
        record.update(exit_code=1, error=str(e))
        raise
    finally:
        record['elapsed_seconds'] = time.perf_counter()-start
        (output/'execution.json').write_text(json.dumps(record, indent=2)+'\n')
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('sample', 'fasta', 'gene-map', 'hmm', 'profiles', 'hmmsearch', 'cath', 'output'):
        p.add_argument('--'+name, required=True)
    p.add_argument('--cpu', type=int, default=2)
    p.add_argument('--runtime-dependency', action='append', default=[])
    p.add_argument('--reuse', action='store_true')
    args = p.parse_args()
    if args.cpu < 1:
        p.error('--cpu must be positive')
    record = run(args)
    print(json.dumps({k:v for k,v in record.items() if k not in ('signature','artifact_hashes')}))


if __name__ == '__main__':
    main()

"""Source-bound presentation scopes, kept separate from primary analyses."""
from copy import deepcopy
import hashlib
import json
import re
from pathlib import Path


def model_digest(result):
    return hashlib.sha256(json.dumps(result,sort_keys=True,ensure_ascii=False,allow_nan=False).encode('utf-8')).hexdigest()


def _assembly_identities(samples):
    """Explicit assembly accessions from sample metadata (never display labels).

    Reuses the pipeline's specimen-identity rules on the rendered model, which
    carries sample metadata but no genome files. Samples without explicit
    assembly metadata stay explicitly unchecked beyond their sample IDs.
    """
    from ...suite_external import cohort_identities
    return cohort_identities({'samples':[{k:v for k,v in s.items() if k not in ('label','display')} for s in samples]})


def type_case_projection(result, external, selected, validate, source_root=None):
    """Build only type rows/sample labels; never merge inventory or overlap."""
    if external is None:external=[]
    if not isinstance(external,list):raise ValueError('external_type_cases must be a list')
    samples=deepcopy(result.get('samples',[]));rows=deepcopy(result.get('type_tables',[]))
    sources=deepcopy(result.get('sources',[]));known={s['sample_id'] for s in samples}
    external_ids=[];bindings=[];source_index={s['id']:s for s in sources}
    if external:
        from ...suite_external import reject_cohort_overlap
        identities=_assembly_identities(samples)
    for envelope in external:
        if not isinstance(envelope,dict) or set(envelope)-{'result','sample_ids','source_refs','result_sha256'}:
            raise ValueError('External type envelope requires result, optional sample_ids and source_refs')
        independent=envelope.get('result');validate(independent)
        refs=envelope.get('source_refs')
        if not isinstance(refs,list) or not refs or any(not isinstance(r,dict) or not isinstance(r.get('path'),str) or not r['path'] or not isinstance(r.get('sha256'),str) or not re.fullmatch(r'[0-9a-f]{64}',r['sha256']) for r in refs):
            raise ValueError('External type cases require source paths and SHA256 bindings')
        digest=model_digest(independent)
        if envelope.get('result_sha256',digest)!=digest:
            raise ValueError('External type result SHA256 does not match its supplied model')
        bound=False
        for ref in refs:
            path=Path(ref['path'])
            if not path.is_absolute():
                if source_root is None:raise ValueError('Relative external type refs require a result root')
                path=Path(source_root)/path
            if not path.is_file():raise ValueError('External type source file is missing')
            content=path.read_bytes()
            if hashlib.sha256(content).hexdigest()!=ref['sha256']:
                raise ValueError('External type source SHA256 mismatch')
            try:bound=bound or model_digest(json.loads(content))==digest
            except (UnicodeDecodeError,ValueError,TypeError):pass
        if not bound:raise ValueError('External type source refs must bind the separately analyzed result JSON')
        available={s['sample_id']:s for s in independent.get('samples',[])};ids=envelope.get('sample_ids',list(available))
        if not isinstance(ids,list) or not ids or any(not isinstance(s,str) for s in ids) or len(set(ids))!=len(ids) or set(ids)-set(available):raise ValueError('Invalid external type sample selection')
        if known&set(ids):raise ValueError('External type sample identity overlaps another primary/external case')
        if any(not any(r['sample_id']==sid for r in independent.get('type_tables',[])) for sid in ids):raise ValueError('External selected sample has no analyzed type rows')
        incoming=_assembly_identities([available[sid] for sid in ids])
        reject_cohort_overlap(incoming,identities)
        unchecked=sorted(r['sample_id'] for r in incoming+identities if not r['assembly_accessions'])
        identity_check=dict(method='explicit_assembly_accession_metadata',
            external_assembly_accessions=sorted(a for r in incoming for a in r['assembly_accessions']),
            status='no_shared_assembly_accession' if not unchecked else 'partly_unchecked_missing_assembly_metadata',
            samples_without_assembly_metadata=unchecked)
        identities.extend(incoming)
        known.update(ids);external_ids.extend(ids)
        for sid in ids:
            sample=deepcopy(available[sid]);sample['cohort']='external';samples.append(sample)
        rows.extend(deepcopy(r) for r in independent['type_tables'] if r['sample_id'] in ids)
        for source in independent.get('sources',[]):
            old=source_index.get(source['id'])
            if old and old.get('kind')!=source.get('kind'):raise ValueError('External type source endpoint conflicts with primary source')
            if not old:sources.append(deepcopy(source));source_index[source['id']]=source
        bindings.append(dict(sample_ids=ids,source_refs=deepcopy(refs),result_sha256=model_digest(independent),scope='independently_analyzed_external_type_cases_only',included_in_primary_aggregates=False,identity_check=identity_check))
    if selected is not None and (len(set(selected))!=len(selected) or set(selected)-known):raise ValueError('type_samples contains unknown or duplicate sample identities')
    if selected is None and external_ids:
        typed=[s['sample_id'] for s in result.get('samples',[]) if any(r['sample_id']==s['sample_id'] for r in rows)]
        highlighted=[s['sample_id'] for s in result.get('samples',[]) if s.get('highlight') and s['sample_id'] in typed]
        selected=(highlighted or typed)+external_ids
    return dict(samples=samples,type_tables=rows,sources=sources),selected,bindings


def fraction_projection(result, scope):
    """Strict I-III counts come from retained CDS identities and explicit roles."""
    if scope=='assigned_I_IV':return None
    if scope!='retained_I_III':raise ValueError('Unknown rm_fraction_scope')
    ledger=result.get('gene_ledger')
    if not isinstance(ledger,list):raise ValueError('retained_I_III fraction scope requires an explicit CDS gene_ledger')
    summaries={r['sample_id']:r for r in result.get('sample_summary',[])};out={};seen=set()
    sample_ids={s['sample_id'] for s in result.get('samples',[])}
    for r in ledger:
        if not isinstance(r,dict) or r.get('sample_id') not in sample_ids or not isinstance(r.get('gene_id'),str) or not r['gene_id']:raise ValueError('Fraction scope needs declared sample and CDS identities')
        key=r['sample_id'],r['gene_id']
        if key in seen:raise ValueError('Duplicate CDS identity in fraction scope')
        seen.add(key)
        if not isinstance(r.get('retained'),bool) or not isinstance(r.get('rm_assignments'),list) or any(not isinstance(a,dict) for a in r.get('rm_assignments',[])):raise ValueError('Strict fraction scope requires explicit retained flags and assignments')
    canonical=('Type I','Type II','Type III','I','II','III')
    for sample in result.get('samples',[]):
        sid=sample['sample_id'];retained=[r for r in ledger if r['sample_id']==sid and r['retained']]
        selected=[r for r in retained if any(a.get('rm_type') in canonical for a in r['rm_assignments'])]
        rr=[r for r in selected if any(a.get('rm_type') in canonical and a.get('role') in ('R','RM','R_component') for a in r['rm_assignments'])]
        unknown=any(not r['rm_assignments'] or any(a.get('rm_type') in (None,'Unresolved') for a in r['rm_assignments']) for r in retained)
        summary=summaries.get(sid,{});state=summary.get('current_retained_counts_status',summary.get('retained_counts_status','unassessed'))
        if unknown and state=='complete':state='partial'
        rstate=state
        if any(a.get('rm_type') in canonical and a.get('role') not in ('R','RM','R_component','M','S','MS','A','C') for r in selected for a in r['rm_assignments']) and rstate=='complete':rstate='partial'
        out[sid]=dict(R_CDS=len(rr),total_CDS=len(selected),R_status=rstate,total_status=state,status='partial' if 'partial' in (state,rstate) else state,scope='retained_I_III',R_member_keys=[[r['sample_id'],r['gene_id']] for r in rr],member_keys=[[r['sample_id'],r['gene_id']] for r in selected])
    return out

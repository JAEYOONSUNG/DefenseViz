"""Automatic source-bound policy evaluation with explicit missing-evidence bounds.

Bounds are conservative over-completions of missing measurements. They do not
create measurements, biological negatives, or new policy thresholds.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, fields
import hashlib
import itertools
import json
import time

from .analysis.evidence import parse_evidence, policy_metadata
from .analysis.policy import Evidence, Policy, evaluate

BOOLEAN_FIELDS = {f.name for f in fields(Evidence) if f.type is bool}
NUMERIC_FIELDS = {f.name for f in fields(Evidence)} - BOOLEAN_FIELDS - {
    'family','role','profiles','raw_nuclease_classes','reference_nuclease_anchor_profile_indices'}
PENALTIES = ('alternative_annotation','conflict','partial')
RM_TYPES = {'Type I','Type II','Type III','Type IV'}
RM_ROLES = {'R','M','S','RM'}


def policy_bounds(values, unknown_axes, policy=None):
    """Evaluate a sound score envelope without treating missing evidence as false.

    Sequence/motif/context/structure support is monotone under the frozen preset.
    Negative ambiguity is reversed; penalty combinations are enumerated because
    alternative-annotation and conflict penalties are mutually suppressing.
    Independent hypothetical completions deliberately overapproximate any source
    correlations, so bounds can be conservative but cannot assert false certainty.
    """
    p = policy or Policy()
    if asdict(p) != asdict(Policy()):
        raise ValueError('Automatic evaluation requires the unchanged selected policy')
    names = {f.name for f in fields(Evidence)}
    if set(values) != names:
        raise ValueError('Automatic evidence must name every policy field')
    unknown = set(unknown_axes)
    if unknown - names or unknown & {'family','role','valid_sequence'}:
        raise ValueError('Known family, role and valid sequence are required for policy bounds')
    for name in BOOLEAN_FIELDS:
        if name not in unknown and type(values[name]) is not bool:
            raise ValueError('Unassessed boolean must be named in unknown_axes: '+name)
    corners=[]
    for positive in (False,True):
        row=deepcopy(values)
        for name in unknown:
            if name in PENALTIES:
                continue
            if name in BOOLEAN_FIELDS:
                row[name] = not positive if name=='hmm_role_ambiguous' else positive
            elif name in NUMERIC_FIELDS:
                if name=='family_compatible':row[name]=positive
                else:row[name]=(0. if 'evalue' in name else 1.) if positive else None
            elif name=='profiles':
                row[name]=deepcopy(values[name] or [])
                if positive:
                    row[name].append(dict(family_compatible=True,role_compatible=True,
                        evalue=0.,query_coverage=1.,model_coverage=1.))
            elif name=='raw_nuclease_classes':
                row[name]=sorted(set(values[name] or []) | ({'PD-ExK'} if positive else set()))
            elif name=='reference_nuclease_anchor_profile_indices':
                row[name]=list(range(len(row.get('profiles') or []))) if positive else []
        # Unknown anchor links may reference the hypothetical profile appended
        # above, independent of iteration order over unknown field names.
        if 'reference_nuclease_anchor_profile_indices' in unknown:
            row['reference_nuclease_anchor_profile_indices']=list(range(len(row['profiles']))) if positive else []
        variable=[n for n in PENALTIES if n in unknown]
        for flags in itertools.product((False,True),repeat=len(variable)):
            candidate=deepcopy(row);candidate.update(zip(variable,flags))
            # Only admissible anchor completions belong to the evidence domain.
            if candidate['reference_nuclease_anchor_support'] is False:
                if 'reference_nuclease_anchor_profile_indices' in unknown:
                    candidate['reference_nuclease_anchor_profile_indices'] = []
            elif not candidate['reference_nuclease_anchor_profile_indices']:
                if 'reference_nuclease_anchor_support' in unknown:
                    candidate['reference_nuclease_anchor_support'] = False
            corners.append(evaluate(parse_evidence(candidate),p))
    scores=[d['score_total'] for d in corners];axes=[d['capped_evidence_axis_count'] for d in corners]
    lo,hi=min(scores),max(scores);amin,amax=min(axes),max(axes)
    keep=(True if lo>=p.retain_cutoff and amin>=p.min_capped_evidence_axes else
          False if hi<p.retain_cutoff or amax<p.min_capped_evidence_axes else None)
    support={d['integrated_R_component_support'] for d in corners}
    return dict(score_min=lo,score_max=hi,score_total=lo if lo==hi else None,
        capped_axes_min=amin,capped_axes_max=amax,keep=keep,
        integrated_R_component_support=next(iter(support)) if len(support)==1 else None,
        assessment='exact_policy_evaluation' if not unknown else 'conservative_policy_bounds',
        unknown_axes=sorted(unknown),hypothetical_evaluations=len(corners),
        thresholds_changed=False,activity_inferred=False)


def evaluation_identity(compiled):
    """Hash content and scientific bindings, independently of local file paths."""
    unordered = {'records','sources','candidate_scope','source_refs','assignment_support',
                 'unknown_axes','reason_codes','gene_ids','source_scope'}
    omitted = {'path','resolved_path','source_path','source_id','elapsed_seconds','created_at'}
    def canonical(value, key=None):
        if isinstance(value, dict):
            return {k:canonical(v,k) for k,v in sorted(value.items()) if k not in omitted}
        if isinstance(value, (list,tuple)):
            items=[canonical(v) for v in value]
            return sorted(items,key=lambda v:json.dumps(v,sort_keys=True,separators=(',',':'))) if key in unordered else items
        return value
    payload=canonical(compiled)
    digest=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    return digest,payload


def refresh_decision_audit(data):
    """Record every changed decision field, including unchanged-state reviews."""
    meta=data['rm_evaluation']
    prior={(r['sample_id'],r['gene_id']):r for r in meta['native_saved_decisions']}
    records={(r['sample_id'],r['gene_id']):r for r in meta['records']}
    transitions=[]
    for row in sorted(data.get('saved_decisions',[]),key=lambda r:(r['sample_id'],r['gene_id'])):
        key=(row['sample_id'],row['gene_id']);before=prior.get(key,{})
        changes=[dict(field=name,before_present=name in before,after_present=name in row,
                      before=deepcopy(before.get(name)),after=deepcopy(row.get(name)))
                 for name in sorted(set(before)|set(row))
                 if (name in before)!=(name in row) or before.get(name)!=row.get(name)]
        if not changes:continue
        record=records.get(key,{})
        transitions.append(dict(sample_id=key[0],gene_id=key[1],previous_state=before.get('state'),
            current_state=row['state'],origin=row.get('origin'),rule_id=record.get('reason'),
            source_scope=deepcopy(record.get('source_scope',[])),
            change='new_integrated_nominee' if key not in prior else 'evaluated_native_nominee',
            state_changed=before.get('state')!=row['state'],field_changes=changes))
    meta['decision_transitions']=transitions
    meta['stored_vs_auto_summary']=dict(
        changed_records=len(transitions),new_records=sum(t['change']=='new_integrated_nominee' for t in transitions),
        state_changed_records=sum(t['state_changed'] for t in transitions),
        field_only_changed_records=sum(not t['state_changed'] for t in transitions),
        changed_field_counts={name:sum(any(c['field']==name for c in t['field_changes']) for t in transitions)
                              for name in sorted({c['field'] for t in transitions for c in t['field_changes']})},
        before_retention_coverage=deepcopy(meta.get('native_saved_retention_coverage',[])),
        after_retention_coverage=deepcopy([c for c in data.get('coverage',[]) if c.get('endpoint') in {'saved_retention','saved_decisions'}]))
    return data


def apply_automatic_rm(data):
    """Evaluate every compiled candidate, retaining all native records for audit."""
    from .automatic_rm import compile_rm_evidence
    started=time.monotonic()
    compiled=compile_rm_evidence(data)
    from .automatic_neighborhood import enrich_neighborhood_evidence
    compiled=enrich_neighborhood_evidence(data, compiled)
    result=deepcopy(data)
    prior={(r['sample_id'],r['gene_id']):deepcopy(r) for r in data.get('saved_decisions',[])}
    if len(prior) != len(data.get('saved_decisions', [])):
        raise ValueError('Duplicate native saved decision')
    decisions=deepcopy(prior);records=[];p=Policy()
    parameters=asdict(p)
    policy_sha=hashlib.sha256(json.dumps(parameters,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    evaluation_id, identity_payload=evaluation_identity(compiled)
    genes={(g['sample_id'],g['gene_id']):g for g in data.get('genes',[])}
    seen=set();evaluations=0
    for source in compiled['records']:
        rec=deepcopy(source);key=(rec['sample_id'],rec['gene_id'])
        if key in seen or key not in genes:raise ValueError('Automatic candidate identity invalid or duplicated')
        seen.add(key)
        old=prior.get(key);e=rec['evidence'];bounds=None
        if old and old['state']!='unassessed':
            state,reason=old['state'],'preserved_existing_assessed_decision'
        elif rec.get('assignment_status')=='out_of_scope':
            state,reason='out_of_scope','no_RM_assignment_under_registered_rules'
        elif (rec.get('assignment_status')!='supported' or e.get('family') not in RM_TYPES
              or e.get('role') not in RM_ROLES):
            state,reason='unassessed','unresolved_RM_type_or_component_assignment'
        elif (e.get('valid_sequence') is not True or not rec.get('protein_sha256')
              or genes[key].get('protein_sha256')!=rec['protein_sha256']):
            state,reason='unassessed','missing_or_conflicting_exact_protein_binding'
        else:
            bounds=policy_bounds(e,rec.get('unknown_axes',[]),p);evaluations+=1
            state='unassessed' if bounds['keep'] is None else 'retained' if bounds['keep'] else 'not_retained'
            reason=('unresolved_policy_bounds_cross_retention_threshold' if bounds['keep'] is None else
                    'automatic_'+bounds['assessment'])
        rec.update(state=state,reason=reason,decision=bounds,previous_saved_decision=old)
        records.append(rec)
        if old and old['state']!='unassessed':continue
        refs=rec.get('source_refs',[])
        if not refs:raise ValueError('Automatic candidate requires exact source references')
        row={**(old or {}),'sample_id':key[0],'gene_id':key[1],
             'protein_sha256':rec.get('protein_sha256'),'state':state,
             'family':e.get('family'),'rm_type':e.get('family') if e.get('family') in RM_TYPES else None,
             'role':e.get('role') if e.get('role') in RM_ROLES else None,
             'source_ref':refs[0],'source_refs':refs,'native_source_ref':(old or {}).get('source_ref'),
             'origin':'automatic_stored_evidence_evaluation','evaluation_id':evaluation_id,
             'source_scope':deepcopy(rec.get('source_scope',[])),
             'policy_sha256':policy_sha,'decision_reason':reason,
             'score_total':bounds['score_total'] if bounds else None,
             'score_min':bounds['score_min'] if bounds else None,
             'score_max':bounds['score_max'] if bounds else None,
             'unknown_evidence_axes':rec.get('unknown_axes',[])}
        for name in ('context_state','orphan','context_source_refs','integrated_R_component_support'):
            row.pop(name,None)
        if bounds and type(bounds['integrated_R_component_support']) is bool:
            row['integrated_R_component_support']=bounds['integrated_R_component_support']
        decisions[key]=row
    from .rm_evaluation import _candidate_integrity_issues
    compiled_scopes = {s['sample_id']: s for s in compiled.get('candidate_scope', [])}
    scopes=[]
    for sample in data['samples']:
        sid=sample['sample_id'];rows=[r for k,r in decisions.items() if k[0]==sid]
        # Exact bound-candidate inventory; unknown retention is explicit. Source
        # conflicts and quarantined identities stay separate in the compiler audit.
        states=[r['state'] for r in rows]
        issues = _candidate_integrity_issues(data, sid, prior)
        source_scope = compiled_scopes.get(sid, {})
        native_coverage=[c.get('status') for c in data.get('coverage',[]) if c.get('sample_id')==sid and c.get('endpoint') in {'saved_retention','saved_decisions'}]
        empty_completed = ('complete' in native_coverage or any(
            c.get('sample_id')==sid and c.get('native_candidates') is True and c.get('candidate_count')==0
            for c in data.get('suite',{}).get('capabilities',[])))
        complete=(bool(rows) or empty_completed) and 'unassessed' not in states and not issues and source_scope.get('status') == 'complete'
        coverage_status = 'complete' if complete else 'partial' if rows or issues else 'missing'
        if not rows and not empty_completed and 'not_run' in native_coverage:
            coverage_status='not_run' 
        scopes.append(dict(sample_id=sid,gene_ids=sorted(r['gene_id'] for r in rows),
            status='complete' if complete else 'partial',coverage_status=coverage_status,
            evaluated_count=sum(r['sample_id']==sid and r['decision'] is not None for r in records),
            reviewed_count=sum(r['sample_id']==sid for r in records),
            assessed_count=sum(r['state'] != 'unassessed' for r in rows),
            source_integrity_issues=issues, compiler_scope=source_scope,
            active_candidate_count=len(rows),native_candidate_count=sum(k[0]==sid for k in prior),
            remaining_unassessed_gene_ids=sorted(r['gene_id'] for r in rows if r['state']=='unassessed'),
            scope='all exact-bound saved nominees and native R-M members; quarantined source rows remain separate'))
    result['saved_decisions']=list(decisions.values())
    definite=sum(r['decision'] is not None and r['state']!='unassessed' for r in records)
    result['rm_evaluation']=dict(schema_version=1,operation=('automatic_evaluation_of_stored_RM_evidence' if definite else 'automatic_policy_bounds_with_unresolved_decisions' if evaluations else 'automatic_features_with_saved_decision_replay'),
        decision_mode='evaluate_unassessed',records=records,candidate_scope=scopes,
        nomination_scope='DNMB-native candidates plus original tool R-M members for integrated review; not an independent fourth detector',
        policy_evaluations=evaluations,policy_decisions=definite,scoring_runs=int(evaluations>0),detector_runs=0,
        preserved_decisions=sum(r.get('reason')=='preserved_existing_assessed_decision' for r in records),
        policy_scope='newly_evaluated_records_only; preserved decisions retain native policy provenance',
        structure_prediction_runs=0,thresholds_changed=False,selected_policy=parameters,
        policy=policy_metadata(p),policy_sha256=policy_sha,evaluation_id=evaluation_id,
        evaluation_identity_payload=identity_payload,
        native_saved_decisions=list(prior.values()),compiler=compiled,
        elapsed_seconds=time.monotonic()-started,activity_inferred=False)
    # Existing imported retention coverage described the saved export, not this
    # new computation. Preserve it in the audit and publish the new scoped status.
    result['rm_evaluation']['native_saved_retention_coverage']=[deepcopy(c) for c in result.get('coverage',[]) if c.get('endpoint') in {'saved_retention','saved_decisions'}]
    result['coverage']=[c for c in result.get('coverage',[]) if c.get('endpoint') not in {'saved_retention','saved_decisions'}]
    result['coverage'].extend(dict(sample_id=s['sample_id'],caller='DNMB',endpoint='saved_retention',status=s['coverage_status'],origin='automatic_stored_evidence_evaluation') for s in scopes)
    return refresh_decision_audit(result)

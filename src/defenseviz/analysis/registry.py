#!/usr/bin/env python3
"""Evidence-layer scanner; no activity inference, scoring, or sequence regex.

The caller supplies qualified, sequence-bound model hits and independently mapped
reference sites. This module validates the evidence contracts and keeps family,
domain, catalytic, module, and structure outcomes separate. It can consume new
genome-wide HMM runs without requiring a native final defense-system call.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path

ORIGINAL_SOURCE = "revision027/pipeline/nuclease_registry.py"
ORIGINAL_SOURCE_SHA256 = "38a155dd32dc685b4a62c266ebbc2ce59c7f0e18ffe09311fe5afed0b56bdc72"

STATES = {'supported', 'incomplete', 'unassessed', 'conflict', 'not_applicable'}
SUBSTRATES = {'DNA', 'RNA', 'dual', 'unknown'}


def finite_number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def valid_interval(start, end, length) -> bool:
    return (isinstance(start, int) and not isinstance(start, bool) and
            isinstance(end, int) and not isinstance(end, bool) and
            isinstance(length, int) and not isinstance(length, bool) and 1 <= start <= end <= length)


def result(state='unassessed', reason='', **details):
    assert state in STATES
    return dict(state=state, reason=reason, **details)


def load_reference(reference: dict, registry_dir: Path) -> str:
    if 'sequence' in reference:
        sequence = reference['sequence']
    else:
        path = Path(reference['sequence_file'])
        if not path.is_absolute(): path = registry_dir / path
        lines = path.read_text(encoding='utf-8').splitlines()
        if sum(x.startswith('>') for x in lines) > 1:
            raise ValueError('Reference must contain exactly one sequence')
        sequence = ''.join(x.strip() for x in lines if not x.startswith('>'))
    if sequence != sequence.upper() or not sequence or not sequence.isalpha():
        raise ValueError('Reference must be uppercase amino-acid sequence without gaps')
    if digest(sequence) != reference['sequence_sha256']:
        raise ValueError('Reference sequence hash mismatch')
    if reference.get('numbering_origin') != '1-based full reference sequence':
        raise ValueError('Explicit full-sequence reference numbering is required')
    positions = []
    for site in reference.get('catalytic_sites', []):
        pos = site['position']
        if not isinstance(pos, int) or isinstance(pos,bool) or not 1 <= pos <= len(sequence):
            raise ValueError('Reference catalytic position out of bounds')
        if sequence[pos - 1] != site['residue']:
            raise ValueError('Reference residue-numbering mismatch')
        if not site.get('evidence_source'):
            raise ValueError('Site must have primary evidence provenance')
        if not site.get('role'):
            raise ValueError('Site must distinguish nuclease from motor or modification')
        positions.append(pos)
    if len(positions) != len(set(positions)):
        raise ValueError('Duplicate reference site')
    return sequence


def _validate_registry(registry: dict, registry_dir: Path = Path('.')) -> None:
    if registry.get('schema_version') != '1.0': raise ValueError('Unknown registry schema')
    seen = set()
    for entry in registry['families']:
        fid = entry['family_id']
        if fid in seen: raise ValueError('Duplicate family identifier')
        seen.add(fid)
        if entry['substrate'] not in SUBSTRATES: raise ValueError('Unknown substrate')
        for model in entry.get('domain_models', []):
            if not (model.get('accession') or model.get('profile')):
                raise ValueError('Model needs an exact accession or profile')
            if model.get('evidence_kind') not in {'family', 'nuclease_domain', 'other_domain'}:
                raise ValueError('Model must separate family and domain evidence')
            if not model.get('source') or not model.get('role'):
                raise ValueError('Model source and biochemical role are required')
            if not model.get('model_sha256'):
                raise ValueError('Exact model checksum is required')
            if model.get('min_model_coverage') is not None and not (finite_number(model['min_model_coverage']) and 0 <= model['min_model_coverage'] <= 1):
                raise ValueError('Invalid model coverage contract')
        for ref in entry.get('references', []): load_reference(ref, registry_dir)
        rule=entry.get('candidate_rule',{}).get('reference_evidence',{})
        if rule.get('enabled') is True:
            for key in ('min_identity','min_query_coverage','min_reference_coverage'):
                if key in rule and not (finite_number(rule[key]) and 0<=rule[key]<=1):
                    raise ValueError('Invalid fraction in reference-evidence contract')
            if 'max_evalue' in rule and not (finite_number(rule['max_evalue']) and rule['max_evalue']>0):
                raise ValueError('Invalid E-value in reference-evidence contract')


def model_match(hit: dict, model: dict) -> bool:
    # Case and punctuation are meaningful. Alternate names must be exact aliases
    # listed in the registry; substring and normalized-name joins are forbidden.
    identities = []
    for field in ('accession', 'profile'):
        if model.get(field): identities.append(hit.get(field) == model[field])
    return (bool(identities) and all(identities) and hit.get('source') == model['source'] and
            hit.get('model_sha256') == model['model_sha256'])


def assess_domains(candidate: dict, entry: dict) -> tuple[dict, list, list]:
    seq = candidate['sequence']; sha = candidate['sequence_sha256']
    usable = []; issues = []; observed = []
    for hit in candidate.get('domain_hits', []):
        models = [m for m in entry.get('domain_models', []) if model_match(hit, m)]
        if not models: continue
        observed.append(hit)
        if hit.get('query_sequence_sha256') != sha:
            issues.append('domain_query_sequence_mismatch'); continue
        if not valid_interval(hit.get('qstart'), hit.get('qend'), len(seq)):
            issues.append('invalid_query_domain_coordinates'); continue
        interval_ok=valid_interval(hit.get('hmm_start'), hit.get('hmm_end'), hit.get('hmm_length'))
        native_family_only=(all(m['evidence_kind']=='family' for m in models) and hit.get('native_family_evidence') is True and
                            finite_number(hit.get('profile_coverage')) and 0<=hit['profile_coverage']<=1)
        if not interval_ok and not native_family_only:
            issues.append('invalid_model_domain_coordinates'); continue
        qualification = hit.get('qualification', {})
        if not (qualification.get('passed') is True and qualification.get('rule') and
                qualification.get('source_sha256') and qualification.get('source')):
            issues.append('domain_qualification_unassessed'); continue
        for model in models:
            coverage = (hit['hmm_end'] - hit['hmm_start'] + 1) / hit['hmm_length'] if interval_ok else hit['profile_coverage']
            required = model.get('min_model_coverage')
            if required is not None and coverage < required:
                issues.append('domain_coverage_incomplete'); continue
            # Methylation and motor evidence cannot provide a nuclease domain.
            if model['evidence_kind'] == 'nuclease_domain' and model['role'] != 'nuclease':
                issues.append('non_nuclease_role_cannot_transfer'); continue
            usable.append(dict(hit=hit, model=model))
    if usable:
        family_qualified=any(u['model']['evidence_kind']=='family' or u['model'].get('family_specific') is True for u in usable)
        state = result('supported' if family_qualified else 'unassessed', 'exact_family_model_evidence' if family_qualified else 'registered_domain_observed_family_unresolved', usable_hits=len(usable), issues=sorted(set(issues)))
    elif observed:
        state = result('incomplete' if any('coverage' in x for x in issues) else 'unassessed',
                       'matching_hit_failed_evidence_contract', issues=sorted(set(issues)))
    else:
        state = result('unassessed', 'no_exact_registered_model_hit')
    return state, usable, issues


def assess_reference_evidence(candidate: dict, entry: dict, ref: dict, registry_dir: Path) -> dict:
    """Separate source-bound reference route. Never fabricates Pfam/domain support."""
    rule = entry.get('candidate_rule', {}).get('reference_evidence', {})
    if rule.get('enabled') is not True:
        return result('unassessed','reference_route_not_registered',reference=ref['accession'])
    reference = load_reference(ref, registry_dir)
    segments = ref.get('catalytic_segments', [])
    sites = [s for s in ref.get('catalytic_sites',[]) if s['role']=='nuclease' and s.get('required_for_reference_subset',True) is not False]
    if not segments or not all(valid_interval(x.get('start'),x.get('end'),len(reference)) for x in segments):
        return result('unassessed','reference_catalytic_segments_not_registered',reference=ref['accession'])
    if not all(any(x['start']<=s['position']<=x['end'] for x in segments) for s in sites):
        return result('conflict','reference_sites_outside_registered_segments',reference=ref['accession'])
    if candidate['sequence_sha256'] == ref['sequence_sha256']:
        return result('supported','exact_full_reference_sequence',reference=ref['accession'],reference_segments=segments,qstart=1,qend=len(candidate['sequence']))
    names = ('min_identity','min_query_coverage','min_reference_coverage','max_evalue')
    if not all(finite_number(rule.get(x)) for x in names) or not rule.get('source') or not rule.get('source_sha256'):
        return result('unassessed','reference_homology_contract_not_pinned',reference=ref['accession'])
    records = [r for r in candidate.get('reference_evidence',[]) if r.get('reference_accession')==ref['accession'] and
               r.get('reference_sequence_sha256')==ref['sequence_sha256'] and r.get('query_sequence_sha256')==candidate['sequence_sha256']]
    for record in records:
        if record.get('method')!='BLASTP' or not record.get('source') or not record.get('source_sha256'): continue
        qf,qt,rf,rt=(record.get(x) for x in ('qstart','qend','rstart','rend'))
        if not valid_interval(qf,qt,len(candidate['sequence'])) or not valid_interval(rf,rt,len(reference)): continue
        ids,alen,ev=(record.get(x) for x in ('identities','alignment_length','evalue'))
        if not all(finite_number(x) for x in (ids,alen,ev)) or not isinstance(ids,int) or not isinstance(alen,int) or alen<=0 or not 0<=ids<=alen or ev<0: continue
        if alen<max(qt-qf+1,rt-rf+1):continue
        identity=ids/alen; qc=(qt-qf+1)/len(candidate['sequence']);rc=(rt-rf+1)/len(reference)
        if not (identity>=rule['min_identity'] and qc>=rule['min_query_coverage'] and rc>=rule['min_reference_coverage'] and ev<=rule['max_evalue']): continue
        if not all(rf<=segment['start']<=segment['end']<=rt for segment in segments): continue
        return result('supported','source_pinned_family_reference_homology',reference=ref['accession'],
                      identity=identity,query_coverage=qc,reference_coverage=rc,evalue=ev,reference_segments=segments,qstart=qf,qend=qt)
    return result('unassessed','no_qualified_reference_alignment',reference=ref['accession'])


def assess_reference(candidate: dict, entry: dict, ref: dict, usable: list, registry_dir: Path, *, reference_route=False) -> dict:
    reference = load_reference(ref, registry_dir)
    allsites=[s for s in ref.get('catalytic_sites',[]) if s['role']=='nuclease']
    sites = [s for s in allsites if s.get('required_for_reference_subset',True) is not False]
    diagnostic_sites=[s for s in allsites if s.get('required_for_reference_subset',True) is False]
    if not sites: return result('unassessed', 'no_registered_nuclease_sites', reference=ref['accession'])
    nucdomains = [u['hit'] for u in usable if u['model']['evidence_kind'] == 'nuclease_domain']
    if reference_route:
        eligibility=assess_reference_evidence(candidate,entry,ref,registry_dir)
        if eligibility['state']!='supported':return eligibility
    elif not nucdomains: return result('unassessed', 'no_qualified_nuclease_domain', reference=ref['accession'])
    elif not any(u['model'].get('family_specific') is True for u in usable if u['model']['evidence_kind']=='nuclease_domain'):
        eligibility=assess_reference_evidence(candidate,entry,ref,registry_dir)
        if eligibility['state']!='supported':
            return result('unassessed','generic_domain_does_not_establish_reference_subfamily',reference=ref['accession'])
    required_methods = ref.get('independent_mapping_methods', [])
    if len(set(required_methods)) < 2:
        return result('unassessed', 'independent_mapping_contract_missing', reference=ref['accession'])
    maps = candidate.get('site_maps', [])
    matched = [m for m in maps if m.get('reference_accession') == ref['accession'] and
               m.get('reference_sequence_sha256') == digest(reference) and
               m.get('query_sequence_sha256') == candidate['sequence_sha256']]
    by_method = {}
    for method in required_methods:
        ms = [m for m in matched if m.get('method') == method and m.get('source_sha256') and m.get('source')]
        if len(ms) != 1:
            return result('unassessed', 'missing_or_duplicate_independent_mapping', reference=ref['accession'], method=method)
        by_method[method] = ms[0]
    mapped = []; failures = []
    for site in sites:
        positions = [m.get('positions', {}).get(str(site['position'])) for m in by_method.values()]
        if any(p is None for p in positions):
            failures.append('required_site_unmapped'); continue
        if len(set(positions)) != 1:
            failures.append('independent_site_mapping_disagreement'); continue
        p = positions[0]
        if not valid_interval(p, p, len(candidate['sequence'])):
            failures.append('query_site_coordinate_invalid'); continue
        if reference_route and not eligibility['qstart']<=p<=eligibility['qend']:
            failures.append('query_site_outside_qualified_reference_alignment'); continue
        if not reference_route and not any(d['qstart'] <= p <= d['qend'] for d in nucdomains):
            failures.append('site_outside_qualified_nuclease_domain'); continue
        allowed = site.get('allowed_residues') or [site['residue']]
        aa = candidate['sequence'][p - 1]
        if aa not in allowed:
            failures.append('required_catalytic_residue_not_conserved')
        mapped.append(dict(reference_position=site['position'], query_position=p,
                           expected=allowed, observed=aa, conserved=aa in allowed))
    if len({r['query_position'] for r in mapped}) != len(mapped):
        failures.append('distinct_reference_sites_collapsed_to_one_query_residue')
    if mapped and not reference_route:
        if ref.get('domain_binding','single')=='single' and not any(all(d['qstart']<=s['query_position']<=d['qend'] for s in mapped) for d in nucdomains):
            failures.append('required_sites_stitched_across_distinct_domain_instances')
        elif ref.get('domain_binding','single')!='single':
            # Explicit split-domain implementations need per-segment, model-bound
            # assignment. This prototype refuses to substitute an arbitrary union.
            failures.append('split_domain_binding_requires_registered_adapter')
    # Site order is a contract, not a universal biological premise; documented
    # circular permutations can set site_order="permuted_supported" explicitly.
    order = [r['query_position'] for r in sorted(mapped, key=lambda r:r['reference_position'])]
    if ref.get('site_order', 'collinear') == 'collinear' and order != sorted(order):
        failures.append('site_mapping_order_conflict')
    diagnostics=[]
    for site in diagnostic_sites:
        positions=[m.get('positions',{}).get(str(site['position'])) for m in by_method.values()]
        good=all(isinstance(p,int) and not isinstance(p,bool) and 1<=p<=len(candidate['sequence']) for p in positions) and len(set(positions))==1
        observed=candidate['sequence'][positions[0]-1] if good else None
        diagnostics.append(dict(reference_position=site['position'],mapped_positions=positions,observed=observed,
                                reference_residue=site['residue'],reference_like=observed in (site.get('allowed_residues') or [site['residue']]) if observed else None,
                                required_for_reference_subset=False))
    if failures:
        state = 'incomplete' if all(x == 'required_site_unmapped' for x in failures) else 'conflict'
        return result(state, 'reference_site_contract_not_met', reference=ref['accession'],
                      failures=sorted(set(failures)), sites=mapped,diagnostic_sites=diagnostics,
                      site_set_completeness=ref.get('site_set_completeness','unspecified_tested_subset'),
                      catalytic_assembly_context=ref.get('catalytic_assembly_context','unassessed'))
    return result('supported', 'qualified_domain_and_independent_reference_site_conservation',
                  reference=ref['accession'], sites=mapped, evidence_route='reference_evidence' if reference_route else 'qualified_domain',
                  site_set_completeness=ref.get('site_set_completeness','unspecified_tested_subset'),
                  catalytic_assembly_context=ref.get('catalytic_assembly_context','unassessed'),diagnostic_sites=diagnostics)


def assess_structure(candidate: dict, entry: dict, supported_ref: dict | None) -> dict:
    if not supported_ref: return result('unassessed', 'no_supported_catalytic_site_mapping')
    model = candidate.get('structural_model')
    if not model: return result('unassessed', 'model_not_available')
    if model.get('role') != 'nuclease': return result('unassessed', 'non_nuclease_model_role')
    if model.get('sequence') != candidate['sequence'] or model.get('sequence_sha256') != candidate['sequence_sha256']:
        return result('conflict', 'exact_full_model_sequence_mismatch')
    if not model.get('source_sha256') or not model.get('source'):
        return result('unassessed', 'model_provenance_missing')
    residues = model.get('residues', {})
    if set(residues) != {str(p) for p in range(1, len(candidate['sequence'])+1)}:
        return result('unassessed', 'model_full_residue_coverage_incomplete')
    for p, row in residues.items():
        if row.get('aa') != candidate['sequence'][int(p)-1]:
            return result('conflict', 'model_residue_identity_mismatch')
        xyz = row.get('CA', [])
        if len(xyz) != 3 or not all(finite_number(x) for x in xyz):
            return result('unassessed', 'missing_or_nonfinite_CA_coordinates')
        if not finite_number(row.get('occupancy')) or row['occupancy'] <= 0:
            return result('unassessed', 'invalid_model_occupancy')
        value = row.get('confidence')
        if not finite_number(value) or not 0 <= value <= 100:
            return result('unassessed', 'invalid_confidence_range')
    if model.get('confidence_kind') != 'pLDDT':
        return result('unassessed', 'confidence_semantics_not_pLDDT')
    site_confidence = [dict(reference_position=s['reference_position'],query_position=s['query_position'],
                            plddt=residues[str(s['query_position'])]['confidence']) for s in supported_ref['sites']]
    vals = [s['plddt'] for s in site_confidence]
    threshold = entry.get('candidate_rule', {}).get('local_plddt_threshold')
    if not finite_number(threshold) or not 0 <= threshold <= 100:
        return result('unassessed', 'reference_local_confidence_contract_not_registered', site_plddt=vals,
                      required_site_confidence=site_confidence)
    mean = sum(vals) / len(vals)
    return result('supported' if mean >= threshold else 'incomplete',
                  'exact_model_catalytic_site_local_confidence', mean_site_plddt=mean, threshold=threshold,
                  required_site_confidence=site_confidence,
                  geometry_state='unassessed', geometry_reason='No template or active-site geometry calculation is inferred from local confidence')


def assess_module(candidate: dict, entry: dict) -> dict:
    requirement = entry.get('candidate_rule', {}).get('module_components')
    if not requirement: return result('not_applicable', 'component_or_domain_registry_entry')
    modules = [m for m in candidate.get('modules', []) if m.get('family_id') == entry['family_id']]
    for m in modules:
        if not m.get('native_system_id') or not m.get('source_sha256') or not m.get('source'):
            continue
        if m.get('sample_id') != candidate.get('sample_id'): continue
        members = m.get('members', [])
        if not members or any(not x.get('protein_sha256') or not x.get('query_id') for x in members): continue
        if not any(x['query_id'] == candidate.get('query_id') and x['protein_sha256'] == candidate['sequence_sha256'] for x in members):
            continue
        present = {x['component'] for x in members}
        if set(requirement) <= present:
            return result('supported', 'complete_native_module_literature_evidence', native_system_id=m['native_system_id'],
                          catalytic_assignment='unassessed_not_transferred_to_members')
    return result('incomplete' if modules else 'unassessed', 'complete_native_module_contract_not_met')


def _assess_candidate(candidate: dict, registry: dict, registry_dir: Path = Path('.')) -> dict:
    seq = candidate['sequence']
    if not seq or seq != seq.upper() or not seq.isalpha() or digest(seq) != candidate['sequence_sha256']:
        raise ValueError('Candidate sequence must be exact, uppercase, ungapped, and hash-bound')
    out = dict(candidate_id=candidate['candidate_id'], sequence_sha256=digest(seq), families=[],
               activity='not_assessed', scoring_changed=False)
    for entry in registry['families']:
        fam, usable, issues = assess_domains(candidate, entry)
        module = assess_module(candidate, entry)
        refelig = [assess_reference_evidence(candidate,entry,ref,registry_dir) for ref in entry.get('references',[])]
        if fam['reason'] == 'no_exact_registered_model_hit' and module['state'] not in {'supported','incomplete'} and not any(r['state']=='supported' for r in refelig):
            continue
        if fam['state']!='supported' and any(r['state']=='supported' for r in refelig):
            fam=result('supported','source_pinned_family_reference_evidence_not_shared_domain_alone',domain_observation=fam)
        nuc = [u for u in usable if u['model']['evidence_kind'] == 'nuclease_domain']
        domain = result('supported', 'qualified_nuclease_domain') if nuc else result('unassessed', 'family_or_other_domain_is_not_nuclease_domain')
        assessments = [assess_reference(candidate, entry, ref, usable, registry_dir) for ref in entry.get('references', [])]
        ref_assessments = [assess_reference(candidate,entry,ref,usable,registry_dir,reference_route=True) for ref in entry.get('references',[])]
        supported = next((x for x in assessments if x['state'] == 'supported'), None)
        if supported:
            cat = result('supported', 'at_least_one_registered_reference_contract_met', references=assessments)
        elif assessments:
            cat = result('conflict' if any(x['state'] == 'conflict' for x in assessments) else
                         'incomplete' if any(x['state'] == 'incomplete' for x in assessments) else 'unassessed',
                         'no_registered_reference_contract_met', references=assessments)
        else: cat = result('unassessed', 'no_registered_catalytic_reference')
        refsite=result('supported' if any(r['state']=='supported' for r in ref_assessments) else
                       'conflict' if any(r['state']=='conflict' for r in ref_assessments) else
                       'incomplete' if any(r['state']=='incomplete' for r in ref_assessments) else 'unassessed',
                       'independent_reference_site_route_without_domain_transfer',references=ref_assessments)
        dna_state = cat['state'] if entry['substrate'] in {'DNA','dual'} else 'not_applicable' if entry['substrate']=='RNA' else 'unassessed'
        structure_by_reference=[dict(reference=r['reference'],assessment=assess_structure(candidate,entry,r)) for r in assessments if r['state']=='supported']
        if len(structure_by_reference)==1:
            structure=structure_by_reference[0]['assessment']
        elif structure_by_reference:
            structure=result('unassessed','multiple_supported_references_no_implicit_best_selection',references=structure_by_reference)
        else:structure=assess_structure(candidate,entry,None)
        out['families'].append(dict(family_id=entry['family_id'], component=entry['component'], substrate=entry['substrate'],
                                    family_evidence=fam, domain_evidence=domain, catalytic_site_evidence=cat,
                                    reference_evidence=refelig,reference_site_evidence=refsite,
                                    DNA_catalytic_site_state=dna_state, module_evidence=module,
                                    structure_evidence=structure,
                                    biological_activity='not_assessed'))
    if not out['families']:
        out['disposition'] = 'unassessed_no_exact_registered_family_or_module_evidence'
    return out


AMINO_ACIDS = frozenset('ACDEFGHIKLMNPQRSTVWYBXZJUO')


def _sha256(value, where):
    if not isinstance(value, str) or len(value) != 64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError(f'{where} must be a lowercase SHA-256 digest')


def validate_registry(registry: dict, registry_dir: Path = Path('.')) -> None:
    """Validate a user-supplied registry; no biological catalogue is hard coded."""
    try:
        if not isinstance(registry, dict) or not isinstance(registry.get('families'), list):
            raise ValueError('Registry must contain a families list')
        for entry in registry['families']:
            if not isinstance(entry.get('family_id'), str) or not entry['family_id'].strip():
                raise ValueError('Registry family_id must be nonempty')
            if not isinstance(entry.get('component'), str) or not entry['component'].strip():
                raise ValueError('Registry component must be explicit')
            names = set()
            for ref in entry.get('references', []):
                if not isinstance(ref.get('accession'), str) or not ref['accession']:
                    raise ValueError('Reference accession must be explicit')
                if ref['accession'] in names:
                    raise ValueError('Duplicate reference accession within family')
                names.add(ref['accession'])
                _sha256(ref.get('sequence_sha256'), 'reference.sequence_sha256')
                seq = load_reference(ref, Path(registry_dir))
                if set(seq) - AMINO_ACIDS:
                    raise ValueError('Reference contains unsupported amino-acid symbols')
                if ref.get('site_order', 'collinear') not in {'collinear', 'permuted_supported'}:
                    raise ValueError('Unknown reference site_order')
                methods = ref.get('independent_mapping_methods', [])
                if not isinstance(methods, list) or any(not isinstance(m, str) or not m for m in methods) or len(methods) != len(set(methods)):
                    raise ValueError('Independent mapping methods must be unique nonempty strings')
                for site in ref.get('catalytic_sites', []):
                    residue = site.get('residue')
                    if not isinstance(residue, str) or len(residue) != 1 or residue not in AMINO_ACIDS:
                        raise ValueError('Reference site residue must be one uppercase amino acid')
                    allowed = site.get('allowed_residues', [residue])
                    if not isinstance(allowed, list) or not allowed or any(not isinstance(a, str) or len(a) != 1 or a not in AMINO_ACIDS for a in allowed):
                        raise ValueError('Allowed reference residues must be an explicit nonempty list of amino acids')
                    if 'required_for_reference_subset' in site and type(site['required_for_reference_subset']) is not bool:
                        raise ValueError('required_for_reference_subset must be boolean')
            for model in entry.get('domain_models', []):
                _sha256(model.get('model_sha256'), 'model.model_sha256')
        _validate_registry(registry, Path(registry_dir))
    except (KeyError, TypeError, AttributeError, OSError) as exc:
        raise ValueError(f'Malformed registry: {exc}') from exc


def _verify_file_binding(source, expected, root, inputs):
    if not isinstance(source, str) or not source:
        raise ValueError('Evidence source must identify a local file')
    _sha256(expected, 'evidence source_sha256')
    if '://' in source:
        raise ValueError('Evidence verification is offline; use a local source file')
    path = Path(source)
    if not path.is_absolute():
        path = root / path
    path = path.resolve()
    if not path.is_file():
        raise ValueError(f'Evidence source file not found: {path}')
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != expected:
        raise ValueError(f'Evidence source hash mismatch: {path}')
    inputs[str(path)] = observed


def assess_candidate(candidate: dict, registry: dict, registry_dir: Path = Path('.'), *,
                     source_root=None, verify_sources=True) -> dict:
    """Assess supplied evidence, checking actual source files by default.

    This evaluates domain/site/module assertions; it does not run HMMER,
    alignments or a structure model. Explicit verify_sources=False is available
    for trusted in-memory fixtures and is prominently recorded in the result.
    """
    validate_registry(registry, registry_dir)
    inputs = {}
    root = Path(source_root) if source_root is not None else Path(registry_dir)
    # Validation reads every external reference, even if no family later matches.
    # Record raw file hashes separately from the validated biological sequence hash.
    for entry in registry['families']:
        for ref in entry.get('references', []):
            if 'sequence' not in ref:
                path = Path(ref['sequence_file'])
                if not path.is_absolute():
                    path = Path(registry_dir) / path
                path = path.resolve()
                inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    try:
        if not isinstance(candidate, dict):
            raise ValueError('Candidate must be an object')
        for k in ('candidate_id', 'sample_id', 'query_id'):
            if not isinstance(candidate.get(k), str) or not candidate[k].strip():
                raise ValueError(f'candidate.{k} must be a nonempty string')
        seq = candidate['sequence']
        if not isinstance(seq, str) or not seq or set(seq) - AMINO_ACIDS:
            raise ValueError('Candidate sequence must be uppercase ungapped amino acids')
        _sha256(candidate.get('sequence_sha256'), 'candidate.sequence_sha256')
        if digest(seq) != candidate['sequence_sha256']:
            raise ValueError('Candidate sequence hash mismatch')
        bound_sources = []
        for key in ('domain_hits', 'reference_evidence', 'site_maps', 'modules'):
            if not isinstance(candidate.get(key, []), list):
                raise ValueError(f'candidate.{key} must be a list')
            if any(not isinstance(r, dict) for r in candidate.get(key, [])):
                raise ValueError(f'candidate.{key} records must be objects')
        for hit in candidate.get('domain_hits', []):
            q = hit.get('qualification', {})
            if q.get('passed') is True:
                bound_sources.append(q)
            if 'target_length' in hit and (type(hit['target_length']) is not int or hit['target_length'] != len(seq)):
                raise ValueError('Domain target_length differs from exact candidate sequence')
        for key in ('reference_evidence', 'site_maps', 'modules'):
            bound_sources.extend(candidate.get(key, []))
        for module in candidate.get('modules', []):
            identities = set()
            for member in module.get('members', []):
                query = member.get('query_id')
                if not isinstance(query, str) or not query:
                    raise ValueError('Module member query_id must be explicit')
                if query in identities:
                    raise ValueError('One module member cannot supply multiple components through duplicate identity')
                identities.add(query)
                _sha256(member.get('protein_sha256'), 'module member protein_sha256')
            if not isinstance(module.get('caller'), str) or not module['caller'].strip():
                raise ValueError('Module requires an explicit native caller, separate from its system identity')
        model = candidate.get('structural_model')
        if model:
            bound_sources.append(model)
        if verify_sources:
            for item in bound_sources:
                if item.get('source') or item.get('source_sha256'):
                    _verify_file_binding(item.get('source'), item.get('source_sha256'), root, inputs)
        out = _assess_candidate(candidate, registry, Path(registry_dir))
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f'Malformed registry candidate: {exc}') from exc
    out.update(sample_id=candidate['sample_id'], query_id=candidate['query_id'],
               source_files_verified=bool(verify_sources), source_sha256=inputs,
               analysis_kind='evaluation_of_supplied_evidence_not_detector_run')
    for family in out['families']:
        module = family['module_evidence']
        if module['state'] == 'supported':
            matching = [m for m in candidate.get('modules', []) if m.get('family_id') == family['family_id'] and m.get('native_system_id') == module['native_system_id'] and m.get('sample_id') == candidate['sample_id']]
            callers = {m['caller'] for m in matching}
            if len(callers) == 1:
                module['native_caller'] = next(iter(callers))
            else:
                family['module_evidence'] = result('unassessed', 'ambiguous_caller_scoped_native_system_identity')
    return out


def assess_registry(registry: dict, candidates: list, *, registry_dir=Path('.'),
                    source_root=None, verify_sources=True) -> list:
    """Batch public API. Candidate identity is a tuple, never a split string."""
    validate_registry(registry, registry_dir)
    if not isinstance(candidates, list):
        raise ValueError('Candidates must be a list')
    seen = set()
    out = []
    for c in candidates:
        if not isinstance(c, dict):
            raise ValueError('Candidate must be an object')
        key = (c.get('sample_id'), c.get('query_id'))
        try:
            if key in seen:
                raise ValueError(f'Duplicate candidate sample/query identity: {key!r}')
            seen.add(key)
        except TypeError as exc:
            raise ValueError('Candidate sample_id/query_id must be strings') from exc
        out.append(assess_candidate(c, registry, registry_dir, source_root=source_root, verify_sources=verify_sources))
    return out

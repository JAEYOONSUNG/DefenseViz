"""Apply the unchanged selected DNMB policy to completed native evidence stages.

Fresh search receipts and exact proteins establish which computations ran. The
native score inputs are explicit booleans, not optimistic completions of unknown
axes. Retention remains an annotation decision, never biological activity.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from dataclasses import asdict
import csv
import hashlib
import json
from pathlib import Path
import re

from .analysis.evidence import parse_evidence, policy_metadata
from .analysis.policy import Policy, Profile, evaluate, from_feature, numeric, truth
from .automatic import evaluation_identity, refresh_decision_audit
from .automatic_context import apply_current_context_review
from .automatic_neighborhood import cognate_subjects, divergent_geometry

_MISSING={'','NA','NaN','nan','None','null','-'}
_FACTORS=('strong','moderate','profile','annotation','motif','context','structure','alternative_annotation','conflict','partial')
_TYPES={'Type I','Type II','Type III','Type IV'}
_ROLES={'R','M','S','RM','Aux','A','C','F'}
# from_feature consumes these native columns. Sequence length, strict optional
# structure credit, and HMM ambiguity are computed and bound explicitly here.
_DEPLOYED_NATIVE_COLUMNS={
    'family_id','enzyme_role','blast_evalue','blast_family_compatible','blast_identity',
    'blast_query_coverage','blast_reference_coverage','blast_role_compatible',
    'curated_blast_evalue','curated_blast_family_compatible','curated_blast_identity',
    'curated_blast_query_coverage','curated_blast_reference_coverage','curated_blast_role_compatible',
    'mtase_motif_verified','rease_motif_hits_raw','rease_motor_verified','rease_nuclease_verified',
    *('factor_'+name for name in _FACTORS)}


def _hash(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def _table(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as f:
        r=csv.DictReader(f,delimiter='\t')
        if not r.fieldnames or len(r.fieldnames)!=len(set(r.fieldnames)):
            raise ValueError(f'{path}: table needs unique column names')
        rows=list(r)
        if any(None in r or None in r.values() for r in rows):
            raise ValueError(f'{path}: malformed row/header widths')
        return rows


def _bool(value, where):
    if str(value).strip().lower() not in {'true','false','t','f','1','0','1.0','0.0'}:
        raise ValueError(f'{where}: an executed boolean observation is required')
    return truth(value)


def _union(intervals):
    if not intervals:return 0
    a,b=sorted(intervals)[0];total=0
    for c,d in sorted(intervals)[1:]:
        if c>b+1:total+=b-a+1;a,b=c,d
        else:b=max(b,d)
    return total+b-a+1


def _fasta(path):
    sequences={};name=None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in sequences:raise ValueError('Duplicate full-search query identity')
                sequences[name]=''
            elif line.strip():
                if name is None:raise ValueError('Invalid full-search FASTA')
                sequences[name]+=''.join(line.split()).upper()
    return {name:(hashlib.sha256(sequence.encode()).hexdigest(),len(sequence)) for name,sequence in sequences.items()}


def _compatible(family, role, type_hints, role_hints):
    types=set(re.split(r'[;,|+/ ]+',type_hints))-_MISSING
    roles=set(re.split(r'[;,|+/ ]+',role_hints))-_MISSING
    native_type=family.replace('Type ','')
    return (native_type in types or native_type=='II' and bool(types&{'IIG','IIM'}),
            role in roles or role=='RM' and {'R','M'}<=roles)


def _profiles(rows, gene, hypothesis, metadata, library_sha):
    groups=defaultdict(list)
    for line,row in rows:
        meta=metadata.get(row['profile_name'])
        if not meta or row['role_hints']!=meta['roles'] or row['type_hints']!=meta['type_hints'] or row['database_sha256']!=meta['library_sha256'] or row['database_sha256']!=library_sha:
            raise ValueError('Native HMM profile metadata/library binding mismatch')
        if row['protein_sha256']!=gene.get('protein_sha256'):
            raise ValueError('Native HMM domain protein hash mismatch')
        if not truth(row['resolved_selected']) or truth(row['off_target']) or row['mapping_status']!='unique':continue
        qlen=int(float(row['sequence_length']));mlen=int(float(row['model_length']))
        if qlen!=gene.get('protein_length') or mlen<=0 or mlen!=int(float(meta['profile_length'])):
            raise ValueError('Native HMM domain length binding mismatch')
        coords=[int(float(row[k])) for k in ('ali_from','ali_to','hmm_from','hmm_to')]
        if not (1<=coords[0]<=coords[1]<=qlen and 1<=coords[2]<=coords[3]<=mlen):
            raise ValueError('Native HMM interval outside its protein/profile')
        groups[(row['profile_name'],mlen)].append((line,row,coords))
    profiles=[];details=[]
    for (name,length),members in sorted(groups.items()):
        fam,role=_compatible(hypothesis['family_id'],hypothesis['enzyme_role'],members[0][1]['type_hints'],members[0][1]['role_hints'])
        query=_union([x[2][:2] for x in members])/gene['protein_length']
        model=_union([x[2][2:] for x in members])/length
        ev=min(float(x[1]['sequence_evalue']) for x in members)
        profile=Profile(fam,role,ev,query,model)
        profiles.append(profile)
        details.append(dict(profile_name=name,source_rows=[x[0] for x in members],**asdict(profile)))
    return profiles,details


def _fresh_hmm_role_ambiguity(rows):
    """Native nomination summary rule, restricted to fresh DNMB domain rows.

    This intentionally has no family/E-value/coverage gate: the frozen native
    nomination summary uses all resolved, non-offtarget, uniquely mapped rows.
    With no native nomination its flag is absent, interpreted as False by the
    original feature adapter. Comparator nomination hints never enter here.
    """
    roles=set()
    for _,row in rows:
        if not truth(row['resolved_selected']) or truth(row['off_target']) or row['mapping_status']!='unique':continue
        hints=row['role_hints'].split(';') if row['role_hints'] else ['unknown']
        if not set(hints)<={'R','M','S','unknown'}:raise ValueError('Invalid fresh native HMM role hint')
        roles.update(hints)
    return bool(roles) and (len(roles-{'unknown'})!=1 or 'unknown' in roles),sorted(roles)


def _divergent(rows, genes):
    positions={};by_contig=defaultdict(list)
    for gid,g in genes.items():
        if all(numeric(g.get(k)) is not None for k in ('start','end')) and g.get('contig') not in _MISSING:
            by_contig[g['contig']].append(g)
    for local in by_contig.values():
        local.sort(key=lambda g:(float(g['start']),float(g['end']),g['locus_tag']))
        for i,g in enumerate(local):positions[g['locus_tag']]=(i,g)
    def eligible(row,role):
        nums=[numeric(row.get('curated_blast_'+s)) for s in ('identity','query_coverage','reference_coverage','evalue')]
        return (row['family_id']=='Type II' and row['enzyme_role']==role and
            row.get('curated_blast_role')==role and row.get('curated_blast_family')=='Type II' and
            truth(row.get('curated_blast_role_compatible')) and truth(row.get('curated_blast_family_compatible')) and
            row.get('partial_status')=='full_length' and not truth(row.get('evidence_conflict')) and
            not truth(row.get('alternative_function_annotation')) and all(n is not None for n in nums) and
            nums[0]>=.5 and min(nums[1:3])>=.8 and nums[3]<=1e-5)
    supported=set();pairs=[]
    for r in [r for r in rows if eligible(r,'R')]:
        for m in [m for m in rows if eligible(m,'M')]:
            if r['query'] not in positions or m['query'] not in positions:raise ValueError('Eligible divergent pair lacks complete coordinates')
            ri,rg=positions[r['query']];mi,mg=positions[m['query']]
            left,right=(rg,mg) if ri<mi else (mg,rg)
            if not divergent_geometry(left,right,abs(ri-mi)-1):continue
            rs=r.get('curated_blast_resolved_subject_id');ms=m.get('curated_blast_resolved_subject_id')
            if not cognate_subjects(rs,ms):continue
            supported.update((r['query'],m['query']))
            pairs.append(dict(r_gene_id=r['query'],m_gene_id=m['query'],r_subject=rs,m_subject=ms,
                              gap_bp=int(right['start'])-int(left['end'])-1,intervening_genes=abs(ri-mi)-1,
                              branch='divergent_cognate_neighborhood_not_same_operon',activity_inferred=False))
    return supported,pairs


def apply_full_evaluation(data, stages):
    """Consume completed stages (sample_id, native_dir, hmm_dir) without mutations.

    Each native directory requires a source-bound stage manifest or completion
    artifact hashes. Original nominees outside the fresh nomination inventory
    are archived as out-of-scope; this is not a biological negative decision.
    """
    result=dict(data);result['suite']=dict(data.get('suite',{}))
    result['derived_artifact_sha256']=dict(data.get('derived_artifact_sha256',{}))
    prior={(r['sample_id'],r['gene_id']):deepcopy(r) for r in data.get('saved_decisions',[])}
    genes={(g['sample_id'],g['gene_id']):g for g in data.get('genes',[])}
    if len(prior)!=len(data.get('saved_decisions',[])) or len(genes)!=len(data.get('genes',[])):
        raise ValueError('Duplicate canonical gene or native decision identity')
    selected={s['sample_id'] for s in stages}
    if len(selected)!=len(stages) or selected!={s['sample_id'] for s in data['samples']}:
        raise ValueError('Completed native stages must cover each selected sample exactly once')
    decisions={};records=[];scopes=[];groups=[];stage_audit=[];fresh_features=[];pairs=[];effectors=[];omitted_contexts={}
    checks={};policy=Policy();parameters=asdict(policy)
    policy_sha=hashlib.sha256(json.dumps(parameters,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    def checked(path,expected=None):
        path=Path(path).resolve();digest=_hash(path)
        if expected is not None and digest!=expected:raise ValueError(f'{path}: completed-stage artifact hash mismatch')
        previous=checks.get(str(path))
        if previous and previous!=digest:raise ValueError(f'{path}: source changed during policy evaluation')
        checks[str(path)]=digest;result['derived_artifact_sha256'][str(path)]=digest
        return dict(path=str(path),sha256=digest)
    for stage in stages:
        sid=stage['sample_id'];native=Path(stage['native_dir']).resolve();hmm=Path(stage['hmm_dir']).resolve()
        cref=checked(native/'completion.json');completion=json.loads((native/'completion.json').read_text())
        if completion.get('passed') is not True or completion.get('full_reference_search') is not True or completion.get('selected_policy_applied') is not False:
            raise ValueError('Native stage is not a completed full-reference evidence computation')
        manifest_path=native/'stage_manifest.json'
        if not manifest_path.exists() and (native.parent/'stage_manifest.json').exists():manifest_path=native.parent/'stage_manifest.json'
        manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else completion
        mref=checked(manifest_path) if manifest_path.exists() else cref
        if manifest.get('exit_code',0)!=0:raise ValueError('Native stage manifest reports failed execution')
        if manifest.get('completion_sha256') and manifest['completion_sha256']!=cref['sha256']:raise ValueError('Completion receipt hash mismatch')
        artifact_hashes=manifest.get('artifacts',manifest.get('artifact_hashes',{}))
        if not artifact_hashes:raise ValueError('Completed native stage needs artifact hashes')
        def artifact(name):
            expected=artifact_hashes.get(name,artifact_hashes.get('native/'+name))
            if isinstance(expected,dict):expected=expected.get('sha256')
            if not expected:raise ValueError(f'Native stage did not bind required artifact {name}')
            return checked(native/name,expected)
        if manifest is not completion:
            expected_completion=manifest.get('completion_sha256') or artifact_hashes.get('completion.json')
            if isinstance(expected_completion,dict):expected_completion=expected_completion.get('sha256')
            if expected_completion!=cref['sha256']:raise ValueError('Completion receipt is not bound by native manifest')
        receipt_path=native.parent/'receipt.json'
        receipt=json.loads(receipt_path.read_text()) if receipt_path.exists() else None
        receipt_ref=None
        if receipt is not None:
            receipt_ref=checked(receipt_path)
            if receipt.get('sample_id')!=sid or receipt.get('exit_code')!=0 or receipt.get('stage_manifest_sha256')!=mref['sha256']:
                raise ValueError('Outer native execution receipt does not bind this completed sample/stage')
            for path,digest in receipt.get('source_inputs',{}).items():checked(path,digest)
            for path,ref in receipt.get('reference_files',{}).items():
                checked(path,ref['sha256']);checked(ref['private_copy'],ref['copy_sha256'])
                if ref['sha256']!=ref['copy_sha256']:raise ValueError('Private reference copy differs from pinned input')
            command=receipt.get('command',[])
            runner=[p for p in command if str(p).endswith('.R')]
            gbff=[p for p in command if str(p).endswith(('.gbff','.gbk','.gb'))]
            if receipt.get('runner_sha256'):
                if len(runner)!=1:raise ValueError('Native runner identity is ambiguous')
                checked(runner[0],receipt['runner_sha256'])
            if receipt.get('gbff_sha256'):
                if len(gbff)!=1:raise ValueError('Native genome identity is ambiguous')
                checked(gbff[0],receipt['gbff_sha256'])
            for filename,field in [('annotations.tsv','input_annotations_sha256'),('hmm_nominees.tsv','input_nominees_sha256')]:
                if manifest.get(field):checked(native.parent/'inputs'/filename,manifest[field])
            if manifest.get('original_inputs_unchanged') is False:raise ValueError('Native inputs changed during the execution')
        elif not manifest.get('input_sha256'):
            raise ValueError('Native stage lacks an input-bound execution receipt')
        href=artifact('native_hypotheses.tsv');gref=artifact('canonical_genes.tsv');lref=artifact('search_ledger.tsv');cmdref=artifact('search_commands.json')
        for path,digest in manifest.get('input_sha256',{}).items():checked(path,digest.get('sha256') if isinstance(digest,dict) else digest)
        binding_path=native.parent/'inputs'/'protein_binding.json'
        protein_binding={}
        if binding_path.is_file():
            expected_binding=manifest.get('input_sha256',{}).get(str(binding_path.resolve()))
            if isinstance(expected_binding,dict):expected_binding=expected_binding.get('sha256')
            if not expected_binding:raise ValueError('Protein normalization receipt is not bound by native stage manifest')
            checked(binding_path,expected_binding);protein_binding=json.loads(binding_path.read_text())
            if protein_binding.get('sample_id')!=sid or protein_binding.get('normalization')!='uppercase whitespace-free; terminal stop removed for HMM only':
                raise ValueError('Unsupported protein normalization identity contract')
            checked(protein_binding['genbank'],protein_binding['genbank_sha256'])
        for path,digest in completion.get('source_hashes',{}).items():checked(path,digest)
        hmm_receipt_path=(hmm/'execution.json').resolve()
        expected_hmm_receipt=manifest.get('input_sha256',{}).get(str(hmm_receipt_path))
        if isinstance(expected_hmm_receipt,dict):expected_hmm_receipt=expected_hmm_receipt.get('sha256')
        if not expected_hmm_receipt:raise ValueError('Native stage manifest does not bind the completed HMM execution receipt')
        hs=checked(hmm_receipt_path,expected_hmm_receipt);hr=json.loads(hmm_receipt_path.read_text())
        if hr.get('exit_code')!=0 or hr.get('inputs_unchanged') is not True or hr.get('hmm',{}).get('exit_code')!=0 or hr.get('resolution',{}).get('exit_code')!=0:
            raise ValueError('Native HMM search/resolution did not complete')
        signature=hr.get('signature',{})
        if signature.get('sample_id')!=sid:raise ValueError('HMM sample identity mismatch')
        for record in signature.get('inputs',{}).values():checked(record['path'],record['sha256'])
        if not {'fasta','gene_map'}<=set(signature.get('inputs',{})):
            raise ValueError('Completed HMM search needs its exact query inventory')
        hmm_sequences=_fasta(signature['inputs']['fasta']['path'])
        hmm_gene_map={}
        for mapping in _table(signature['inputs']['gene_map']['path']):
            if mapping.get('sample_id')!=sid or mapping['locus_tag'] in hmm_gene_map:
                raise ValueError('HMM query inventory has conflicting sample/gene identity')
            expected=(mapping['protein_sha256'],int(mapping['aa_length']))
            if hmm_sequences.get(mapping['query_id'])!=expected:
                raise ValueError('HMM query map does not bind its actual protein FASTA')
            hmm_gene_map[mapping['locus_tag']]=expected
        dref=checked(hmm/'all_domains.tsv',hr['artifact_hashes']['all_domains.tsv'])
        metadata_record=signature['inputs']['profiles'];pref=checked(metadata_record['path'],metadata_record['sha256'])
        metadata={r['profile_name']:r for r in _table(metadata_record['path'])}
        canonical=_table(native/'canonical_genes.tsv');cg={r['locus_tag']:r for r in canonical}
        if len(cg)!=len(canonical):raise ValueError('Duplicate canonical native gene')
        complete_proteome={}
        for row in canonical:
            sequence=''.join(row.get('translation','').split()).upper().rstrip('*')
            if sequence not in {'','NA'}:
                complete_proteome[row['locus_tag']]=(hashlib.sha256(sequence.encode()).hexdigest(),len(sequence))
        if complete_proteome!=hmm_gene_map:
            raise ValueError('Completed HMM query inventory is not the entire exact canonical translated proteome')
        hypotheses=_table(native/'native_hypotheses.tsv');hx={r['query']:r for r in hypotheses}
        if len(hx)!=len(hypotheses) or len(hypotheses)!=completion.get('final_candidates'):
            raise ValueError('Native candidate count/identity mismatch')
        domains=defaultdict(list)
        for n,row in enumerate(_table(hmm/'all_domains.tsv'),2):
            if row.get('sample_id')!=sid:raise ValueError('HMM domain sample mismatch')
            domains[row['locus_tag']].append((n,row))
        searched={r['locus_tag'] for r in _table(native/'search_ledger.tsv') if _bool(r['completed'],'search ledger')}
        commands=json.loads((native/'search_commands.json').read_text())
        if searched and not commands:raise ValueError('Native search ledger lacks BLAST command receipts')
        blast_sequences={}
        for command in commands:
            argv=command['argv']
            if command['exit_code']!=0 or '-db' not in argv or '-subject' in argv:raise ValueError('Full native BLAST search did not complete')
            checked(argv[argv.index('-query')+1],command['query_sha256'])
            checked(argv[argv.index('-out')+1],command['output_sha256'])
            for query,value in _fasta(argv[argv.index('-query')+1]).items():
                if query in blast_sequences and blast_sequences[query]!=value:raise ValueError('BLAST query identity differs across completed stages')
                blast_sequences[query]=value
        for row in hypotheses:
            missing_features=_DEPLOYED_NATIVE_COLUMNS-set(row)
            if missing_features:raise ValueError('Native deployment feature columns missing: '+', '.join(sorted(missing_features)))
            gid=row['query'];gene=genes.get((sid,gid));native_gene=cg.get(gid)
            if gene is None or native_gene is None:raise ValueError('Fresh native nominee is not a canonical gene')
            sequence=''.join(native_gene.get('translation','').split()).upper()
            if sequence=='NA':sequence=''
            digest=hashlib.sha256(sequence.encode()).hexdigest() if sequence else None
            clean_sequence=sequence.rstrip('*')
            clean_digest=hashlib.sha256(clean_sequence.encode()).hexdigest() if clean_sequence else None
            if digest!=gene.get('protein_sha256'):
                if not (clean_sequence and protein_binding.get('source_protein_sha256',{}).get(gid)==gene.get('protein_sha256') and
                        hmm_gene_map.get(gid)==(clean_digest,len(clean_sequence))):
                    raise ValueError(f'{sid}/{gid}: fresh native protein hash mismatch')
            if sequence and gid not in searched:raise ValueError(f'{sid}/{gid}: translated native nominee lacks completed BLAST search')
            if sequence and blast_sequences.get(gid) not in {(digest,len(sequence)),(clean_digest,len(clean_sequence))}:
                raise ValueError(f'{sid}/{gid}: native nominee absent from exact completed BLAST FASTA')
            if sequence and hmm_gene_map.get(gid)!=(clean_digest,len(clean_sequence)):
                raise ValueError(f'{sid}/{gid}: native nominee absent from exact completed HMM FASTA')
            for name in _FACTORS:_bool(row.get('factor_'+name),'factor_'+name)
            row['family_id']=row.get('family_id') if row.get('family_id') not in _MISSING else 'Unresolved'
            row['enzyme_role']=row.get('enzyme_role') if row.get('enzyme_role') not in _MISSING else 'unknown'
            row['sequence_length_aa']=len(clean_sequence)
            row['source_protein_sha256']=gene.get('protein_sha256')
            row['canonical_protein_sha256']=digest
            row['hmm_protein_sha256']=clean_digest
        support,local_pairs=_divergent(hypotheses,cg)
        pairs.extend(dict(sample_id=sid,source_refs=[href,gref],**p) for p in local_pairs)
        row_numbers={r['query']:n for n,r in enumerate(hypotheses,2)}
        for row in hypotheses:
            gid=row['query'];gene=genes[sid,gid];sequence_length=row['sequence_length_aa']
            hmm_gene=dict(gene,protein_sha256=row['hmm_protein_sha256'],protein_length=sequence_length)
            ps,profile_details=_profiles(domains[gid],hmm_gene,row,metadata,signature['inputs']['hmm']['sha256']) if sequence_length else ([],[])
            union_ambiguity=row.get('hmm_role_ambiguous')
            native_ambiguity,native_roles=_fresh_hmm_role_ambiguity(domains[gid]) if sequence_length else (False,[])
            row['hmm_role_ambiguous_nomination_union']=union_ambiguity
            row['hmm_role_ambiguous']=native_ambiguity
            native_factor={k:_bool(row['factor_'+k],k) for k in _FACTORS}
            native_total=max(0,max(35*native_factor['strong'],25*native_factor['moderate'],25*native_factor['profile'],10*native_factor['annotation'])+
                20*native_factor['motif']+20*native_factor['context']+25*native_factor['structure']-
                15*native_factor['alternative_annotation']-15*(native_factor['conflict'] and not native_factor['alternative_annotation'])-10*native_factor['partial'])
            if native_total!=numeric(row.get('score_total')):raise ValueError(f'{sid}/{gid}: captured native factor/score mismatch')
            if native_factor['structure'] and not truth(row.get('numeric_structure_admissible_structural_evidence')):
                raise ValueError('Native structure credit lacks an exact numerical geometry assessment')
            # No optional structure model is asserted when the full core did not
            # execute that stage. This is zero *credit*, not measured absence.
            row['numeric_structure_admissible_structural_evidence']=native_factor['structure']
            ev=parse_evidence(asdict(from_feature(row,ps,divergent=gid in support)))
            decision=evaluate(ev,policy)
            decision.update(score_min=decision['score_total'],score_max=decision['score_total'],
                            assessment='exact_policy_evaluation',unknown_axes=[],hypothetical_evaluations=0)
            if not ev.valid_sequence and decision['keep']:raise ValueError('Untranslated candidate cannot be retained')
            state='retained' if decision['keep'] else 'not_retained'
            sources=[{**href,'row':row_numbers[gid]},gref,lref,cmdref,dref,hs,pref,mref]
            record=dict(sample_id=sid,gene_id=gid,protein_sha256=gene.get('protein_sha256'),evidence=asdict(ev),
                decision=decision,state=state,reason='completed_native_evidence_selected_policy',unknown_axes=[],assignment_status='evaluated',
                source_refs=sources,source_scope=['current_full_native_nomination'],previous_saved_decision=deepcopy(prior.get((sid,gid))),
                profile_details=profile_details,optional_structure_status=row.get('score_structure_state','not_assessed'),
                native_factor_score=native_total,genome_search_scope='complete_native_nominees_not_genomewide_blast')
            record['translation_status']='translated' if sequence_length else 'no_translated_protein'
            record['hmm_role_ambiguity']=dict(nomination_union_flag=union_ambiguity,scoring_native_flag=native_ambiguity,
                native_role_hints=native_roles,source_refs=[dref,hs,pref],scope='fresh_native_resolved_nonofftarget_unique_domains_only')
            record['protein_binding']=dict(source_sha256=gene.get('protein_sha256'),canonical_sha256=row['canonical_protein_sha256'],
                hmm_sha256=row['hmm_protein_sha256'],normalization='terminal_stop_removal_only',
                source_receipt=protein_binding.get('normalization') if gene.get('protein_sha256')!=row['canonical_protein_sha256'] else 'exact_canonical_translation')
            record['evidence']['raw_nuclease_classes']=sorted(record['evidence']['raw_nuclease_classes'])
            records.append(record)
            if decision['keep'] and ev.role in {'R','RM'} and decision['integrated_R_component_support'] is True:
                effectors.append(dict(sample_id=sid,gene_id=gid,protein_sha256=gene['protein_sha256'],
                    rule_id='full_native_selected_policy_nuclease_R',assignment_kind='individual_nuclease_effector',
                    status='supported_putative_effector_annotation',evidence_grade='selected_policy_nuclease_R_component_support',
                    family=ev.family,role=ev.role,policy_sha256=policy_sha,policy_version=policy.version,
                    source_refs=deepcopy(sources),origin='full_native_evidence_evaluation',automatic=True,
                    motif_state='unassessed',activity='unassessed',activity_inferred=False,
                    support_predicate='retained_and_R_or_RM_and_integrated_R_component_support',
                    validated_motif_site_inferred=False))
            decisions[sid,gid]=dict(sample_id=sid,gene_id=gid,protein_sha256=gene.get('protein_sha256'),state=state,
                family=ev.family,rm_type=ev.family if ev.family in _TYPES else None,role=ev.role if ev.role in _ROLES else None,
                native_evidence_role=ev.role,score_total=decision['score_total'],score_min=decision['score_total'],score_max=decision['score_total'],
                source_ref=sources[0],source_refs=sources,origin='full_native_evidence_evaluation',decision_reason=record['reason'],
                integrated_R_component_support=decision['integrated_R_component_support'],policy_sha256=policy_sha,
                unknown_evidence_axes=[],source_scope=record['source_scope'])
            fresh_features.append(dict(sample_id=sid,gene_id=gid,protein_sha256=gene.get('protein_sha256'),source_ref=sources[0],native_fields=dict(row)))
        for row in hypotheses:
            gid=row['query']
            for prefix,typ,required in [('typei','Type I',['R','M','S']),('typeii','Type II',['R','M']),('typeiii','Type III',['R','M'])]:
                text=row.get(prefix+'_context_partners','')
                no_partners=not text or text in _MISSING
                # The native Type III scan explicitly exports an M-only
                # boundary with no partner text. Missing partner metadata by
                # itself is not a completed singleton observation.
                singleton_mod=(prefix=='typeiii' and row['family_id']==typ and row['enzyme_role']=='M' and
                    row.get('typeiii_context_status')=='missing_res_check_neighbors' and
                    row.get('typeiii_context_roles')=='M' and row.get('typeiii_context_anchor')==gid)
                if no_partners and not singleton_mod:continue
                members={gid}
                for partner in ([] if no_partners else text.split(' | ')):
                    choices=[g for g in cg if partner==g or partner.startswith(g+':')]
                    if len(choices)!=1:raise ValueError(f'{sid}/{gid}: native context partner identity unresolved')
                    members.add(choices[0])
                if any(m not in hx for m in members):raise ValueError(f'{sid}/{gid}: native context partner not evaluated')
                group_key=[sid,typ,sorted(members)]
                group_id='full_native_context:'+hashlib.sha256(json.dumps(group_key,sort_keys=True).encode()).hexdigest()
                missing_protein=sorted(m for m in members if not genes[sid,m].get('protein_sha256'))
                if missing_protein:
                    omitted=omitted_contexts.setdefault(group_id,dict(group_id=group_id,sample_id=sid,family=typ,
                        member_keys=[[sid,m] for m in sorted(members)],
                        translated_member_keys=[[sid,m] for m in sorted(members) if m not in missing_protein],
                        untranslated_member_keys=[[sid,m] for m in missing_protein],
                        native_context_state=row.get(prefix+'_context_status'),required_roles=required,
                        reason_code='untranslated_context_member',source_refs=[gref],orphan_inferred=False))
                    ref={**href,'row':row_numbers[gid]}
                    if ref not in omitted['source_refs']:omitted['source_refs'].append(ref)
                    continue
                groups.append(dict(group_id=group_id,sample_id=sid,provider='DNMB',run_id=href['sha256'],membership_kind='saved_context',
                    family=typ,context_state=row.get(prefix+'_context_status'),required_roles=required,
                    members=[dict(gene_id=m,protein_sha256=genes[sid,m]['protein_sha256']) for m in sorted(members)],
                    source_refs=[{**href,'row':row_numbers[gid]}],origin='fresh_native_context_scan'))
        extra=sorted(g for s,g in prior if s==sid and g not in hx)
        for gid in extra:
            old=prior[sid,gid]
            decisions[sid,gid]=dict(old,state='out_of_scope',origin='fresh_nomination_scope_reconciliation',
                decision_reason='outside_current_native_nomination_inventory_not_a_biological_negative',source_ref=href,
                source_refs=[href,gref,lref],source_scope=['outside_current_full_native_nomination'])
            for field in ('context_state','orphan','integrated_R_component_support','score_total','score_min','score_max'):
                decisions[sid,gid].pop(field,None)
        scopes.append(dict(sample_id=sid,gene_ids=sorted(hx),status='complete',coverage_status='complete',
            active_candidate_count=len(hx),native_candidate_count=sum(s==sid for s,g in prior),evaluated_count=len(hx),reviewed_count=len(hx),
            assessed_count=len(hx),remaining_unassessed_gene_ids=[],scope='current_full_native_nomination',
            translated_count=sum(bool(r['sequence_length_aa']) for r in hypotheses),
            untranslated_gene_ids=sorted(r['query'] for r in hypotheses if not r['sequence_length_aa']),
            outside_current_nomination_gene_ids=extra,source_refs=[href,mref]))
        stage_audit.append(dict(sample_id=sid,native_stage_source=cref,hmm_stage_source=hs,manifest_source=mref,execution_receipt_source=receipt_ref,
            native_run_kind=completion.get('run_kind'),new_blast_commands=sum(not c.get('reuse_scope') for c in commands),
            reused_exact_blast_commands=sum(bool(c.get('reuse_scope')) for c in commands),
            canonical_genes=len(canonical),complete_hmm_proteins=len(complete_proteome),
            full_reference_search=True,genomewide_blast=False,final_candidates=len(hx),translated_searched=sum(bool(r['sequence_length_aa']) for r in hypotheses),
            optional_structure_status='not_executed_by_full_core'))
    compiled=dict(records=[{k:r[k] for k in ('sample_id','gene_id','protein_sha256','evidence')} for r in records],selected_policy=parameters)
    eval_id,identity=evaluation_identity(compiled)
    for row in decisions.values():row['evaluation_id']=eval_id
    result['saved_decisions']=list(decisions.values())
    old_generated_effectors=[e for e in data.get('effector_evidence',[]) if e.get('origin')=='full_native_evidence_evaluation']
    result['effector_evidence']=[e for e in data.get('effector_evidence',[]) if e not in old_generated_effectors]+effectors
    for effector in effectors:effector['evaluation_id']=eval_id
    old_groups=[g for g in data.get('groups',[]) if g.get('provider') in {'DNMB','Review'} and g.get('membership_kind')=='saved_context']
    unique_groups={g['group_id']:g for g in groups}
    for g in groups:
        existing=unique_groups[g['group_id']]
        for ref in g['source_refs']:
            if ref not in existing['source_refs']:existing['source_refs'].append(ref)
    result['groups']=[g for g in data.get('groups',[]) if g not in old_groups]+list(unique_groups.values())
    result['suite']['archived_native_context_groups']=deepcopy(old_groups)
    result['suite']['archived_native_features']=data.get('suite',{}).get('features',[])
    result['suite']['features']=fresh_features
    original_coverage=[c for c in data.get('coverage',[]) if c.get('endpoint') in {'saved_retention','saved_decisions'}]
    result['coverage']=[c for c in data.get('coverage',[]) if c not in original_coverage]
    result['coverage'].extend(dict(sample_id=s['sample_id'],caller='DNMB',endpoint='saved_retention',status='complete',origin='full_native_evidence_evaluation') for s in scopes)
    result['rm_evaluation']=dict(schema_version=1,operation='fresh_native_evidence_selected_policy_evaluation',decision_mode='evaluate_current_full_nomination',
        records=records,candidate_scope=scopes,policy_evaluations=len(records),policy_decisions=len(records),scoring_runs=int(bool(records)),
        detector_runs=0,homology_search_runs=0,hmm_search_runs=0,structure_prediction_runs=0,thresholds_changed=False,
        completed_native_stages_consumed=len(stages),completed_hmm_stages_consumed=len(stages),
        selected_policy=parameters,policy=policy_metadata(policy),policy_sha256=policy_sha,evaluation_id=eval_id,evaluation_identity_payload=identity,
        native_saved_decisions=list(prior.values()),native_saved_retention_coverage=deepcopy(original_coverage),
        generated_policy_effector_evidence=len(effectors),archived_generated_effector_evidence=deepcopy(old_generated_effectors),
        archived_native_context_groups=deepcopy(old_groups),full_evidence_stages=stage_audit,divergent_context_pairs=pairs,
        preserved_decisions=0,policy_scope='all current completed native nominees; prior-only records outside current nomination',activity_inferred=False)
    for path,digest in checks.items():
        if _hash(path)!=digest:raise ValueError(f'{path}: source changed during full evaluation')
    result=apply_current_context_review(result)
    result['rm_evaluation']['current_context_review']['omitted_context_groups']=list(omitted_contexts.values())
    partial_context={(r['sample_id'],r['family']) for r in omitted_contexts.values()}
    result['rm_evaluation']['current_context_review']['observed_inventory_coverage']=[dict(sample_id=s['sample_id'],family=typ,
        status='partial' if (s['sample_id'],typ) in partial_context else 'complete',
        reason_codes=['untranslated_context_member'] if (s['sample_id'],typ) in partial_context else [],
        omitted_group_ids=sorted(r['group_id'] for r in omitted_contexts.values() if r['sample_id']==s['sample_id'] and r['family']==typ),
        scope='current_retained_groups_from_completed_context_scan',source_refs=s['source_refs']) for s in scopes for typ in ('Type I','Type III')]
    return refresh_decision_audit(result)

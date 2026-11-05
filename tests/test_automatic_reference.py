"""Fresh reference measurements cannot borrow stale metadata or mutate calls."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

import pytest

from defenseviz import automatic_reference as module

SEQ="MAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQAMQKIQHELGVPQQMQGRVGEAALSKRWEK"


def fixture(tmp_path, *, subject="R.Example", extra=None):
    cache=tmp_path/'cache';cache.mkdir()
    (cache/'rebase_db.fasta').write_text('>'+subject+'\n'+SEQ+'\n>irrelevant\n'+SEQ*2+'\n')
    metadata=dict(subject_id=subject,sequence_sha256=sha256(SEQ.encode()).hexdigest(),reference_family="Type II",reference_role="R")
    rows=[metadata]+(extra or [])
    (cache/'reference_metadata.json').write_text(json.dumps(dict(schema_version=1,records=rows)))
    digest=sha256(SEQ.encode()).hexdigest();ref=dict(path="original.tsv",sha256="a"*64,row=2)
    data=dict(samples=[dict(sample_id="sample")], genes=[dict(sample_id="sample",gene_id="gene",protein_sha256=digest,protein_length=len(SEQ))],
        calls=[dict(sample_id="sample",gene_id="gene",caller="PADLOC",family="RM")],suite=dict(
            reference_matches=[dict(sample_id="sample",gene_id="gene",reference_id=subject,source_ref=ref)],
            annotations=[dict(sample_id="sample",gene_id="gene",source_ref=ref,native_fields=dict(translation=SEQ))]))
    return data,cache


def mock_blast(monkeypatch, *, no_hit=False, bad_length=False):
    def run(command, commands, *, stdout_path=None):
        commands.append(dict(command=command,exit_code=0,elapsed_seconds=0))
        if "-version" in command:return "synthetic BLAST 1"
        outfile=Path(command[command.index('-out')+1])
        values=['q0','r0','100',str(len(SEQ)),'0','0','1',str(len(SEQ)),'1',str(len(SEQ)), '1e-25','120',str(len(SEQ)+(1 if bad_length else 0)),str(len(SEQ))]
        outfile.write_text('' if no_hit else '\t'.join(values)+'\n')
        if stdout_path is not None:Path(stdout_path).write_text('')
        return ''
    monkeypatch.setattr(module,'_run',run)


def test_exact_proteins_current_metadata_single_alignment_and_no_native_mutation(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path);before=deepcopy(data);mock_blast(monkeypatch)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    assert data==before and result['calls']==before['calls']
    row=result['suite']['reference_alignments'][0]
    assert row['native_fields']['identity']==1 and row['native_fields']['reference_length']==len(SEQ)
    assert row['protein_sha256']==sha256(SEQ.encode()).hexdigest()
    assert row['reference_sequence_sha256']==sha256(SEQ.encode()).hexdigest()
    assert row['operation']=='fresh_reference_alignment' and row['command_exit_code']==0
    receipt=result['suite']['automatic_reference']
    assert receipt['detector_runs']==0 and receipt['homology_search_runs']==1
    assert receipt['reference_residues']==len(SEQ)*3
    command=receipt['commands'][-1]['command']
    assert command[command.index('-dbsize')+1]==str(len(SEQ)*3)
    assert '-max_hsps' not in command
    assert row['search_scope']=='saved_nominated_subjects' and row['search_is_exhaustive'] is False
    assert set(result['inputs'])=={str(cache/'rebase_db.fasta'),str(cache/'reference_metadata.json')}


def test_completed_no_alignment_is_explicit_not_missing_measurement(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path);mock_blast(monkeypatch,no_hit=True)
    row=module.enrich_reference_evidence(data,cache,tmp_path/'work')['suite']['reference_alignments'][0]
    assert row['alignment_status']=='no_alignment' and 'identity' not in row['native_fields']
    assert row['command_exit_code']==0 and row['reference_metadata_source_refs']


@pytest.mark.parametrize('text,expected',[('putative Type II restriction enzyme','Type II'),('putative Type III restriction enzyme','Type III'),
    ('Type IIG restriction enzyme/methyltransferase','Type II'),('Type IV restriction enzyme','Type IV'),('Type I restriction enzyme','Type I')])
def test_explicit_reference_description_has_boundary_aware_type(text,expected):
    value=module._metadata_assignment(dict(enz_type=text,rm_type='Type_II'))
    assert value['reference_family']==expected
    assert value['family_metadata_disagreement']==(expected!='Type II')


def test_reference_assignment_unknown_and_ambiguity_not_inferred_from_name():
    assert module._metadata_assignment(dict(enzyme_name='M.Example'))['reference_role'] is None
    value=module._metadata_assignment(dict(enz_type='Type I / Type III restriction enzyme',rm_type='Type_II'))
    assert value['reference_family'] is None and value['ambiguous_description']
    assert module._metadata_assignment(dict(enz_type='Type IIG restriction enzyme/methyltransferase'))['reference_role']=='RM'


def test_duplicate_reference_name_is_unresolved_even_if_one_is_closer(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path,extra=[dict(subject_id='R.Example',sequence=SEQ[:-1],reference_family='Type II',reference_role='R')])
    mock_blast(monkeypatch)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    assert result['suite']['reference_alignments']==[]
    receipt=result['suite']['automatic_reference'];assert receipt['homology_search_runs']==0
    assert receipt['unresolved'][0]['reason']=='unresolved_current_reference'
    mapping=json.loads((tmp_path/'work/reference_mapping.json').read_text())
    assert mapping[0]['reason']=='ambiguous_reference_name'


def test_fasta_metadata_sequence_mismatch_cannot_supply_subject_length(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path);metadata=json.loads((cache/'reference_metadata.json').read_text())
    metadata['records'][0]['sequence_sha256']='b'*64
    (cache/'reference_metadata.json').write_text(json.dumps(metadata));mock_blast(monkeypatch)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    assert result['suite']['reference_alignments']==[]
    assert json.loads((tmp_path/'work/reference_mapping.json').read_text())[0]['reason']=='reference_metadata_fasta_sequence_mismatch'


def test_query_hash_conflict_cannot_receive_current_reference_evidence(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path);data['suite']['annotations'][0]['native_fields']['translation']=SEQ+'A';mock_blast(monkeypatch)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    assert not result['suite']['reference_alignments']
    assert result['suite']['automatic_reference']['unresolved'][0]['reason']=='missing_or_conflicting_exact_protein_binding'


def test_context_labels_are_never_reference_proteins(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path);data['suite']['reference_matches'][0]['reference_id']='typeI_context:R:gene';mock_blast(monkeypatch)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    assert not result['suite']['reference_alignments']


def test_blast_reported_protein_lengths_must_match_actual_bound_inputs(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path);mock_blast(monkeypatch,bad_length=True)
    with pytest.raises(ValueError,match='length binding mismatch'):
        module.enrich_reference_evidence(data,cache,tmp_path/'work')
    assert json.loads((tmp_path/'work/reference_evaluation.json').read_text())['status']=='failed'


def test_reference_cache_cannot_be_an_output_directory(tmp_path):
    data,cache=fixture(tmp_path)
    with pytest.raises(ValueError,match='read-only cache'):
        module.enrich_reference_evidence(data,cache,cache/'work')


def test_only_terminal_rebase_separator_is_normalized_with_receipt(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path)
    (cache/'rebase_db.fasta').write_text('>R.Example\n'+SEQ+'<>\n')
    meta=dict(subject_id='R.Example',sequence=SEQ+'<>',seq_length=len(SEQ)+2,reference_family='Type II',reference_role='R')
    (cache/'reference_metadata.json').write_text(json.dumps(dict(schema_version=1,records=[meta])))
    mock_blast(monkeypatch)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    row=result['suite']['reference_alignments'][0]
    assert row['native_fields']['reference_length']==len(SEQ)
    mapping=json.loads((tmp_path/'work/reference_mapping.json').read_text())[0]
    assert mapping['sequence_normalizations']==['terminal_REBASE_record_separator_removed']
    assert mapping['reference_source_sequence_sha256']==sha256((SEQ+'<>').encode()).hexdigest()
    assert result['suite']['automatic_reference']['reference_residues']==len(SEQ)
    assert module._reference_sequence(SEQ[:10]+'<>'+SEQ[10:])[0] is None


def test_terminal_stop_query_preserves_raw_binding_and_records_normalization(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path)
    raw=SEQ+'*';raw_digest=sha256(raw.encode()).hexdigest()
    data['genes'][0].update(protein_sha256=raw_digest,protein_length=len(raw))
    data['suite']['annotations'][0]['native_fields']['translation']=raw
    mock_blast(monkeypatch)
    before=deepcopy(data)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    row=result['suite']['reference_alignments'][0]
    assert data==before
    assert row['protein_sha256']==row['query_source_sequence_sha256']==raw_digest
    assert row['query_sequence_sha256']==sha256(SEQ.encode()).hexdigest()
    assert row['query_sequence_length']==len(SEQ)
    assert row['query_sequence_normalizations']==['native_nonletter_removal']
    query=tmp_path/'work/queries.faa'
    assert query.read_text()=='>q0\n'+SEQ+'\n'
    assert row['query_sequence_source_ref']['sha256']==sha256(query.read_bytes()).hexdigest()
    assert result['suite']['automatic_reference']['aligned_pairs']==1


def test_normalized_query_hash_cannot_substitute_for_mismatched_raw_source(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path)
    # Raw source contains a terminal marker, but the canonical hash does not.
    data['suite']['annotations'][0]['native_fields']['translation']=SEQ+'*'
    mock_blast(monkeypatch)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    assert result['suite']['reference_alignments']==[]
    assert result['suite']['automatic_reference']['homology_search_runs']==0
    assert result['suite']['automatic_reference']['unresolved'][0]['reason']=='missing_or_conflicting_exact_protein_binding'


@pytest.mark.parametrize('text', [
    'anti-restriction protein', 'antirestriction protein', 'non-restriction protein',
    'not a restriction enzyme', 'lacks restriction activity',
])
def test_negative_restriction_description_cannot_assign_catalytic_role(text):
    value=module._metadata_assignment(dict(enz_type=text,subunit='R'))
    assert value['reference_role'] is None
    assert value['role_description_conflict'] is True
    assert value['role_metadata_disagreement'] is True


@pytest.mark.parametrize('text', [
    'Type II restriction-modification system', 'Type IIG restriction-modification system',
    'restriction-related protein', 'DNA restriction protein',
])
def test_generic_system_description_does_not_assign_a_restriction_component(text):
    assert module._metadata_assignment(dict(enz_type=text))['reference_role'] is None


@pytest.mark.parametrize('text,expected', [
    ('Type II restriction enzyme', 'R'), ('Type IV restriction endonuclease', 'R'),
    ('Type I restriction-modification system subunit R', 'R'),
    ('Type IIG restriction enzyme/methyltransferase', 'RM'),
    ('Type I specificity subunit', 'S'), ('Type III methyltransferase', 'M'),
])
def test_explicit_component_descriptions_preserve_reference_role(text,expected):
    assert module._metadata_assignment(dict(enz_type=text))['reference_role']==expected


def test_hsp_merge_uses_union_coverage_and_native_length_weighted_identity():
    a=dict(pident=90.,length=30.,bitscore=100.,evalue=1e-20,qstart=1.,qend=30.,sstart=1.,send=30.)
    b=dict(pident=50.,length=20.,bitscore=60.,evalue=1e-10,qstart=25.,qend=44.,sstart=41.,send=60.)
    c=dict(pident=70.,length=5.,bitscore=10.,evalue=1e-3,qstart=60.,qend=64.,sstart=61.,send=65.)
    merged=module._collapse_hsps([a,b,c])
    assert merged['dnmb_aligned_query_length']==49  # overlap once; gap excluded
    assert merged['dnmb_aligned_subject_length']==55
    assert merged['length']==49 and merged['dnmb_hsp_count']==3
    assert merged['pident']==pytest.approx((90*30+50*20+70*5)/55)
    assert merged['bitscore']==170 and merged['evalue']==1e-20
    assert (merged['qstart'],merged['qend'])==(1,64)


def test_native_incomplete_hsp_intervals_fall_back_to_best_bitscore_hsp():
    a=dict(pident=90.,length=30.,bitscore=100.,evalue=1e-20,qstart=1.,qend=30.,sstart=1.,send=30.)
    b=dict(pident=50.,length=20.,bitscore=60.,evalue=1e-10,qstart=None,qend=None,sstart=41.,send=60.)
    merged=module._collapse_hsps([a,b])
    assert merged['pident']==90 and merged['length']==30
    assert merged['dnmb_aligned_query_length']==30 and merged['dnmb_hsp_count']==1


def test_all_pair_hsps_are_kept_and_aggregated_with_line_provenance(tmp_path,monkeypatch):
    data,cache=fixture(tmp_path)
    def run(command,commands,*,stdout_path=None):
        commands.append(dict(command=command,exit_code=0,elapsed_seconds=0))
        if '-version' in command:return 'synthetic BLAST 1'
        outfile=Path(command[command.index('-out')+1])
        out=[]
        for start,end,pident,bits in [(1,30,90,100),(41,len(SEQ),80,80)]:
            out.append('\t'.join(map(str,['q0','r0',pident,end-start+1,0,0,start,end,start,end,1e-30,bits,len(SEQ),len(SEQ)])))
        outfile.write_text('\n'.join(out)+'\n')
        if stdout_path is not None:Path(stdout_path).write_text('')
        return ''
    monkeypatch.setattr(module,'_run',run)
    result=module.enrich_reference_evidence(data,cache,tmp_path/'work')
    row=result['suite']['reference_alignments'][0]; n=row['native_fields']
    assert n['dnmb_hsp_count']==2
    assert n['query_coverage']==n['reference_coverage']==pytest.approx((len(SEQ)-10)/len(SEQ))
    assert n['identity']==pytest.approx((90*30+80*(len(SEQ)-40))/(len(SEQ)-10)/100)
    assert row['source_ref']['rows']==[1,2]
    assert len((tmp_path/'work/reference_alignments.tsv').read_text().splitlines())==2

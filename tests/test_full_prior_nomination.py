from copy import deepcopy
import csv,hashlib,importlib.util,json
from pathlib import Path
import pytest

from defenseviz import native_pipeline as module


def fixture(tmp_path):
    source=tmp_path/'saved_hits.tsv';source.write_text('query\tfamily_id\tenzyme_role\tscore_total\nx\tType IV\tR\t99\n')
    digest=hashlib.sha256(source.read_bytes()).hexdigest();seq='MAC*';sha=hashlib.sha256(seq.encode()).hexdigest()
    data=dict(genes=[dict(sample_id='s',gene_id='x',protein_sha256=sha)],input_sha256={str(source):digest},
       suite=dict(sources=[dict(sample_id='s',path='native/saved_hits.tsv',resolved_path=str(source),sha256=digest)],
       annotations=[dict(sample_id='s',gene_id='x',native_fields=dict(locus_tag='x',translation=seq,product='protein',REBASEfinder_old='stale'))]),
       saved_decisions=[dict(sample_id='s',gene_id='x',protein_sha256=sha,family='Type IV',role='R',score_total=99,
          source_ref=dict(path='native/saved_hits.tsv',sha256=digest,row=2))])
    directory=tmp_path/'input';directory.mkdir();(directory/'genome.gbff').write_text('original genome')
    out=tmp_path/'work';out.mkdir()
    return data,directory,out,source


def test_prior_candidate_identity_only_no_inherited_roles_or_points(tmp_path):
    data,directory,out,_=fixture(tmp_path);before=deepcopy(data)
    module._prepare_sample(data,'s',directory,out)
    rows=list(csv.DictReader((out/'saved_rm_nominees.tsv').open(),delimiter='\t'))
    assert len(rows)==1 and rows[0]['locus_tag']=='x'
    assert set(rows[0])=={'sample_id','locus_tag','source_protein_sha256','normalized_protein_sha256','aa_length','source_refs'}
    assert rows[0]['source_protein_sha256']==hashlib.sha256(b'MAC*').hexdigest()
    assert rows[0]['normalized_protein_sha256']==hashlib.sha256(b'MAC').hexdigest() and rows[0]['aa_length']=='3'
    assert json.loads(rows[0]['source_refs'])[0]['row']==2
    assert data==before


def test_saved_source_mutation_is_not_a_valid_prior_nomination(tmp_path):
    data,directory,out,source=fixture(tmp_path);source.write_text(source.read_text()+'changed')
    with pytest.raises(ValueError,match='source hash mismatch'):module._prepare_sample(data,'s',directory,out)


def test_prior_protein_mismatch_never_reuses_locus_name(tmp_path):
    data,directory,out,_=fixture(tmp_path);data['saved_decisions'][0]['protein_sha256']='a'*64
    with pytest.raises(ValueError,match='protein binding mismatch'):module._prepare_sample(data,'s',directory,out)


def test_no_translation_remains_explicit_unsearchable_audit(tmp_path):
    data,directory,out,_=fixture(tmp_path)
    data['genes'].append(dict(sample_id='s',gene_id='pseudo',protein_sha256=None))
    data['saved_decisions'].append(dict(sample_id='s',gene_id='pseudo',protein_sha256=None))
    module._prepare_sample(data,'s',directory,out)
    audit=json.loads((out/'saved_rm_nomination.json').read_text())
    assert audit['excluded']==[dict(gene_id='pseudo',reason='no_translated_protein')]
    assert audit['role_type_score_inherited'] is False


def test_mismatched_sample_does_not_nominate_cross_genome_locus(tmp_path):
    data,directory,out,_=fixture(tmp_path);data['saved_decisions'][0]['sample_id']='other'
    module._prepare_sample(data,'s',directory,out)
    assert list(csv.DictReader((out/'saved_rm_nominees.tsv').open(),delimiter='\t'))==[]

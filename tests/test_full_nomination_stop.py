import csv
import hashlib
from pathlib import Path
import pytest
from defenseviz.native_pipeline import _nomination_union


def case(sequence='MAAA*', reported=4, model='REase_II_00013', evalue='1e-80', model_len=4, aligned=4):
    digest=hashlib.sha256(sequence.encode()).hexdigest()
    gene=dict(sample_id='s',gene_id='g',protein_sha256=digest,protein_length=len(sequence))
    profile=dict(sample_id='s',gene_id='g',namespace='PADLOC',
        source_ref=dict(path='native.domtblout',sha256='0'*64,row=2),
        native_fields=dict(query_name=model,target_name='g',query_length=str(model_len),
            target_length=str(reported),ali_from='1',ali_to=str(aligned),hmm_from='1',
            hmm_to=str(aligned),evalue=evalue))
    return dict(genes=[gene],suite=dict(profiles=[profile],annotations=[dict(sample_id='s',gene_id='g',
        native_fields=dict(locus_tag='g',translation=sequence))]))


def execute(tmp_path,data):
    source=tmp_path/'fresh.tsv';out=tmp_path/'union.tsv'
    source.write_text('locus_tag\tprofile_name\trole_hint\tprotein_sha256\tsequence_length\n')
    stats=_nomination_union(data,'s',source,out)
    with out.open()as f:rows=list(csv.DictReader(f,delimiter='\t'))
    return stats,rows


@pytest.mark.parametrize('sequence,reported,coverage',[('MAAA*',4,'1.0'),('MAAA*',5,'0.8'),('MAAA',4,'1.0')])
def test_exact_terminal_stop_representations_bind_clean_native_query(tmp_path,sequence,reported,coverage):
    stats,rows=execute(tmp_path,case(sequence,reported))
    assert stats['saved_qualified_nominees']==1 and len(rows)==1
    assert rows[0]['protein_sha256']==hashlib.sha256(b'MAAA').hexdigest()
    assert rows[0]['sequence_length']=='4'
    assert rows[0]['source_protein_sha256']==hashlib.sha256(sequence.encode()).hexdigest()
    assert rows[0]['query_coverage']==coverage


def test_real_sequence_length_difference_stays_rejected(tmp_path):
    stats,rows=execute(tmp_path,case(reported=3,aligned=3))
    assert stats['saved_qualified_nominees']==0 and rows==[]


def test_changed_protein_hash_cannot_authorize_normalization(tmp_path):
    data=case();data['genes'][0]['protein_sha256']=hashlib.sha256(b'MCCC*').hexdigest()
    with pytest.raises(ValueError,match='binding mismatch'):execute(tmp_path,data)


def test_conflicting_donor_cannot_authorize_normalization(tmp_path):
    data=case();data['suite']['annotations'].append(dict(sample_id='s',gene_id='g',native_fields=dict(locus_tag='g',translation='MCCC*')))
    with pytest.raises(ValueError,match='Conflicting'):execute(tmp_path,data)


def test_no_profile_maxima_mixed_across_models(tmp_path):
    data=case(model_len=40)
    weak=case(model='NEase_II_00007',evalue='1.0')['suite']['profiles'][0]
    data['suite']['profiles'].append(weak)
    stats,rows=execute(tmp_path,data)
    assert stats['saved_qualified_nominees']==0 and rows==[]

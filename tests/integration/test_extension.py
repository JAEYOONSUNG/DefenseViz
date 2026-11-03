import json
import shutil
import pytest
from tests.analysis.test_registry import data
from tests.analysis.test_integration import evidence
from defenseviz.analysis.admission import entry_sha256
from defenseviz.pipeline import run
from defenseviz.assessment import assess


def project(data):
    reg,c,rule,root=data
    reference=reg['families'][0]['references'][0]
    sequence=reference.pop('sequence');reference['sequence_file']='reference.faa'
    (root/'reference.faa').write_text('>ref\n'+sequence+'\n')
    rule['entry_sha256']=entry_sha256(reg['families'][0])
    (root/'samples.tsv').write_text('sample_id\tlabel\narbitrary\tSynthetic\n')
    records=[dict(sample_id='arbitrary',gene_id='unseen_gene',evidence=evidence(old_strong=True,old_context=True,old_motif=False,curated_identity=.9,curated_query_coverage=.9,curated_reference_coverage=.9,curated_evalue=1e-20,curated_role_compatible=True,curated_family_compatible=True),provenance={'sequence_sha256':c['sequence_sha256']})]
    for name,value in [('registry.json',reg),('candidates.json',[c]),('rules.json',[rule]),('evidence.json',records)]:
        (root/name).write_text(json.dumps(value))
    config={'schema_version':1,'samples':'samples.tsv','calls':[],'evidence':'evidence.json','evidence_extension':{'registry':'registry.json','candidates':'candidates.json','admission_rules':'rules.json'}}
    p=root/'project.json';p.write_text(json.dumps(config));return p


def test_extension_integration_hashes_and_relocation(data,tmp_path):
    config=project(data)
    run(config,tmp_path/'one')
    result=json.loads((tmp_path/'one/analysis.json').read_text())
    assert result['reference_anchor_admission']['records'][0]['status']=='admitted'
    assert result['decisions'][0]['decision']['score_motif']==20
    receipt=json.loads((tmp_path/'one/receipt.json').read_text())
    assert {'reference.faa','synthetic_evidence.txt'} <= {x['path'] for x in receipt['inputs']}
    moved=tmp_path/'moved';moved.mkdir()
    for p in config.parent.iterdir():
        if p.is_file():shutil.copy2(p,moved/p.name)
    run(moved/'project.json',tmp_path/'two')
    assert result==json.loads((tmp_path/'two/analysis.json').read_text())


def test_assessment_command_records_reference_file(data,tmp_path):
    project(data)
    assess(tmp_path/'registry.json',tmp_path/'candidates.json',tmp_path/'assessment')
    receipt=json.loads((tmp_path/'assessment/receipt.json').read_text())
    assert 'reference.faa' in {x['path'] for x in receipt['inputs']}
    assert receipt['detector_runs']==0


def test_malformed_empty_extension_never_succeeds(data,tmp_path):
    config=project(data)
    (tmp_path/'candidates.json').write_text('[]')
    (tmp_path/'registry.json').write_text('{"invalid":true}')
    (tmp_path/'rules.json').write_text('{"invalid":true}')
    with pytest.raises(ValueError):run(config,tmp_path/'bad')
    assert not (tmp_path/'bad').exists()

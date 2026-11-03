import copy
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import pytest
from defenseviz.demo import write_project
from defenseviz.pipeline import run


def test_end_to_end(tmp_path):
    config=write_project(tmp_path/'input')
    before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in config.parent.iterdir()}
    output=tmp_path/'result'
    summary=run(config,output,argv=['test'])
    result=json.loads((output/'analysis.json').read_text())
    assert summary['samples']==3
    by={r['sample_id']:r for r in result['sample_summary']}
    assert [by[s]['union_cds'] for s in ('isolate-a','isolate-b','isolate-c')]==[5,3,1]
    assert by['isolate-b']['retained_dnmb_only_cds']==1
    assert by['isolate-c']['caller_counts']['PADLOC'] is None
    assert all(r['family']!='domain-only' for r in result['cells'])
    receipt=json.loads((output/'receipt.json').read_text())
    assert receipt['argv']==['test'] and receipt['exit_code']==0 and receipt['detector_runs']==0
    for name,sha in receipt['outputs'].items():
        assert hashlib.sha256((output/name).read_bytes()).hexdigest()==sha
    assert before=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in config.parent.iterdir()}
    with pytest.raises(FileExistsError):run(config,output)


def test_identifier_invariance(tmp_path):
    from defenseviz.io import load_project
    from defenseviz.analysis import integrate
    config=write_project(tmp_path/'input');data=load_project(config)
    a=integrate(data['samples'],data['calls'],data['evidence_records'],coverage=data['coverage'])
    renamed=copy.deepcopy(data);ids={s['sample_id']:f'new__sample/{i}:renamed' for i,s in enumerate(data['samples'])}
    for key in ('samples','calls','evidence_records','coverage'):
        for row in renamed[key]:
            row['sample_id']=ids[row['sample_id']]
            if 'gene_id' in row:row['gene_id']='new|'+row['gene_id']
    b=integrate(renamed['samples'],renamed['calls'],renamed['evidence_records'],coverage=renamed['coverage'])
    def without_ids(rows):return [{k:v for k,v in r.items() if k not in ('sample_id','gene_id')} for r in rows]
    assert without_ids(a['sample_summary'])==without_ids(b['sample_summary'])
    assert without_ids(a['decisions'])==without_ids(b['decisions'])


def test_bad_input_no_partial_result(tmp_path):
    config=write_project(tmp_path/'input');doc=json.loads(config.read_text());doc['options']['max_rows_per_page']=0;config.write_text(json.dumps(doc))
    with pytest.raises((ValueError,TypeError)):run(config,tmp_path/'output')
    assert not (tmp_path/'output').exists()
    assert not list(tmp_path.glob('.output.staging-*'))


def test_cli_actionable_missing_file(tmp_path):
    p=subprocess.run([sys.executable,'-m','defenseviz','validate',str(tmp_path/'missing.json')],capture_output=True,text=True)
    assert p.returncode==2 and 'defenseviz:' in p.stderr and 'Traceback' not in p.stderr

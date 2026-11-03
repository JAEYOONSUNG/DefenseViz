"""Small entirely synthetic dataset distributed with the program."""
from __future__ import annotations
import csv
import json
import shutil
import tempfile
from pathlib import Path


def synthetic_evidence(family='Type II', role='M', **overrides):
    value = dict(family=family, role=role, valid_sequence=True, identity=.9, query_coverage=.95,
                 reference_coverage=.95, evalue=1e-30, role_compatible=True, family_compatible=True,
                 hmm_role_ambiguous=False, profiles=[], raw_nuclease_classes=[], legacy_nuclease=False,
                 legacy_motor=False, legacy_mtase_motif=False, old_strong=True, old_moderate=False,
                 old_profile=False, annotation=True, old_motif=True, old_context=True,
                 strict_structure=False, alternative_annotation=False, conflict=False, partial=False)
    value.update(overrides)
    return value


def write_project(directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=False)
    def table(name, rows, fields=None):
        with (directory / name).open('w', encoding='utf-8', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=fields or list(rows[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(rows)
    table('samples.tsv', [
        dict(sample_id='isolate-a', label='Example isolate A', genus='Examplegenus', species='example A', highlight='true'),
        dict(sample_id='isolate-b', label='Example isolate B', genus='Examplegenus', species='example B', highlight='false'),
        dict(sample_id='isolate-c', label='Example isolate C', genus='Othergenus', species='example C', highlight='false')])
    calls=[]
    def call(s,g,caller,family,category,system,role='',rm_type='',kind='final_member'):
        calls.append(dict(sample_id=s,gene_id=g,caller=caller,family=family,category=category,system_id=system,role=role,rm_type=rm_type,evidence_kind=kind))
    for caller in ('DefenseFinder','PADLOC'):
        call('isolate-a','gene-r',caller,'RM_TypeI','R-M','rm1','R','Type I')
        call('isolate-a','gene-m',caller,'RM_TypeI','R-M','rm1','M','Type I')
        call('isolate-a','gene-s',caller,'RM_TypeI','R-M','rm1','S','Type I')
    call('isolate-a','gene-g', 'PADLOC','Gabija','Nuclease','gab1')
    call('isolate-a','gene-x', 'DefenseFinder','Uncatalogued family','Unresolved','novel1')
    call('isolate-b','gene-c', 'DefenseFinder','Cas-example','CRISPR-Cas','cas1')
    call('isolate-b','gene-b', 'PADLOC','BREX','Modification','brex1')
    call('isolate-c','gene-z', 'DefenseFinder','Uncatalogued family','Unresolved','novel2')
    call('isolate-a','domain-only', 'DefenseFinder','Gabija','Nuclease','','','',kind='hmm_hit')
    table('calls.tsv',calls)
    table('coverage.tsv',[dict(sample_id=s,caller=c,status='missing' if s=='isolate-c' and c=='PADLOC' else 'complete') for s in ('isolate-a','isolate-b','isolate-c') for c in ('DefenseFinder','PADLOC')])
    evidence=[dict(sample_id='isolate-a',gene_id=g,evidence=synthetic_evidence('Type I',role),provenance={'kind':'synthetic_demo_only'}) for g,role in [('gene-r','R'),('gene-m','M'),('gene-s','S')]]
    evidence.append(dict(sample_id='isolate-b',gene_id='dnmb-only',evidence=synthetic_evidence(),provenance={'kind':'synthetic_demo_only'}))
    (directory/'evidence.json').write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf-8')
    project={'schema_version':1,'samples':'samples.tsv','calls':[{'format':'normalized','path':'calls.tsv'}], 'coverage':'coverage.tsv','evidence':'evidence.json','options':{'title':'DefenseViz · synthetic demonstration'}}
    path=directory/'project.json';path.write_text(json.dumps(project,indent=2)+'\n',encoding='utf-8')
    return path


def demo(output_dir: Path) -> dict:
    from .pipeline import run
    output_dir=Path(output_dir).absolute()
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(f'Demo output already exists: {output_dir}')
    output_dir.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix=f'.{output_dir.name}.demo-',dir=output_dir.parent))
    try:
        config=write_project(stage/'inputs')
        result=run(config,stage/'results',argv=['defenseviz','demo','--out',str(output_dir)])
        if output_dir.exists():raise FileExistsError(f'Output appeared during demo: {output_dir}')
        stage.rename(output_dir)
        return {'output':str(output_dir),'config':str(output_dir/'inputs/project.json'), 'figure':str(output_dir/'results/figure2_ab.pdf'),'samples':result['samples'],'data_kind':'synthetic'}
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise

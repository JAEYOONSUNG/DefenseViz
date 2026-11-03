"""Source-bound registry evaluation; no sequence searches or structure prediction."""
from __future__ import annotations
import json
import os
from pathlib import Path
import shutil
import tempfile
import time
from .pipeline import digest, software_state, write_json, portable_provenance


def _strict_json(path):
    def pairs(values):
        result={}
        for k,v in values:
            if k in result:raise ValueError(f'Duplicate JSON key {k!r} in {path}')
            result[k]=v
        return result
    def bad(token):raise ValueError(f'Nonfinite JSON number {token!r} in {path}')
    return json.loads(path.read_text(encoding='utf-8'),object_pairs_hook=pairs,parse_constant=bad)


def assess(registry_path, candidates_path, output_dir, *, source_root=None, argv=None):
    from .analysis.registry import assess_registry
    start=time.monotonic();rp=Path(registry_path).resolve();cp=Path(candidates_path).resolve();out=Path(output_dir).absolute()
    if out.exists() or out.is_symlink():raise FileExistsError(f'Output already exists: {out}')
    inputs={str(rp):digest(rp),str(cp):digest(cp)}
    registry=_strict_json(rp);candidates=_strict_json(cp)
    if not isinstance(registry,dict) or not isinstance(candidates,list):raise ValueError('Registry must be an object and candidates a list')
    for family in registry.get('families',[]):
        for reference in family.get('references',[]):
            if reference.get('sequence_file'):
                path=Path(reference['sequence_file']);path=(rp.parent/path).resolve() if not path.is_absolute() else path
                inputs[str(path)]=digest(path)
    results=assess_registry(registry,candidates,registry_dir=rp.parent,source_root=source_root or rp.parent,verify_sources=True)
    for row in results:
        for name,sha in row.get('source_sha256',{}).items():
            if name in inputs and inputs[name]!=sha:raise ValueError(f'Input changed: {name}')
            inputs[name]=sha
    for name,sha in inputs.items():
        if digest(Path(name))!=sha:raise ValueError(f'Input changed during assessment: {name}')
    out.parent.mkdir(parents=True,exist_ok=True);stage=Path(tempfile.mkdtemp(prefix=f'.{out.name}.staging-',dir=out.parent))
    try:
        write_json(stage/'assessment.json',portable_provenance(results,rp.parent))
        write_json(stage/'receipt.json',dict(operation='saved_domain_and_reference_evidence_assessment',detector_runs=0,structure_prediction_runs=0,
            input_path_base='registry directory',inputs=[dict(path=os.path.relpath(name,rp.parent),sha256=sha) for name,sha in sorted(inputs.items())],
            software=software_state(),argv=argv,exit_code=0,elapsed_seconds=time.monotonic()-start,
            outputs={'assessment.json':digest(stage/'assessment.json')},activity='not_assessed'))
        if out.exists():raise FileExistsError(f'Output appeared during assessment: {out}')
        stage.rename(out)
        return dict(output=str(out),candidates=len(results),activity='not_assessed')
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True);raise

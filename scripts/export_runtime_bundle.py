#!/usr/bin/env python3
"""Export local DNMB resources without copying sample data or changing originals.

Requires an installed DefenseViz (or its source checkout on PYTHONPATH).
System R, BLAST and HMMER remain explicit platform dependencies.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')


def copy_file(source,target,bindings):
    source=Path(source).resolve();target=Path(target)
    before=sha(source);target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,target)
    if sha(source)!=before or sha(target)!=before:raise ValueError('Changed resource: '+str(source))
    bindings[str(source)]=before


def copy_tree(source,target,bindings,include=None,omitted=None):
    source=Path(source).resolve();target=Path(target);target.mkdir(parents=True,exist_ok=False)
    for item in sorted(source.rglob('*')):
        rel=item.relative_to(source)
        if any(p in {'.git','.DS_Store','__pycache__'} for p in rel.parts):continue
        if item.is_symlink() and item.is_dir():raise ValueError('Directory symlink needs explicit packaging: '+str(item))
        if item.is_file():
            if include is not None and not include(rel):
                if omitted is not None:omitted.append(dict(source=str(item),relative=rel.as_posix(),sha256=sha(item),reason='outside runtime package source/assets allowlist'))
                continue
            copy_file(item,target/rel,bindings)


def package_asset_allowed(relative,name,source=False):
    """Keep installable code and audited references, never historical samples."""
    parts=Path(relative).parts
    if source:
        if len(parts)==1:return parts[0] in {'DESCRIPTION','NAMESPACE','LICENSE','LICENSE.md','LICENSE.txt','NEWS','NEWS.md','configure','configure.win','cleanup','cleanup.win'}
        if parts[0] in {'R','man','src'}:return True
        if parts[0]!='inst':return False
        parts=parts[1:]
    if parts[0]=='validation':return False
    if parts[0]=='scripts':return Path(parts[-1]).suffix in {'.R','.r','.py','.sh'}
    if parts[0]=='extdata':
        return name=='DNMB' and len(parts)>=2 and parts[1] in {'rebasefinder_structure_refs.tsv','rebasefinder_structure_templates'}
    return not source


def inventory(root):
    root=Path(root)
    return {p.relative_to(root).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size)
            for p in sorted(root.rglob('*')) if p.is_file() and p!=root/'bundle_manifest.json'}


def dependency_closure(cfg,work,env):
    """Resolve declared imports plus every namespace measured by the full probe."""
    from defenseviz.native_pipeline import _capture_r_runtime
    runtime,resources=_capture_r_runtime(cfg,work,env)
    script=Path(work)/'dependency_closure.R';output=Path(work)/'dependency_closure.json'
    script.write_text('''a<-commandArgs(TRUE)
.libPaths(c(a[-c(1,2)],.libPaths()))
runtime<-jsonlite::fromJSON(a[1])
db<-installed.packages(lib.loc=.libPaths(),noCache=TRUE)
db<-db[!duplicated(db[,"Package"]),,drop=FALSE]
seeds<-runtime$packages$name
deps<-tools::package_dependencies(seeds,db=db,which=c("Depends","Imports","LinkingTo"),recursive=TRUE)
names<-sort(unique(c(seeds,unlist(deps))))
names<-setdiff(names,"R")
rows<-lapply(names,function(n){p<-normalizePath(find.package(n));d<-read.dcf(file.path(p,"DESCRIPTION"));priority<-if("Priority"%in%colnames(d))d[1,"Priority"]else "";list(name=n,path=p,version=as.character(packageVersion(n,lib.loc=dirname(p))),license=unname(d[1,"License"]),system=priority%in%c("base","recommended")&&startsWith(p,paste0(normalizePath(R.home("library")),"/")))})
jsonlite::write_json(list(R_version=R.version.string,R_home=R.home(),packages=rows),a[2],auto_unbox=TRUE,pretty=TRUE)
''')
    command=[cfg['rscript'],'--vanilla',str(script),str(Path(work)/'r_runtime.json'),str(output),*cfg['r_libraries']]
    start=time.monotonic();p=subprocess.run(command,capture_output=True,text=True,env=env)
    write(Path(work)/'dependency_closure_execution.json',dict(argv=command,exit_code=p.returncode,
        elapsed_seconds=time.monotonic()-start,stdout=p.stdout,stderr=p.stderr))
    if p.returncode:raise RuntimeError('R dependency closure failed: '+p.stderr)
    return runtime,json.loads(output.read_text()),resources


def cath_wrapper():
    return '''#!/bin/sh
set -eu
bundle_bin=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
export DYLD_LIBRARY_PATH="$bundle_bin/../lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="$bundle_bin/../lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$bundle_bin/cath-resolve-hits.bin" "$@"
'''


REBUILD_SCRIPT=r'''#!/usr/bin/env python3
"""Rebuild the two exact preserved source packages into a new local R library."""
import argparse,hashlib,json,os,pathlib,shutil,subprocess,time
p=argparse.ArgumentParser();p.add_argument('--library',required=True);p.add_argument('--rscript');a=p.parse_args()
b=pathlib.Path(__file__).resolve().parents[1]
m=json.loads((b/'bundle_manifest.json').read_text())
def sha(x):
 h=hashlib.sha256()
 with x.open('rb') as f:
  for v in iter(lambda:f.read(1048576),b''):h.update(v)
 return h.hexdigest()
for name,r in m['artifacts'].items():
 if sha(b/name)!=r['sha256']:raise SystemExit('Bundle changed: '+name)
cfg=json.loads((b/'engine.json').read_text());rscript=a.rscript or cfg['rscript']
lib=pathlib.Path(a.library).expanduser().resolve()
if lib==b or b in lib.parents:raise SystemExit('Use a new library outside the preserved bundle')
lib.mkdir(parents=True,exist_ok=False)
for source in (b/'r-library').iterdir():
 if source.name not in ('DNMB','DefenseViz'):shutil.copytree(source,lib/source.name)
env=dict(os.environ);env['R_LIBS']=str(lib);env['R_LIBS_USER']=str(lib)
records=[]
for name in ('DefenseViz','DNMB'):
 # The R executable alongside Rscript is the same configured installation.
 R=pathlib.Path(rscript).resolve().with_name('R')
 cmd=[str(R),'CMD','INSTALL','--library='+str(lib),str(b/'sources'/name)]
 start=time.monotonic();proc=subprocess.run(cmd,capture_output=True,text=True,env=env)
 records.append(dict(argv=cmd,exit_code=proc.returncode,elapsed_seconds=time.monotonic()-start,stdout=proc.stdout,stderr=proc.stderr))
 (lib/'rebuild_receipt.json').write_text(json.dumps(records,indent=2)+'\n')
 if proc.returncode:raise SystemExit('R installation failed; see '+str(lib/'rebuild_receipt.json'))
print(str(lib))
'''


def export_bundle(config,out,dnmb_source,legacy_source,evidence=()):
    from defenseviz.native_pipeline import read_engine_config,_runtime_environment,_resolve_blastdbcmd
    cfg=read_engine_config(config);out=Path(out).expanduser().resolve()
    sources={'DNMB':Path(dnmb_source).resolve(),'DefenseViz':Path(legacy_source).resolve()}
    for name,path in sources.items():
        desc=(path/'DESCRIPTION').read_text()
        if not any(line.strip()=='Package: '+name for line in desc.splitlines()):
            raise ValueError('Wrong source package: '+str(path))
    for source in [*sources.values(),Path(cfg['reference_cache']),*map(Path,cfg['r_libraries'])]:
        if out==source or source in out.parents:raise ValueError('Bundle must be outside its source directories')
    out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();bindings={str(Path(config).resolve()):sha(config)}
    probe=out/'provenance'/'original-runtime';probe.mkdir(parents=True)
    runtime,closure,runtime_hashes=dependency_closure(cfg,probe,_runtime_environment(cfg,out/'validation-cache'))
    write(out/'provenance'/'measured_runtime_assets.json',runtime_hashes)
    omitted=[]
    for package in closure['packages']:
        if not package['system']:
            name=package['name'];allow=(lambda rel,n=name:package_asset_allowed(rel,n)) if name in sources else None
            copy_tree(package['path'],out/'r-library'/name,bindings,allow,omitted)
    for name,path in sources.items():
        copy_tree(path,out/'sources'/name,bindings,lambda rel,n=name:package_asset_allowed(rel,n,source=True),omitted)
    write(out/'provenance/omitted_nonruntime_package_assets.json',omitted)
    cache=Path(cfg['reference_cache'])
    required=[cache/'rebase_data.rds',cache/'rebase_db.fasta',*sorted(cache.glob('rebase_db.fasta.*'))]
    for p in required:copy_file(p,out/'databases/rebase'/p.name,bindings)
    for name in ('REBASE_protein_seqs.txt','rebase_reference_state.rds'):
        if (cache/name).is_file():copy_file(cache/name,out/'databases/rebase'/name,bindings)
    copy_file(cfg['hmm_database'],out/'databases/rm-hmm/RM_HMMs.hmm',bindings)
    copy_file(cfg['profile_metadata'],out/'databases/rm-hmm/profile_roles.tsv',bindings)
    # The configured binary wrapper may contain source-workspace paths. Preserve
    # it only as provenance and generate the executable wrapper from relative paths.
    original_cath=Path(cfg['cath_resolve_hits'])
    copy_file(original_cath,out/'provenance/original-cath-wrapper',bindings)
    # read_engine_config resolves symlinks; loader aliases such as libicu*.68
    # must retain their configured names when copied to a different directory.
    original_cfg=json.loads(Path(config).read_text())
    candidates=[Path(p).expanduser() for p in original_cfg.get('runtime_dependencies',[])]
    candidates=[p if p.is_absolute() else Path(config).resolve().parent/p for p in candidates]
    binary=next((p for p in candidates if p.name=='cath-resolve-hits.bin'),None)
    if binary is None:
        if original_cath.read_bytes().startswith(b'#!'):raise ValueError('CATH wrapper needs its binary in runtime_dependencies')
        binary=original_cath
    copy_file(binary,out/'bin/cath-resolve-hits.bin',bindings)
    runtime_dependencies=['bin/cath-resolve-hits.bin']
    targets={}
    for dep in candidates:
        if dep.resolve()==binary.resolve():continue
        if dep.name in targets and targets[dep.name]!=dep.resolve():raise ValueError('Dependency filename collision: '+dep.name)
        targets[dep.name]=dep.resolve()
        copy_file(dep,out/'lib'/dep.name,bindings)
        runtime_dependencies.append('lib/'+dep.name)
    (out/'bin/cath-resolve-hits').write_text(cath_wrapper());(out/'bin/cath-resolve-hits').chmod(0o755)
    copied=dict(schema_version=1,rscript=cfg['rscript'],r_libraries=['r-library'],
        reference_cache='databases/rebase',hmm_database='databases/rm-hmm/RM_HMMs.hmm',
        profile_metadata='databases/rm-hmm/profile_roles.tsv',blastp=cfg['blastp'],
        blastdbcmd=_resolve_blastdbcmd(cfg),hmmsearch=cfg['hmmsearch'],
        cath_resolve_hits='bin/cath-resolve-hits',runtime_dependencies=runtime_dependencies,cpu=cfg['cpu'])
    write(out/'engine.json',copied)
    (out/'setup').mkdir();(out/'setup/rebuild_r.py').write_text(REBUILD_SCRIPT);(out/'setup/rebuild_r.py').chmod(0o755)
    for index,p in enumerate(map(Path,evidence)):
        target=out/'provenance'/'source-evidence'/f'{index:02d}-{p.name}'
        if p.is_dir():copy_tree(p,target,bindings)
        else:copy_file(p,target,bindings)
    write(out/'provenance'/'source_bindings.json',bindings)
    source_state={}
    for name,path in sources.items():
        source_state[name]={}
        for label,args in [('commit',['rev-parse','HEAD']),('status',['status','--porcelain=v1']),('remote',['remote','-v'])]:
            proc=subprocess.run(['git','--no-optional-locks','-C',str(path),*args],capture_output=True,text=True)
            source_state[name][label]=dict(exit_code=proc.returncode,stdout=proc.stdout,stderr=proc.stderr)
    write(out/'provenance/source_git.json',source_state)
    write(out/'provenance/package_licenses.json',[{k:p[k] for k in ('name','version','license','system')} for p in closure['packages']])
    (out/'README.txt').write_text('DNMB local runtime bundle\n\nengine.json contains relative bundled-resource paths. System R, BLAST and HMMER remain external. Copied R packages and CATH/ICU are platform-bound.\nRebuild exact DNMB/legacy source packages: python setup/rebuild_r.py --library /new/r-library\nUse the rebuilt library in a new engine configuration and run the full validation.\nReference/HMM data and source evidence are local preservation assets, not a grant to redistribute. See provenance/package_licenses.json and provided source-evidence. No sample datasets are included.\n')
    for path,expected in bindings.items():
        if sha(path)!=expected:raise ValueError('Original resource changed during export: '+path)
    manifest=dict(schema_version=1,created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        platform=dict(system=platform.system(),machine=platform.machine(),python=platform.python_version()),
        R_version=runtime['R_version'],system_executables={k:dict(path=copied[k],sha256=sha(copied[k])) for k in ('rscript','blastp','blastdbcmd','hmmsearch')},
        packages=len(closure['packages']),bundled_packages=sum(not p['system'] for p in closure['packages']),
        source_inputs_unchanged=True,elapsed_seconds=time.monotonic()-started,artifacts=inventory(out))
    write(out/'bundle_manifest.json',manifest)
    return manifest


def verify_bundle(bundle,report):
    from defenseviz.native_pipeline import read_engine_config,_runtime_environment,_capture_r_runtime,_validate_reference_index
    bundle=Path(bundle).resolve();report=Path(report).resolve()
    if report==bundle or bundle in report.parents:raise ValueError('Verification report must be outside the preserved bundle')
    manifest=json.loads((bundle/'bundle_manifest.json').read_text())
    actual=inventory(bundle)
    if actual!=manifest['artifacts']:raise ValueError('Bundle file inventory or content changed')
    for record in manifest['system_executables'].values():
        if sha(record['path'])!=record['sha256']:raise ValueError('Configured system executable changed: '+record['path'])
    cfg=read_engine_config(bundle/'engine.json');report.mkdir(parents=True,exist_ok=False)
    env=_runtime_environment(cfg,report/'cache')
    runtime,hashes=_capture_r_runtime(cfg,report,env)
    if runtime['R_version']!=manifest['R_version']:raise ValueError('R version differs from the preserved runtime')
    r_home=(Path(runtime['R_home'])/'library').resolve();library=(bundle/'r-library').resolve()
    for package in runtime['packages']:
        p=Path(package['path']).resolve()
        if not (p.is_relative_to(library) or p.is_relative_to(r_home)):
            raise ValueError('Runtime escaped preserved R libraries: '+str(p))
    command=[cfg['cath_resolve_hits'],'--version'];start=time.monotonic()
    p=subprocess.run(command,capture_output=True,text=True,env=env)
    write(report/'cath_execution.json',dict(argv=command,exit_code=p.returncode,elapsed_seconds=time.monotonic()-start,stdout=p.stdout,stderr=p.stderr))
    if p.returncode:raise RuntimeError('Relocated CATH executable failed: '+p.stderr)
    db=_validate_reference_index(Path(cfg['reference_cache'])/'rebase_db.fasta',cfg['blastdbcmd'],report/'reference_validation.json',env)
    result=dict(passed=True,bundle=str(bundle),bundle_manifest_sha256=sha(bundle/'bundle_manifest.json'),
        measured_R_namespaces=len(runtime['packages']),R_packages_bound_to_bundle_or_R_installation=True,
        reference_entries=db['entries'],CATH_exit_code=0,full_analysis_executed=False)
    write(report/'verification.json',result);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('export');p.add_argument('--engine-config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--dnmb-source',required=True);p.add_argument('--legacy-source',required=True)
    p.add_argument('--evidence',action='append',default=[])
    v=sub.add_parser('verify');v.add_argument('--bundle',required=True);v.add_argument('--report',required=True)
    a=parser.parse_args()
    if a.command=='export':result=export_bundle(a.engine_config,a.out,a.dnmb_source,a.legacy_source,a.evidence)
    else:result=verify_bundle(a.bundle,a.report)
    print(json.dumps({k:v for k,v in result.items() if k!='artifacts'},indent=2))


if __name__=='__main__':main()

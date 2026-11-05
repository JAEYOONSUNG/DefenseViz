#!/usr/bin/env python3
"""Install pinned public DNMB R namespaces into a fresh private directory."""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

SOURCES={
    'DNMB':dict(url='https://github.com/JAEYOONSUNG/DNMB.git',commit='976489e9967e9e4dfacc61f7e83e66b1da0f4f2b'),
    'DefenseViz':dict(url='https://github.com/JAEYOONSUNG/DefenseViz.git',commit='53a1dc4b8c6595945762e21354fbdc86d2b3731a'),
}
PATTERNS=['/DESCRIPTION','/NAMESPACE','/LICENSE','/LICENSE.md','/LICENSE.txt','/NEWS','/NEWS.md','/R/','/man/','/src/']
HELPERS=['Genbank_organizer','.dnmb_write_protein_fasta','.dnmb_module_clean_annotation_key',
    '.dnmb_na_vector_like','.dnmb_normalize_translation','%||%','arrange.vars','.dnmb_db_cache_root',
    '.dnmb_module_empty_optional_long_table','.dnmb_download_asset','.dnmb_nonempty_file',
    '.dnmb_remote_asset_metadata','.dnmb_transactional_replace','.dnmb_copy_local_asset',
    '.dnmb_external_timeout_seconds','dnmb_detect_binary','dnmb_run_external','.dnmb_parse_http_headers','.dnmb_db_validate_key']

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def write(path,obj):Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')

def run(argv,records,receipt,cwd=None,env=None):
    start=time.monotonic();p=subprocess.run(list(map(str,argv)),cwd=cwd,env=env,text=True,capture_output=True)
    record=dict(argv=list(map(str,argv)),cwd=str(cwd) if cwd else None,exit_code=p.returncode,
                elapsed_seconds=time.monotonic()-start,stdout=p.stdout,stderr=p.stderr)
    records.append(record);write(receipt,records)
    if p.returncode:raise RuntimeError('Command failed; inspect '+str(receipt)+'\n'+p.stderr[-3000:])
    return p.stdout

def fetch_source(name,spec,destination,records,receipt):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=False)
    run(['git','init','--quiet',str(destination)],records,receipt)
    run(['git','-C',str(destination),'remote','add','origin',spec['url']],records,receipt)
    run(['git','-C',str(destination),'config','core.sparseCheckout','true'],records,receipt)
    patterns=list(PATTERNS)
    if name=='DNMB':patterns+=['/inst/scripts/','/inst/extdata/rebasefinder_structure_refs.tsv','/inst/extdata/rebasefinder_structure_templates/']
    (destination/'.git/info/sparse-checkout').write_text('\n'.join(patterns)+'\n')
    run(['git','-C',str(destination),'fetch','--filter=blob:none','--depth','1','origin',spec['commit']],records,receipt)
    run(['git','-C',str(destination),'checkout','--detach',spec['commit']],records,receipt)
    observed=run(['git','--no-optional-locks','-C',str(destination),'rev-parse','HEAD'],records,receipt).strip()
    if observed!=spec['commit']:raise ValueError('Downloaded source commit does not match pin: '+name)
    status=run(['git','--no-optional-locks','-C',str(destination),'status','--porcelain=v1'],records,receipt)
    if status.strip():raise ValueError('Downloaded source working tree is modified: '+name)
    description=(destination/'DESCRIPTION').read_text()
    if 'Package: '+name not in description:raise ValueError('Wrong source package: '+name)
    license=(destination/'LICENSE').read_text()
    if not ('MIT License' in license or 'Permission is hereby granted, free of charge' in license):
        raise ValueError('Expected MIT license missing: '+name)
    artifacts={p.relative_to(destination).as_posix():sha(p) for p in destination.rglob('*')
               if p.is_file() and '.git' not in p.relative_to(destination).parts}
    if any(p.endswith(('.gbk','.gbff','.xlsx')) for p in artifacts):raise ValueError('Unexpected genome/sample artifact in selected source')
    return dict(**spec,observed_commit=observed,working_tree_clean=True,selected_source_sha256=artifacts)


DEPENDENCY_PROBE='''a<-commandArgs(TRUE)
.libPaths(c(a[-c(1,2,3)],R.home("library")))
req<-character()
for(p in a[1:2]){d<-read.dcf(file.path(p,"DESCRIPTION"));for(k in intersect(c("Depends","Imports","LinkingTo"),colnames(d)))req<-c(req,trimws(gsub("\\\\([^)]*\\\\)","",strsplit(d[1,k],",")[[1]])))}
req<-setdiff(unique(c(req,"jsonlite","digest","withr","Peptides","gtools","openxlsx","plyr","reshape2","stringr","tidyr")),c("R","DNMB","DefenseViz",""))
missing<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)){text<-paste("Missing R dependencies:",paste(missing,collapse=", "));writeLines(text,a[3]);message(text);quit(status=2)}
rows<-lapply(sort(loadedNamespaces()),function(n){p<-normalizePath(if(n=="base")file.path(R.home("library"),"base")else getNamespaceInfo(asNamespace(n),"path"));list(name=n,path=p,version=as.character(packageVersion(n,lib.loc=dirname(p))))})
jsonlite::write_json(list(R_version=R.version.string,R_home=R.home(),packages=rows),a[3],auto_unbox=TRUE,pretty=TRUE)
'''

HELPER_PROBE='''a<-commandArgs(TRUE);.libPaths(c(a[-c(1,2,3)],R.home("library")))
loadNamespace("DNMB");loadNamespace("DefenseViz")
names<-jsonlite::fromJSON(a[2]);expected<-setNames(vector("list",length(names)),names)
for(p in list.files(file.path(a[1],"R"),pattern="[.]R$",full.names=TRUE))for(e in parse(p))if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%names){k<-as.character(e[[2]]);expected[[k]]<-c(expected[[k]],list(e[[3]]))}
rows<-lapply(names,function(n){f<-get(n,envir=asNamespace("DNMB"),inherits=FALSE);matched<-any(vapply(expected[[n]],function(e)identical(formals(f),e[[2]],ignore.srcref=TRUE)&&identical(body(f),e[[3]],ignore.srcref=TRUE),logical(1)));actual<-paste(deparse(as.call(list(as.name("function"),formals(f),body(f))),width.cutoff=500L),collapse="\\n");list(name=n,source_body_matches=matched,body_sha256=digest::digest(actual,algo="sha256",serialize=FALSE))})
packages<-lapply(c("DNMB","DefenseViz"),function(n)list(name=n,path=normalizePath(getNamespaceInfo(asNamespace(n),"path")),version=as.character(packageVersion(n))))
jsonlite::write_json(list(helpers=rows,packages=packages),a[3],auto_unbox=TRUE,pretty=TRUE)
if(!all(vapply(rows,function(r)r$source_body_matches,logical(1))))stop("Installed helper bodies differ from pinned source")
'''

def install_runtime(out,rscript,dependency_libraries):
    out=Path(out).expanduser().resolve()
    executable=shutil.which(rscript) if '/' not in rscript else rscript
    if not executable or not Path(executable).is_file():raise FileNotFoundError('Rscript not found: '+rscript)
    executable=Path(executable).resolve();R=executable.with_name('R')
    if not R.is_file():raise FileNotFoundError('R executable must be beside Rscript: '+str(R))
    libraries=[Path(p).expanduser().resolve() for p in dependency_libraries]
    if not libraries:raise ValueError('Provide at least one existing --dependency-library from a compatible DNMBsuite R environment')
    for p in libraries:
        if not p.is_dir():raise FileNotFoundError('Dependency library missing: '+str(p))
        if out==p or p in out.parents:raise ValueError('Destination must be outside dependency libraries')
    out.mkdir(parents=True,exist_ok=False);records=[];receipt=out/'commands.json';sources={};start=time.monotonic()
    env=dict(os.environ);env['GIT_OPTIONAL_LOCKS']='0'
    library=out/'library';library.mkdir()
    env['R_LIBS']=os.pathsep.join(map(str,[library,*libraries]));env['R_LIBS_USER']=env['R_LIBS'];env['R_LIBS_SITE']=''
    try:
        for name,spec in SOURCES.items():sources[name]=fetch_source(name,spec,out/'sources'/name,records,receipt)
        write(out/'source_receipt.json',sources)
        probe=out/'dependency_probe.R';probe.write_text(DEPENDENCY_PROBE)
        run([executable,'--vanilla',probe,out/'sources/DNMB',out/'sources/DefenseViz',out/'dependencies.json',*libraries],records,receipt,env=env)
        dependencies=json.loads((out/'dependencies.json').read_text())
        allowed=[*libraries,Path(dependencies['R_home']).resolve()/'library']
        escapes=[p for p in dependencies['packages'] if not any(Path(p['path']).resolve().is_relative_to(d) for d in allowed)]
        if escapes:raise ValueError('Dependencies loaded outside explicit libraries: '+', '.join(p['path'] for p in escapes))
        for name in ('DefenseViz','DNMB'):
            run([R,'CMD','INSTALL','--library='+str(library),out/'sources'/name],records,receipt,env=env)
        probe=out/'helper_probe.R';probe.write_text(HELPER_PROBE);write(out/'helpers.json',HELPERS)
        run([executable,'--vanilla',probe,out/'sources/DNMB',out/'helpers.json',out/'helper_verification.json',library,*libraries],records,receipt,env=env)
        helpers=json.loads((out/'helper_verification.json').read_text())
        if any(not Path(p['path']).is_relative_to(library) for p in helpers['packages']):raise ValueError('Old DNMB namespace was loaded instead of new installation')
        for name,source in sources.items():
            for rel,digest in source['selected_source_sha256'].items():
                if sha(out/'sources'/name/rel)!=digest:raise ValueError('Pinned source changed during installation: '+name+'/'+rel)
        result=dict(schema_version=1,installed=True,source='pinned_public_GitHub_commits',sources=sources,
            Rscript=dict(path=str(executable),sha256=sha(executable)),R_version=dependencies['R_version'],
            r_libraries=[str(library),*map(str,libraries)],automatic_dependency_installations=0,
            helper_source_matches=len(helpers['helpers']),full_analysis_executed=False,
            created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start)
        write(out/'runtime_installation.json',result);return result
    except Exception as error:
        write(out/'failure.json',dict(installed=False,error=str(error),elapsed_seconds=time.monotonic()-start))
        raise RuntimeError('Native runtime installation failed; inspect '+str(out)+': '+str(error)) from error

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True);p.add_argument('--rscript',default='Rscript')
    p.add_argument('--dependency-library',action='append',required=True)
    a=p.parse_args();result=install_runtime(a.out,a.rscript,a.dependency_library)
    print(json.dumps({k:v for k,v in result.items() if k!='sources'},indent=2))

if __name__=='__main__':main()

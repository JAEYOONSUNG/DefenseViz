# Full fixed native annotation/HMM nomination + reference search + context core.
# Source libraries/data are read-only. Every generated artifact is under out.
a <- commandArgs(TRUE); stopifnot(length(a)%in%c(8L,9L));reuse_primary<-if(length(a)==9L&&nzchar(a[9]))normalizePath(a[9])else NULL
runtime <- jsonlite::fromJSON(normalizePath(a[1])); policy <- normalizePath(a[2]); gb <- normalizePath(a[3]); annfile <- normalizePath(a[4]); hmmfile <- normalizePath(a[5]); cache <- normalizePath(a[6]); out <- normalizePath(a[7]); cpu <- as.integer(a[8]); started <- proc.time()[['elapsed']]
stopifnot(dir.exists(out), !file.exists(file.path(out,'completion.json')))
.libPaths(c(unlist(runtime$r_libraries), .libPaths()))
blast_executable <- normalizePath(runtime$blastp)
stopifnot(nzchar(Sys.which("blastp")), identical(normalizePath(Sys.which("blastp")), blast_executable))
metadata_source <- file.path(policy,"native_metadata.R"); sys.source(metadata_source,envir=.GlobalEnv)
Sys.setenv(DNMB_CACHE_ROOT=cache, DEFENSEVIZ_CACHE_DIR=file.path(cache,'db_modules/rebasefinder/cache'), XDG_CACHE_HOME=cache, TMPDIR=file.path(out,'tmp'))
dir.create(file.path(out,'tmp'),showWarnings=FALSE)
suppressPackageStartupMessages(library(dplyr))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
wt <- function(x,n){x<-as.data.frame(x);for(k in names(x))if(is.list(x[[k]]))x[[k]]<-vapply(x[[k]],function(z)paste(z,collapse=';'),character(1));write.table(x,file.path(out,n),sep='\t',row.names=FALSE,quote=TRUE,na='NA')}
json <- function(x,n)jsonlite::write_json(x,file.path(out,n),auto_unbox=TRUE,pretty=TRUE,na='null')
sha <- function(p) digest::digest(file=p,algo='sha256')
source_paths <- metadata_source
dv <- new.env(parent=asNamespace('DefenseViz'))
for(n in c('01_data_loading.R','01_hmm_nomination.R','02_mtase_detection.R','03_rease_detection.R','04_operon_analysis.R','06_rebase_comparison.R','07_scoring.R','08_pipeline.R')){p<-file.path(policy,'source/DefenseViz',n);sys.source(p,envir=dv);source_paths<-c(source_paths,p)}
dn <- new.env(parent=asNamespace('DNMB'))
for(n in c('db_module_rebasefinder.R','rebasefinder_motifs.R','rebasefinder_metadata.R','rebasefinder_hmm_nomination.R')){p<-file.path(policy,'source/DNMB',n);sys.source(p,envir=dn);source_paths<-c(source_paths,p)}
for(n in c('rebasefinder_multifactor_score.R','rebasefinder_curation.R')){p<-file.path(policy,'source',n);sys.source(p,envir=dn);source_paths<-c(source_paths,p)}
source_hashes<-setNames(as.list(vapply(source_paths,sha,character(1))),source_paths)
# Export the exact unchanged native score inputs for downstream selected-policy evaluation.
original_score <- dn$.dnmb_rebasefinder_multifactor_score
dn$.dnmb_rebasefinder_multifactor_score <- function(strong,moderate,profile,annotation,motif,context,structure,alternative_annotation,conflict,partial,policy=dn$.dnmb_rebasefinder_score_policy()) {
 z <- original_score(strong,moderate,profile,annotation,motif,context,structure,alternative_annotation,conflict,partial,policy)
 values <- list(strong=strong,moderate=moderate,profile=profile,annotation=annotation,motif=motif,context=context,structure=structure,alternative_annotation=alternative_annotation,conflict=conflict,partial=partial)
 for (k in names(values)) z[[paste0('factor_',k)]] <- values[[k]] %in% TRUE
 z
}

# Canonical coordinates/products/protein from original GBFF, independent donor annotations.
DNMB::Genbank_organizer(gb,save_output=FALSE)
g <- as.data.frame(get('genbank_table',envir=.GlobalEnv),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(g$locus_tag)); wt(g,'canonical_genes.tsv')
ann<-read.delim(annfile,check.names=FALSE,stringsAsFactors=FALSE,quote='"',comment.char='')
normalize_protein<-function(x){z<-toupper(gsub('[[:space:]]','',as.character(x)));sub('\\*+$','',z)}
extras<-setdiff(names(ann),names(g));extras<-extras[!grepl('^REBASEfinder_|^native_annotation_|^sequence_',extras)]
join<-vector('list',nrow(g));g$sequence_sha256<-NA_character_
for(i in seq_len(nrow(g))){
 ix<-which(ann$locus_tag==g$locus_tag[i]);seq<-normalize_protein(g$translation[i]);same<-ix[normalize_protein(ann$translation[ix])==seq];same<-same[!is.na(same)]
 for(k in extras){vals<-unique(as.character(ann[[k]][same]));vals<-vals[!is.na(vals)&nzchar(vals)&!vals %in% c('NA','-','None')];g[i,k]<-if(length(vals))paste(vals,collapse='; ')else NA_character_}
 if(!is.na(seq)&&nzchar(seq)&&seq!='NA')g$sequence_sha256[i]<-digest::digest(seq,algo='sha256',serialize=FALSE)
 join[[i]]<-data.frame(locus_tag=g$locus_tag[i],donor_rows=paste(same,collapse=';'),exact_protein_match=length(same)>0)
}
wt(do.call(rbind,join),'annotation_join.tsv');wt(g,'genes_enriched.tsv')
# Only safe serialized input is sent to the native collector.
input<-dn$.dnmb_rebasefinder_prepare_input(g,out)
genes<-dv$load_dnmb(input);stopifnot(!anyDuplicated(genes$locus_tag))
# Adapter-only input schema selection. Domain matching/mapping stays native.
.defenseviz_select_pfam_input <- function(data) {
  columns <- names(data)
  schemas <- list(
    legacy_pfam = "PFAMs",
    eggnog_pfam = "EggNOG_PFAMs",
    interpro_pfam_accession = c("Signature accession_Pfam", "Signature.accession_Pfam")
  )
  # read.delim repairs colliding source names with .1/.2 suffixes; do not
  # silently treat one repaired copy as the authoritative annotation field.
  repaired <- grepl("^(PFAMs|EggNOG_PFAMs|Signature[ .]accession_Pfam)\\.[0-9]+$", columns)
  if (any(repaired)) stop("Ambiguous repaired Pfam domain column: ", paste(columns[repaired], collapse=", "), call.=FALSE)
  for (schema in names(schemas)) {
    matches <- columns[columns %in% schemas[[schema]]]
    if (length(matches)>1L) stop("Ambiguous Pfam domain columns for ", schema, ": ", paste(matches, collapse=", "), call.=FALSE)
  }
  supported <- columns[columns %in% unlist(schemas, use.names=FALSE)]
  selected <- if ("PFAMs" %in% supported) "PFAMs" else if (length(supported)) supported[1L] else NA_character_
  source <- if (is.na(selected)) NA_character_ else names(schemas)[vapply(schemas, function(x) selected %in% x, logical(1))]
  values <- if (is.na(selected)) character() else as.character(data[[selected]])
  nonplaceholder <- !is.na(values) & nzchar(trimws(values)) & !trimws(values) %in% c("NA", "-", "None")
  list(column=selected, source=source, supported_columns_in_input_order=supported,
       ignored_pfam_columns=columns[grepl("pfam", columns, ignore.case=TRUE) & !columns %in% supported],
       selected_nonplaceholder_cells=sum(nonplaceholder),
       rule="Existing exact PFAMs preference; otherwise first supported domain column in native input order. Values do not determine selection or fallback.")
}

pfam_input<-.defenseviz_select_pfam_input(genes);json(pfam_input,'pfam_input_selection.json')
m<-dv$detect_methyltransferases(genes,pfam_col=pfam_input$column,product_col='product');m<-dv$.defenseviz_dnmb_candidate_compat(m)
r<-dv$detect_restriction_enzymes(genes,pfam_col=pfam_input$column,product_col='product',seq_col='translation',search_motifs=FALSE)
h<-dv$.defenseviz_read_hmm_candidates(hmmfile,genes,'locus_tag')
# Previously exported candidate identities remain search nominees only. Their
# saved type/role/curation values are deliberately absent from this interface.
priorfile<-file.path(dirname(annfile),'saved_rm_nominees.tsv');stopifnot(file.exists(priorfile))
prior<-read.delim(priorfile,check.names=FALSE,stringsAsFactors=FALSE,quote='"',comment.char='')
required<-c('sample_id','locus_tag','source_protein_sha256','normalized_protein_sha256','aa_length','source_refs')
stopifnot(identical(names(prior),required),!anyDuplicated(prior$locus_tag))
if(nrow(prior)){
 stopifnot(length(unique(prior$sample_id))==1L)
 gi<-match(prior$locus_tag,genes$locus_tag);stopifnot(!anyNA(gi))
 for(i in seq_len(nrow(prior))){
  seq<-normalize_protein(genes$translation[gi[i]])
  stopifnot(nchar(seq)==prior$aa_length[i],digest::digest(seq,algo='sha256',serialize=FALSE)==prior$normalized_protein_sha256[i])
  donors<-which(ann$locus_tag==prior$locus_tag[i]);raw<-toupper(gsub('[[:space:]]','',as.character(ann$translation[donors])))
  stopifnot(any(vapply(raw,function(s)digest::digest(s,algo='sha256',serialize=FALSE),character(1))==prior$source_protein_sha256[i]))
 }
}
wt(prior,'saved_nomination_input.tsv')
legacy<-unique(c(m$locus_tag,r$locus_tag));ids<-unique(c(legacy,h$genes$locus_tag,prior$locus_tag))
hnew<-h$genes[!h$genes$locus_tag %in% legacy,,drop=FALSE]
pnew<-genes[genes$locus_tag %in% prior$locus_tag & !genes$locus_tag %in% c(legacy,h$genes$locus_tag),,drop=FALSE]
if(nrow(h$summary))h$summary$hmm_new_candidate<-!h$summary$locus_tag %in% legacy
wt(m,'mtase_candidates.tsv');wt(r,'rease_candidates.tsv');wt(h$summary,'hmm_summary.tsv')
wt(data.frame(locus_tag=genes$locus_tag,annotation_nominee=genes$locus_tag %in% legacy,hmm_nominee=genes$locus_tag %in% h$genes$locus_tag,prior_nominee=genes$locus_tag %in% prior$locus_tag,native_nominee=genes$locus_tag %in% ids),'nomination_ledger.tsv')
refdir<-file.path(cache,'db_modules/rebasefinder/cache');refpath<-file.path(refdir,'rebase_data.rds');db<-file.path(refdir,'rebase_db.fasta')
stopifnot(all(file.exists(c(refpath,db,paste0(db,c('.phr','.pin','.psq'))))))
ref<-readRDS(refpath)
# Project only typed metadata; exact FASTA IDs/proteins bind all BLAST joins.
projection <- .dnmb_project_reference(ref); ref <- projection$reference
idxmap <- .dnmb_reference_index(ref,db)
changed <- with(projection$audit, is.na(old_family)!=is.na(family) | is.na(old_role)!=is.na(role) | (!is.na(old_family)&!is.na(family)&old_family!=family) | (!is.na(old_role)&!is.na(role)&old_role!=role) | is.na(previous_sequence_length)!=is.na(biological_sequence_length) | (!is.na(previous_sequence_length)&!is.na(biological_sequence_length)&previous_sequence_length!=biological_sequence_length))
wt(projection$audit[changed | projection$audit$type_conflict | projection$audit$role_conflict,,drop=FALSE],'metadata_projection_audit.tsv')
json(list(reference_rows=nrow(ref),encoded_rows=length(idxmap),projected_rows=sum(changed),type_conflicts=sum(projection$audit$type_conflict),role_conflicts=sum(projection$audit$role_conflict),repair='explicit typed description and numbered prefix; ambiguous metadata stays unknown',fasta_identity_verified=TRUE,terminal_markers_removed=sum(projection$audit$terminal_marker_removed),reference_lengths='biological amino acids; whitespace and terminal <> excluded',references_changed_on_disk=FALSE),'metadata_projection.json')
original_reference_rows<-dn$.dnmb_rebasefinder_reference_rows
# Native fallback remains for non-BLAST labels. Exact encoded subjects use exact resource rows.
dn$.dnmb_rebasefinder_reference_rows<-function(subject_ids,rebase_data,alignment_length=NULL,preferred_family=NULL,preferred_role=NULL){
 z<-original_reference_rows(subject_ids,rebase_data,alignment_length,preferred_family,preferred_role);ii<-unname(idxmap[as.character(subject_ids)]);ok<-!is.na(ii)
 z$reference_index[ok]<-ii[ok];z$reference_name[ok]<-ref$enzyme_name[ii[ok]];z$reference_aa_len[ok]<-ref$seq_length[ii[ok]];z$reference_family[ok]<-gsub('_',' ',ref$rm_type[ii[ok]],fixed=TRUE);z$reference_role[ok]<-ref$subunit[ii[ok]];z$reference_enz_type[ok]<-ref$enz_type[ii[ok]];z$reference_is_gold_standard[ok]<-ref$is_gold_standard[ii[ok]] %in% TRUE;z$reference_is_putative[ok]<-grepl('putative',ref$enz_type[ii[ok]],ignore.case=TRUE);z$reference_resolved_subject_id[ok]<-as.character(subject_ids)[ok];z
}
dn$.dnmb_rebasefinder_best_from_blast<-function(blast_tbl,rebase_data=NULL,min_identity=.10).dnmb_exact_best_blast(blast_tbl,ref,idxmap,min_identity)
dn$.dnmb_rebasefinder_cached_rebase_data<-function(cache_root=NULL)ref
dn$.dnmb_rebasefinder_rebase_fasta<-function(cache_root=NULL)db
# Fully local search; explicit exit code distinguishes zero HSPs from failure.
commands<-list();search_ledger<-list()
run_blast<-function(query_fasta,rebase_fasta,output_file,evalue=1e-5,num_threads=cpu,max_target_seqs=5,verbose=TRUE){
 before<-proc.time()[['elapsed']];args<-c('-query',query_fasta,'-db',db,'-out',output_file,'-evalue',as.character(evalue),'-num_threads',as.character(num_threads),'-max_target_seqs',as.character(max_target_seqs),'-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore')
 code<-system2(blast_executable,vapply(args,shQuote,character(1)),stdout=paste0(output_file,'.stdout'),stderr=paste0(output_file,'.stderr'));commands[[length(commands)+1L]]<<-list(argv=c(blast_executable,args),exit_code=code,elapsed_seconds=proc.time()[['elapsed']]-before,query_sha256=sha(query_fasta),output_sha256=if(file.exists(output_file))sha(output_file)else NULL)
 json(commands,'search_commands.json');if(code!=0L||!file.exists(output_file))stop('BLAST failed')
 search_ledger[[length(search_ledger)+1L]]<<-data.frame(locus_tag=sub('^>','',grep('^>',readLines(query_fasta),value=TRUE)),search_stage=if(length(commands)==1L)'primary'else 'supplemental',completed=TRUE,query_fasta=query_fasta)
 output_file
}
dn$.dnmb_rebasefinder_run_blastp<-run_blast
q<-genes[genes$locus_tag %in% ids & !is.na(genes$translation)&nzchar(genes$translation),,drop=FALSE]
qfile<-file.path(out,'query_sequences.fasta');get('.dnmb_write_protein_fasta',envir=dn,inherits=TRUE)(data.frame(protein_label=q$locus_tag,protein_seq=q$translation),qfile)
rawfile<-file.path(out,'blast_results.txt')
if(!is.null(reuse_primary)){
 oldq<-file.path(reuse_primary,'query_sequences.fasta');oldraw<-file.path(reuse_primary,'blast_results.txt');oldcmd<-jsonlite::read_json(file.path(reuse_primary,'search_commands.json'),simplifyVector=FALSE)
 stopifnot(file.exists(oldq),file.exists(oldraw),sha(oldq)==sha(qfile),length(oldcmd)>0,oldcmd[[1]]$exit_code==0)
 # Launcher verifies identical full reference hashes. Exact query hash binds the search scope.
 stopifnot(file.copy(oldraw,rawfile,overwrite=FALSE));commands[[1]]<-oldcmd[[1]];commands[[1]]$reuse_scope<-'exact_query_and_reference_full_database_search';commands[[1]]$source_commands_sha256<-sha(file.path(reuse_primary,'search_commands.json'));commands[[1]]$reused_output_sha256<-sha(oldraw);json(commands,'search_commands.json')
 search_ledger[[1]]<-data.frame(locus_tag=sub('^>','',grep('^>',readLines(qfile),value=TRUE)),search_stage='primary_reused_exact',completed=TRUE,query_fasta=qfile)
}else if(nrow(q))run_blast(qfile,db,rawfile)else file.create(rawfile)
parse<-function(path){if(!file.exists(path)||file.info(path)$size==0)return(data.frame());z<-dv$parse_blast_results(path,verbose=FALSE);if(nrow(z)&&anyNA(unname(idxmap[as.character(z$rebase_enzyme)])))stop('BLAST subject absent from verified reference FASTA');z}
raw<-parse(rawfile);filtered<-if(nrow(raw))dv$filter_blast_results(raw,min_identity=.10,min_length=50,max_evalue=.001,verbose=FALSE)else raw
# Single row per query only at the legacy comprehensive-table entrypoint; full HSP retained separately.
details<-NULL
if(nrow(filtered)){
 b<-filtered[order(filtered$query_id,-filtered$pct_identity,-filtered$bitscore,filtered$evalue),,drop=FALSE];b<-b[!duplicated(b$query_id),,drop=FALSE];rr<-dn$.dnmb_rebasefinder_reference_rows(b$rebase_enzyme,ref,alignment_length=b$length)
 details<-data.frame(query_id=b$query_id,blast_match=b$rebase_enzyme,blast_identity=b$pct_identity,blast_length=b$length,blast_evalue=b$evalue,blast_bitscore=b$bitscore,rm_type=rr$reference_family,subunit=rr$reference_role,rec_seq=ref$rec_seq[rr$reference_index],blast_qstart=b$qstart,blast_qend=b$qend,blast_sstart=b$sstart,blast_send=b$send)
}
operons<-dv$identify_rm_operons(m,r,genes,max_intervening=1L,max_distance=5000);if(nrow(operons))operons<-dv$predict_rm_type_from_operon(operons)
oi<-if(nrow(operons))dplyr::distinct(dplyr::select(operons,dplyr::all_of('locus_tag'),dplyr::any_of(c('operon_id','partner_locus_tag','partner_product'))))else NULL
ma<-dv$create_annotated_table(m,blast_details=details,operon_info=oi,dnmb_data=genes,id_col='locus_tag',component_type='MTase');ra<-dv$create_annotated_table(r,blast_details=details,operon_info=oi,dnmb_data=genes,id_col='locus_tag',component_type='REase');ha<-dv$create_annotated_table(hnew,blast_details=details,operon_info=oi,dnmb_data=genes,id_col='locus_tag',component_type='HMM nominee')
pa<-dv$create_annotated_table(pnew,blast_details=details,operon_info=oi,dnmb_data=genes,id_col='locus_tag',component_type='Saved R-M nominee')
cs<-genes[genes$locus_tag %in% ids,c('locus_tag','translation')];cm<-dv$detect_catalytic_motifs(cs,seq_col='translation',id_col='locus_tag',component_type='all');cms<-if(nrow(cm))dv$summarize_catalytic_motifs(cm,id_col='locus_tag')else data.frame()
comp<-dv$create_comprehensive_rm_table_v2(mtase_annotated=ma,rease_annotated=ra,trd_data=NULL,catalytic_summary=cms,filtered_ids=if(nrow(filtered))unique(filtered$query_id)else character(),id_col='locus_tag',nomination_annotated=dplyr::bind_rows(ha,pa))
x<-dn$.dnmb_rebasefinder_normalize_hits(comp);stopifnot(!anyDuplicated(x$query))
select<-function(x)dn$.dnmb_rebasefinder_select_primary_blast(x,genes,blast_tbl=filtered,rebase_data=ref,min_identity=.10,min_length=50,max_evalue=.001)
context<-function(x){x<-dn$.dnmb_rebasefinder_add_typeiii_context(x,genes,max_operon_gap=5000,max_intervening=1L,max_neighbors=3L);x<-dn$.dnmb_rebasefinder_add_typei_context(x,genes,max_operon_gap=5000,max_intervening=1L,max_neighbors=4L);dn$.dnmb_rebasefinder_add_typeii_context(x,genes,max_operon_gap=5000,max_intervening=1L,max_neighbors=3L)}
x<-select(x);x<-context(x);x<-dn$.dnmb_rebasefinder_add_typeiv_candidates(x,genes);x<-dn$.dnmb_rebasefinder_add_motif_evidence(x,genes)
supp<-dn$.dnmb_rebasefinder_supplemental_rebase_blast(x,genes,output_dir=out,blast_min_identity=.10,blast_min_length=50,cache_root=cache,cpu=cpu,download_rebase_if_missing=FALSE);x<-supp$hits
sraw<-parse(file.path(out,'supplemental_rebase_blast/supplemental_blast_results.txt'));if(nrow(sraw)){sf<-dv$filter_blast_results(sraw,min_identity=.10,min_length=50,max_evalue=.001,verbose=FALSE);filtered<-dplyr::bind_rows(filtered,sf)}
x<-select(x);x<-context(x)
# A later context pass can nominate additional translated neighbors. Close the
# searched inventory before curation; a successful but empty search is recorded.
translated <- !is.na(genes$translation)&nzchar(normalize_protein(genes$translation))&normalize_protein(genes$translation)!='NA'
closed <- .dnmb_close_context_search(x,
 candidate_ids=function(z)z$query[z$query %in% genes$locus_tag[translated]],
 searched_ids=function()if(length(search_ledger))unique(unlist(lapply(search_ledger,function(z)as.character(z$locus_tag))))else character(),
 search=function(missing,round){
  stage<-file.path(out,sprintf('context_search_%03d',round));dir.create(stage)
  qq<-genes[match(missing,genes$locus_tag),,drop=FALSE]
  qpath<-file.path(stage,'queries.faa'); rawpath<-file.path(stage,'blast_results.txt')
  get('.dnmb_write_protein_fasta',envir=dn,inherits=TRUE)(data.frame(protein_label=qq$locus_tag,protein_seq=qq$translation),qpath)
  written<-sub('^>','',grep('^>',readLines(qpath),value=TRUE))
  if(!setequal(written,missing))stop('Context search FASTA omitted a translated candidate')
  run_blast(qpath,db,rawpath,max_target_seqs=10)
  cr<-parse(rawpath);sraw<<-dplyr::bind_rows(sraw,cr)
  if(nrow(cr))filtered<<-dplyr::bind_rows(filtered,dv$filter_blast_results(cr,min_identity=.10,min_length=50,max_evalue=.001,verbose=FALSE))
 },refresh=function(z)context(select(z)),max_rounds=nrow(genes))
x<-closed$hits;json(closed$inventory,'search_inventory_closure.json')
x<-dn$.dnmb_rebasefinder_add_split_r_evidence(x,genes);x<-dn$.dnmb_rebasefinder_add_motif_evidence(x,genes);x<-dn$.dnmb_rebasefinder_add_sequence_partial_status(x,genes)
if(length(commands)&&any(vapply(commands,function(z)z$exit_code!=0L,logical(1))))stop('A native search failed; full execution cannot be claimed')
pre<-x;x<-dn$.dnmb_rebasefinder_curate_hits(x,genes);x<-dn$.dnmb_rebasefinder_attach_hmm_nomination(x,h$summary);stopifnot(!anyDuplicated(x$query))
# Curation scores above are native preparatory evidence. Python applies fixed selected policy afterward.
wt(comp,'native_comprehensive.tsv');wt(operons,'initial_native_operons.tsv');wt(cm,'native_collector_catalytic_motifs.tsv');wt(cms,'native_collector_catalytic_summary.tsv');wt(raw,'raw_blast.tsv');wt(sraw,'supplemental_raw_blast.tsv');wt(filtered,'filtered_blast.tsv');if(!is.null(details))wt(details,'native_blast_details.tsv');wt(pre,'before_curation_hypotheses.tsv');wt(x,'native_hypotheses.tsv');wt(x[x$curation_keep %in% TRUE,,drop=FALSE],'native_retained.tsv');wt(supp$status,'supplemental_status.tsv');if(length(search_ledger))wt(do.call(rbind,search_ledger),'search_ledger.tsv')
json(commands,'search_commands.json')
if(!length(search_ledger))wt(data.frame(locus_tag=character(),search_stage=character(),completed=logical(),query_fasta=character()),'search_ledger.tsv')
# Actual namespaces/DLLs must match the launcher-measured runtime, not only the
# configured library search path. Capture before the success marker is written.
runtime_packages <- lapply(sort(loadedNamespaces()),function(n)list(name=n,
 path=normalizePath(if(n=='base')file.path(R.home('library'),'base')else getNamespaceInfo(asNamespace(n),'path')),
 version=as.character(utils::packageVersion(n))))
runtime_dlls <- lapply(getLoadedDLLs(),function(x)as.character(x[['path']]))
json(list(R_version=R.version.string,R_home=R.home(),packages=runtime_packages,loaded_DLLs=runtime_dlls),'runtime_packages.json')
json(list(passed=TRUE,run_kind=if(is.null(reuse_primary))'fresh_native_RM_core_with_complete_nominated_query_search'else 'fresh_native_RM_core_using_exact_prior_full_database_search',rows=nrow(genes),native_nominees=length(ids),final_candidates=nrow(x),native_preparatory_retained=sum(x$curation_keep %in% TRUE),primary_queries=nrow(q),primary_HSP=nrow(raw),supplemental_HSP=nrow(sraw),full_reference_search=TRUE,genomewide_search=FALSE,selected_policy_applied=FALSE,structure_modeling='not_executed_by_core; supply exact existing structure evidence downstream',source_hashes=source_hashes,elapsed_seconds=proc.time()[['elapsed']]-started),'completion.json')

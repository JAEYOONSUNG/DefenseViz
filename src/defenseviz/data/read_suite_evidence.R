# Read-only projection of saved DNMB evidence. No packages or project code loaded.
# Each stdout line is a small JSON object; source row indices are one-based.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L || !file.exists(args[[1L]])) stop('one existing stage cache is required')
quote_json <- function(x) {
  chars <- utf8ToInt(enc2utf8(x))
  escaped <- vapply(chars,function(n) {
    if(n==34L) return('\\"')
    if(n==92L) return('\\\\')
    if(n<32L) return(sprintf('\\u%04x',n))
    intToUtf8(n)
  },character(1))
  paste0('"',paste0(escaped,collapse=''),'"')
}
json <- function(x) {
  if(is.null(x)) return('null')
  if(is.data.frame(x)) return(paste0('[',paste(vapply(seq_len(nrow(x)),function(i)json(as.list(x[i,,drop=FALSE])),character(1)),collapse=','),']'))
  if(inherits(x,'POSIXt')||inherits(x,'Date')) x<-as.character(x)
  if(is.list(x)) {
    if(!is.null(names(x))) return(paste0('{',paste(paste0(vapply(names(x),quote_json,character(1)),':',vapply(x,json,character(1))),collapse=','),'}'))
    return(paste0('[',paste(vapply(x,json,character(1)),collapse=','),']'))
  }
  if(length(x)!=1L) return(paste0('[',paste(vapply(seq_along(x),function(i)json(x[[i]]),character(1)),collapse=','),']'))
  if(is.na(x)) return('null')
  if(is.logical(x)) return(if(x)'true' else 'false')
  if(is.numeric(x)) return(if(is.finite(x))format(x,scientific=FALSE,trim=TRUE,digits=17) else 'null')
  quote_json(as.character(x))
}
emit <- function(x) cat(json(x),'\n',sep='')
x <- readRDS(args[[1L]])
if(!is.list(x)) stop('stage cache must be a list')
r <- x$module_results$REBASEfinder$module_result
if(is.null(r)) {emit(list(kind='metadata',schema_version=1L,status='module_absent'));quit(status=0L)}
if(!is.list(r)) stop('saved REBASEfinder module must be a list')
emit(list(kind='metadata',schema_version=1L,status='available',R_version=R.version.string,
          written_at=x$written_at,genbank_files=x$signature$genbank_signature$files,
          module_status=r$status,module_ok=r$ok))
p <- r$pipeline_result
if(is.null(p)) p<-list()
selected <- c('rm_comprehensive','blast_raw','blast_filtered','rebase_results','operon_data','catalytic_motifs','catalytic_summary')
tables <- list(hits=r$hits)
for(n in selected) tables[[paste0('pipeline_result/',n)]] <- p[[n]]
# Annotation text is already imported through the gene table. Keep fields that
# identify the exact protein and native R-M observations, avoiding huge repeats.
for(n in names(tables)) {
  v<-tables[[n]]
  if(!is.data.frame(v)) next
  emit(list(kind='table',table=n,row_count=nrow(v),columns=names(v)))
  keep<-names(v)[!grepl('^(Signature|Score)|nt_seq|^note$|^amino_acid$|^Mw$|^pI$|^anticodon$|^EC_number$|^gene_synonym$|^old_locus_tag$',names(v))]
  for(i in seq_len(nrow(v))) emit(list(kind='row',table=n,row=i,fields=as.list(v[i,keep,drop=FALSE])))
}

# Fixed, read-only REBASE metadata projection. No DNMB code is sourced.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("expected metadata RDS, requested subject IDs, output TSV")
x <- readRDS(args[[1]])
if (!is.data.frame(x) || !all(c("enzyme_name","sequence") %in% names(x)))
  stop("REBASE metadata must contain enzyme_name and sequence")
requested <- readLines(args[[2]],warn=FALSE)
seqs <- toupper(gsub("[[:space:]]", "", as.character(x$sequence)))
valid <- !is.na(x$sequence) & nchar(seqs)>50L
x <- x[valid,,drop=FALSE]; seqs <- seqs[valid]
names_original <- as.character(x$enzyme_name)
# Reproduce only the native FASTA identifier encoding, not biological inference.
ids <- gsub("\\s+","_",names_original,perl=TRUE)
ids <- gsub("[|;,()\\[\\]\\{\\}]","_",ids,perl=TRUE)
ids <- gsub("[^A-Za-z0-9_.-]","",ids)
dup <- duplicated(ids) | duplicated(ids,fromLast=TRUE)
ids[dup] <- paste0(ids[dup],"_",seq_len(nrow(x))[dup])
selected <- which(ids %in% requested | names_original %in% requested)
col <- function(name) if(name %in% names(x)) as.character(x[[name]][selected]) else rep(NA_character_,length(selected))
out <- data.frame(fasta_id=ids[selected], enzyme_name=names_original[selected],
  sequence=seqs[selected], seq_length=nchar(seqs[selected]),
  rm_type=col("rm_type"), subunit=col("subunit"), enz_type=col("enz_type"),
  metadata_row=which(valid)[selected],stringsAsFactors=FALSE)
write.table(out,args[[3]],sep="\t",quote=TRUE,row.names=FALSE,na="")

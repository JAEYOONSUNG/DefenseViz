# In-memory reference projection. Never changes the reference cache or thresholds.
.dnmb_project_reference <- function(ref) {
  canonical <- c('Type I', 'Type II', 'Type III', 'Type IV')
  parse_type <- function(text) {
    if (is.na(text)) return(NA_character_)
    if(grepl('Type[ _-]*[IVX]+[ ]*[/+-][ ]*(Type[ _-]*)?[IVX]+',text,perl=TRUE,ignore.case=TRUE)) return('conflict')
    tokens <- regmatches(text, gregexpr('\\bType[ _-]*[A-Za-z]+', text, perl=TRUE, ignore.case=TRUE))[[1]]
    if (!length(tokens)) return(NA_character_)
    roman <- sub('^Type[ _-]*', '', tokens, ignore.case=TRUE)
    roman <- sub('^([IVX]+).*$', '\\1', toupper(roman))
    if (any(!roman %in% c('I','II','III','IV')) || length(unique(roman)) != 1L) return('conflict')
    paste('Type', roman[1])
  }
  description <- as.character(ref$enz_type)
  kinds <- unique(description)
  explicit <- vapply(kinds, parse_type, character(1))[match(description,kinds)]
  family <- gsub('_',' ',as.character(ref$rm_type),fixed=TRUE)
  family[!family %in% canonical] <- NA_character_
  family[!is.na(explicit)] <- explicit[!is.na(explicit)]
  family[family == 'conflict' & !is.na(family)] <- NA_character_
  prefix <- rep(NA_character_,nrow(ref))
  named <- grepl('^[MRS][0-9]*[.]',ref$enzyme_name,ignore.case=TRUE)
  prefix[named] <- toupper(substr(ref$enzyme_name[named],1,1))
  is_m <- grepl('methyltransferase|methylase',description,ignore.case=TRUE)
  is_r <- grepl('restriction enzyme|endonuclease',description,ignore.case=TRUE)
  is_s <- grepl('specificity',description,ignore.case=TRUE)
  fused <- grepl('restriction enzyme[ /-]*(and[ ]*)?methyltransferase|restriction[ /-]methylation enzyme',description,ignore.case=TRUE)
  described <- rep(NA_character_,nrow(ref))
  described[is_m & !is_r & !is_s] <- 'M'
  described[is_r & !is_m & !is_s] <- 'R'
  described[is_s & !is_m & !is_r] <- 'S'
  conflicting <- ((is_m+is_r+is_s)>1L & !fused) | (!is.na(prefix) & !is.na(described) & prefix != described & !fused)
  role <- as.character(ref$subunit)
  role[!role %in% c('M','R','S','RM')] <- NA_character_
  role[!is.na(described)] <- described[!is.na(described)]
  role[!is.na(prefix)] <- prefix[!is.na(prefix)]
  role[fused] <- 'RM'
  role[conflicting] <- NA_character_
  audit <- data.frame(reference_index=seq_len(nrow(ref)), enzyme_name=ref$enzyme_name,
    old_family=ref$rm_type, family=gsub(' ','_',family,fixed=TRUE),
    old_role=ref$subunit, role=role, type_conflict=explicit %in% 'conflict', role_conflict=conflicting)
  ref$rm_type <- audit$family
  ref$subunit <- role
  clean <- sub('<>$','',toupper(gsub('[[:space:]]','',as.character(ref$sequence))))
  audit$raw_sequence_length <- nchar(ref$sequence)
  audit$previous_sequence_length <- ref$seq_length
  audit$biological_sequence_length <- ifelse(is.na(ref$sequence),NA_integer_,nchar(clean))
  audit$terminal_marker_removed <- !is.na(ref$sequence) & grepl('<>$',gsub('[[:space:]]','',as.character(ref$sequence)))
  ref$seq_length <- audit$biological_sequence_length
  list(reference=ref,audit=audit)
}

.dnmb_reference_index <- function(ref, fasta) {
  # Reproduce the pinned FASTA writer's RAW length filter before ID suffixing.
  valid <- !is.na(ref$sequence) & nchar(ref$sequence)>50L
  rows <- which(valid)
  encoded <- gsub('\\s+','_',as.character(ref$enzyme_name[rows]),perl=TRUE)
  encoded <- gsub('[|;,()\\[\\]\\{\\}]','_',encoded,perl=TRUE)
  encoded <- gsub('[^A-Za-z0-9_.-]','',encoded)
  duplicate <- duplicated(encoded)|duplicated(encoded,fromLast=TRUE)
  encoded[duplicate] <- paste0(encoded[duplicate],'_',seq_along(encoded)[duplicate])
  if (anyDuplicated(encoded) || any(!nzchar(encoded))) stop('Reference FASTA encoding is not unique')
  lines <- readLines(fasta,warn=FALSE)
  heads <- which(startsWith(lines,'>'))
  if (!length(heads) || heads[1]!=1L) stop('Reference FASTA has no valid first header')
  actual_ids <- substring(lines[heads],2L)
  if (!identical(actual_ids,encoded)) stop('Reference FASTA IDs do not match reference rows')
  boundaries <- c(heads,length(lines)+1L)
  sequences <- if(all(diff(boundaries)==2L)) lines[heads+1L] else vapply(seq_along(heads),function(i) {
    if(boundaries[i+1L]<=heads[i]+1L) return('')
    paste(lines[seq.int(heads[i]+1L,boundaries[i+1L]-1L)],collapse='')
  },character(1))
  normalize <- function(x) toupper(gsub('[[:space:]]','',as.character(x)))
  if (!identical(normalize(sequences),normalize(ref$sequence[rows]))) stop('Reference FASTA proteins do not match reference rows')
  setNames(rows,encoded)
}

.dnmb_exact_best_blast <- function(blast_tbl, ref, idxmap, min_identity=.10) {
  if (is.null(blast_tbl) || !nrow(blast_tbl)) return(data.frame())
  if (!'pct_identity' %in% names(blast_tbl) && 'pident' %in% names(blast_tbl)) blast_tbl$pct_identity <- blast_tbl$pident/100
  b <- as.data.frame(blast_tbl[!is.na(blast_tbl$pct_identity)&blast_tbl$pct_identity>=min_identity,,drop=FALSE])
  if(!nrow(b)) return(data.frame())
  subjects <- if('rebase_enzyme' %in% names(b)) as.character(b$rebase_enzyme) else as.character(b$subject_id)
  ii <- unname(idxmap[subjects])
  if(anyNA(ii)) stop('BLAST subject absent from verified reference FASTA')
  gold <- ref$is_gold_standard[ii] %in% TRUE
  rec <- as.character(ref$rec_seq[ii]); rec[is.na(rec)] <- ''
  # Preserve native identity/gold/recognition-length ranking, with exact subject rows.
  order <- order(b$query_id,-b$pct_identity,-as.integer(gold),-nchar(rec),seq_len(nrow(b)))
  keep <- order[!duplicated(b$query_id[order])]
  result <- data.frame(query_id=b$query_id[keep],best_rebase_match=subjects[keep],best_match_identity=b$pct_identity[keep],
    rm_type=gsub('_',' ',ref$rm_type[ii[keep]],fixed=TRUE),subunit=ref$subunit[ii[keep]],
    rec_seq=ref$rec_seq[ii[keep]],is_gold_standard=gold[keep])
  for(n in intersect(c('evalue','bitscore','length'),names(b))) result[[n]] <- b[[n]][keep]
  result
}

.dnmb_close_context_search <- function(hits,candidate_ids,searched_ids,search,refresh,max_rounds) {
  rounds <- 0L
  repeat {
    candidates <- unique(candidate_ids(hits))
    searched <- unique(searched_ids())
    missing <- setdiff(candidates,searched)
    if(!length(missing)) break
    rounds <- rounds+1L
    if(rounds>max_rounds) stop('Context search inventory did not converge')
    search(missing,rounds)
    if(length(setdiff(missing,searched_ids()))) stop('Context search did not complete every translated query')
    hits <- refresh(hits)
  }
  list(hits=hits,inventory=list(completed=TRUE,rounds=rounds,translated_candidates=length(candidates),
    searched_candidates=length(intersect(candidates,searched)),unsearched_candidates=missing))
}

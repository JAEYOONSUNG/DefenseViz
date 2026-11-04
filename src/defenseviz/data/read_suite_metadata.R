# Fixed data reader. No project code is sourced, no expressions are evaluated.
# stdout is JSON; only the named native RDS metadata files are opened.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[[1]])) stop("one existing suite directory is required")
quote_json <- function(x) {
  code <- utf8ToInt(enc2utf8(x))
  escaped <- vapply(code, function(n) {
    if (n == 34L) return('\\"')
    if (n == 92L) return('\\\\')
    if (n < 32L) return(sprintf('\\u%04x', n))
    intToUtf8(n)
  }, character(1))
  paste0('"', paste0(escaped, collapse = ''), '"')
}
json <- function(x) {
  if (is.null(x)) return('null')
  if (is.data.frame(x)) return(paste0('[', paste(vapply(seq_len(nrow(x)), function(i) json(as.list(x[i, , drop=FALSE])), character(1)), collapse=','), ']'))
  if (inherits(x, 'POSIXt') || inherits(x, 'Date')) x <- as.character(x)
  if (is.list(x)) {
    if (!is.null(names(x))) return(paste0('{', paste(paste0(vapply(names(x), quote_json, character(1)), ':', vapply(x, json, character(1))), collapse=','), '}'))
    return(paste0('[', paste(vapply(x, json, character(1)), collapse=','), ']'))
  }
  if (length(x) != 1L) return(paste0('[', paste(vapply(seq_along(x), function(i) json(x[[i]]), character(1)), collapse=','), ']'))
  if (is.na(x)) return('null')
  if (is.logical(x)) return(if (x) 'true' else 'false')
  if (is.numeric(x)) return(if (is.finite(x)) format(x, scientific=FALSE, trim=TRUE, digits=17) else 'null')
  quote_json(as.character(x))
}
read_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error=function(e) list(read_error=conditionMessage(e)))
}
dir <- args[[1]]
mods <- c('DefenseFinder','PADLOC','DefensePredictor','REBASEfinder')
out <- list(reader_schema_version=1L, R_version=R.version.string, sentinels=list(), cache=NULL, signatures=list())
for (m in mods) {
  md <- file.path(dir, paste0('dnmb_module_', tolower(m)))
  out$sentinels[m] <- list(read_safe(file.path(md, paste0('.dnmb_',tolower(m),'_complete'))))
  out$signatures[m] <- list(read_safe(file.path(md, '.input_signature.rds')))
}
c <- read_safe(file.path(dir, '.dnmb_module_stage_cache.rds'))
if (!is.null(c)) {
  if (!is.null(c$read_error)) out$cache <- c else {
    out$cache <- list(written_at=c$written_at, signature=c$signature, modules=list())
    for (m in mods) {
      b <- c$module_results[[m]]$module_result
      if (is.null(b)) next
      counts <- list()
      for (n in names(b)) if (is.data.frame(b[[n]])) counts[[n]] <- nrow(b[[n]])
      out$cache$modules[[m]] <- list(ok=b$ok, status=b$status, files=b$files, table_counts=counts)
    }
  }
}
cat(json(out), '\n', sep='')

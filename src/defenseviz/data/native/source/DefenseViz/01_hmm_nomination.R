# Explicit HMM nomination input. These fields nominate queries; they do not
# set family, subunit, confidence, BLAST thresholds, or curation outcomes.
.defenseviz_read_hmm_candidates <- function(hmm_candidates, dnmb_data, id_col = "locus_tag") {
  if (is.null(hmm_candidates)) return(list(domains = data.frame(), genes = data.frame(), summary = data.frame()))
  if (is.character(hmm_candidates) && length(hmm_candidates) == 1L) {
    if (!file.exists(hmm_candidates)) stop("HMM nomination TSV does not exist: ", hmm_candidates)
    nominees <- utils::read.delim(hmm_candidates, stringsAsFactors = FALSE,
      check.names = FALSE, quote = "", comment.char = "", colClasses = "character")
  } else if (is.data.frame(hmm_candidates)) {
    nominees <- as.data.frame(hmm_candidates, stringsAsFactors = FALSE)
  } else stop("hmm_candidates must be NULL, a TSV path, or a data frame")
  required <- c("locus_tag", "profile_name", "role_hint")
  missing <- setdiff(required, names(nominees))
  if (length(missing)) stop("HMM nomination missing required columns: ", paste(missing, collapse = ", "))
  if (!id_col %in% names(dnmb_data) || !"translation" %in% names(dnmb_data)) stop("Input genes lack locus identifier or translation")
  # Empty, header-only nominations are a valid no-hit result.
  if (!nrow(nominees)) return(list(domains = nominees, genes = dnmb_data[FALSE, , drop = FALSE], summary = data.frame()))
  identity_required <- c("protein_sha256", "sequence_length")
  missing_identity <- setdiff(identity_required, names(nominees))
  if (length(missing_identity)) stop("HMM nomination missing identity columns: ", paste(missing_identity, collapse = ", "))
  for (col in c(required, identity_required)) {
    nominees[[col]] <- trimws(as.character(nominees[[col]]))
    if (anyNA(nominees[[col]]) || any(!nzchar(nominees[[col]]))) stop("Empty HMM nomination field: ", col)
  }
  invalid_roles <- setdiff(nominees$role_hint, c("R", "M", "S", "unknown"))
  if (length(invalid_roles)) stop("Invalid HMM role_hint: ", paste(invalid_roles, collapse = ", "))
  ids <- as.character(dnmb_data[[id_col]])
  if (anyDuplicated(ids[!is.na(ids) & nzchar(ids)])) stop("Ambiguous duplicate input locus identifiers; explicit unique mapping required")
  idx <- match(nominees$locus_tag, ids)
  if (anyNA(idx)) stop("Unknown HMM nomination loci: ", paste(unique(nominees$locus_tag[is.na(idx)]), collapse = ", "))
  protein <- toupper(gsub("[[:space:]]", "", as.character(dnmb_data$translation[idx])))
  protein <- sub("\\*+$", "", protein)
  if (anyNA(protein) || any(!nzchar(protein))) stop("HMM nominee has no translated protein")
  if ("protein_sha256" %in% names(nominees)) {
    supplied <- tolower(as.character(nominees$protein_sha256))
    if (any(!grepl("^[0-9a-f]{64}$", supplied))) stop("Invalid HMM nominee protein SHA256 format")
    check <- rep(TRUE, length(supplied))
    actual <- vapply(protein, digest::digest, character(1), algo = "sha256", serialize = FALSE)
    if (any(check & (supplied != actual))) stop("HMM nominee protein SHA256 mismatch at: ", paste(unique(nominees$locus_tag[check & (supplied != actual)]), collapse = ", "))
  }
  if ("sequence_length" %in% names(nominees)) {
    raw <- as.character(nominees$sequence_length)
    check <- !is.na(raw) & nzchar(raw)
    len <- suppressWarnings(as.numeric(raw))
    if (any(check & (is.na(len) | len != nchar(protein)))) stop("HMM nominee sequence_length mismatch")
  }
  if (!"source" %in% names(nominees)) nominees$source <- "hmmsearch"
  # Retain every supplied evidence row (including duplicates) for provenance.
  # One query per native locus; sorting makes domain-order immaterial.
  loci <- unique(nominees$locus_tag)
  collapse <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")
  summary <- do.call(rbind, lapply(loci, function(locus) {
    x <- nominees[nominees$locus_tag == locus, , drop = FALSE]
    roles <- sort(unique(x$role_hint[x$role_hint != "unknown"]))
    data.frame(locus_tag = locus, hmm_nominated = TRUE,
      hmm_profiles = collapse(x$profile_name), hmm_role_hints = collapse(x$role_hint),
      hmm_role_hint = if (length(roles) == 1L && !"unknown" %in% x$role_hint) roles else NA_character_,
      hmm_role_ambiguous = length(roles) != 1L || "unknown" %in% x$role_hint,
      hmm_evidence_rows = nrow(x), hmm_nomination_sources = collapse(x$source),
      stringsAsFactors = FALSE)
  }))
  gene_rows <- dnmb_data[match(loci, ids), , drop = FALSE]
  gene_rows <- dplyr::left_join(gene_rows, summary, by = setNames("locus_tag", id_col))
  list(domains = nominees, genes = gene_rows, summary = summary)
}

# Same existing DNMB duplicate/cap behavior, expressed as an explicit local
# option rather than temporary mutation of the DefenseViz namespace.
.defenseviz_dnmb_candidate_compat <- function(result) {
  id_col <- intersect(c("locus_tag", "protein_id"), names(result))[1]
  if (!is.na(id_col) && nrow(result) > 0) result <- result[!duplicated(result[[id_col]]), , drop = FALSE]
  if ("detection_method" %in% names(result) && nrow(result) > 500L) {
    keep <- result$detection_method != "keyword_only" | result$mtase_confidence %in% c("high", "medium")
    if (sum(keep) > 0) result <- result[keep, , drop = FALSE]
  }
  result
}

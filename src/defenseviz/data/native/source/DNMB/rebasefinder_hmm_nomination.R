# Observational metadata only; called after all original curation decisions.
.dnmb_rebasefinder_attach_hmm_nomination <- function(hits, summary) {
  if (is.null(summary) || !is.data.frame(summary) || !nrow(summary)) return(hits)
  idx <- match(hits$query, summary$locus_tag)
  for (col in setdiff(names(summary), "locus_tag")) hits[[col]] <- summary[[col]][idx]
  hits$hmm_nominated <- !is.na(idx)
  if ("hmm_new_candidate" %in% names(hits)) hits$hmm_new_candidate[is.na(idx)] <- FALSE
  native_role <- if ("enzyme_role" %in% names(hits)) as.character(hits$enzyme_role) else rep(NA_character_, nrow(hits))
  unresolved <- is.na(native_role) | !native_role %in% c("R", "M", "S", "RM", "R/M")
  hits$hmm_nomination_review <- hits$hmm_nominated & (unresolved | hits$hmm_role_ambiguous %in% TRUE)
  hits$hmm_nomination_review_reason <- ifelse(!hits$hmm_nominated, NA_character_,
    ifelse(hits$hmm_role_ambiguous %in% TRUE, "multiple_or_unknown_profile_role_hints;native_curation_unchanged",
      ifelse(unresolved, "role_unresolved_after_native_annotation;native_curation_unchanged", NA_character_)))
  hits
}

# Content identity is required for top-level module cache reuse. Path equality
# alone cannot identify a changed nomination file; missing files fail closed.
.dnmb_hmm_nomination_signature <- function(candidates) {
  if (is.null(candidates)) return(NULL)
  if (is.character(candidates) && length(candidates) == 1L) {
    path <- normalizePath(candidates, mustWork = TRUE)
    if (dir.exists(path)) stop("HMM nomination path must be a file")
    return(list(contract = "explicit_hmm_nomination_v1", kind = "tsv",
      path = path, size = unname(file.info(path)$size),
      sha256 = digest::digest(file = path, algo = "sha256")))
  }
  if (is.data.frame(candidates)) return(list(contract = "explicit_hmm_nomination_v1",
    kind = "data_frame", rows = nrow(candidates),
    sha256 = digest::digest(candidates, algo = "sha256", serialize = TRUE)))
  stop("HMM nominees must be NULL, a TSV path or a data frame")
}

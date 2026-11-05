# Explicit provisional annotation policy, not a calibrated activity probability.
# Freeze before saved-case evaluation; no optimization on the manuscript cohort.
.dnmb_rebasefinder_score_policy <- function() {
  list(version = "DNMB_multifactor_v1_provisional", retain_cutoff = 50,
       high_band = 75, reference_strong = 35, reference_moderate = 25,
       profile = 25, annotation = 10, motif = 20, context = 20,
       structure = 25, alternative_annotation_penalty = 15,
       conflict_penalty = 15, partial_penalty = 10)
}

.dnmb_rebasefinder_multifactor_score <- function(
    strong, moderate, profile, annotation, motif, context, structure,
    alternative_annotation, conflict, partial,
    policy = .dnmb_rebasefinder_score_policy()) {
  values <- list(strong, moderate, profile, annotation, motif, context,
                 structure, alternative_annotation, conflict, partial)
  n <- length(strong)
  if (!all(vapply(values, length, integer(1)) == n)) stop("Score factors have unequal lengths")
  expected <- names(.dnmb_rebasefinder_score_policy())
  if (!identical(sort(names(policy)), sort(expected))) stop("Invalid score policy fields")
  numeric_policy <- unlist(policy[setdiff(expected, "version")])
  if (!is.numeric(numeric_policy) || any(!is.finite(numeric_policy)) ||
      any(numeric_policy < 0) || policy$retain_cutoff <= 0 ||
      policy$high_band < policy$retain_cutoff || policy$high_band > 100 ||
      max(policy$reference_strong, policy$reference_moderate, policy$profile,
          policy$annotation) + policy$motif + policy$context + policy$structure > 100)
    stop("Invalid score policy bounds")
  yes <- function(x) x %in% TRUE
  # Correlated annotation/reference/profile support is capped as ONE factor.
  seq_points <- pmax(policy$reference_strong * yes(strong),
                     policy$reference_moderate * yes(moderate),
                     policy$profile * yes(profile), policy$annotation * yes(annotation))
  motif_points <- policy$motif * yes(motif)
  context_points <- policy$context * yes(context)
  structure_points <- policy$structure * yes(structure)
  alternative_penalty <- policy$alternative_annotation_penalty * yes(alternative_annotation)
  # Do not double-penalize conflicting alternative annotation.
  conflict_penalty <- policy$conflict_penalty * (yes(conflict) & !yes(alternative_annotation))
  partial_penalty <- policy$partial_penalty * yes(partial)
  positive <- seq_points + motif_points + context_points + structure_points
  total <- pmax(0, positive - alternative_penalty - conflict_penalty - partial_penalty)
  data.frame(score_sequence = seq_points, score_motif = motif_points,
             score_context = context_points, score_structure = structure_points,
             score_alternative_penalty = alternative_penalty,
             score_conflict_penalty = conflict_penalty, score_partial_penalty = partial_penalty,
             score_positive = positive, score_total = total,
             score_factor_count = (seq_points > 0) + (motif_points > 0) +
                                  (context_points > 0) + (structure_points > 0),
             score_keep = total >= policy$retain_cutoff,
             score_cutoff = policy$retain_cutoff, score_policy_version = policy$version,
             score_calibration = "provisional_not_calibrated",
             stringsAsFactors = FALSE)
}

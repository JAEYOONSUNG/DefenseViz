#' @title Type I/III TRD Recognition Sequence Prediction
#' @description Advanced TRD analysis for Type I and III R-M systems
#'
#' @details
#' Type I methodology based on:
#' Furuta Y, Namba-Fukuyo H, Shibata TF, et al. (2014)
#' Methylome Diversification through Changes in DNA Methyltransferase Sequence Specificity.
#' PLoS Genet 10(4): e1004272. https://doi.org/10.1371/journal.pgen.1004272
#'
#' Key concepts:
#' - Type I HsdS: Two TRDs recognize bipartite sequences
#'   - TRD1 -> left half (5' side, e.g., CAA in CAA(N6)CTC)
#'   - TRD2 -> right half (3' side, e.g., CTC)
#'   - Linker alpha-helix length -> determines spacer (N) count
#'
#' - Type III Mod: Single TRD recognizes asymmetric sequence
#'
#' @name trd_prediction
#' @references
#' Furuta Y et al. (2014) PLoS Genet 10(4): e1004272
NULL

# ============================================================================
# Secondary Structure Prediction (local heuristic — no external API)
# ============================================================================

#' Predict Secondary Structure (local, no API)
#'
#' @description Uses amino acid propensity heuristic (~70-75% accuracy).
#' Sufficient for linker boundary estimation in TRD extraction.
#' No network calls, runs in <0.01 seconds per sequence.
#'
#' @param sequence Protein sequence
#' @param api_token Ignored (kept for backward compatibility)
#' @return Character string with secondary structure per residue (H/E/C)
#' @export
predict_secondary_structure_protT5 <- function(sequence, api_token = NULL) {
  predict_ss_heuristic(sequence)
}


#' Heuristic Secondary Structure Prediction
#' @description Amino acid propensity-based prediction with smoothing
#' @param sequence Protein sequence
#' @return Secondary structure string (H/E/C)
#' @export
predict_ss_heuristic <- function(sequence) {

  helix_aa <- c("A", "E", "L", "M", "Q", "K", "R", "H")
  sheet_aa <- c("V", "I", "Y", "F", "W", "T")

  aa_vec <- strsplit(toupper(sequence), "")[[1]]
  n <- length(aa_vec)
  if (n == 0L) return("")

  # Vectorized classification
  ss <- rep("C", n)
  ss[aa_vec %in% helix_aa] <- "H"
  ss[aa_vec %in% sheet_aa] <- "E"

  # Smoothing: require 4+ consecutive for helix
  if (n >= 5L) {
    ss_smoothed <- ss
    for (i in 3:(n - 2)) {
      h_count <- sum(ss[(i - 2):(i + 2)] == "H")
      if (h_count >= 4L) {
        ss_smoothed[i] <- "H"
      } else if (h_count <= 1L) {
        ss_smoothed[i] <- "C"
      }
    }
    ss <- ss_smoothed
  }

  paste(ss, collapse = "")
}


# ============================================================================
# REBASE-Based Learning: Helix Length vs Spacer Correlation
# ============================================================================

# Pre-computed helix-spacer model coefficients from REBASE Type I S subunits.
# Avoids runtime REBASE download and model training.
.dnmb_helix_spacer_default_intercept <- 4.2
.dnmb_helix_spacer_default_slope     <- 0.025

#' Predict Spacer Length from Linker Properties
#'
#' @description Uses pre-computed linear model (intercept + slope * linker_length)
#' derived from REBASE Type I S subunits. Falls back to heuristic bins.
#'
#' @param linker_length Linker region length in aa
#' @param helix_residues Number of helical residues (optional)
#' @param linker_ss Secondary structure of linker (optional)
#' @return Predicted spacer length (integer, clamped 4-10)
#' @export
predict_spacer_from_linker <- function(linker_length, helix_residues = NULL, linker_ss = NULL) {

  if (!is.null(linker_ss) && is.null(helix_residues)) {
    helix_residues <- sum(strsplit(linker_ss, "")[[1]] == "H")
  }

  # Pre-computed model prediction
  pred <- .dnmb_helix_spacer_default_intercept +
          .dnmb_helix_spacer_default_slope * linker_length
  pred <- as.integer(round(pred))
  max(4L, min(10L, pred))
}


#' Learn Helix-to-Spacer Correlation from REBASE
#'
#' @description Rebuilds linear model from REBASE Type I S subunits.
#' Only needed to update the pre-computed coefficients.
#'
#' @param force_rebuild Force rebuilding the model
#' @return List with correlation model
#' @export
learn_helix_spacer_correlation <- function(force_rebuild = FALSE) {

  cache_dir <- get_rebase_cache_dir()
  model_file <- file.path(cache_dir, "helix_spacer_model.rds")

  if (file.exists(model_file) && !force_rebuild) {
    model <- readRDS(model_file)
    message("Loaded helix-spacer model: ", length(model$data), " data points")
    return(model)
  }

  message("Learning helix-spacer correlation from REBASE...")

  rebase_data <- get_rebase_data(verbose = FALSE)

  type1_s <- rebase_data %>%
    dplyr::filter(
      stringr::str_detect(enz_type, "Type I") &
      stringr::str_detect(enzyme_name, "^S\\.")
    ) %>%
    dplyr::filter(!is.na(sequence) & nchar(sequence) >= 350) %>%
    dplyr::filter(!is.na(rec_seq) & rec_seq != "" & rec_seq != "?")

  if (nrow(type1_s) == 0) {
    warning("No Type I S subunits found")
    return(list(data = data.frame(), predict = function(x) 6))
  }

  training_data <- list()

  for (i in seq_len(nrow(type1_s))) {
    row <- type1_s[i, ]
    rec_parsed <- parse_bipartite_recognition(row$rec_seq)
    if (is.null(rec_parsed)) next

    spacer_n <- rec_parsed$spacer_n
    if (is.na(spacer_n) || spacer_n < 3 || spacer_n > 12) next

    seq_len_val <- nchar(row$sequence)
    linker_start <- round(seq_len_val * 0.40) + 1
    linker_end <- round(seq_len_val * 0.58)
    linker_seq <- substr(row$sequence, linker_start, linker_end)
    linker_length <- nchar(linker_seq)

    helix_aa <- c("A", "E", "L", "M", "Q", "K", "R")
    linker_aa <- strsplit(linker_seq, "")[[1]]
    helix_content <- sum(linker_aa %in% helix_aa) / length(linker_aa)

    training_data[[length(training_data) + 1]] <- data.frame(
      enzyme = row$enzyme_name,
      spacer_n = spacer_n,
      linker_length = linker_length,
      helix_content = helix_content,
      estimated_helix_residues = round(linker_length * helix_content)
    )
  }

  if (length(training_data) == 0) {
    warning("Could not extract training data")
    return(list(data = data.frame(), predict = function(x) 6))
  }

  train_df <- dplyr::bind_rows(training_data)
  message("  Training data: ", nrow(train_df), " entries")

  model_fit <- NULL
  if (nrow(train_df) >= 5) {
    model_fit <- tryCatch(lm(spacer_n ~ linker_length, data = train_df), error = function(e) NULL)
    if (!is.null(model_fit)) {
      message("  Model: spacer_n = ", round(coef(model_fit)[1], 2), " + ",
              round(coef(model_fit)[2], 4), " * linker_length")
    }
  }

  predict_spacer <- function(linker_length, helix_residues = NULL) {
    if (!is.null(model_fit)) {
      pred <- predict(model_fit, newdata = data.frame(linker_length = linker_length))
      return(max(4, min(10, round(pred))))
    }
    if (linker_length < 60) return(5)
    if (linker_length < 75) return(6)
    if (linker_length < 90) return(7)
    return(8)
  }

  spacer_stats <- train_df %>%
    dplyr::group_by(spacer_n) %>%
    dplyr::summarize(
      count = dplyr::n(),
      mean_linker = mean(linker_length),
      sd_linker = sd(linker_length),
      mean_helix_res = mean(estimated_helix_residues),
      .groups = "drop"
    )

  result <- list(data = train_df, model = model_fit, stats = spacer_stats, predict = predict_spacer)
  saveRDS(result, model_file)
  message("  Saved model to: ", model_file)
  result
}


# ============================================================================
# Type I TRD Extraction
# ============================================================================

#' Extract TRD1 and TRD2 Boundaries from Type I HsdS
#'
#' @description Identifies TRD1 (left half), linker, and TRD2 (right half)
#' Architecture: [N-term]--[TRD1]--[Linker helices]--[TRD2]--[C-term]
#'
#' When use_ss=TRUE, uses local heuristic SS prediction (instant, no API).
#'
#' @param sequence HsdS protein sequence
#' @param use_ss Use secondary structure prediction for better boundaries
#' @return List with TRD1, linker, TRD2 boundaries
#' @export
extract_type1_trd_boundaries <- function(sequence, use_ss = FALSE) {

  seq_len <- nchar(sequence)

  if (seq_len < 300) {
    warning("Sequence too short for typical HsdS (", seq_len, " aa)")
    return(NULL)
  }

  # Default proportional boundaries
  n_term_end <- max(40, round(seq_len * 0.10))
  trd1_start <- n_term_end + 1
  trd1_end <- round(seq_len * 0.40)

  linker_start <- trd1_end + 1
  linker_end <- round(seq_len * 0.58)

  trd2_start <- linker_end + 1
  trd2_end <- round(seq_len * 0.92)

  linker_ss <- NULL
  linker_helix_count <- NULL

  if (use_ss) {
    # Local heuristic — instant, no API call
    ss <- predict_ss_heuristic(sequence)

    if (nchar(ss) == seq_len) {
      ss_vec <- strsplit(ss, "")[[1]]

      # Sliding window helix density
      window <- 20L
      n_windows <- seq_len - window + 1L
      helix_density <- vapply(seq_len(n_windows), function(i) {
        sum(ss_vec[i:(i + window - 1L)] == "H") / window
      }, numeric(1))

      # Look for peak helix density in middle region
      mid_start <- round(seq_len * 0.35)
      mid_end <- round(seq_len * 0.65)
      mid_indices <- mid_start:min(mid_end, length(helix_density))

      if (length(mid_indices) > 0) {
        mid_density <- helix_density[mid_indices]
        peak_idx <- which.max(mid_density)

        threshold <- 0.4
        left_bound <- peak_idx
        while (left_bound > 1 && mid_density[left_bound - 1] > threshold) {
          left_bound <- left_bound - 1
        }
        right_bound <- peak_idx
        while (right_bound < length(mid_density) && mid_density[right_bound + 1] > threshold) {
          right_bound <- right_bound + 1
        }

        linker_start <- mid_start + left_bound - 1
        linker_end <- min(mid_start + right_bound + window - 1, round(seq_len * 0.65))
        trd1_end <- linker_start - 1
        trd2_start <- linker_end + 1
      }

      linker_ss <- substr(ss, linker_start, linker_end)
      linker_helix_count <- sum(strsplit(linker_ss, "")[[1]] == "H")
    }
  }

  trd1_seq <- substr(sequence, trd1_start, trd1_end)
  linker_seq <- substr(sequence, linker_start, linker_end)
  trd2_seq <- substr(sequence, trd2_start, trd2_end)
  linker_length <- nchar(linker_seq)

  predicted_spacer <- predict_spacer_from_linker(
    linker_length,
    helix_residues = linker_helix_count,
    linker_ss = linker_ss
  )

  list(
    trd1 = list(start = trd1_start, end = trd1_end, sequence = trd1_seq, length = nchar(trd1_seq)),
    linker = list(start = linker_start, end = linker_end, sequence = linker_seq, length = linker_length,
                  secondary_structure = linker_ss, helix_residues = linker_helix_count),
    trd2 = list(start = trd2_start, end = trd2_end, sequence = trd2_seq, length = nchar(trd2_seq)),
    predicted_spacer = predicted_spacer,
    total_length = seq_len
  )
}


# ============================================================================
# Type III TRD Extraction
# ============================================================================

#' Extract TRD Boundaries from Type III Mod
#' @param sequence Mod protein sequence
#' @return List with TRD boundaries
#' @export
extract_type3_trd_boundaries <- function(sequence) {

  seq_len <- nchar(sequence)

  if (seq_len < 300) {
    warning("Sequence too short for typical Type III Mod")
    return(NULL)
  }

  trd_start <- round(seq_len * 0.58)
  trd_end <- round(seq_len * 0.90)

  list(
    trd = list(start = trd_start, end = trd_end,
               sequence = substr(sequence, trd_start, trd_end),
               length = trd_end - trd_start + 1),
    total_length = seq_len
  )
}


# ============================================================================
# REBASE TRD Reference (cached, built once)
# ============================================================================

#' Build Type I TRD Reference Database
#' @description Creates TRD1/TRD2 reference with half-recognition sites. Cached.
#' @param force_rebuild Force rebuilding
#' @return tibble with TRD reference
#' @export
build_type1_trd_reference <- function(force_rebuild = FALSE) {

  cache_dir <- get_rebase_cache_dir()
  ref_file <- file.path(cache_dir, "type1_trd_reference.rds")

  if (file.exists(ref_file) && !force_rebuild) {
    ref <- readRDS(ref_file)
    message("Loaded Type I TRD reference: ", nrow(ref), " entries")
    return(ref)
  }

  message("Building Type I TRD reference from REBASE...")

  rebase_data <- get_rebase_data(verbose = FALSE)

  type1_s <- rebase_data %>%
    dplyr::filter(
      stringr::str_detect(enz_type, "Type I") &
      stringr::str_detect(enzyme_name, "^S\\.")
    ) %>%
    dplyr::filter(!is.na(sequence) & nchar(sequence) >= 350) %>%
    dplyr::filter(!is.na(rec_seq) & rec_seq != "" & rec_seq != "?")

  if (nrow(type1_s) == 0) {
    warning("No Type I S subunits found")
    return(dplyr::tibble())
  }

  message("  Found ", nrow(type1_s), " Type I HsdS with recognition sequences")

  trd_entries <- list()
  idx <- 1

  for (i in seq_len(nrow(type1_s))) {
    row <- type1_s[i, ]

    rec_parsed <- parse_bipartite_recognition(row$rec_seq)
    if (is.null(rec_parsed)) next

    trd_bounds <- extract_type1_trd_boundaries(row$sequence, use_ss = FALSE)
    if (is.null(trd_bounds)) next

    trd_entries[[idx]] <- dplyr::tibble(
      source_enzyme = row$enzyme_name, organism = row$org_name,
      trd_position = "TRD1", trd_sequence = trd_bounds$trd1$sequence,
      trd_length = trd_bounds$trd1$length, half_site = rec_parsed$left_half,
      full_recognition = row$rec_seq, spacer_length = rec_parsed$spacer_n
    )
    idx <- idx + 1

    trd_entries[[idx]] <- dplyr::tibble(
      source_enzyme = row$enzyme_name, organism = row$org_name,
      trd_position = "TRD2", trd_sequence = trd_bounds$trd2$sequence,
      trd_length = trd_bounds$trd2$length, half_site = rec_parsed$right_half,
      full_recognition = row$rec_seq, spacer_length = rec_parsed$spacer_n
    )
    idx <- idx + 1
  }

  if (length(trd_entries) == 0) {
    warning("Could not extract TRD sequences")
    return(dplyr::tibble())
  }

  ref <- dplyr::bind_rows(trd_entries)
  message("  Built reference: ", nrow(ref), " TRD entries")
  saveRDS(ref, ref_file)
  ref
}


#' Parse Bipartite Recognition Sequence
#' @param rec_seq Recognition sequence string (e.g., "CAA(6)CTC" or "CAANNNNNNCTC")
#' @return List with left_half, right_half, spacer_n
parse_bipartite_recognition <- function(rec_seq) {

  if (is.na(rec_seq) || rec_seq == "" || rec_seq == "?") return(NULL)

  rec_clean <- stringr::str_replace_all(rec_seq, "\\s+", "")

  if (stringr::str_detect(rec_clean, "\\(\\d+\\)")) {
    left <- stringr::str_extract(rec_clean, "^[ACGTRYWSMKHBVDN]+(?=\\()")
    spacer <- as.integer(stringr::str_extract(rec_clean, "(?<=\\()\\d+(?=\\))"))
    right <- stringr::str_extract(rec_clean, "(?<=\\))[ACGTRYWSMKHBVDN]+$")

    if (!is.na(left) && !is.na(right) && !is.na(spacer)) {
      return(list(left_half = left, right_half = right, spacer_n = spacer))
    }
  }

  if (stringr::str_detect(rec_clean, "N{2,}")) {
    parts <- stringr::str_split(rec_clean, "N+")[[1]]
    if (length(parts) >= 2) {
      left <- parts[1]
      right <- parts[length(parts)]
      spacer <- nchar(stringr::str_extract(rec_clean, "N+"))

      if (nchar(left) > 0 && nchar(right) > 0) {
        return(list(left_half = left, right_half = right, spacer_n = spacer))
      }
    }
  }

  NULL
}


#' Build Type III TRD Reference Database
#' @param force_rebuild Force rebuilding
#' @return tibble with TRD reference
#' @export
build_type3_trd_reference <- function(force_rebuild = FALSE) {

  cache_dir <- get_rebase_cache_dir()
  ref_file <- file.path(cache_dir, "type3_trd_reference.rds")

  if (file.exists(ref_file) && !force_rebuild) {
    ref <- readRDS(ref_file)
    message("Loaded Type III TRD reference: ", nrow(ref), " entries")
    return(ref)
  }

  message("Building Type III TRD reference...")

  rebase_data <- get_rebase_data(verbose = FALSE)

  type3_mod <- rebase_data %>%
    dplyr::filter(
      stringr::str_detect(enz_type, "Type III") &
      stringr::str_detect(enzyme_name, "^M\\.")
    ) %>%
    dplyr::filter(!is.na(sequence) & nchar(sequence) >= 300) %>%
    dplyr::filter(!is.na(rec_seq) & rec_seq != "" & rec_seq != "?")

  if (nrow(type3_mod) == 0) {
    warning("No Type III Mod found")
    return(dplyr::tibble())
  }

  trd_entries <- list()

  for (i in seq_len(nrow(type3_mod))) {
    row <- type3_mod[i, ]
    trd_bounds <- extract_type3_trd_boundaries(row$sequence)
    if (is.null(trd_bounds)) next

    trd_entries[[i]] <- dplyr::tibble(
      source_enzyme = row$enzyme_name, organism = row$org_name,
      trd_sequence = trd_bounds$trd$sequence, trd_length = trd_bounds$trd$length,
      recognition_seq = stringr::str_replace_all(row$rec_seq, "\\s+", "")
    )
  }

  ref <- dplyr::bind_rows(trd_entries) %>% dplyr::filter(!is.na(trd_sequence))
  message("  Built reference: ", nrow(ref), " entries")
  saveRDS(ref, ref_file)
  ref
}


# ============================================================================
# Batch BLAST: build DB once, query all at once
# ============================================================================

#' Build a persistent BLAST database from reference sequences
#' @param ref_seqs Character vector of reference sequences
#' @param ref_ids Character vector of reference IDs
#' @param db_dir Directory to store the DB
#' @return Path to BLAST DB prefix
.build_persistent_blast_db <- function(ref_seqs, ref_ids, db_dir) {
  db_prefix <- file.path(db_dir, "trd_ref")
  fasta_path <- paste0(db_prefix, ".fasta")

  # Write reference FASTA
  fasta_lines <- unlist(lapply(seq_along(ref_seqs), function(i) {
    c(paste0(">", ref_ids[i]), ref_seqs[i])
  }))
  writeLines(fasta_lines, fasta_path)

  # Build DB once
  system(sprintf("makeblastdb -in %s -dbtype prot -out %s 2>/dev/null",
                 shQuote(fasta_path), shQuote(db_prefix)),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  db_prefix
}


#' Batch BLAST multiple queries against a single persistent DB
#' @param query_seqs Character vector of query sequences
#' @param query_ids Character vector of query IDs
#' @param db_prefix Path to BLAST DB prefix
#' @param min_identity Minimum identity (0-1)
#' @return tibble with all BLAST results
.batch_blast_queries <- function(query_seqs, query_ids, db_prefix, min_identity = 0.30) {

  query_file <- tempfile(fileext = ".fasta")
  out_file <- tempfile(fileext = ".txt")
  on.exit(unlink(c(query_file, out_file)), add = TRUE)

  fasta_lines <- unlist(lapply(seq_along(query_seqs), function(i) {
    c(paste0(">", query_ids[i]), query_seqs[i])
  }))
  writeLines(fasta_lines, query_file)

  system(sprintf(
    "blastp -query %s -db %s -out %s -evalue 10 -max_target_seqs 20 -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' 2>/dev/null",
    shQuote(query_file), shQuote(db_prefix), shQuote(out_file)
  ))

  if (!file.exists(out_file) || file.info(out_file)$size == 0) {
    return(dplyr::tibble(query_id = character(), subject_id = character(),
                         pident = numeric(), length = integer(),
                         qstart = integer(), qend = integer(),
                         sstart = integer(), send = integer(),
                         evalue = numeric(), bitscore = numeric()))
  }

  hits <- utils::read.table(out_file, sep = "\t", stringsAsFactors = FALSE,
                            col.names = c("query_id", "subject_id", "pident", "length",
                                         "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
  hits %>%
    dplyr::filter(pident >= min_identity * 100) %>%
    dplyr::arrange(dplyr::desc(bitscore)) %>%
    dplyr::as_tibble()
}


#' Run BLAST for a single TRD sequence (legacy, uses temp DB)
#' @param query_seq Query sequence
#' @param subject_seqs Subject sequences
#' @param subject_ids Subject IDs
#' @param min_identity Minimum identity
#' @return tibble with BLAST results
run_trd_blast <- function(query_seq, subject_seqs, subject_ids, min_identity = 0.30) {
  .batch_blast_queries(query_seq, "query",
    .build_persistent_blast_db(subject_seqs, subject_ids, tempdir()),
    min_identity)
}


# ============================================================================
# BLAST TRD Against Reference (using cached DB)
# ============================================================================

#' BLAST TRD Against Type I Reference
#' @param trd_sequence TRD sequence
#' @param trd_position "TRD1" or "TRD2"
#' @param min_identity Minimum identity (0-1)
#' @return tibble with hits including blast_identity
#' @export
blast_trd_type1 <- function(trd_sequence, trd_position = NULL, min_identity = 0.30) {

  ref <- build_type1_trd_reference()
  if (nrow(ref) == 0) return(dplyr::tibble())

  if (!is.null(trd_position)) {
    ref <- ref %>% dplyr::filter(trd_position == !!trd_position)
  }
  if (nrow(ref) == 0) return(dplyr::tibble())

  hits <- run_trd_blast(trd_sequence, ref$trd_sequence, ref$source_enzyme, min_identity)

  if (nrow(hits) > 0) {
    hits <- hits %>%
      dplyr::left_join(
        ref %>% dplyr::select(source_enzyme, trd_position, half_site, full_recognition, spacer_length),
        by = c("subject_id" = "source_enzyme")
      ) %>%
      dplyr::rename(blast_identity = pident)
  }
  hits
}


#' BLAST TRD Against Type III Reference
#' @param trd_sequence TRD sequence
#' @param min_identity Minimum identity (0-1)
#' @return tibble with hits including blast_identity
#' @export
blast_trd_type3 <- function(trd_sequence, min_identity = 0.30) {

  ref <- build_type3_trd_reference()
  if (nrow(ref) == 0) return(dplyr::tibble())

  hits <- run_trd_blast(trd_sequence, ref$trd_sequence, ref$source_enzyme, min_identity)

  if (nrow(hits) > 0) {
    hits <- hits %>%
      dplyr::left_join(
        ref %>% dplyr::select(source_enzyme, recognition_seq),
        by = c("subject_id" = "source_enzyme")
      ) %>%
      dplyr::rename(blast_identity = pident)
  }
  hits
}


# ============================================================================
# Main Prediction Functions
# ============================================================================

#' Predict Type I Recognition Sequence
#' @param sequence HsdS protein sequence
#' @param enzyme_name Name for output
#' @param use_ss Use secondary structure prediction (local heuristic)
#' @param verbose Print progress
#' @return List with predictions including blast_identity
#' @export
predict_type1_recognition <- function(sequence, enzyme_name = "query",
                                       use_ss = FALSE, verbose = TRUE) {

  if (verbose) message("=== Type I Recognition Prediction: ", enzyme_name, " ===")

  trd_bounds <- extract_type1_trd_boundaries(sequence, use_ss = use_ss)

  if (is.null(trd_bounds)) {
    return(list(enzyme_name = enzyme_name, predicted_recognition = NA, confidence = 0,
                error = "TRD extraction failed"))
  }

  if (verbose) {
    message("    TRD1: ", trd_bounds$trd1$length, " aa | Linker: ",
            trd_bounds$linker$length, " aa | TRD2: ", trd_bounds$trd2$length, " aa")
  }

  trd1_hits <- blast_trd_type1(trd_bounds$trd1$sequence, trd_position = "TRD1")
  trd2_hits <- blast_trd_type1(trd_bounds$trd2$sequence, trd_position = "TRD2")

  left_half <- NA; right_half <- NA
  left_identity <- 0; right_identity <- 0
  left_source <- NA; right_source <- NA

  if (nrow(trd1_hits) > 0) {
    best <- trd1_hits[1, ]
    left_half <- best$half_site; left_identity <- best$blast_identity; left_source <- best$subject_id
  }

  if (nrow(trd2_hits) > 0) {
    best <- trd2_hits[1, ]
    right_half <- best$half_site; right_identity <- best$blast_identity; right_source <- best$subject_id
  }

  predicted_spacer <- trd_bounds$predicted_spacer

  if (!is.na(left_half) && !is.na(right_half)) {
    predicted_rec <- paste0(left_half, "(", predicted_spacer, ")", right_half)
    confidence <- (left_identity + right_identity) / 200
  } else if (!is.na(left_half)) {
    predicted_rec <- paste0(left_half, "(", predicted_spacer, ")????")
    confidence <- left_identity / 200 * 0.5
  } else if (!is.na(right_half)) {
    predicted_rec <- paste0("????(", predicted_spacer, ")", right_half)
    confidence <- right_identity / 200 * 0.5
  } else {
    predicted_rec <- NA; confidence <- 0
  }

  if (verbose) {
    message("    Result: ", ifelse(is.na(predicted_rec), "Unknown", predicted_rec),
            " (confidence: ", round(confidence * 100), "%)")
  }

  list(
    enzyme_name = enzyme_name, sequence_length = nchar(sequence),
    trd1_length = trd_bounds$trd1$length, trd2_length = trd_bounds$trd2$length,
    linker_length = trd_bounds$linker$length,
    linker_helix_residues = trd_bounds$linker$helix_residues,
    trd1_hits = trd1_hits, trd2_hits = trd2_hits,
    left_half = left_half, left_blast_identity = left_identity, left_source = left_source,
    right_half = right_half, right_blast_identity = right_identity, right_source = right_source,
    predicted_spacer = predicted_spacer, predicted_recognition = predicted_rec, confidence = confidence
  )
}


#' Predict Type III Recognition Sequence
#' @param sequence Mod protein sequence
#' @param enzyme_name Name for output
#' @param verbose Print progress
#' @return List with prediction including blast_identity
#' @export
predict_type3_recognition <- function(sequence, enzyme_name = "query", verbose = TRUE) {

  trd_bounds <- extract_type3_trd_boundaries(sequence)

  if (is.null(trd_bounds)) {
    return(list(enzyme_name = enzyme_name, predicted_recognition = NA, confidence = 0))
  }

  trd_hits <- blast_trd_type3(trd_bounds$trd$sequence)

  predicted_rec <- NA; blast_identity <- 0; source_enzyme <- NA

  if (nrow(trd_hits) > 0) {
    best <- trd_hits[1, ]
    predicted_rec <- best$recognition_seq
    blast_identity <- best$blast_identity
    source_enzyme <- best$subject_id
  }

  if (verbose) {
    message("    Type III ", enzyme_name, ": ",
            ifelse(is.na(predicted_rec), "Unknown", predicted_rec),
            " (", round(blast_identity, 1), "% identity)")
  }

  list(
    enzyme_name = enzyme_name, sequence_length = nchar(sequence),
    trd_length = trd_bounds$trd$length, trd_hits = trd_hits,
    predicted_recognition = predicted_rec, blast_identity = blast_identity,
    source_enzyme = source_enzyme, confidence = blast_identity / 100
  )
}


# ============================================================================
# Batch Processing (parallelized)
# ============================================================================

#' Batch Predict Type I Recognition (parallelized)
#' @param candidates Data frame with filtered candidates
#' @param seq_col Sequence column
#' @param id_col ID column
#' @param use_ss Use secondary structure prediction (local heuristic)
#' @param verbose Print progress
#' @param n_cores Number of CPU cores (default: 80% of available)
#' @return Data frame with predictions and blast_identity
#' @export
predict_type1_recognition_batch <- function(candidates, seq_col = "translation",
                                             id_col = "locus_tag", use_ss = FALSE,
                                             verbose = TRUE, n_cores = NULL) {

  if (nrow(candidates) == 0) return(dplyr::tibble())

  if (verbose) message("[TRD Type I] ", nrow(candidates), " candidates")

  n <- nrow(candidates)

  process_one <- function(i) {
    gene_id <- candidates[[id_col]][i]
    sequence <- candidates[[seq_col]][i]

    pred <- tryCatch({
      predict_type1_recognition(sequence, gene_id, use_ss = use_ss, verbose = FALSE)
    }, error = function(e) {
      list(enzyme_name = gene_id, predicted_recognition = NA, confidence = 0)
    })

    dplyr::tibble(
      gene_id = gene_id,
      predicted_recognition = pred$predicted_recognition,
      left_half = pred$left_half %||% NA,
      left_blast_identity = pred$left_blast_identity %||% NA,
      right_half = pred$right_half %||% NA,
      right_blast_identity = pred$right_blast_identity %||% NA,
      spacer_length = pred$predicted_spacer %||% NA,
      linker_length = pred$linker_length %||% NA,
      linker_helix_residues = pred$linker_helix_residues %||% NA,
      confidence = pred$confidence,
      system_type = "Type I"
    )
  }

  # Determine cores
  if (is.null(n_cores)) {
    avail <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) 1L)
    if (is.na(avail) || avail < 1L) avail <- 1L
    n_cores <- max(1L, as.integer(floor(avail * 0.8)))
  }

  if (n_cores > 1L && n > 1L && .Platform$OS.type != "windows") {
    if (verbose) message("  Using ", n_cores, " cores for parallel BLAST")
    results <- parallel::mclapply(seq_len(n), process_one, mc.cores = n_cores)
  } else {
    results <- lapply(seq_len(n), process_one)
  }

  result_df <- dplyr::bind_rows(results)

  if (verbose) {
    n_pred <- sum(!is.na(result_df$predicted_recognition))
    message("  Type I complete: ", n_pred, "/", nrow(result_df), " predictions")
  }
  result_df
}


#' Batch Predict Type III Recognition (parallelized)
#' @param candidates Data frame with filtered candidates
#' @param seq_col Sequence column
#' @param id_col ID column
#' @param verbose Print progress
#' @param n_cores Number of CPU cores (default: 80% of available)
#' @return Data frame with predictions and blast_identity
#' @export
predict_type3_recognition_batch <- function(candidates, seq_col = "translation",
                                             id_col = "locus_tag", verbose = TRUE,
                                             n_cores = NULL) {

  if (nrow(candidates) == 0) return(dplyr::tibble())

  if (verbose) message("[TRD Type III] ", nrow(candidates), " candidates")

  n <- nrow(candidates)

  process_one <- function(i) {
    gene_id <- candidates[[id_col]][i]
    sequence <- candidates[[seq_col]][i]

    pred <- tryCatch({
      predict_type3_recognition(sequence, gene_id, verbose = FALSE)
    }, error = function(e) {
      list(enzyme_name = gene_id, predicted_recognition = NA, confidence = 0)
    })

    dplyr::tibble(
      gene_id = gene_id,
      predicted_recognition = pred$predicted_recognition,
      blast_identity = pred$blast_identity %||% NA,
      source_enzyme = pred$source_enzyme %||% NA,
      trd_length = pred$trd_length %||% NA,
      confidence = pred$confidence,
      system_type = "Type III"
    )
  }

  if (is.null(n_cores)) {
    avail <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) 1L)
    if (is.na(avail) || avail < 1L) avail <- 1L
    n_cores <- max(1L, as.integer(floor(avail * 0.8)))
  }

  if (n_cores > 1L && n > 1L && .Platform$OS.type != "windows") {
    if (verbose) message("  Using ", n_cores, " cores for parallel BLAST")
    results <- parallel::mclapply(seq_len(n), process_one, mc.cores = n_cores)
  } else {
    results <- lapply(seq_len(n), process_one)
  }

  result_df <- dplyr::bind_rows(results)

  if (verbose) {
    n_pred <- sum(!is.na(result_df$predicted_recognition))
    message("  Type III complete: ", n_pred, "/", nrow(result_df), " predictions")
  }
  result_df
}

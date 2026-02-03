#' @title Structural Analysis Module
#' @description Functions for analyzing protein secondary structure (ProstT5) to estimate TRD parameters
#' @name structural_analysis
#'
#' @details
#' Type I S subunits have a characteristic domain architecture:
#' - Two TRD (Target Recognition Domain) regions at N and C termini
#' - Central coiled-coil alpha-helical spacer between TRDs
#' - Conserved N-repeat regions (NNNN patterns) flanking TRDs
#'
#' The alpha-helix spacer length determines the distance between two half-sites
#' in the bipartite recognition sequence (e.g., AACNNNNNNGTGC where N is the spacer)
#'
#' This module uses ProstT5 3-state secondary structure predictions (H=helix, E=strand, C=coil)
#' to estimate the alpha-helical spacer length and predict the N-repeat distance in recognition sequences.
NULL

#' Export Sequences for ProstT5 Prediction
#'
#' @description Exports S subunit sequences to FASTA for ProstT5 prediction
#'
#' @param s_subunits Data frame containing S subunit sequences (must have id_col and 'translation')
#' @param output_file Output FASTA path
#' @param id_col ID column name
#'
#' @return Path to output file
#' @export
export_sequences_for_prostt5 <- function(s_subunits, output_file, id_col = "locus_tag") {

  if (nrow(s_subunits) == 0) {
    warning("No sequences to export")
    return(NULL)
  }

  if (!"translation" %in% names(s_subunits)) {
    stop("Column 'translation' not found in data")
  }

  # Filter valid sequences
  valid_seqs <- s_subunits %>%
    dplyr::filter(!is.na(translation) & nchar(translation) > 50) %>%
    # Remove stop characters if present
    dplyr::mutate(translation = gsub("\\*", "", translation))

  if (nrow(valid_seqs) == 0) {
    warning("No valid protein sequences found")
    return(NULL)
  }

  # Write FASTA
  lines <- character()
  for (i in seq_len(nrow(valid_seqs))) {
    lines <- c(lines, paste0(">", valid_seqs[[id_col]][i]))
    lines <- c(lines, valid_seqs$translation[i])
  }

  writeLines(lines, output_file)
  message(sprintf("Exported %d sequences to %s", nrow(valid_seqs), output_file))
  return(output_file)
}


#' Parse ProstT5 Prediction Output
#'
#' @description Parses ProstT5 3-state output (FASTA-like format or CSV)
#' Format assumption: FASTA where sequence lines are the structure strings (H, E, C)
#'
#' @param prediction_file Path to ProstT5 output file
#'
#' @return Tibble with structure predictions
#' @export
parse_prostt5_output <- function(prediction_file) {

  if (!file.exists(prediction_file)) {
    stop("Prediction file not found: ", prediction_file)
  }

  lines <- readLines(prediction_file)
  header_idx <- grep("^>", lines)

  if (length(header_idx) == 0) {
    stop("Invalid ProstT5 output format (no headers)")
  }

  ids <- sub("^>", "", lines[header_idx])
  # Clean potential descriptions from IDs
  ids <- sub("\\s.*", "", ids)

  structures <- character(length(header_idx))

  for (i in seq_along(header_idx)) {
    start <- header_idx[i] + 1
    end <- if (i < length(header_idx)) header_idx[i+1] - 1 else length(lines)
    structures[i] <- paste(lines[start:end], collapse = "")
  }

  dplyr::tibble(
    locus_tag = ids,
    structure_seq = structures,
    seq_length = nchar(structures)
  )
}


#' Analyze S Subunit Structure for Alpha-Helix Spacer
#'
#' @description Analyzes S subunit structure to find the central alpha-helical spacer
#' based on long stretches of 'H' (Helix) prediction.
#'
#' Type I S subunits have characteristic architecture:
#' - N-terminal TRD (coil/strand dominant)
#' - Central alpha-helical coiled-coil spacer (long H stretches)
#' - C-terminal TRD (coil/strand dominant)
#'
#' @param structure_data Data from parse_prostt5_output()
#' @param min_helix_len Minimum length to consider as part of the spacer (default 15)
#'
#' @return Data with estimated spacer analysis
#' @export
analyze_s_subunit_structure <- function(structure_data, min_helix_len = 15) {

  # Function to analyze single sequence structure
  analyze_single_seq <- function(struct_seq, seq_id) {
    if (is.na(struct_seq) || nchar(struct_seq) == 0) {
      return(list(
        spacer_len_aa = 0,
        n_helices = 0,
        max_helix_aa = 0,
        central_helix_start = NA,
        central_helix_end = NA,
        trd1_region = NA,
        trd2_region = NA
      ))
    }

    # Run Length Encoding to find segments
    chars <- strsplit(struct_seq, "")[[1]]
    rle_res <- rle(chars)

    # Find all Helix segments ('H')
    helix_indices <- which(rle_res$values == "H")
    helix_lengths <- rle_res$lengths[helix_indices]

    # Calculate positions
    cumsum_lengths <- cumsum(rle_res$lengths)
    helix_starts <- c(1, cumsum_lengths[-length(cumsum_lengths)] + 1)[helix_indices]
    helix_ends <- cumsum_lengths[helix_indices]

    # Filter by minimum length
    valid_idx <- which(helix_lengths >= min_helix_len)

    if (length(valid_idx) == 0) {
      return(list(
        spacer_len_aa = 0,
        n_helices = 0,
        max_helix_aa = 0,
        central_helix_start = NA,
        central_helix_end = NA,
        trd1_region = NA,
        trd2_region = NA
      ))
    }

    valid_helices <- helix_lengths[valid_idx]
    valid_starts <- helix_starts[valid_idx]
    valid_ends <- helix_ends[valid_idx]

    # Find the central spacer region
    # The central alpha-helical spacer is typically the longest continuous helix region
    # in the middle of the protein (not at the very ends)
    seq_len <- nchar(struct_seq)
    middle_region_start <- seq_len * 0.2  # 20% from start
    middle_region_end <- seq_len * 0.8    # 80% from start

    # Filter helices in the central region
    central_idx <- which(valid_starts >= middle_region_start & valid_ends <= middle_region_end)

    if (length(central_idx) > 0) {
      # Sum up central helices as they might form a coiled-coil spacer
      central_helices <- valid_helices[central_idx]
      central_starts <- valid_starts[central_idx]
      central_ends <- valid_ends[central_idx]

      # The spacer might be fragmented by small loops
      # Calculate total helical content in central region
      total_spacer_aa <- sum(central_helices)

      # Find the dominant central helix
      max_central_idx <- which.max(central_helices)
      max_helix_aa <- central_helices[max_central_idx]
      central_helix_start <- central_starts[max_central_idx]
      central_helix_end <- central_ends[max_central_idx]

      # Define TRD regions based on central helix position
      trd1_region <- paste0("1-", central_helix_start - 1)
      trd2_region <- paste0(central_helix_end + 1, "-", seq_len)

    } else {
      # Fallback: use the longest helix overall
      total_spacer_aa <- sum(valid_helices)
      max_idx <- which.max(valid_helices)
      max_helix_aa <- valid_helices[max_idx]
      central_helix_start <- valid_starts[max_idx]
      central_helix_end <- valid_ends[max_idx]
      trd1_region <- NA
      trd2_region <- NA
    }

    list(
      spacer_len_aa = total_spacer_aa,
      n_helices = length(valid_idx),
      max_helix_aa = max_helix_aa,
      central_helix_start = central_helix_start,
      central_helix_end = central_helix_end,
      trd1_region = trd1_region,
      trd2_region = trd2_region
    )
  }

  # Apply analysis to each sequence
  results <- structure_data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      analysis = list(analyze_single_seq(structure_seq, locus_tag))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      spacer_len_aa = sapply(analysis, `[[`, "spacer_len_aa"),
      n_helices = sapply(analysis, `[[`, "n_helices"),
      max_helix_aa = sapply(analysis, `[[`, "max_helix_aa"),
      central_helix_start = sapply(analysis, `[[`, "central_helix_start"),
      central_helix_end = sapply(analysis, `[[`, "central_helix_end"),
      trd1_region = sapply(analysis, `[[`, "trd1_region"),
      trd2_region = sapply(analysis, `[[`, "trd2_region")
    ) %>%
    dplyr::select(-analysis)

  # Calculate physical distances and N-repeat prediction
  results <- estimate_trd_distance(results)

  return(results)
}


#' Estimate TRD Distance and N-repeat Length
#'
#' @description Converts alpha-helix length to estimated physical distance and
#' predicted N-repeat length in the recognition sequence.
#'
#' Assumptions:
#' - Alpha helix rise: ~1.5 Angstrom per residue
#' - DNA B-form rise: ~3.4 Angstrom per base pair
#' - Each amino acid in the spacer contributes to DNA-binding distance
#'
#' The N-repeat (spacer) in Type I recognition sequences corresponds to the
#' physical distance spanned by the alpha-helical coiled-coil domain.
#'
#' @param struct_data Data frame with `spacer_len_aa`
#'
#' @return Annotated data with distance estimates
#' @export
estimate_trd_distance <- function(struct_data) {

  # Conversion factors based on structural biology
  HELIX_RISE_PER_AA <- 1.5   # Angstroms per amino acid in alpha helix
  DNA_RISE_PER_BP <- 3.4     # Angstroms per base pair in B-DNA

  # Coiled-coil correction factor
  # Coiled-coils have ~1.5A rise per residue, similar to regular alpha helix
  # but the supercoil geometry affects the effective length
  COILED_COIL_FACTOR <- 0.9  # Slight reduction due to supercoiling

  struct_data %>%
    dplyr::mutate(
      # Calculate physical length of the spacer
      estimated_length_angstrom = spacer_len_aa * HELIX_RISE_PER_AA * COILED_COIL_FACTOR,

      # Estimate the DNA span in base pairs
      estimated_bp_span = estimated_length_angstrom / DNA_RISE_PER_BP,

      # Round to nearest integer for N-repeat prediction
      predicted_n_repeat = round(estimated_bp_span),

      # Confidence based on spacer length
      spacer_confidence = dplyr::case_when(
        spacer_len_aa >= 50 ~ "high",     # Clear central spacer
        spacer_len_aa >= 30 ~ "medium",   # Moderate spacer
        spacer_len_aa >= 15 ~ "low",      # Short spacer
        TRUE ~ "very_low"                  # Minimal helix
      ),

      # Round for readability
      estimated_bp_span = round(estimated_bp_span, 1)
    )
}


#' Analyze S Subunit with Pre-computed ProstT5 Results
#'
#' @description Complete analysis of S subunits using pre-computed ProstT5 predictions
#'
#' @param s_subunits Data frame with S subunit sequences
#' @param prostt5_file Path to ProstT5 output file
#' @param id_col ID column name
#'
#' @return Combined analysis results
#' @export
analyze_s_subunit_with_prostt5 <- function(s_subunits, prostt5_file, id_col = "locus_tag") {

  # Parse ProstT5 output
  structure_data <- parse_prostt5_output(prostt5_file)

  # Analyze structure
  analysis_results <- analyze_s_subunit_structure(structure_data)

  # Join with original S subunit data
  results <- s_subunits %>%
    dplyr::left_join(
      analysis_results,
      by = setNames("locus_tag", id_col)
    )

  # Add summary statistics
  message(sprintf("Analyzed %d S subunits", nrow(results)))
  message(sprintf("  Mean spacer length: %.1f aa", mean(results$spacer_len_aa, na.rm = TRUE)))
  message(sprintf("  Predicted N-repeat range: %d-%d bp",
                  min(results$predicted_n_repeat, na.rm = TRUE),
                  max(results$predicted_n_repeat, na.rm = TRUE)))

  return(results)
}


#' Predict Recognition Sequence Pattern from S Subunit Structure
#'
#' @description Generates a predicted recognition sequence pattern based on
#' TRD analysis and N-repeat prediction
#'
#' @param s_analysis S subunit analysis results
#' @param rebase_data REBASE data for TRD comparison (optional)
#'
#' @return Data with predicted recognition patterns
#' @export
predict_recognition_pattern <- function(s_analysis, rebase_data = NULL) {

  if (nrow(s_analysis) == 0) {
    return(dplyr::tibble())
  }

  s_analysis %>%
    dplyr::mutate(
      # Generate predicted pattern based on N-repeat length
      predicted_pattern = dplyr::case_when(
        is.na(predicted_n_repeat) ~ NA_character_,
        predicted_n_repeat <= 3 ~ paste0("(3-4bp)N(", predicted_n_repeat, ")(3-4bp)"),
        predicted_n_repeat <= 6 ~ paste0("(3-5bp)N(", predicted_n_repeat, ")(3-5bp)"),
        predicted_n_repeat <= 10 ~ paste0("(3-6bp)N(", predicted_n_repeat, ")(3-6bp)"),
        TRUE ~ paste0("(5-7bp)N(", predicted_n_repeat, ")(5-7bp)")
      ),

      # Pattern description
      pattern_description = dplyr::case_when(
        predicted_n_repeat <= 3 ~ "Short spacer - likely EcoRI-like specificity",
        predicted_n_repeat <= 6 ~ "Medium spacer - typical Type I pattern",
        predicted_n_repeat <= 10 ~ "Long spacer - extended recognition",
        TRUE ~ "Very long spacer - unusual pattern"
      )
    )
}


#' Visualize S Subunit Structure
#'
#' @description Creates a visualization of S subunit structure predictions
#'
#' @param s_analysis S subunit analysis results
#' @param output_file Output file path (optional)
#'
#' @return ggplot object
#' @export
plot_s_subunit_structure <- function(s_analysis, output_file = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  if (nrow(s_analysis) == 0) {
    message("No S subunit data to plot")
    return(NULL)
  }

  # Prepare plot data
  plot_data <- s_analysis %>%
    dplyr::filter(!is.na(spacer_len_aa)) %>%
    dplyr::mutate(
      label = paste0(locus_tag, "\n(N", predicted_n_repeat, ")")
    )

  # Create structure diagram
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(label, -spacer_len_aa))) +
    # TRD1 region (if available)
    ggplot2::geom_rect(
      ggplot2::aes(xmin = as.numeric(factor(label)) - 0.4,
                    xmax = as.numeric(factor(label)) - 0.1,
                    ymin = 0, ymax = 100),
      fill = "#E41A1C", alpha = 0.7
    ) +
    # Central spacer
    ggplot2::geom_rect(
      ggplot2::aes(xmin = as.numeric(factor(label)) - 0.1,
                    xmax = as.numeric(factor(label)) + 0.1,
                    ymin = 0, ymax = spacer_len_aa),
      fill = "#377EB8", alpha = 0.7
    ) +
    # TRD2 region (if available)
    ggplot2::geom_rect(
      ggplot2::aes(xmin = as.numeric(factor(label)) + 0.1,
                    xmax = as.numeric(factor(label)) + 0.4,
                    ymin = 0, ymax = 100),
      fill = "#4DAF4A", alpha = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = spacer_len_aa + 10, label = paste0(predicted_n_repeat, "N")),
      size = 3
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      x = "S Subunit (locus_tag)",
      y = "Alpha-helix spacer length (aa)",
      title = "Type I S Subunit Structure Analysis",
      subtitle = "Blue = Central alpha-helical spacer, Red/Green = TRD regions"
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 10, height = 6, dpi = 300)
    message("S subunit structure plot saved to: ", output_file)
  }

  return(p)
}


#' Generate ProstT5 Command for External Execution
#'
#' @description Generates the command to run ProstT5 prediction externally
#'
#' @param input_fasta Path to input FASTA file
#' @param output_file Path for output predictions
#'
#' @return Command string for running ProstT5
#' @export
generate_prostt5_command <- function(input_fasta, output_file) {

  # ProstT5 can be run via Python or command line
  # This generates example commands for both methods

  message("To run ProstT5 secondary structure prediction:")
  message("\n=== Option 1: Using Python API ===")
  message(sprintf('
from transformers import T5EncoderModel, T5Tokenizer
import torch
from Bio import SeqIO

# Load model
tokenizer = T5Tokenizer.from_pretrained("Rostlab/ProstT5", do_lower_case=False)
model = T5EncoderModel.from_pretrained("Rostlab/ProstT5")

# Read sequences
sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse("%s", "fasta")}

# Predict (simplified - see ProstT5 documentation for full code)
# Output format: FASTA with structure string (H/E/C) as sequence
', input_fasta))

  message("\n=== Option 2: Using ColabFold or web server ===")
  message("Upload your FASTA to: https://predictprotein.org/ or similar service")
  message(sprintf("Input file: %s", input_fasta))
  message(sprintf("Save 3-state predictions (H/E/C) to: %s", output_file))

  return(list(
    input = input_fasta,
    output = output_file
  ))
}

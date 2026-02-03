#' @title Restriction Enzyme Detection Module
#' @description Functions for identifying restriction endonucleases based on domains and catalytic motifs
#' @name rease_detection
NULL

#' Define Restriction Enzyme PFAM Domains
#'
#' @description Returns PFAM domains associated with restriction endonucleases
#'
#' @param type RE type: "all", "Type_I", "Type_II", "Type_III", "Type_IV"
#'
#' @return A tibble with PFAM IDs, names, and RE types
#' @export
get_rease_pfam_domains <- function(type = "all") {

  rease_domains <- dplyr::tibble(
    pfam_id = c(
      # Type I REase
      "PF04313", "PF13155", "PF00270",
      # Type II REase - various families
      "PF09019", "PF02806", "PF01420", "PF02234", "PF05869",
      "PF10544", "PF12950", "PF13333", "PF13392",
      # Type III REase
      "PF00176",
      # Type IV REase (McrBC, Mrr)
      "PF04471", "PF01844",
      # Common nuclease domains
      "PF09820", "PF00012", "PF02463", "PF13540",
      # HNH nuclease
      "PF01844", "PF13392", "PF14279",
      # PD-(D/E)XK nuclease superfamily
      "PF12705", "PF09491"
    ),
    pfam_name = c(
      "Helicase_RecQ_C", "HsdR_N", "DEAD",
      "Endonuc_TypeII", "McrB_C", "Endonuclease_NS", "PspGI", "MTS_N6_adenine",
      "Endonuclease_1", "Endonuclease_7", "Endonuc-II-G1", "Endonuclease_VI",
      "SNF2_N",
      "McrC", "HNH",
      "Endonuclease_1_2", "HSP70", "RecN_CSD", "Endonuc_Zn",
      "HNH", "Endonuclease_VI", "HNH_3",
      "PD_D_E_XK", "HSDR_N"
    ),
    re_type = c(
      "Type_I", "Type_I", "Type_I",
      "Type_II", "Type_II", "Type_II", "Type_II", "Type_II",
      "Type_II", "Type_II", "Type_II", "Type_II",
      "Type_III",
      "Type_IV", "HNH",
      "nuclease", "other", "other", "nuclease",
      "HNH", "HNH", "HNH",
      "PD-D/E-XK", "Type_I"
    ),
    description = c(
      "RecQ helicase C-terminal", "Type I R HsdR N-terminal", "DEAD/DEAH box helicase",
      "Type II endonuclease", "McrB nuclease", "Endonuclease NS", "PspGI endonuclease", "Associated MTase domain",
      "Endonuclease 1", "Endonuclease 7", "Type II endonuclease group 1", "Endonuclease VI",
      "SNF2 helicase N-terminal",
      "McrC family endonuclease", "HNH nuclease",
      "Endonuclease 1 type 2", "HSP70 (not RE)", "RecN coiled-coil", "Zinc finger nuclease",
      "HNH endonuclease", "Endonuclease VI", "HNH nuclease 3",
      "PD-(D/E)XK nuclease", "Type I HsdR N-terminal"
    )
  )

  if (type != "all") {
    rease_domains <- dplyr::filter(rease_domains, re_type == type)
  }

  return(rease_domains)
}


#' Define Catalytic Motifs for Nucleases
#'
#' @description Returns regex patterns for known nuclease catalytic motifs
#'
#' @return A tibble with motif patterns and descriptions
#' @export
get_nuclease_catalytic_motifs <- function() {

  dplyr::tibble(
    motif_name = c(
      "PD-D/E-XK",
      "HNH",
      "GIY-YIG",
      "PLD",
      "HALF_PIPE",
      "Type_II_metal_A",
      "Type_II_metal_B"
    ),
    pattern = c(
      # PD-(D/E)XK superfamily - most common REase motif
      "P[DE].{10,30}[DE]x{0,2}K",
      # HNH nuclease motif
      "H[A-Z]{1,2}N[A-Z]{10,40}H",
      # GIY-YIG nuclease
      "G[IV]Y.{10,30}Y[IV]G",
      # Phospholipase D-like
      "H[A-Z]K[A-Z]{3,6}D",
      # HALF_PIPE
      "HxD.{15,25}H.{10,20}D",
      # Type II metal binding site A
      "[DE].{0,2}[DE].{20,40}[DEQ].{0,2}K",
      # Type II metal binding site B
      "EX[A-Z]K"
    ),
    description = c(
      "PD-(D/E)XK nuclease superfamily",
      "HNH/EndoVII nuclease",
      "GIY-YIG nuclease family",
      "Phospholipase D family nuclease",
      "HALF_PIPE nuclease",
      "Type II restriction enzyme metal binding A",
      "Type II restriction enzyme EXK motif"
    ),
    confidence = c(
      "high", "high", "high", "medium", "medium", "medium", "medium"
    )
  )
}


#' Detect Restriction Enzymes by PFAM Domain
#'
#' @description Identifies restriction endonucleases using PFAM domains
#'
#' @param dnmb_data DNMB data frame
#' @param pfam_col Column containing PFAM annotations
#'
#' @return Tibble with identified REases
#' @export
detect_rease_by_pfam <- function(dnmb_data, pfam_col = NULL) {

  # Auto-detect PFAM column
  if (is.null(pfam_col)) {
    pfam_cols <- names(dnmb_data)[stringr::str_detect(names(dnmb_data), "(?i)pfam")]
    if ("PFAMs" %in% pfam_cols) {
      pfam_col <- "PFAMs"
    } else if (length(pfam_cols) > 0) {
      pfam_col <- pfam_cols[1]
    } else {
      stop("No PFAM column found")
    }
    message("Using PFAM column: ", pfam_col)
  }

  rease_domains <- get_rease_pfam_domains("all")

  # Build pattern
  pfam_pattern <- paste(rease_domains$pfam_id, collapse = "|")
  name_pattern <- paste(rease_domains$pfam_name, collapse = "|")
  combined_pattern <- paste0("(?i)(", pfam_pattern, "|", name_pattern, ")")

  rease_hits <- dnmb_data %>%
    dplyr::filter(!is.na(.data[[pfam_col]])) %>%
    dplyr::filter(stringr::str_detect(.data[[pfam_col]], combined_pattern))

  if (nrow(rease_hits) == 0) {
    message("No restriction enzymes found by PFAM")
    return(dplyr::tibble())
  }

  # Extract matched domains
  rease_hits <- rease_hits %>%
    dplyr::mutate(
      matched_re_pfam = stringr::str_extract_all(.data[[pfam_col]], combined_pattern),
      matched_re_pfam_str = sapply(matched_re_pfam, paste, collapse = ";")
    ) %>%
    dplyr::select(-matched_re_pfam) %>%
    # Classify RE type
    dplyr::mutate(
      re_type_pfam = dplyr::case_when(
        stringr::str_detect(matched_re_pfam_str, "(?i)HsdR|Helicase_RecQ|HSDR") ~ "Type_I",
        stringr::str_detect(matched_re_pfam_str, "(?i)Endonuc_TypeII|PspGI|Endonuc-II") ~ "Type_II",
        stringr::str_detect(matched_re_pfam_str, "(?i)SNF2|Type.?III") ~ "Type_III",
        stringr::str_detect(matched_re_pfam_str, "(?i)McrB|McrC|Mrr") ~ "Type_IV",
        stringr::str_detect(matched_re_pfam_str, "(?i)HNH") ~ "HNH_nuclease",
        stringr::str_detect(matched_re_pfam_str, "(?i)PD_D_E|Endonuclease") ~ "putative_nuclease",
        TRUE ~ "unknown"
      )
    )

  message(sprintf("Identified %d potential REases by PFAM", nrow(rease_hits)))
  return(rease_hits)
}


#' Detect Restriction Enzymes by Product Keywords
#'
#' @description Identifies REases using product annotation keywords
#'
#' @param dnmb_data DNMB data frame
#' @param product_col Product description column
#'
#' @return Tibble with identified REases
#' @export
detect_rease_by_keyword <- function(dnmb_data, product_col = "product") {

  rease_keywords <- c(
    # General
    "restriction", "endonuclease", "REase",
    # Type-specific subunits
    "hsdR", "hsdS", "res subunit", "restriction subunit",
    # Type specific
    "type I restriction", "type II restriction",
    "type III restriction", "type IV restriction",
    # Modified base restriction
    "McrB", "McrC", "Mrr", "modified cytosine",
    # Nuclease activity
    "DNA cleavage", "DNA cutting"
  )

  # Exclusion patterns - NOT R-M system (Cas, CRISPR, etc.)
  rease_exclude <- c(
    "CRISPR", "Cas[0-9]", "Cas protein", "Csm", "Cmr", "Csx", "Cpf1",
    "transposase", "integrase", "recombinase",
    "toxin", "antitoxin",
    "competence", "ComK"
  )

  pattern <- paste0("(?i)(", paste(rease_keywords, collapse = "|"), ")")
  exclude_pattern <- paste0("(?i)(", paste(rease_exclude, collapse = "|"), ")")

  rease_hits <- dnmb_data %>%
    dplyr::filter(!is.na(.data[[product_col]])) %>%
    dplyr::filter(stringr::str_detect(.data[[product_col]], pattern)) %>%
    # Exclude non-RM nucleases
    dplyr::filter(!stringr::str_detect(.data[[product_col]], exclude_pattern)) %>%
    # Classify by keyword
    dplyr::mutate(
      re_type_keyword = dplyr::case_when(
        stringr::str_detect(.data[[product_col]], "(?i)type.?I[^IV]|hsdR") ~ "Type_I",
        stringr::str_detect(.data[[product_col]], "(?i)type.?II") ~ "Type_II",
        stringr::str_detect(.data[[product_col]], "(?i)type.?III") ~ "Type_III",
        stringr::str_detect(.data[[product_col]], "(?i)type.?IV|Mcr|Mrr") ~ "Type_IV",
        TRUE ~ "putative"
      )
    )

  message(sprintf("Identified %d potential REases by keyword (after exclusion)", nrow(rease_hits)))
  return(rease_hits)
}


#' Detect Nuclease by Catalytic Motif
#'
#' @description Searches protein sequences for known nuclease catalytic motifs
#'
#' @param dnmb_data DNMB data frame
#' @param seq_col Column containing protein sequences (default: "translation")
#'
#' @return Tibble with catalytic motif hits
#' @export
detect_nuclease_by_motif <- function(dnmb_data, seq_col = "translation") {

  if (!seq_col %in% names(dnmb_data)) {
    stop("Sequence column '", seq_col, "' not found")
  }

  motifs <- get_nuclease_catalytic_motifs()

  # Initialize result
  motif_hits <- dnmb_data %>%
    dplyr::filter(!is.na(.data[[seq_col]])) %>%
    dplyr::mutate(
      catalytic_motifs = NA_character_,
      motif_count = 0L
    )

  # Search each motif
  for (i in seq_len(nrow(motifs))) {
    motif_name <- motifs$motif_name[i]
    pattern <- motifs$pattern[i]

    # Simplified motif search (basic regex)
    # For more accurate detection, would need Biostrings
    has_motif <- stringr::str_detect(
      motif_hits[[seq_col]],
      stringr::regex(pattern, ignore_case = FALSE)
    )

    # Update results
    motif_hits <- motif_hits %>%
      dplyr::mutate(
        catalytic_motifs = dplyr::if_else(
          has_motif,
          dplyr::if_else(
            is.na(catalytic_motifs),
            motif_name,
            paste(catalytic_motifs, motif_name, sep = ";")
          ),
          catalytic_motifs
        ),
        motif_count = motif_count + as.integer(has_motif)
      )
  }

  # Filter to only genes with motifs
  motif_hits <- motif_hits %>%
    dplyr::filter(!is.na(catalytic_motifs))

  message(sprintf("Found %d proteins with nuclease catalytic motifs", nrow(motif_hits)))
  return(motif_hits)
}


#' Combine REase Detection Results
#'
#' @description Integrates all REase detection methods
#'
#' @param dnmb_data DNMB data frame
#' @param pfam_col PFAM column (auto-detect if NULL)
#' @param product_col Product column
#' @param seq_col Sequence column for motif search
#' @param search_motifs Whether to search catalytic motifs (slower)
#'
#' @return Combined tibble with all REase candidates
#' @export
detect_restriction_enzymes <- function(dnmb_data,
                                        pfam_col = NULL,
                                        product_col = "product",
                                        seq_col = "translation",
                                        search_motifs = FALSE) {

  # PFAM-based detection
  pfam_results <- tryCatch(
    detect_rease_by_pfam(dnmb_data, pfam_col),
    error = function(e) {
      message("PFAM REase detection skipped: ", e$message)
      dplyr::tibble()
    }
  )

  # Keyword-based detection
  keyword_results <- tryCatch(
    detect_rease_by_keyword(dnmb_data, product_col),
    error = function(e) {
      message("Keyword REase detection skipped: ", e$message)
      dplyr::tibble()
    }
  )

  # Motif-based detection (optional)
  motif_results <- dplyr::tibble()
  if (search_motifs && seq_col %in% names(dnmb_data)) {
    motif_results <- tryCatch(
      detect_nuclease_by_motif(dnmb_data, seq_col),
      error = function(e) {
        message("Motif detection skipped: ", e$message)
        dplyr::tibble()
      }
    )
  }

  # Get ID column
  id_col <- intersect(c("locus_tag", "protein_id"), names(dnmb_data))[1]

  # Combine all results
  all_ids <- unique(c(
    if (nrow(pfam_results) > 0) pfam_results[[id_col]] else character(),
    if (nrow(keyword_results) > 0) keyword_results[[id_col]] else character(),
    if (nrow(motif_results) > 0) motif_results[[id_col]] else character()
  ))

  if (length(all_ids) == 0) {
    message("No restriction enzymes identified")
    return(dplyr::tibble())
  }

  # Build combined table
  combined <- dnmb_data %>%
    dplyr::filter(.data[[id_col]] %in% all_ids)

  # Add PFAM info
  if (nrow(pfam_results) > 0) {
    combined <- combined %>%
      dplyr::left_join(
        pfam_results %>% dplyr::select(dplyr::all_of(id_col),
                                        matched_re_pfam_str, re_type_pfam),
        by = id_col
      )
  }

  # Add keyword info
  if (nrow(keyword_results) > 0) {
    combined <- combined %>%
      dplyr::left_join(
        keyword_results %>% dplyr::select(dplyr::all_of(id_col), re_type_keyword),
        by = id_col
      )
  }

  # Add motif info
  if (nrow(motif_results) > 0) {
    combined <- combined %>%
      dplyr::left_join(
        motif_results %>% dplyr::select(dplyr::all_of(id_col),
                                         catalytic_motifs, motif_count),
        by = id_col
      )
  }

  # Consensus RE type and confidence
  if (!"catalytic_motifs" %in% names(combined)) {
    combined$catalytic_motifs <- NA_character_
  }

  combined <- combined %>%
    dplyr::mutate(
      re_type = dplyr::coalesce(re_type_pfam, re_type_keyword, "putative"),
      detection_evidence = paste(
        if ("re_type_pfam" %in% names(.)) dplyr::if_else(!is.na(re_type_pfam), "pfam", "") else "",
        if ("re_type_keyword" %in% names(.)) dplyr::if_else(!is.na(re_type_keyword), "keyword", "") else "",
        if ("catalytic_motifs" %in% names(.)) dplyr::if_else(!is.na(catalytic_motifs), "motif", "") else "",
        sep = "+"
      ),
      rease_confidence = dplyr::case_when(
        stringr::str_count(detection_evidence, "\\+") >= 2 ~ "high",
        !is.na(re_type_pfam) | !is.na(catalytic_motifs) ~ "medium",
        TRUE ~ "low"
      )
    ) %>%
    dplyr::mutate(
      detection_evidence = stringr::str_replace_all(detection_evidence, "^\\+|\\+$|\\+\\+", "+"),
      detection_evidence = stringr::str_replace(detection_evidence, "^\\+", "")
    )

  # Relocate key columns to front for visibility
  combined <- combined %>%
    dplyr::relocate(
      dplyr::any_of(c("locus_tag", "protein_id", "product",
                       "re_type", "rease_confidence", "detection_evidence",
                       "matched_re_pfam_str", "re_type_pfam", "re_type_keyword",
                       "catalytic_motifs", "motif_count"))
    )

  # Sort by confidence (high first), then by locus_tag
  id_col <- intersect(c("locus_tag", "protein_id"), names(combined))[1]
  combined <- combined %>%
    dplyr::mutate(
      .conf_order = factor(rease_confidence, levels = c("high", "medium", "low"))
    ) %>%
    dplyr::arrange(.conf_order, .data[[id_col]]) %>%
    dplyr::select(-.conf_order)

  message(sprintf("Total: %d restriction enzyme candidates identified", nrow(combined)))

  return(combined)
}

#' @title Methyltransferase Detection Module
#' @description Functions for identifying DNA methyltransferases based on domains and annotations
#' @name mtase_detection
NULL

#' Define DNA Methyltransferase PFAM Domains
#'
#' @description Returns a list of PFAM domains associated with DNA methyltransferases
#'
#' @param type Methyltransferase type: "all", "m6A", "m4C", "m5C"
#'
#' @return A tibble with PFAM IDs, names, and methylation types
#' @export
get_mtase_pfam_domains <- function(type = "all") {
  
  # Comprehensive MTase PFAM domains based on REBASE classification
  mtase_domains <- dplyr::tibble(
    pfam_id = c(
      # N6-adenine methyltransferases (m6A)
      "PF02086", "PF05869", "PF13659", "PF02384", 
      # N4-cytosine methyltransferases (m4C)
      "PF01555", 
      # C5-cytosine methyltransferases (m5C)
      "PF00145", "PF12564",
      # Type I MTase
      "PF02605",
      # Type II MTase
      "PF00440",
      # Dam methylase
      "PF05971",
      # Other DNA MTases
      "PF01035", "PF13651", "PF13847", "PF01021"
    ),
    pfam_name = c(
      "N6_Mtase", "MTS_N6_adenine", "Methyltransf_26", "N6_N4_Mtase",
      "N4_Cytosine_Mtase",
      "DNA_methylase", "Cytosine_DNA_MT",
      "Meth_Mtase_mod",
      "DNA_methylase",
      "Dam",
      "HNH_3", "Methyltransf_25", "Methyltransf_31", "S-AdoMet_synt_short"
    ),
    mtase_type = c(
      "m6A", "m6A", "m6A", "m6A",
      "m4C",
      "m5C", "m5C",
      "Type_I",
      "Type_II",
      "Dam",
      "other", "other", "other", "cofactor"
    ),
    description = c(
      "N6-adenine DNA methyltransferase", 
      "N-terminal domain of N6-adenine MTase",
      "Methyltransferase class 26",
      "N6 adenine and N4 cytosine MTase",
      "N4-cytosine DNA methyltransferase",
      "C-5 cytosine DNA methyltransferase",
      "Cytosine-specific DNA methyltransferase",
      "Type I restriction enzyme M subunit",
      "Type II DNA methyltransferase",
      "Dam methyltransferase",
      "HNH endonuclease (sometimes with MTase)",
      "Methyltransferase class 25",
      "Methyltransferase class 31",
      "S-adenosylmethionine synthetase"
    )
  )
  
  # Filter by type
  if (type != "all") {
    mtase_domains <- dplyr::filter(mtase_domains, mtase_type == type)
  }
  
  return(mtase_domains)
}


#' Detect Methyltransferases by PFAM Domain
#'
#' @description Identifies methyltransferases in DNMB data using PFAM annotations
#'
#' @param dnmb_data DNMB data frame
#' @param pfam_col Column containing PFAM annotations (auto-detect if NULL)
#' @param min_score Minimum domain score to consider (if score column available)
#'
#' @return Tibble with identified methyltransferases and their classification
#' @export
detect_mtase_by_pfam <- function(dnmb_data, pfam_col = NULL, min_score = NULL) {
  
  # Auto-detect PFAM column
  if (is.null(pfam_col)) {
    pfam_cols <- names(dnmb_data)[stringr::str_detect(names(dnmb_data), "(?i)pfam")]
    
    # Prefer eggNOG PFAM column
    if ("PFAMs" %in% pfam_cols) {
      pfam_col <- "PFAMs"
    } else if (length(pfam_cols) > 0) {
      pfam_col <- pfam_cols[1]
    } else {
      stop("No PFAM column found in data. Please specify pfam_col parameter.")
    }
    message("Using PFAM column: ", pfam_col)
  }
  
  # Get MTase PFAM domains
  mtase_domains <- get_mtase_pfam_domains("all")
  
  # Create search pattern
  pfam_pattern <- paste(mtase_domains$pfam_id, collapse = "|")
  pfam_name_pattern <- paste(mtase_domains$pfam_name, collapse = "|")
  combined_pattern <- paste0("(?i)(", pfam_pattern, "|", pfam_name_pattern, ")")
  
  # Filter genes with MTase domains
  mtase_hits <- dnmb_data %>%
    dplyr::filter(!is.na(.data[[pfam_col]])) %>%
    dplyr::filter(stringr::str_detect(.data[[pfam_col]], combined_pattern))
  
  if (nrow(mtase_hits) == 0) {
    message("No methyltransferases found by PFAM domain search")
    return(dplyr::tibble())
  }
  
  # Extract matched PFAM domains and annotate
  mtase_hits <- mtase_hits %>%
    dplyr::mutate(
      matched_pfam = stringr::str_extract_all(.data[[pfam_col]], combined_pattern),
      matched_pfam_str = sapply(matched_pfam, paste, collapse = ";")
    ) %>%
    dplyr::select(-matched_pfam)
  
  # Classify MTase type based on domains
  mtase_hits <- classify_mtase_type(mtase_hits, mtase_domains)
  
  message(sprintf("Identified %d potential methyltransferases by PFAM", nrow(mtase_hits)))
  
  return(mtase_hits)
}


#' Detect Methyltransferases by Product Keywords
#'
#' @description Identifies methyltransferases using product annotation keywords
#'
#' @param dnmb_data DNMB data frame
#' @param product_col Column containing product description
#'
#' @return Tibble with identified methyltransferases
#' @export
detect_mtase_by_keyword <- function(dnmb_data, product_col = "product") {
  
  # MTase keywords
  mtase_keywords <- c(
    # General
    "methyltransferase", "methylase", "MTase",
    # Specific types
    "N-6 adenine", "N6-adenine", "N6 adenine", "m6A",
    "N-4 cytosine", "N4-cytosine", "N4 cytosine", "m4C",
    "C-5 cytosine", "C5-cytosine", "cytosine-5", "m5C",
    # Type-specific
    "hsdM", "mod gene", "modification methylase",
    "dam methylase", "dcm methylase",
    # DNA modification
    "DNA modification", "DNA methylation"
  )
  
  pattern <- paste0("(?i)(", paste(mtase_keywords, collapse = "|"), ")")
  
  mtase_hits <- dnmb_data %>%
    dplyr::filter(!is.na(.data[[product_col]])) %>%
    dplyr::filter(stringr::str_detect(.data[[product_col]], pattern))
  
  # Infer methylation type from keywords
  mtase_hits <- mtase_hits %>%
    dplyr::mutate(
      mtase_type_keyword = dplyr::case_when(
        stringr::str_detect(.data[[product_col]], "(?i)N-?6|m6A|adenine") ~ "m6A",
        stringr::str_detect(.data[[product_col]], "(?i)N-?4|m4C") ~ "m4C",
        stringr::str_detect(.data[[product_col]], "(?i)C-?5|m5C|cytosine-5") ~ "m5C",
        stringr::str_detect(.data[[product_col]], "(?i)dam") ~ "Dam",
        stringr::str_detect(.data[[product_col]], "(?i)dcm") ~ "Dcm",
        TRUE ~ "unknown"
      )
    )
  
  message(sprintf("Identified %d potential methyltransferases by keyword", nrow(mtase_hits)))
  
  return(mtase_hits)
}


#' Classify Methyltransferase Type
#'
#' @description Internal function to classify MTase type based on PFAM domains
#'
#' @param mtase_data Data with MTase hits
#' @param mtase_domains Reference PFAM domain table
#'
#' @return Data with added mtase_type column
#' @keywords internal
classify_mtase_type <- function(mtase_data, mtase_domains) {
  
  mtase_data %>%
    dplyr::mutate(
      mtase_type_pfam = dplyr::case_when(
        # m6A
        stringr::str_detect(matched_pfam_str, "(?i)PF02086|PF05869|PF13659|PF02384|N6") ~ "m6A",
        # m4C
        stringr::str_detect(matched_pfam_str, "(?i)PF01555|N4") ~ "m4C",
        # m5C
        stringr::str_detect(matched_pfam_str, "(?i)PF00145|PF12564|DNA_methylase|Cytosine") ~ "m5C",
        # Dam
        stringr::str_detect(matched_pfam_str, "(?i)PF05971|Dam") ~ "Dam",
        # Type I
        stringr::str_detect(matched_pfam_str, "(?i)PF02605|Meth_Mtase") ~ "Type_I",
        TRUE ~ "unknown"
      )
    )
}


#' Combine MTase Detection Results
#'
#' @description Combines PFAM-based and keyword-based MTase detection
#'
#' @param dnmb_data DNMB data frame
#' @param pfam_col PFAM column name (auto-detect if NULL)
#' @param product_col Product description column name
#'
#' @return Combined tibble with all identified MTases
#' @export
detect_methyltransferases <- function(dnmb_data, 
                                       pfam_col = NULL, 
                                       product_col = "product") {
  
  # Detect by PFAM
  pfam_results <- tryCatch(
    detect_mtase_by_pfam(dnmb_data, pfam_col),
    error = function(e) {
      message("PFAM detection skipped: ", e$message)
      dplyr::tibble()
    }
  )
  
  # Detect by keyword
  keyword_results <- tryCatch(
    detect_mtase_by_keyword(dnmb_data, product_col),
    error = function(e) {
      message("Keyword detection skipped: ", e$message)
      dplyr::tibble()
    }
  )
  
  # Get ID column
  id_col <- intersect(c("locus_tag", "protein_id"), names(dnmb_data))[1]

  # Collect all unique IDs from both detection methods
  all_ids <- unique(c(
    if (nrow(pfam_results) > 0) pfam_results[[id_col]] else character(),
    if (nrow(keyword_results) > 0) keyword_results[[id_col]] else character()
  ))

  if (length(all_ids) == 0) {
    message("No methyltransferases identified")
    return(dplyr::tibble())
  }

  # Start from original data to preserve all columns
  combined <- dnmb_data %>%
    dplyr::filter(.data[[id_col]] %in% all_ids)

  # Add PFAM detection info
  if (nrow(pfam_results) > 0) {
    combined <- combined %>%
      dplyr::left_join(
        pfam_results %>% dplyr::select(dplyr::all_of(id_col), matched_pfam_str, mtase_type_pfam),
        by = id_col
      )
  } else {
    combined$matched_pfam_str <- NA_character_
    combined$mtase_type_pfam <- NA_character_
  }

  # Add keyword detection info
  if (nrow(keyword_results) > 0) {
    combined <- combined %>%
      dplyr::left_join(
        keyword_results %>% dplyr::select(dplyr::all_of(id_col), mtase_type_keyword),
        by = id_col
      )
  } else {
    combined$mtase_type_keyword <- NA_character_
  }

  # Consensus type and detection method
  combined <- combined %>%
    dplyr::mutate(
      mtase_type = dplyr::case_when(
        !is.na(mtase_type_pfam) & mtase_type_pfam != "unknown" ~ mtase_type_pfam,
        !is.na(mtase_type_keyword) & mtase_type_keyword != "unknown" ~ mtase_type_keyword,
        TRUE ~ "putative"
      ),
      detection_method = dplyr::case_when(
        !is.na(mtase_type_pfam) & !is.na(mtase_type_keyword) ~ "pfam+keyword",
        !is.na(mtase_type_pfam) ~ "pfam_only",
        !is.na(mtase_type_keyword) ~ "keyword_only",
        TRUE ~ "unknown"
      )
    )
  
  # Add confidence score based on detection method
  combined <- combined %>%
    dplyr::mutate(
      mtase_confidence = dplyr::case_when(
        detection_method == "pfam+keyword" ~ "high",
        detection_method == "pfam_only" ~ "medium",
        detection_method == "keyword_only" ~ "low",
        TRUE ~ "very_low"
      )
    )

  # Relocate key columns to front for visibility
  combined <- combined %>%
    dplyr::relocate(
      dplyr::any_of(c("locus_tag", "protein_id", "product",
                       "mtase_type", "mtase_confidence", "detection_method",
                       "matched_pfam_str", "mtase_type_pfam", "mtase_type_keyword"))
    )

  # Sort by confidence (high first), then by locus_tag
  id_col <- intersect(c("locus_tag", "protein_id"), names(combined))[1]
  combined <- combined %>%
    dplyr::mutate(
      .conf_order = factor(mtase_confidence, levels = c("high", "medium", "low", "very_low"))
    ) %>%
    dplyr::arrange(.conf_order, .data[[id_col]]) %>%
    dplyr::select(-.conf_order)

  message(sprintf("Total: %d methyltransferases identified", nrow(combined)))

  return(combined)
}

#' @title TRD Extraction Module
#' @description Functions for extracting and analyzing Target Recognition Domains (TRDs)
#' @name trd_extraction
NULL

#' Define TRD-Associated PFAM Domains
#'
#' @description Returns PFAM domains associated with Target Recognition Domains
#'
#' @return Tibble with TRD-related PFAM information
#' @export
get_trd_pfam_domains <- function() {
  
  dplyr::tibble(
    pfam_id = c(
      # Type I specificity (S) subunit TRDs
      "PF01420",  # Specificity subunit central region
      "PF04851",  # Type I restriction enzyme S subunit
      "PF18759",  # Type I RM specificity subunit
      # Type III TRD regions
      "PF01555",  # N4 cytosine MTase (has TRD region)
      "PF02086",  # N6 adenine MTase
      # HsdS-like
      "PF13156",  # S subunit N-terminal
      "PF08463"   # S subunit C-terminal
    ),
    pfam_name = c(
      "HsdS_N",
      "TER_restrict",
      "TypeI_RM_S",
      "N4_cytosine_Mtase",
      "N6_Mtase",
      "HsdS_N2",
      "HsdS_C"
    ),
    trd_type = c(
      "Type_I_S",
      "Type_I_S",
      "Type_I_S",
      "Type_III",
      "Type_III",
      "Type_I_S",
      "Type_I_S"
    ),
    description = c(
      "Type I specificity subunit N-terminal TRD",
      "Type I restriction enzyme S subunit",
      "Type I R-M specificity subunit",
      "N4-cytosine methyltransferase (contains TRD)",
      "N6-adenine methyltransferase (contains TRD)",
      "HsdS N-terminal domain",
      "HsdS C-terminal domain"
    )
  )
}


#' Define TRD Sequence Motifs
#'
#' @description Returns conserved sequence motifs found in TRD regions
#'
#' @return Tibble with TRD motif patterns
#' @export
get_trd_sequence_motifs <- function() {
  
  dplyr::tibble(
    motif_name = c(
      "TRD1_conserved",
      "TRD2_conserved",
      "TypeI_central",
      "TypeIII_trd"
    ),
    pattern = c(
      # Type I TRD1 conserved region
      "T[A-Z]{2,4}G.{10,30}[DE].{5,15}[RK]",
      # Type I TRD2 conserved region
      "[FY].{5,10}D.{2,5}N.{10,30}[WY]",
      # Type I central conserved region between TRDs
      "TGEL.{5,20}FQI",
      # Type III TRD-like region
      "[DE]P[A-Z]{3,6}F.{10,25}[RK].{3,8}[DE]"
    ),
    description = c(
      "Type I TRD1 conserved region",
      "Type I TRD2 conserved region",
      "Type I central conserved sequence",
      "Type III TRD region"
    ),
    system_type = c("Type_I", "Type_I", "Type_I", "Type_III")
  )
}


#' Extract TRD Regions by PFAM
#'
#' @description Identifies genes with TRD-associated PFAM domains
#'
#' @param dnmb_data DNMB data frame
#' @param pfam_col PFAM column name
#'
#' @return Tibble with TRD-containing genes
#' @export
extract_trd_by_pfam <- function(dnmb_data, pfam_col = NULL) {
  
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
  }
  
  trd_domains <- get_trd_pfam_domains()
  
  # Build pattern
  pattern <- paste0("(?i)(", 
                    paste(c(trd_domains$pfam_id, trd_domains$pfam_name), collapse = "|"),
                    ")")
  
  trd_hits <- dnmb_data %>%
    dplyr::filter(!is.na(.data[[pfam_col]])) %>%
    dplyr::filter(stringr::str_detect(.data[[pfam_col]], pattern))
  
  if (nrow(trd_hits) == 0) {
    message("No TRD domains found by PFAM")
    return(dplyr::tibble())
  }
  
  # Classify TRD type
  trd_hits <- trd_hits %>%
    dplyr::mutate(
      trd_pfam_match = stringr::str_extract_all(.data[[pfam_col]], pattern),
      trd_pfam_str = sapply(trd_pfam_match, paste, collapse = ";"),
      trd_type = dplyr::case_when(
        stringr::str_detect(.data[[pfam_col]], "(?i)HsdS|TypeI|TER_restrict") ~ "Type_I_specificity",
        stringr::str_detect(.data[[pfam_col]], "(?i)N4_cytosine|N6_Mtase") ~ "Type_III_TRD",
        TRUE ~ "putative_TRD"
      )
    ) %>%
    dplyr::select(-trd_pfam_match)
  
  message(sprintf("Found %d genes with TRD-associated PFAM domains", nrow(trd_hits)))
  
  return(trd_hits)
}


#' Extract TRD Regions by Keyword
#'
#' @description Identifies TRD-containing genes by product keywords
#'
#' @param dnmb_data DNMB data frame
#' @param product_col Product column name
#'
#' @return Tibble with TRD-containing genes
#' @export
extract_trd_by_keyword <- function(dnmb_data, product_col = "product") {

  # Strict TRD keywords - avoid generic terms like "specificity" alone

  trd_keywords <- c(
    "HsdS", "S subunit", "specificity subunit",
    "target recognition domain", "TRD",
    "type I.*specificity", "type III.*mod",
    "restriction.*specificity"
  )

  # Exclusion patterns to avoid false positives (competence, transport, etc.)
  exclude_keywords <- c(
    "competence", "ComK", "ComG", "ComE",
    "ABC transporter", "substrate.?specificity",
    "sigma factor", "transcription",
    "binding protein"
  )
  
  pattern <- paste0("(?i)(", paste(trd_keywords, collapse = "|"), ")")
  exclude_pattern <- paste0("(?i)(", paste(exclude_keywords, collapse = "|"), ")")

  trd_hits <- dnmb_data %>%
    dplyr::filter(!is.na(.data[[product_col]])) %>%
    dplyr::filter(stringr::str_detect(.data[[product_col]], pattern)) %>%
    # Exclude false positives
    dplyr::filter(!stringr::str_detect(.data[[product_col]], exclude_pattern))

  if (nrow(trd_hits) > 0) {
    trd_hits <- trd_hits %>%
      dplyr::mutate(
        trd_type_keyword = dplyr::case_when(
          stringr::str_detect(.data[[product_col]], "(?i)HsdS|type.?I.*S|specificity.?subunit") ~ "Type_I_specificity",
          stringr::str_detect(.data[[product_col]], "(?i)type.?III|mod") ~ "Type_III_TRD",
          TRUE ~ "putative_TRD"
        )
      )
  }

  message(sprintf("Found %d genes with TRD keywords (after exclusion filter)", nrow(trd_hits)))
  
  return(trd_hits)
}


#' Search TRD Sequence Motifs
#'
#' @description Searches protein sequences for conserved TRD motifs
#'
#' @param dnmb_data DNMB data with sequences
#' @param seq_col Sequence column name
#'
#' @return Tibble with motif search results
#' @export
search_trd_motifs <- function(dnmb_data, seq_col = "translation") {
  
  if (!seq_col %in% names(dnmb_data)) {
    stop("Sequence column not found: ", seq_col)
  }
  
  motifs <- get_trd_sequence_motifs()
  
  # Initialize
  trd_motif_hits <- dnmb_data %>%
    dplyr::filter(!is.na(.data[[seq_col]])) %>%
    dplyr::mutate(
      trd_motifs_found = NA_character_,
      trd_motif_count = 0L
    )
  
  # Search each motif
  for (i in seq_len(nrow(motifs))) {
    motif_name <- motifs$motif_name[i]
    pattern <- motifs$pattern[i]
    
    has_motif <- stringr::str_detect(
      trd_motif_hits[[seq_col]],
      stringr::regex(pattern, ignore_case = FALSE)
    )
    
    trd_motif_hits <- trd_motif_hits %>%
      dplyr::mutate(
        trd_motifs_found = dplyr::if_else(
          has_motif,
          dplyr::if_else(
            is.na(trd_motifs_found),
            motif_name,
            paste(trd_motifs_found, motif_name, sep = ";")
          ),
          trd_motifs_found
        ),
        trd_motif_count = trd_motif_count + as.integer(has_motif)
      )
  }
  
  # Filter to hits only
  trd_motif_hits <- trd_motif_hits %>%
    dplyr::filter(trd_motif_count > 0)
  
  message(sprintf("Found %d proteins with TRD sequence motifs", nrow(trd_motif_hits)))
  
  return(trd_motif_hits)
}


#' Extract TRD Sequences for Analysis
#'
#' @description Extracts TRD sequences for downstream comparison with REBASE
#'
#' @param trd_data Data with TRD genes
#' @param seq_col Sequence column
#' @param extract_region Whether to try extracting just TRD region (requires domain boundaries)
#'
#' @return Tibble with extracted TRD sequences ready for REBASE comparison
#' @export
extract_trd_sequences <- function(trd_data, 
                                   seq_col = "translation",
                                   extract_region = FALSE) {
  
  id_col <- intersect(c("locus_tag", "protein_id"), names(trd_data))[1]
  
  if (!seq_col %in% names(trd_data)) {
    stop("Sequence column not found: ", seq_col)
  }
  
  trd_seqs <- trd_data %>%
    dplyr::filter(!is.na(.data[[seq_col]])) %>%
    dplyr::select(
      dplyr::all_of(id_col),
      sequence = dplyr::all_of(seq_col),
      dplyr::any_of(c("product", "trd_type", "trd_pfam_str", "trd_type_keyword"))
    ) %>%
    dplyr::mutate(
      seq_length = nchar(sequence),
      # Create FASTA header
      fasta_header = paste0(">", .data[[id_col]], " | ",
                            dplyr::coalesce(trd_type, "putative_TRD"))
    )
  
  message(sprintf("Extracted %d TRD sequences for analysis", nrow(trd_seqs)))
  
  return(trd_seqs)
}


#' Combine TRD Extraction Results
#'
#' @description Integrates all TRD extraction methods
#'
#' @param dnmb_data DNMB data frame
#' @param pfam_col PFAM column
#' @param product_col Product column
#' @param seq_col Sequence column
#' @param search_motif_sequences Whether to search sequence motifs
#'
#' @return Combined tibble with all TRD candidates
#' @export
extract_trd_regions <- function(dnmb_data,
                                 pfam_col = NULL,
                                 product_col = "product",
                                 seq_col = "translation",
                                 search_motif_sequences = FALSE) {
  
  # PFAM-based extraction
  pfam_results <- tryCatch(
    extract_trd_by_pfam(dnmb_data, pfam_col),
    error = function(e) {
      message("PFAM TRD extraction skipped: ", e$message)
      dplyr::tibble()
    }
  )
  
  # Keyword-based extraction
  keyword_results <- tryCatch(
    extract_trd_by_keyword(dnmb_data, product_col),
    error = function(e) {
      message("Keyword TRD extraction skipped: ", e$message)
      dplyr::tibble()
    }
  )
  
  # Motif-based search
  motif_results <- dplyr::tibble()
  if (search_motif_sequences && seq_col %in% names(dnmb_data)) {
    motif_results <- tryCatch(
      search_trd_motifs(dnmb_data, seq_col),
      error = function(e) {
        message("Motif TRD search skipped: ", e$message)
        dplyr::tibble()
      }
    )
  }
  
  # Get ID column
  id_col <- intersect(c("locus_tag", "protein_id"), names(dnmb_data))[1]
  
  # Combine all IDs
  all_ids <- unique(c(
    if (nrow(pfam_results) > 0) pfam_results[[id_col]] else character(),
    if (nrow(keyword_results) > 0) keyword_results[[id_col]] else character(),
    if (nrow(motif_results) > 0) motif_results[[id_col]] else character()
  ))
  
  if (length(all_ids) == 0) {
    message("No TRD regions identified")
    return(dplyr::tibble())
  }
  
  # Build combined table
  combined <- dnmb_data %>%
    dplyr::filter(.data[[id_col]] %in% all_ids)
  
  # Add PFAM TRD info
  if (nrow(pfam_results) > 0) {
    combined <- combined %>%
      dplyr::left_join(
        pfam_results %>% dplyr::select(dplyr::all_of(id_col), 
                                        trd_pfam_str, trd_type),
        by = id_col
      )
  }
  
  # Add keyword TRD info
  if (nrow(keyword_results) > 0) {
    combined <- combined %>%
      dplyr::left_join(
        keyword_results %>% dplyr::select(dplyr::all_of(id_col), trd_type_keyword),
        by = id_col
      )
  }
  
  # Add motif info
  if (nrow(motif_results) > 0) {
    combined <- combined %>%
      dplyr::left_join(
        motif_results %>% dplyr::select(dplyr::all_of(id_col), 
                                         trd_motifs_found, trd_motif_count),
        by = id_col
      )
  }
  
  # Ensure columns exist even if joins didn't happen
  if (!"trd_type" %in% names(combined)) combined$trd_type <- NA_character_
  if (!"trd_type_keyword" %in% names(combined)) combined$trd_type_keyword <- NA_character_
  if (!"trd_motifs_found" %in% names(combined)) combined$trd_motifs_found <- NA_character_
  
  # Consensus TRD classification
  combined <- combined %>%
    dplyr::mutate(
      trd_classification = dplyr::coalesce(trd_type, trd_type_keyword, "putative_TRD"),
      trd_detection_method = paste(
        dplyr::if_else(!is.na(trd_type), "pfam", ""),
        dplyr::if_else(!is.na(trd_type_keyword), "keyword", ""),
        dplyr::if_else(!is.na(trd_motifs_found), "motif", ""),
        sep = "+"
      ),
      trd_confidence = dplyr::case_when(
        stringr::str_count(trd_detection_method, "\\+") >= 2 ~ "high",
        stringr::str_detect(trd_detection_method, "pfam") ~ "medium",
        TRUE ~ "low"
      )
    ) %>%
    dplyr::mutate(
      trd_detection_method = stringr::str_replace_all(trd_detection_method, "NA\\+?|\\+NA", ""),
      trd_detection_method = stringr::str_replace_all(trd_detection_method, "^\\+|\\+$", "")
    )
  
  message(sprintf("Total: %d TRD regions identified", nrow(combined)))
  
  return(combined)
}


#' Write TRD Sequences to FASTA
#'
#' @description Exports TRD sequences to FASTA format for BLAST analysis
#'
#' @param trd_data TRD data with sequences
#' @param output_file Output FASTA file path
#' @param seq_col Sequence column name
#'
#' @return Invisibly returns the output file path
#' @export
write_trd_fasta <- function(trd_data, output_file, seq_col = "translation") {
  
  id_col <- intersect(c("locus_tag", "protein_id"), names(trd_data))[1]
  
  if (!seq_col %in% names(trd_data)) {
    stop("Sequence column not found: ", seq_col)
  }
  
  # Build FASTA content
  fasta_lines <- character()
  
  for (i in seq_len(nrow(trd_data))) {
    header <- paste0(">", trd_data[[id_col]][i])
    
    # Add classification if available
    if ("trd_classification" %in% names(trd_data)) {
      header <- paste(header, trd_data$trd_classification[i], sep = " | ")
    }
    if ("product" %in% names(trd_data)) {
      header <- paste(header, trd_data$product[i], sep = " | ")
    }
    
    seq <- trd_data[[seq_col]][i]
    
    if (!is.na(seq) && nchar(seq) > 0) {
      # Wrap sequence at 60 characters
      seq_wrapped <- paste(
        stringr::str_sub(seq, seq(1, nchar(seq), 60), seq(60, nchar(seq) + 59, 60)),
        collapse = "\n"
      )
      
      fasta_lines <- c(fasta_lines, header, seq_wrapped)
    }
  }
  
  # Write to file
  writeLines(fasta_lines, output_file)
  
  message(sprintf("Wrote %d TRD sequences to %s", nrow(trd_data), output_file))
  
  invisible(output_file)
}

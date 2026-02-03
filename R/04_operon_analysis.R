#' @title Operon Analysis Module
#' @description Functions for analyzing R-M system operons with type-specific completeness
#' @name operon_analysis
NULL

#' Identify R-M System Operons
#'
#' @description Identifies R-M operons by pairing MTase and REase genes.
#' Allows maximum 1 intervening gene between R-M components.
#' Requires same system type for pairing.
#'
#' @param mtase_data MTase candidates from detect_methyltransferases()
#' @param rease_data REase candidates from detect_restriction_enzymes()
#' @param all_data Full DNMB data for context
#' @param max_intervening Maximum intervening genes allowed (default: 1)
#' @param max_distance Maximum bp distance between components (default: 5000)
#'
#' @return Tibble with R-M operon assignments and completeness
#' @export
identify_rm_operons <- function(mtase_data,
                                 rease_data,
                                 all_data,
                                 max_intervening = 1,
                                 max_distance = 5000) {

  id_col <- intersect(c("locus_tag", "protein_id"), names(all_data))[1]

  # Ensure numeric positions
  all_data <- all_data %>%
    dplyr::mutate(
      start = as.numeric(start),
      end = as.numeric(end)
    )

  # Sort all genes by position
  all_sorted <- all_data %>%
    dplyr::arrange(contig, start) %>%
    dplyr::mutate(.gene_idx = dplyr::row_number())

  # Get MTase and REase info with positions

  mtase_ids <- if (nrow(mtase_data) > 0) mtase_data[[id_col]] else character()
  rease_ids <- if (nrow(rease_data) > 0) rease_data[[id_col]] else character()

  if (length(mtase_ids) == 0 && length(rease_ids) == 0) {
    message("No R-M genes to analyze")
    return(dplyr::tibble())
  }

  # Mark R-M genes in sorted data
  rm_genes <- all_sorted %>%
    dplyr::filter(.data[[id_col]] %in% c(mtase_ids, rease_ids)) %>%
    dplyr::mutate(
      is_mtase = .data[[id_col]] %in% mtase_ids,
      is_rease = .data[[id_col]] %in% rease_ids
    )

  # Add type information from detection results (including REBASE type if available)
  if (nrow(mtase_data) > 0) {
    mtase_cols <- intersect(c("mtase_type", "rebase_type"), names(mtase_data))
    if (length(mtase_cols) > 0) {
      rm_genes <- rm_genes %>%
        dplyr::left_join(
          mtase_data %>% dplyr::select(dplyr::all_of(c(id_col, mtase_cols))),
          by = id_col
        )
    }
  }
  if (!"mtase_type" %in% names(rm_genes)) rm_genes$mtase_type <- NA_character_

  if (nrow(rease_data) > 0) {
    rease_cols <- intersect(c("re_type", "rebase_type"), names(rease_data))
    if (length(rease_cols) > 0) {
      # Avoid duplicate rebase_type column
      if ("rebase_type" %in% names(rm_genes) && "rebase_type" %in% rease_cols) {
        rm_genes <- rm_genes %>%
          dplyr::left_join(
            rease_data %>%
              dplyr::select(dplyr::all_of(id_col), re_type,
                            rebase_type_r = dplyr::any_of("rebase_type")),
            by = id_col
          ) %>%
          dplyr::mutate(
            rebase_type = dplyr::coalesce(rebase_type, rebase_type_r)
          ) %>%
          dplyr::select(-dplyr::any_of("rebase_type_r"))
      } else {
        rm_genes <- rm_genes %>%
          dplyr::left_join(
            rease_data %>% dplyr::select(dplyr::all_of(c(id_col, rease_cols))),
            by = id_col
          )
      }
    }
  }
  if (!"re_type" %in% names(rm_genes)) rm_genes$re_type <- NA_character_
  if (!"rebase_type" %in% names(rm_genes)) rm_genes$rebase_type <- NA_character_

  # Infer system type for each gene
  # PRIORITY: REBASE type > Product annotation > Detection type

  rm_genes <- rm_genes %>%
    dplyr::mutate(
      inferred_type = dplyr::case_when(
        # HIGHEST PRIORITY: REBASE-derived type (most reliable)
        !is.na(rebase_type) & rebase_type != "unknown" ~ rebase_type,

        # Type I: MUST have Hsd keyword in product (HsdM, HsdR, HsdS)
        stringr::str_detect(product, "(?i)hsd[MRS]|type.?I.*(restriction|modification|subunit)") ~ "Type_I",
        stringr::str_detect(re_type, "(?i)^Type_I$") ~ "Type_I",
        stringr::str_detect(mtase_type, "(?i)^Type_I$") ~ "Type_I",

        # Type III: Mod/Res system
        stringr::str_detect(product, "(?i)type.?III|\\bmod\\b|\\bres\\b") ~ "Type_III",
        stringr::str_detect(re_type, "(?i)Type_III") ~ "Type_III",

        # Type IV: McrBC, Mrr (modified base restriction)
        stringr::str_detect(product, "(?i)type.?IV|mcr[ABC]|\\bmrr\\b") ~ "Type_IV",
        stringr::str_detect(re_type, "(?i)Type_IV") ~ "Type_IV",

        # Type II: explicit or default for M+R pairs
        stringr::str_detect(product, "(?i)type.?II") ~ "Type_II",
        stringr::str_detect(re_type, "(?i)Type_II") ~ "Type_II",

        # Orphan Dam/Dcm
        stringr::str_detect(mtase_type, "(?i)Dam|Dcm") ~ "orphan",

        # Default: unknown (will be resolved at operon level)
        TRUE ~ "unknown"
      )
    )

  # Pair R-M genes into operons
  operons <- pair_rm_genes(rm_genes, all_sorted, id_col, max_intervening, max_distance)

  if (nrow(operons) == 0) {
    message("No R-M operons identified")
    return(dplyr::tibble())
  }

  # Calculate completeness for each operon
  operons <- calculate_operon_completeness(operons, rm_genes, id_col)

  # Clean output - keep essential columns
  operons <- operons %>%
    dplyr::relocate(
      dplyr::any_of(c(id_col, "product", "operon_id", "rm_role", "system_type",
                       "completeness", "completeness_detail", "n_components"))
    )

  # Summary
  n_operons <- dplyr::n_distinct(operons$operon_id)
  complete_count <- operons %>%
    dplyr::filter(completeness == "complete") %>%
    dplyr::pull(operon_id) %>%
    dplyr::n_distinct()

  message(sprintf("Identified %d R-M operons (%d complete)", n_operons, complete_count))

  return(operons)
}


#' Pair R-M Genes into Operons
#'
#' @description Internal function to pair MTase and REase genes
#'
#' @keywords internal
pair_rm_genes <- function(rm_genes, all_sorted, id_col, max_intervening, max_distance) {

  if (nrow(rm_genes) == 0) return(dplyr::tibble())

  # Group by contig
  contigs <- unique(rm_genes$contig)
  operon_list <- list()
  operon_counter <- 0


  for (ctg in contigs) {
    ctg_genes <- rm_genes %>%
      dplyr::filter(contig == ctg) %>%
      dplyr::arrange(.gene_idx)

    if (nrow(ctg_genes) == 0) next

    # Track which genes are already assigned
    assigned <- rep(FALSE, nrow(ctg_genes))

    for (i in seq_len(nrow(ctg_genes))) {
      if (assigned[i]) next

      current_gene <- ctg_genes[i, ]
      operon_members <- list(current_gene)
      assigned[i] <- TRUE

      # Look for nearby R-M genes to form operon
      if (i < nrow(ctg_genes)) {
        for (j in (i + 1):nrow(ctg_genes)) {
          if (assigned[j]) next

          candidate <- ctg_genes[j, ]

          # Check distance (gene index difference)
          gene_gap <- candidate$.gene_idx - current_gene$.gene_idx - 1

          if (gene_gap > max_intervening) break  # Too far, stop looking

          # Check bp distance
          bp_dist <- abs(candidate$start - current_gene$end)
          if (bp_dist > max_distance) next

          # Check type compatibility
          type1 <- current_gene$inferred_type
          type2 <- candidate$inferred_type
          if (type1 != "unknown" && type2 != "unknown" && type1 != type2) next

          # Valid pairing - add to operon
          operon_members[[length(operon_members) + 1]] <- candidate
          assigned[j] <- TRUE

          # Update current_gene for next iteration
          current_gene <- candidate
        }
      }

      # Create operon if we have members
      if (length(operon_members) > 0) {
        operon_counter <- operon_counter + 1
        operon_df <- dplyr::bind_rows(operon_members)
        operon_df$operon_id <- operon_counter
        operon_list[[operon_counter]] <- operon_df
      }
    }
  }

  if (length(operon_list) == 0) return(dplyr::tibble())

  result <- dplyr::bind_rows(operon_list)

  # Assign rm_role
  result <- result %>%
    dplyr::mutate(
      rm_role = dplyr::case_when(
        is_mtase & is_rease ~ "M+R",
        is_mtase ~ "M",
        is_rease ~ "R",
        TRUE ~ "other"
      )
    )

  return(result)
}


#' Calculate Operon Completeness
#'
#' @description Calculates completeness based on system type
#' - Type I: M + R + S (3 components)
#' - Type II: M + R (2 components)
#' - Type III: Mod + Res (2 components)
#' - Type IV: typically single R component
#'
#' @keywords internal
calculate_operon_completeness <- function(operons, rm_genes, id_col) {

  # Check for Type I specific subunits (HsdM, HsdR, HsdS) in each operon
  operons <- operons %>%
    dplyr::mutate(
      # Strict S subunit detection - only HsdS or explicit "specificity subunit"
      is_s_subunit = stringr::str_detect(product, "(?i)hsdS|specificity.?subunit|type.?I.*S.?subunit"),
      # Check for explicit Hsd subunits
      has_hsd_keyword = stringr::str_detect(product, "(?i)hsd[MRS]")
    )

  # Determine system type per operon
  operon_types <- operons %>%
    dplyr::group_by(operon_id) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_mtase = sum(is_mtase, na.rm = TRUE),
      n_rease = sum(is_rease, na.rm = TRUE),
      n_s_subunit = sum(is_s_subunit, na.rm = TRUE),
      has_any_hsd = any(has_hsd_keyword, na.rm = TRUE),
      types_found = paste(unique(inferred_type[inferred_type != "unknown"]), collapse = ";"),
      products_concat = paste(product, collapse = " | "),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Type I requires HsdM/R/S keywords
      is_type_I = has_any_hsd | stringr::str_detect(types_found, "Type_I"),
      system_type = dplyr::case_when(
        # Type I: must have Hsd keywords
        is_type_I ~ "Type_I",
        # Type III
        stringr::str_detect(types_found, "Type_III") ~ "Type_III",
        # Type IV
        stringr::str_detect(types_found, "Type_IV") ~ "Type_IV",
        # Type II (explicit or default for M+R)
        stringr::str_detect(types_found, "Type_II") ~ "Type_II",
        n_mtase > 0 & n_rease > 0 ~ "Type_II",  # Default for M+R pairs
        # Orphans
        n_mtase > 0 & n_rease == 0 ~ "orphan_M",
        n_rease > 0 & n_mtase == 0 ~ "orphan_R",
        TRUE ~ "unknown"
      )
    )

  # Calculate completeness per operon (use operon_types which already has subunit info)
  completeness_info <- operon_types %>%
    dplyr::mutate(
      has_M = n_mtase > 0,
      has_R = n_rease > 0,
      has_S = n_s_subunit > 0,
      n_components = n_mtase + n_rease
    ) %>%
    dplyr::mutate(
      completeness = dplyr::case_when(
        # Type I: needs M + R + S
        system_type == "Type_I" & has_M & has_R & has_S ~ "complete",
        system_type == "Type_I" & has_M & has_R ~ "partial_no_S",
        system_type == "Type_I" & has_M & has_S ~ "partial_no_R",
        system_type == "Type_I" & has_R & has_S ~ "partial_no_M",
        system_type == "Type_I" ~ "incomplete",

        # Type II: needs M + R
        system_type == "Type_II" & has_M & has_R ~ "complete",
        system_type == "Type_II" & has_M ~ "M_only",
        system_type == "Type_II" & has_R ~ "R_only",

        # Type III: needs Mod + Res
        system_type == "Type_III" & has_M & has_R ~ "complete",
        system_type == "Type_III" & has_M ~ "Mod_only",
        system_type == "Type_III" & has_R ~ "Res_only",

        # Type IV: typically R only
        system_type == "Type_IV" & has_R ~ "complete",

        # Orphans
        system_type == "orphan_M" ~ "orphan_MTase",
        system_type == "orphan_R" ~ "orphan_REase",

        TRUE ~ "unknown"
      ),
      completeness_detail = dplyr::case_when(
        system_type == "Type_I" ~ paste0("M:", has_M, " R:", has_R, " S:", has_S),
        TRUE ~ paste0("M:", has_M, " R:", has_R)
      )
    )

  # Join back to operons
  operons <- operons %>%
    dplyr::left_join(
      completeness_info %>% dplyr::select(operon_id, system_type, completeness,
                                           completeness_detail, n_components),
      by = "operon_id"
    )

  return(operons)
}


#' Summarize R-M Operons
#'
#' @description Creates a summary table of R-M operons
#'
#' @param operon_data Data from identify_rm_operons()
#'
#' @return Summary tibble with one row per operon
#' @export
summarize_rm_operons <- function(operon_data) {

  id_col <- intersect(c("locus_tag", "protein_id"), names(operon_data))[1]

  operon_data %>%
    dplyr::group_by(operon_id) %>%
    dplyr::summarise(
      system_type = dplyr::first(system_type),
      completeness = dplyr::first(completeness),
      n_genes = dplyr::n(),
      genes = paste(.data[[id_col]], collapse = ", "),
      products = paste(product, collapse = " | "),
      contig = dplyr::first(contig),
      start_pos = min(start, na.rm = TRUE),
      end_pos = max(end, na.rm = TRUE),
      span_bp = end_pos - start_pos,
      .groups = "drop"
    ) %>%
    dplyr::arrange(system_type, dplyr::desc(completeness == "complete"), operon_id)
}


#' Predict R-M Type from Operon (Legacy wrapper)
#'
#' @description Legacy function for compatibility - now integrated into identify_rm_operons
#'
#' @param operon_data Operon data
#'
#' @return Data with predicted type (passthrough if already present)
#' @export
predict_rm_type_from_operon <- function(operon_data) {

  if ("system_type" %in% names(operon_data)) {
    # Already has type info, just rename for compatibility
    if (!"predicted_type" %in% names(operon_data)) {
      operon_data$predicted_type <- operon_data$system_type
    }
    if (!"operon_type" %in% names(operon_data)) {
      operon_data$operon_type <- operon_data$completeness
    }
    return(operon_data)
  }

  # Fallback for old data format
  operon_data %>%
    dplyr::mutate(
      predicted_type = "unknown",
      operon_type = "unknown"
    )
}

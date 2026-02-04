#' @title RMscan Main Pipeline
#' @description Main wrapper function for complete R-M system analysis
#' @name rmscan_pipeline
NULL

#' Run Complete R-M System Analysis Pipeline
#'
#' @description Executes the full RMscan analysis pipeline on DNMB data
#'
#' Pipeline order:
#' 1. Load data & annotation analysis
#' 2. Detect MTase (PFAM + keyword) - initial candidates
#' 3. Detect REase (PFAM + keyword) - initial candidates
#' 4. Operon (genome context) analysis
#' 5. REBASE BLAST comparison
#' 6. **FILTER** - Remove low-quality BLAST matches (identity < 10%, length < 50)
#' 7. Create annotated tables (BLAST results + operon info, important cols first)
#' 7b. Detect catalytic motifs (PDxDxK, DPPY, PxGxG, Walker, DEAD/DEAH, etc.)
#' 8. TRD extraction (for filtered candidates)
#' 9. Extract recognition sequences from REBASE
#' 10. Create comprehensive summary table
#' 11. Generate R-M type heatmap (blast_identity >= 50%)
#' 12. Export results (R-M_REBASE_analysis.xlsx)
#'
#' @param input_file Path to DNMB annotation file (Excel or CSV)
#' @param output_dir Output directory for results
#' @param compare_rebase Whether to run REBASE comparison
#' @param search_motifs Whether to search catalytic/TRD motifs
#' @param blast_min_identity Minimum BLAST identity for filtering (default 0.10 = 10%)
#' @param blast_min_length Minimum alignment length for filtering (default 50)
#' @param max_operon_gap Maximum gap for operon detection
#' @param max_intervening Maximum intervening genes for operon
#' @param verbose Print progress messages
#' @param save_intermediates Save intermediate data to global environment for debugging
#'
#' @return List containing all analysis results
#' @export
rmscan_pipeline <- function(input_file,
                             output_dir = ".",
                             compare_rebase = TRUE,
                             search_motifs = TRUE,
                             blast_min_identity = 0.10,
                             blast_min_length = 50,
                             max_operon_gap = 5000,
                             max_intervening = 1,
                             verbose = TRUE,
                             save_intermediates = TRUE) {

  # Start timing
  start_time <- Sys.time()

  if (verbose) message("\n", paste(rep("=", 60), collapse = ""))
  if (verbose) message("RMscan - Restriction-Modification System Analysis Pipeline")
  if (verbose) message(paste(rep("=", 60), collapse = ""))
  if (verbose) message("Started: ", format(start_time, "%Y-%m-%d %H:%M:%S"))

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }

  # Helper function to save to global environment
  save_to_env <- function(name, value) {
    if (save_intermediates) {
      assign(paste0("rmscan_", name), value, envir = .GlobalEnv)
      if (verbose) message("  -> Saved to global env: rmscan_", name)
    }
  }

  # ============================================
  # Step 1: Load Data & Annotation Analysis
  # ============================================
  if (verbose) message("\n[Step 1/12] Loading DNMB data & analyzing annotations...")

  dnmb_data <- load_dnmb(input_file)
  annotation_sources <- detect_annotation_sources(dnmb_data)
  id_col <- intersect(c("locus_tag", "protein_id"), names(dnmb_data))[1]

  save_to_env("dnmb", dnmb_data)
  save_to_env("annotation_sources", annotation_sources)

  # ============================================
  # Step 2: Detect Methyltransferases (Initial Candidates)
  # ============================================
  if (verbose) message("\n[Step 2/12] Detecting methyltransferases (PFAM + keyword)...")

  mtase_data <- tryCatch({
    detect_methyltransferases(
      dnmb_data,
      pfam_col = if (annotation_sources$pfam) NULL else NA,
      product_col = "product"
    )
  }, error = function(e) {
    warning("MTase detection failed: ", e$message)
    dplyr::tibble()
  })

  save_to_env("mtase", mtase_data)
  if (verbose) message("  Found ", nrow(mtase_data), " initial MTase candidates")

  # ============================================
  # Step 3: Detect Restriction Enzymes (Initial Candidates)
  # ============================================
  if (verbose) message("\n[Step 3/12] Detecting restriction enzymes (PFAM + keyword)...")

  rease_data <- tryCatch({
    detect_restriction_enzymes(
      dnmb_data,
      pfam_col = if (annotation_sources$pfam) NULL else NA,
      product_col = "product",
      seq_col = "translation",
      search_motifs = FALSE
    )
  }, error = function(e) {
    warning("REase detection failed: ", e$message)
    dplyr::tibble()
  })

  save_to_env("rease", rease_data)
  if (verbose) message("  Found ", nrow(rease_data), " initial REase candidates")

  # ============================================
  # Step 4: Operon (Genome Context) Analysis
  # ============================================
  if (verbose) message("\n[Step 4/12] Analyzing operon/genome context...")

  operon_data <- tryCatch({
    if (nrow(mtase_data) > 0 || nrow(rease_data) > 0) {
      rm_operons <- identify_rm_operons(
        mtase_data,
        rease_data,
        dnmb_data,
        max_intervening = max_intervening,
        max_distance = max_operon_gap
      )
      if (nrow(rm_operons) > 0) {
        rm_operons <- predict_rm_type_from_operon(rm_operons)
      }
      rm_operons
    } else {
      dplyr::tibble()
    }
  }, error = function(e) {
    warning("Operon analysis failed: ", e$message)
    dplyr::tibble()
  })

  save_to_env("operon", operon_data)
  if (verbose && nrow(operon_data) > 0) {
    message("  Found ", dplyr::n_distinct(operon_data$operon_id), " operons")
  }

  # ============================================
  # Step 5: REBASE BLAST Comparison
  # ============================================
  rebase_results <- dplyr::tibble()
  blast_raw <- dplyr::tibble()
  blast_filtered <- dplyr::tibble()
  rebase_db <- dplyr::tibble()

  if (compare_rebase) {
    if (verbose) message("\n[Step 5/12] REBASE comparison (download, BLAST)...")

    # Load REBASE database
    tryCatch({
      rebase_db <- get_rebase_data(verbose = FALSE)
      save_to_env("rebase_db", rebase_db)
    }, error = function(e) {
      message("  Note: Could not load REBASE database: ", e$message)
    })

    # Run BLAST
    rebase_results <- tryCatch({
      query_ids <- unique(c(
        if (nrow(mtase_data) > 0) mtase_data[[id_col]] else character(),
        if (nrow(rease_data) > 0) rease_data[[id_col]] else character()
      ))

      if (length(query_ids) > 0) {
        query_seqs <- dnmb_data %>%
          dplyr::filter(.data[[id_col]] %in% query_ids) %>%
          dplyr::select(dplyr::all_of(c(id_col, "translation")))
        save_to_env("query_seqs", query_seqs)

        run_rebase_comparison(
          query_data = query_seqs,
          output_dir = output_dir,
          use_local_blast = TRUE
        )
      } else {
        dplyr::tibble()
      }
    }, error = function(e) {
      warning("REBASE comparison failed: ", e$message)
      dplyr::tibble()
    })

    # Parse raw BLAST results
    blast_file <- file.path(output_dir, "blast_results.txt")
    if (file.exists(blast_file) && file.info(blast_file)$size > 0) {
      blast_raw <- parse_blast_results(blast_file)
      save_to_env("blast_raw", blast_raw)
      if (verbose) message("  Parsed ", nrow(blast_raw), " raw BLAST hits")
    }
  }

  # ============================================
  # Step 6: FILTER BLAST Results
  # ============================================
  if (verbose) message("\n[Step 6/12] Filtering BLAST results (identity >= ",
                       blast_min_identity * 100, "%, length >= ", blast_min_length, ")...")

  if (nrow(blast_raw) > 0) {
    blast_filtered <- filter_blast_results(
      blast_raw,
      min_identity = blast_min_identity,
      min_length = blast_min_length,
      max_evalue = 1e-3,
      verbose = verbose
    )
    save_to_env("blast_filtered", blast_filtered)

    # Update rebase_results with filtered data
    if (nrow(blast_filtered) > 0 && nrow(rebase_db) > 0) {
      rebase_results <- get_best_rebase_match(blast_filtered, rebase_db, min_identity = blast_min_identity)
      save_to_env("rebase", rebase_results)
    }
  }

  # Get filtered candidate IDs (only candidates with quality BLAST hits proceed)
  filtered_query_ids <- if (nrow(blast_filtered) > 0) {
    unique(blast_filtered$query_id)
  } else {
    character()
  }

  if (verbose) message("  Candidates passing filter: ", length(filtered_query_ids))

  # ============================================
  # Step 7: Create Annotated Tables (BLAST details + operon info)
  # ============================================
  if (verbose) message("\n[Step 7/12] Creating annotated tables with BLAST details...")

  # Prepare operon info for joining (operon_id만)
  operon_info <- if (nrow(operon_data) > 0) {
    operon_data %>%
      dplyr::select(
        dplyr::all_of(id_col),
        dplyr::any_of(c("operon_id", "partner_locus_tag", "partner_product"))
      ) %>%
      dplyr::distinct()
  } else {
    NULL
  }

  # Prepare BLAST details for joining (include ALL BLAST columns)
  blast_details <- if (nrow(blast_filtered) > 0 && nrow(rebase_db) > 0) {
    # Get best hit per query with full details
    best_hits <- blast_filtered %>%
      dplyr::group_by(query_id) %>%
      dplyr::slice_max(order_by = pct_identity, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      # ★ Strip _\d+ suffix from rebase_enzyme for matching
      dplyr::mutate(.rebase_enzyme_clean = sub("_\\d+$", "", rebase_enzyme))

    # REBASE lookup - enz_type에서 rm_type 추출, enzyme_name에서 subunit 추출
    rebase_lookup <- rebase_db %>%
      dplyr::select(enzyme_name, enz_type, rec_seq) %>%
      dplyr::group_by(enzyme_name) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        # ★ enz_type에서 Type I~IV만 추출 (Type IIG 등 서브타입 제외)
        rm_type = dplyr::case_when(
          grepl("Type I[^IV]|Type I$", enz_type) ~ "Type I",
          grepl("Type II[^I]|Type II$", enz_type) ~ "Type II",
          grepl("Type III", enz_type) ~ "Type III",
          grepl("Type IV", enz_type) ~ "Type IV",
          TRUE ~ NA_character_
        ),
        # ★ enzyme_name 앞의 접두사에서 subunit 추출 (M., R., S. -> M, R, S, 없으면 R)
        subunit = dplyr::case_when(
          grepl("^M\\.", enzyme_name) ~ "M",
          grepl("^R\\.", enzyme_name) ~ "R",
          grepl("^S\\.", enzyme_name) ~ "S",
          TRUE ~ "R"  # 접두사 없으면 R
        )
      )

    # Join using cleaned enzyme name
    best_hits %>%
      dplyr::left_join(
        rebase_lookup,
        by = c(".rebase_enzyme_clean" = "enzyme_name")
      ) %>%
      dplyr::select(
        query_id,
        # BLAST results
        blast_match = rebase_enzyme,
        blast_identity = pct_identity,
        blast_length = length,
        blast_evalue = evalue,
        blast_bitscore = bitscore,
        # REBASE info - enz_type에서 추출한 rm_type, enzyme_name에서 추출한 subunit
        rm_type,
        subunit,
        rec_seq,
        # BLAST alignment details
        blast_qstart = qstart,
        blast_qend = qend,
        blast_sstart = sstart,
        blast_send = send
      )
  } else {
    NULL
  }

  save_to_env("blast_details", blast_details)

  # ★ DEBUG: Check blast_details content
  if (!is.null(blast_details) && nrow(blast_details) > 0) {
    if (verbose) {
      message("  BLAST details ready: ", nrow(blast_details), " hits")
      message("    Columns: ", paste(names(blast_details)[1:min(8, ncol(blast_details))], collapse = ", "), "...")
      message("    Sample query_ids: ", paste(head(blast_details$query_id, 3), collapse = ", "))
      n_with_rec_seq <- sum(!is.na(blast_details$rec_seq) & blast_details$rec_seq != "" & blast_details$rec_seq != "?")
      message("    Hits with rec_seq: ", n_with_rec_seq)
    }
  } else {
    if (verbose) message("  WARNING: No BLAST details available for annotation!")
  }

  # Create annotated MTase table
  mtase_annotated <- create_annotated_table(
    mtase_data,
    blast_details = blast_details,
    operon_info = operon_info,
    dnmb_data = dnmb_data,
    id_col = id_col,
    component_type = "MTase"
  )
  save_to_env("mtase_annotated", mtase_annotated)
  if (verbose) message("  Created mtase_annotated: ", nrow(mtase_annotated), " rows, ", ncol(mtase_annotated), " cols")

  # Create annotated REase table
  rease_annotated <- create_annotated_table(
    rease_data,
    blast_details = blast_details,
    operon_info = operon_info,
    dnmb_data = dnmb_data,
    id_col = id_col,
    component_type = "REase"
  )
  save_to_env("rease_annotated", rease_annotated)
  if (verbose) message("  Created rease_annotated: ", nrow(rease_annotated), " rows, ", ncol(rease_annotated), " cols")

  # ============================================
  # Step 7b: Detect Catalytic Motifs (for all candidates)
  # ============================================
  if (verbose) message("\n[Step 7b] Detecting catalytic motifs (PDxDxK, DPPY, PxGxG)...")

  catalytic_motifs <- dplyr::tibble()
  catalytic_summary <- dplyr::tibble()

  tryCatch({
    # Get sequences for all candidates
    all_candidate_ids <- unique(c(
      if (nrow(mtase_data) > 0) mtase_data[[id_col]] else character(),
      if (nrow(rease_data) > 0) rease_data[[id_col]] else character()
    ))

    if (length(all_candidate_ids) > 0) {
      candidate_seqs <- dnmb_data %>%
        dplyr::filter(.data[[id_col]] %in% all_candidate_ids) %>%
        dplyr::select(dplyr::all_of(c(id_col, "translation")))

      # Detect catalytic motifs
      catalytic_motifs <- detect_catalytic_motifs(
        candidate_seqs,
        seq_col = "translation",
        id_col = id_col,
        component_type = "all"
      )
      save_to_env("catalytic_motifs", catalytic_motifs)

      if (nrow(catalytic_motifs) > 0) {
        # Summarize motifs per gene
        catalytic_summary <- summarize_catalytic_motifs(catalytic_motifs, id_col = id_col)
        save_to_env("catalytic_summary", catalytic_summary)

        if (verbose) {
          n_mtase_complete <- sum(catalytic_summary$mtase_complete, na.rm = TRUE)
          n_rease_motif <- sum(catalytic_summary$rease_nuclease_present | catalytic_summary$rease_helicase_present, na.rm = TRUE)
          message(sprintf("  MTase with DPPY+PxGxG (complete): %d", n_mtase_complete))
          message(sprintf("  REase with PDxDxK motif: %d", n_rease_motif))
        }
      }
    }
  }, error = function(e) {
    warning("Catalytic motif detection failed: ", e$message)
  })

  # ============================================
  # Step 8: TRD Extraction (for filtered candidates)
  # ============================================
  if (verbose) message("\n[Step 8/12] TRD extraction...")

  trd_data <- dplyr::tibble()

  if (length(filtered_query_ids) > 0) {
    trd_data <- tryCatch({
      extract_trd_regions(
        dnmb_data %>% dplyr::filter(.data[[id_col]] %in% filtered_query_ids),
        pfam_col = if (annotation_sources$pfam) NULL else NA,
        product_col = "product",
        seq_col = "translation",
        search_motif_sequences = search_motifs
      )
    }, error = function(e) {
      warning("TRD extraction failed: ", e$message)
      dplyr::tibble()
    })
    save_to_env("trd", trd_data)
    if (verbose) message("  TRD data extracted for ", nrow(trd_data), " candidates")
  } else {
    if (verbose) message("  No filtered candidates for TRD extraction")
    save_to_env("trd", dplyr::tibble())
  }

  # NOTE: Scoring & Classification removed (not needed)

  # ============================================
  # Step 9: Recognition Sequence Extraction from REBASE
  # ============================================
  if (verbose) message("\n[Step 9/12] Extracting recognition sequences from REBASE BLAST results...")

  recognition_predictions <- dplyr::tibble()

  if (nrow(rebase_results) > 0) {
    recognition_predictions <- tryCatch({
      # Get recognition sequences from REBASE BLAST results
      rebase_results %>%
        dplyr::filter(!is.na(rec_seq) & rec_seq != "" & rec_seq != "?") %>%
        dplyr::select(
          gene_id = query_id,
          rm_type,
          predicted_recognition_seq = rec_seq
        ) %>%
        dplyr::mutate(
          prediction_source = "REBASE_BLAST",
          trd_applicable = dplyr::case_when(
            # Type I: S subunit is key
            grepl("Type.?I[^I]|Type_I$", rm_type) ~ "Type_I_check_S_subunit",
            # Type II: M subunit
            grepl("Type.?II", rm_type) ~ "Type_II_M_subunit",
            # Type III: Mod subunit
            grepl("Type.?III", rm_type) ~ "Type_III_Mod_subunit",
            TRUE ~ "check_manually"
          )
        )
    }, error = function(e) {
      warning("Recognition prediction failed: ", e$message)
      dplyr::tibble()
    })

    save_to_env("recognition", recognition_predictions)
    if (verbose) message("  Predicted recognition sequences: ", nrow(recognition_predictions))
  }

  # ============================================
  # Step 10: Create Comprehensive Summary Table
  # ============================================
  if (verbose) message("\n[Step 10/12] Creating comprehensive R-M summary table...")

  # Get trd_data - use local variable if exists, otherwise try global env
  trd_data_for_comprehensive <- if (exists("trd_data") && !is.null(trd_data) && nrow(trd_data) > 0) {
    trd_data
  } else if (exists("rmscan_trd", envir = .GlobalEnv)) {
    get("rmscan_trd", envir = .GlobalEnv)
  } else {
    dplyr::tibble()
  }

  # catalytic_summary는 Step 7b에서 이미 생성됨 - 직접 사용
  # (catalytic_summary 변수가 이 스코프에 이미 존재함)

  if (verbose && nrow(catalytic_summary) > 0) {
    message("  Catalytic motif data: ", nrow(catalytic_summary), " genes with motifs")
  }

  rm_comprehensive <- create_comprehensive_rm_table_v2(
    mtase_annotated = mtase_annotated,
    rease_annotated = rease_annotated,
    trd_data = trd_data_for_comprehensive,
    catalytic_summary = catalytic_summary,
    filtered_ids = filtered_query_ids,
    id_col = id_col
  )

  # NOTE: R-M System Type Classification removed (rm_type info from REBASE is sufficient)

  save_to_env("rm_comprehensive", rm_comprehensive)
  if (verbose) message("  Created comprehensive table: ", nrow(rm_comprehensive), " candidates")

  # ============================================
  # Step 11: Generate Heatmap (blast_identity >= 0.5)
  # ============================================
  if (verbose) message("\n[Step 11/12] Generating R-M type heatmap...")

  heatmap_file <- NULL
  tryCatch({
    heatmap_file <- generate_rm_type_heatmap(
      rm_comprehensive,
      output_dir = output_dir,
      min_identity = 0.5,
      id_col = id_col,
      verbose = verbose
    )
    if (!is.null(heatmap_file)) {
      save_to_env("heatmap_file", heatmap_file)
    }
  }, error = function(e) {
    warning("Heatmap generation failed: ", e$message)
  })

  # ============================================
  # Step 12: Export Results (xlsx only)
  # ============================================
  if (verbose) message("\n[Step 12/12] Exporting results...")

  # Main xlsx output
  xlsx_file <- file.path(output_dir, "R-M_REBASE_analysis.xlsx")
  export_rm_xlsx(rm_comprehensive, xlsx_file, verbose = verbose)

  # ============================================
  # Finish
  # ============================================
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")

  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("Analysis Complete!")
    message(sprintf("Time elapsed: %.1f seconds", as.numeric(elapsed)))
    message("Results saved to: ", output_dir)
    message(paste(rep("=", 60), collapse = ""))
  }

  # Return all results
  return(list(
    # Main results
    rm_comprehensive = rm_comprehensive,
    mtase_annotated = mtase_annotated,
    rease_annotated = rease_annotated,

    # Intermediate data
    dnmb_data = dnmb_data,
    mtase_data = mtase_data,
    rease_data = rease_data,
    operon_data = operon_data,
    blast_raw = blast_raw,
    blast_filtered = blast_filtered,
    rebase_results = rebase_results,
    recognition_predictions = recognition_predictions,

    # Catalytic motif data
    catalytic_motifs = catalytic_motifs,
    catalytic_summary = catalytic_summary,

    # Output files
    xlsx_file = xlsx_file,
    heatmap_file = heatmap_file,

    # Metadata
    filtered_candidate_ids = filtered_query_ids,
    output_dir = output_dir,
    elapsed_time = elapsed
  ))
}


#' Create Annotated Table with BLAST Details
#'
#' @description Creates annotated table with BLAST results and operon info,
#' with important columns at the front
#'
#' @param rm_data MTase or REase data
#' @param blast_details BLAST details table
#' @param operon_info Operon information
#' @param dnmb_data Original DNMB data
#' @param id_col ID column name
#' @param component_type "MTase" or "REase"
#'
#' @return Annotated tibble with columns in priority order
#' @export
create_annotated_table <- function(rm_data,
                                    blast_details = NULL,
                                    operon_info = NULL,
                                    dnmb_data = NULL,
                                    id_col = "locus_tag",
                                    component_type = "MTase") {

  if (is.null(rm_data) || nrow(rm_data) == 0) {
    return(dplyr::tibble())
  }

  result <- rm_data
  # NOTE: rm_component 제거됨

  # Add BLAST details
  if (!is.null(blast_details) && nrow(blast_details) > 0) {
    result <- result %>%
      dplyr::left_join(
        blast_details %>% dplyr::rename(!!id_col := query_id),
        by = id_col
      )
  } else {
    # Add empty BLAST columns
    result$blast_match <- NA_character_
    result$blast_identity <- NA_real_
    result$blast_length <- NA_integer_
    result$blast_evalue <- NA_real_
    result$blast_bitscore <- NA_real_
    result$rm_type <- NA_character_
    result$subunit <- NA_character_
    result$rec_seq <- NA_character_
  }

  # Add operon info (operon_id만 가져옴)
  if (!is.null(operon_info) && nrow(operon_info) > 0) {
    # operon_id와 partner 정보만 선택
    operon_cols_to_join <- operon_info %>%
      dplyr::select(
        dplyr::all_of(id_col),
        dplyr::any_of(c("operon_id", "partner_locus_tag", "partner_product"))
      ) %>%
      dplyr::distinct()

    result <- result %>% dplyr::left_join(operon_cols_to_join, by = id_col)
  }

  # Add DNMB info (select key columns only)
  if (!is.null(dnmb_data) && nrow(dnmb_data) > 0) {
    dnmb_cols <- dnmb_data %>%
      dplyr::select(
        dplyr::all_of(id_col),
        dplyr::any_of(c("protein_id", "gene", "product", "old_locus_tag",
                        "seqid", "start", "end", "strand"))
      ) %>%
      dplyr::distinct()

    # Avoid duplicate columns
    existing_cols <- intersect(names(result), names(dnmb_cols))
    existing_cols <- setdiff(existing_cols, id_col)
    if (length(existing_cols) > 0) {
      dnmb_cols <- dnmb_cols %>% dplyr::select(-dplyr::any_of(existing_cols))
    }
    result <- result %>% dplyr::left_join(dnmb_cols, by = id_col)
  }

  # Reorder columns
  priority_cols <- c(
    # ===== 1. ID & BASIC INFO =====
    id_col,
    "product",
    "gene",

    # ===== 2. BLAST results =====
    "rm_type",        # ★ REBASE rm_type 그대로
    "subunit",        # ★ REBASE subunit 그대로
    "blast_match",
    "rec_seq",
    "blast_identity", "blast_length", "blast_evalue", "blast_bitscore",
    "blast_qstart", "blast_qend", "blast_sstart", "blast_send",

    # ===== 3. Operon info (operon_id만) =====
    "operon_id",
    "partner_locus_tag", "partner_product",

    # ===== 4. Location =====
    "seqid", "start", "end", "strand"
  )

  existing_priority <- intersect(priority_cols, names(result))

  result <- result %>%
    dplyr::select(dplyr::all_of(existing_priority), dplyr::everything())

  return(result)
}


#' Create Comprehensive R-M Table V2
#'
#' @description Creates comprehensive summary combining MTase + REase with all info.
#'
#' Sorting: passed_blast_filter (TRUE first) -> blast_identity -> blast_evalue -> operon_id -> locus_tag
#' Column order: ID, product, gene, passed_blast_filter, BLAST/REBASE info, catalytic motifs, operon, location
#'
#' @param mtase_annotated Annotated MTase data
#' @param rease_annotated Annotated REase data
#' @param trd_data TRD data (optional)
#' @param catalytic_summary Catalytic motif summary (optional)
#' @param filtered_ids IDs that passed BLAST filtering
#' @param id_col ID column name
#'
#' @return Comprehensive tibble sorted by passed_blast_filter, blast_identity
#' @export
create_comprehensive_rm_table_v2 <- function(mtase_annotated = NULL,
                                              rease_annotated = NULL,
                                              trd_data = NULL,
                                              catalytic_summary = NULL,
                                              filtered_ids = character(),
                                              id_col = "locus_tag") {

  # Combine MTase and REase
  rm_combined <- dplyr::bind_rows(mtase_annotated, rease_annotated)

  if (nrow(rm_combined) == 0) {
    message("No R-M candidates to combine")
    return(dplyr::tibble())
  }

  # Add "passed_filter" flag
  rm_combined <- rm_combined %>%
    dplyr::mutate(
      passed_blast_filter = .data[[id_col]] %in% filtered_ids
    )

  # NOTE: scored_data and classified_data joining removed (not needed)

  # Add TRD info if available
  if (!is.null(trd_data) && nrow(trd_data) > 0) {
    trd_cols <- trd_data %>%
      dplyr::select(
        dplyr::any_of(c(id_col, "locus_tag", "protein_id")),
        dplyr::any_of(c("trd_count", "trd_regions", "trd_motifs", "trd_pfam",
                        "has_trd", "trd_start", "trd_end", "trd_sequence"))
      ) %>%
      dplyr::distinct()

    # Get the ID column that exists in trd_cols
    trd_id_col <- intersect(c(id_col, "locus_tag", "protein_id"), names(trd_cols))[1]

    if (!is.na(trd_id_col)) {
      # Rename to match if needed
      if (trd_id_col != id_col && id_col %in% names(rm_combined)) {
        trd_cols <- trd_cols %>% dplyr::rename(!!id_col := !!trd_id_col)
      }

      # Avoid duplicate columns
      existing <- intersect(names(rm_combined), setdiff(names(trd_cols), id_col))
      if (length(existing) > 0) {
        trd_cols <- trd_cols %>% dplyr::select(-dplyr::any_of(existing))
      }

      rm_combined <- rm_combined %>%
        dplyr::left_join(trd_cols, by = id_col)
    }
  }

  # Add catalytic motif info if available
  if (!is.null(catalytic_summary) && is.data.frame(catalytic_summary) && nrow(catalytic_summary) > 0) {
    message("  Joining catalytic motif data (", nrow(catalytic_summary), " rows)...")
    message("    catalytic_summary columns: ", paste(names(catalytic_summary), collapse = ", "))

    # Get the ID column that exists in catalytic_summary
    motif_id_col <- intersect(c(id_col, "locus_tag", "protein_id", "gene_id"), names(catalytic_summary))[1]

    if (!is.na(motif_id_col)) {
      message("    Using ID column: ", motif_id_col)

      # Select columns to join - SIMPLIFIED (positions에 motif 정보 포함)
      motif_cols <- catalytic_summary %>%
        dplyr::select(
          dplyr::all_of(motif_id_col),
          dplyr::any_of(c(
            # ===== REase =====
            "rease_motif_positions",      # ★ Positions with motif names (REQUIRED)
            "rease_nuclease_present",     # Has nuclease motif
            "rease_helicase_present",     # Has Walker + DEAD core
            "rease_helicase_complete",    # Walker + DEAD + (SAT or QxxR)

            # ===== MTase =====
            "mtase_motif_positions",      # ★ Positions with motif names (REQUIRED)
            "mtase_complete",             # Both catalytic + SAM present

            # ===== Confidence =====
            "motif_confidence"            # high/medium/low/none
          ))
        ) %>%
        dplyr::distinct()

      message("    Selected ", ncol(motif_cols), " columns for joining")

      # Rename ID column if needed
      if (motif_id_col != id_col) {
        motif_cols <- motif_cols %>%
          dplyr::rename(!!id_col := !!rlang::sym(motif_id_col))
      }

      # Avoid duplicate columns
      existing <- intersect(names(rm_combined), setdiff(names(motif_cols), id_col))
      if (length(existing) > 0) {
        message("    Removing duplicate columns: ", paste(existing, collapse = ", "))
        motif_cols <- motif_cols %>% dplyr::select(-dplyr::any_of(existing))
      }

      # Join
      n_before <- nrow(rm_combined)
      rm_combined <- rm_combined %>%
        dplyr::left_join(motif_cols, by = id_col)

      # Check what was joined
      n_with_motif <- sum(!is.na(rm_combined$mtase_complete) | !is.na(rm_combined$rease_nuclease_present), na.rm = TRUE)
      message("    Joined: ", n_with_motif, " rows have motif data")

      # Add columns as NA if they don't exist after join (for consistency)
      motif_cols_expected <- c(
        # REase (positions에 motif 정보 포함)
        "rease_motif_positions",
        "rease_nuclease_present", "rease_helicase_present", "rease_helicase_complete",
        # MTase (positions에 motif 정보 포함)
        "mtase_motif_positions",
        "mtase_complete",
        # Confidence
        "motif_confidence"
      )
      for (col in motif_cols_expected) {
        if (!col %in% names(rm_combined)) {
          rm_combined[[col]] <- NA
        }
      }
    } else {
      message("    WARNING: No matching ID column found in catalytic_summary")
    }
  } else {
    message("  No catalytic motif data to join")
    # Add empty motif columns for consistency
    rm_combined$rease_motif_positions <- NA_character_
    rm_combined$rease_nuclease_present <- NA
    rm_combined$rease_helicase_present <- NA
    rm_combined$rease_helicase_complete <- NA
    rm_combined$mtase_motif_positions <- NA_character_
    rm_combined$mtase_complete <- NA
    rm_combined$motif_confidence <- NA_character_
  }

  # ===== REMOVE UNWANTED COLUMNS =====
  # Remove columns that are not needed in final output
  cols_to_remove <- c(
    # From annotated/detection tables - 불필요
    "operon_type", "system_type", "predicted_type",
    "rm_component", "rm_system_type", "rm_type_evidence", "rm_type_scores",
    "re_type", "rease_confidence", "detection_evidence",
    "matched_re_pfam_str", "re_type_pfam", "re_type_keyword",
    "catalytic_motifs", "rm_classification",
    # Detection source columns
    "pfam_detected", "keyword_detected", "motif_detected",
    "pfam_match", "keyword_match", "detection_method",
    "is_type_i_hsd", "is_type_iii", "is_mrr", "is_type_iv",
    # Score columns - 불필요
    "total_score", "mtase_score", "rease_score", "operon_score",
    "trd_score", "rebase_score", "normalized_score", "confidence_level",
    # Old motif column names
    "rease_nuclease_motifs", "rease_helicase_motifs", "rease_motifs",
    "mtase_catalytic_motifs", "mtase_sam_motifs", "mtase_motifs",
    # Redundant REBASE columns (이전 이름들)
    "rebase_type", "rebase_rm_type", "rebase_subunit",
    "rebase_enz_type", "rebase_org",
    # Other
    "PFAM", "pfam", "COG", "EC_number", "old_locus_tag", "protein_id",
    ".rebase_enzyme_clean"
  )
  rm_combined <- rm_combined %>%
    dplyr::select(-dplyr::any_of(cols_to_remove))

  # ===== REORDER COLUMNS =====
  priority_cols <- c(
    # ===== 1. ID & BASIC INFO =====
    id_col,              # locus_tag
    "product",
    "gene",
    "passed_blast_filter",

    # ===== 2. BLAST/REBASE RESULTS =====
    "rm_type",           # ★ REBASE rm_type 원본 (Type I, Type II 등)
    "subunit",           # ★ REBASE subunit 원본 (M, R, S 등)
    "blast_match",       # REBASE enzyme name
    "rec_seq",           # Recognition sequence
    "blast_identity", "blast_length", "blast_evalue", "blast_bitscore",

    # ===== 3. CATALYTIC MOTIFS =====
    "motif_confidence",
    "rease_motif_positions",
    "rease_nuclease_present",
    "rease_helicase_present",
    "rease_helicase_complete",
    "mtase_motif_positions",
    "mtase_complete",

    # ===== 4. OPERON INFO =====
    "operon_id",
    "partner_locus_tag", "partner_product",

    # ===== 5. TRD info =====
    "trd_count", "has_trd", "trd_regions", "trd_motifs",

    # ===== 6. LOCATION =====
    "seqid", "start", "end", "strand"
  )

  existing_priority <- intersect(priority_cols, names(rm_combined))
  other_cols <- setdiff(names(rm_combined), existing_priority)

  rm_combined <- rm_combined %>%
    dplyr::select(dplyr::all_of(existing_priority), dplyr::everything())

  # ===== SORTING =====
  # 1. passed_blast_filter (TRUE first)
  # 2. blast_identity (high first)
  # 3. blast_evalue (low first)
  # 4. operon_id (grouped, NA at the end)
  # 5. locus_tag (alphabetical)

  # Ensure sorting columns exist
  if (!"passed_blast_filter" %in% names(rm_combined)) {
    rm_combined$passed_blast_filter <- FALSE
  }
  if (!"operon_id" %in% names(rm_combined)) {
    rm_combined$operon_id <- NA_character_
  }
  if (!"blast_identity" %in% names(rm_combined)) {
    rm_combined$blast_identity <- NA_real_
  }
  if (!"blast_evalue" %in% names(rm_combined)) {
    rm_combined$blast_evalue <- NA_real_
  }

  # Sort: operon_id로 그룹핑 (같은 operon끼리 붙음) > passed_blast_filter > blast_identity > blast_evalue
  # operon별 max identity 계산해서 operon 그룹 순서 결정
  rm_combined <- rm_combined %>%
    dplyr::group_by(operon_id) %>%
    dplyr::mutate(
      .operon_max_identity = max(blast_identity, na.rm = TRUE),
      .operon_has_passed = any(passed_blast_filter, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      dplyr::desc(.operon_has_passed),     # Operon with passed filter first
      dplyr::desc(.operon_max_identity),   # Operon with higher max identity first
      operon_id,                            # ★ Same operon_id grouped together
      dplyr::desc(blast_identity),          # Within operon: high identity first
      blast_evalue,                         # Then low evalue
      .data[[id_col]]                       # Then by locus_tag
    ) %>%
    dplyr::select(-.operon_max_identity, -.operon_has_passed)

  message(sprintf("Created comprehensive R-M table: %d candidates (%d passed filter), %d columns",
                  nrow(rm_combined), sum(rm_combined$passed_blast_filter, na.rm = TRUE), ncol(rm_combined)))

  # Report operon grouping
  if ("operon_id" %in% names(rm_combined)) {
    n_with_operon <- sum(!is.na(rm_combined$operon_id))
    n_operons <- dplyr::n_distinct(rm_combined$operon_id, na.rm = TRUE)
    message(sprintf("  %d candidates in %d operons (sorted by operon groups)", n_with_operon, n_operons))
  }

  # ★ VERIFY BLAST results are present
  if ("blast_match" %in% names(rm_combined)) {
    n_with_blast <- sum(!is.na(rm_combined$blast_match))
    message(sprintf("  ★ BLAST results: %d candidates have BLAST matches", n_with_blast))

    if ("rm_type" %in% names(rm_combined)) {
      n_with_rmtype <- sum(!is.na(rm_combined$rm_type))
      message(sprintf("  ★ rm_type from REBASE: %d candidates", n_with_rmtype))
    }

    if ("rec_seq" %in% names(rm_combined)) {
      n_with_rec_seq <- sum(!is.na(rm_combined$rec_seq) & rm_combined$rec_seq != "" & rm_combined$rec_seq != "?", na.rm = TRUE)
      message(sprintf("  ★ Recognition sequences: %d candidates have rec_seq", n_with_rec_seq))
    }
  } else {
    message("  WARNING: blast_match column missing from comprehensive table!")
  }

  # ★ VERIFY motif columns
  if ("motif_confidence" %in% names(rm_combined)) {
    message(sprintf("  ★ Motif confidence: high=%d, medium=%d, low=%d, none=%d",
                    sum(rm_combined$motif_confidence == "high", na.rm = TRUE),
                    sum(rm_combined$motif_confidence == "medium", na.rm = TRUE),
                    sum(rm_combined$motif_confidence == "low", na.rm = TRUE),
                    sum(rm_combined$motif_confidence == "none", na.rm = TRUE)))
  }

  # ★ VERIFY motif positions exist
  if ("rease_motif_positions" %in% names(rm_combined)) {
    n_rease_pos <- sum(!is.na(rm_combined$rease_motif_positions), na.rm = TRUE)
    message(sprintf("  ★ REase with motif positions: %d", n_rease_pos))
  }
  if ("mtase_motif_positions" %in% names(rm_combined)) {
    n_mtase_pos <- sum(!is.na(rm_combined$mtase_motif_positions), na.rm = TRUE)
    message(sprintf("  ★ MTase with motif positions: %d", n_mtase_pos))
  }

  return(rm_combined)
}


#' Export Results V2
#'
#' @description Exports all results with better organization
#'
#' @param mtase_annotated Annotated MTase data
#' @param rease_annotated Annotated REase data
#' @param rm_comprehensive Comprehensive table
#' @param operon_data Operon data
#' @param blast_raw Raw BLAST results
#' @param blast_filtered Filtered BLAST results
#' @param recognition_predictions Recognition predictions
#' @param output_dir Output directory
#' @param formats Export formats
#'
#' @return List of exported file paths
#' @export
export_results_v2 <- function(mtase_annotated,
                               rease_annotated,
                               rm_comprehensive,
                               operon_data,
                               blast_raw,
                               blast_filtered,
                               recognition_predictions,
                               output_dir,
                               formats = c("csv", "xlsx")) {

  exported_files <- character()

  # Helper to safely check nrow
  safe_nrow <- function(x) if (!is.null(x) && is.data.frame(x)) nrow(x) else 0

  # === CSV Export ===
  if ("csv" %in% formats) {
    # Main result
    if (safe_nrow(rm_comprehensive) > 0) {
      f <- file.path(output_dir, "rm_comprehensive.csv")
      utils::write.csv(rm_comprehensive, f, row.names = FALSE)
      exported_files <- c(exported_files, f)
      message("  Exported: rm_comprehensive.csv")
    }

    # Annotated tables
    if (safe_nrow(mtase_annotated) > 0) {
      f <- file.path(output_dir, "mtase_annotated.csv")
      utils::write.csv(mtase_annotated, f, row.names = FALSE)
      exported_files <- c(exported_files, f)
    }

    if (safe_nrow(rease_annotated) > 0) {
      f <- file.path(output_dir, "rease_annotated.csv")
      utils::write.csv(rease_annotated, f, row.names = FALSE)
      exported_files <- c(exported_files, f)
    }

    # BLAST results
    if (safe_nrow(blast_raw) > 0) {
      f <- file.path(output_dir, "blast_raw.csv")
      utils::write.csv(blast_raw, f, row.names = FALSE)
      exported_files <- c(exported_files, f)
    }

    if (safe_nrow(blast_filtered) > 0) {
      f <- file.path(output_dir, "blast_filtered.csv")
      utils::write.csv(blast_filtered, f, row.names = FALSE)
      exported_files <- c(exported_files, f)
    }

    # Other results
    if (safe_nrow(operon_data) > 0) {
      f <- file.path(output_dir, "operons.csv")
      utils::write.csv(operon_data, f, row.names = FALSE)
      exported_files <- c(exported_files, f)
    }

    if (safe_nrow(recognition_predictions) > 0) {
      f <- file.path(output_dir, "recognition_sequences.csv")
      utils::write.csv(recognition_predictions, f, row.names = FALSE)
      exported_files <- c(exported_files, f)
    }
  }

  # === Excel Export ===
  if ("xlsx" %in% formats) {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      message("Package 'writexl' needed for Excel export")
    } else {
      sheets <- list()

      # Main sheets
      if (safe_nrow(rm_comprehensive) > 0) sheets$RM_Comprehensive <- rm_comprehensive
      if (safe_nrow(mtase_annotated) > 0) sheets$MTase_Annotated <- mtase_annotated
      if (safe_nrow(rease_annotated) > 0) sheets$REase_Annotated <- rease_annotated
      if (safe_nrow(blast_filtered) > 0) sheets$BLAST_Filtered <- blast_filtered
      if (safe_nrow(operon_data) > 0) sheets$Operons <- operon_data
      if (safe_nrow(recognition_predictions) > 0) sheets$Recognition_Seqs <- recognition_predictions

      # Summary sheet
      sheets$Summary <- dplyr::tibble(
        Category = c("Total MTase candidates", "Total REase candidates",
                     "Passed BLAST filter", "Total operons",
                     "With recognition sequences"),
        Count = c(
          safe_nrow(mtase_annotated),
          safe_nrow(rease_annotated),
          if (safe_nrow(rm_comprehensive) > 0 && "passed_blast_filter" %in% names(rm_comprehensive)) {
            sum(rm_comprehensive$passed_blast_filter, na.rm = TRUE)
          } else 0,
          if (safe_nrow(operon_data) > 0 && "operon_id" %in% names(operon_data)) {
            dplyr::n_distinct(operon_data$operon_id)
          } else 0,
          safe_nrow(recognition_predictions)
        )
      )

      if (length(sheets) > 0) {
        xlsx_file <- file.path(output_dir, "rm_analysis_results.xlsx")
        writexl::write_xlsx(sheets, xlsx_file)
        exported_files <- c(exported_files, xlsx_file)
        message("  Exported: rm_analysis_results.xlsx")
      }
    }
  }

  return(exported_files)
}


#' Print RMscan Results Summary
#'
#' @description Prints formatted summary of analysis results
#'
#' @param results Results from rmscan_pipeline()
#'
#' @export
print_rmscan_summary <- function(results) {

  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("RMscan Analysis Summary\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  safe_nrow <- function(x) if (!is.null(x) && is.data.frame(x)) nrow(x) else 0

  cat("Input Data:\n")
  cat(sprintf("  Total genes: %d\n", safe_nrow(results$dnmb_data)))
  cat("\n")

  cat("Initial Candidates:\n")
  cat(sprintf("  MTases detected: %d\n", safe_nrow(results$mtase_data)))
  cat(sprintf("  REases detected: %d\n", safe_nrow(results$rease_data)))
  cat("\n")

  cat("BLAST Results:\n")
  cat(sprintf("  Raw BLAST hits: %d\n", safe_nrow(results$blast_raw)))
  cat(sprintf("  Filtered BLAST hits: %d\n", safe_nrow(results$blast_filtered)))
  cat(sprintf("  Candidates passing filter: %d\n", length(results$filtered_candidate_ids)))
  cat("\n")

  if (safe_nrow(results$rm_comprehensive) > 0) {
    cat("Comprehensive Table:\n")
    cat(sprintf("  Total candidates: %d\n", nrow(results$rm_comprehensive)))
    if ("passed_blast_filter" %in% names(results$rm_comprehensive)) {
      cat(sprintf("  Passed BLAST filter: %d\n",
                  sum(results$rm_comprehensive$passed_blast_filter, na.rm = TRUE)))
    }
    if ("rec_seq" %in% names(results$rm_comprehensive)) {
      n_with_rec <- sum(!is.na(results$rm_comprehensive$rec_seq) &
                        results$rm_comprehensive$rec_seq != "" &
                        results$rm_comprehensive$rec_seq != "?", na.rm = TRUE)
      cat(sprintf("  With recognition seq: %d\n", n_with_rec))
    }
    cat("\n")
  }

  cat(sprintf("Results saved to: %s\n", results$output_dir))
  if (!is.null(results$elapsed_time)) {
    cat(sprintf("Analysis time: %.1f seconds\n", as.numeric(results$elapsed_time)))
  }
  cat(paste(rep("=", 60), collapse = ""), "\n")

  # Print available data
  cat("\nData in R environment (rmscan_*):\n")
  rmscan_vars <- ls(pattern = "^rmscan_", envir = .GlobalEnv)
  if (length(rmscan_vars) > 0) {
    for (v in rmscan_vars) {
      obj <- get(v, envir = .GlobalEnv)
      if (is.data.frame(obj)) {
        cat(sprintf("  %-30s: %d x %d\n", v, nrow(obj), ncol(obj)))
      }
    }
  }

  cat("\nKey tables:\n")
  cat("  rmscan_mtase_annotated   - MTase + BLAST + operon (중요 열 앞쪽)\n")
  cat("  rmscan_rease_annotated   - REase + BLAST + operon (중요 열 앞쪽)\n")
  cat("  rmscan_rm_comprehensive  - 종합 테이블 (passed_blast_filter 플래그 포함)\n")
  cat("  rmscan_blast_filtered    - 필터링된 BLAST 결과\n")
}


#' Quick R-M Scan (Minimal Analysis)
#'
#' @description Fast analysis using only keyword-based detection
#'
#' @param input_file Path to DNMB file
#' @param verbose Print messages
#'
#' @return Tibble with R-M candidates
#' @export
rmscan_quick <- function(input_file, verbose = TRUE) {

  if (verbose) message("Running quick R-M scan...")

  dnmb_data <- load_dnmb(input_file)
  rm_candidates <- extract_rm_candidates(dnmb_data)

  if (verbose) message(sprintf("Found %d potential R-M genes", nrow(rm_candidates)))

  return(rm_candidates)
}


#' Generate R-M Type Dot Plot
#'
#' @description Creates dot plot for high-identity BLAST hits grouped by rm_type and subunit
#' Style: Figure S2b template with background colors by Type
#'
#' @param rm_data Comprehensive R-M data
#' @param output_dir Output directory
#' @param min_identity Minimum blast_identity threshold (default 0.5)
#' @param id_col ID column name
#' @param verbose Print messages
#'
#' @return Path to generated plot file
#' @export
generate_rm_type_heatmap <- function(rm_data,
                                      output_dir = ".",
                                      min_identity = 0.5,
                                      id_col = "locus_tag",
                                      verbose = TRUE) {

  if (is.null(rm_data) || nrow(rm_data) == 0) {
    if (verbose) message("  No data for dotplot")
    return(NULL)
  }

  # Filter for high identity hits
  high_identity <- rm_data %>%
    dplyr::filter(!is.na(blast_identity) & blast_identity >= min_identity)

  if (nrow(high_identity) == 0) {
    if (verbose) message("  No candidates with blast_identity >= ", min_identity)
    return(NULL)
  }

  # Process subunit and rm_type (subunit은 이미 enzyme_name에서 추출되어 R이 기본값)
  high_identity <- high_identity %>%
    dplyr::mutate(
      subunit_clean = dplyr::case_when(
        is.na(subunit) | subunit == "" ~ "R",
        TRUE ~ toupper(subunit)
      ),
      rm_type_clean = dplyr::case_when(
        is.na(rm_type) | rm_type == "" ~ "Unknown",
        TRUE ~ rm_type
      )
    )

  # Create System column (Type + Subunit)
  high_identity <- high_identity %>%
    dplyr::mutate(
      System = stringr::str_c(rm_type_clean, " (", subunit_clean, ")")
    )

  # ★ 전체 시스템 목록 고정 (항상 이 순서로 표시)
  all_systems <- c(
    "Type I (M)", "Type I (R)", "Type I (S)",
    "Type II (M)", "Type II (R)",
    "Type III (M)", "Type III (R)",
    "Type IV (R)"
  )

  # Count by System
  count_data <- high_identity %>%
    dplyr::count(rm_type_clean, subunit_clean, System, name = "Count")

  # 전체 시스템에 대해 0 카운트 추가 (존재하지 않는 시스템도 표시)
  all_systems_df <- dplyr::tibble(
    System = all_systems,
    rm_type_clean = stringr::str_extract(all_systems, "Type I{1,3}|Type IV"),
    subunit_clean = stringr::str_extract(all_systems, "(?<=\\()[MRS](?=\\))")
  )

  # 기존 데이터와 merge (없는 시스템은 Count = 0)
  count_data <- all_systems_df %>%
    dplyr::left_join(
      count_data %>% dplyr::select(System, Count),
      by = "System"
    ) %>%
    dplyr::mutate(
      Count = dplyr::if_else(is.na(Count), 0L, Count),
      System = factor(System, levels = all_systems),
      subunit_clean = factor(subunit_clean, levels = c("M", "R", "S"))
    )

  if (verbose) {
    n_with_count <- sum(count_data$Count > 0)
    message(sprintf("  Dotplot data: %d systems with data (of %d total)",
                    n_with_count, length(all_systems)))
    message("    Systems: ", paste(all_systems, collapse = ", "))
  }

  # ★ 고정된 색상 (Subunit별 - Type IV도 R subunit으로 취급)
  subunit_colors <- c(
    "M" = "#6A5ACD",
    "R" = "#E8B923",
    "S" = "#87CEEB"
  )

  # ★ 고정된 배경 위치 (Type별)
  # Type I: 1-3, Type II: 4-5, Type III: 6-7, Type IV: 8
  bg_rects <- dplyr::tibble(
    rm_type = c("Type I", "Type II", "Type III", "Type IV"),
    xmin = c(0.5, 3.5, 5.5, 7.5),
    xmax = c(3.5, 5.5, 7.5, 8.5),
    fill = c("#FF9800", "#4CAF50", "#2196F3", "#9C27B0")
  )

  n_systems <- length(all_systems)

  # ★ Count > 0인 것만 점으로 표시
  plot_data <- count_data %>% dplyr::filter(Count > 0)

  # Build plot
  p <- ggplot2::ggplot(count_data, ggplot2::aes(x = System, y = 1))

  # ★ 고정된 배경 사각형 추가
  for (i in seq_len(nrow(bg_rects))) {
    p <- p + ggplot2::annotate(
      "rect",
      xmin = bg_rects$xmin[i], xmax = bg_rects$xmax[i],
      ymin = 0.5, ymax = 1.5,
      fill = bg_rects$fill[i], alpha = 0.15
    )
  }

  # ★ Count > 0인 것만 점 추가
  if (nrow(plot_data) > 0) {
    max_count <- max(plot_data$Count, na.rm = TRUE)
    p <- p +
      ggplot2::geom_point(
        data = plot_data,
        ggplot2::aes(size = Count, fill = subunit_clean),
        shape = 21,
        stroke = 0.8,
        color = "black"
      ) +
      ggplot2::scale_size_continuous(
        name = "Count",
        range = c(4, 14),
        breaks = seq(1, max_count, by = max(1, floor(max_count/4))),
        limits = c(1, max_count)
      ) +
      ggplot2::scale_fill_manual(
        name = "Subunit",
        values = subunit_colors,
        labels = c(
          "M" = "M (Methyltransferase)",
          "R" = "R (Restriction enzyme)",
          "S" = "S (Specificity)"
        ),
        guide = ggplot2::guide_legend(override.aes = list(size = 6, stroke = 0.8))
      )
  }

  p <- p +
    ggplot2::scale_x_discrete(drop = FALSE) +  # ★ 빈 레벨도 표시
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::labs(
      title = paste0("R-M System Distribution (BLAST identity >= ", min_identity * 100, "%)")
    ) +
    ggplot2::coord_cartesian(ylim = c(0.5, 1.5))

  # Save plot (8 시스템 고정, 정사각형에 가까운 비율)
  plot_file_pdf <- file.path(output_dir, "RM_system_dotplot.pdf")
  plot_file_png <- file.path(output_dir, "RM_system_dotplot.png")

  tryCatch({
    ggplot2::ggsave(plot_file_pdf, p, width = 8, height = 1.5)
    ggplot2::ggsave(plot_file_png, p, width = 8, height = 1.5, dpi = 300)

    if (verbose) {
      message("  Saved dotplot: ", plot_file_pdf)
      message("  Saved dotplot: ", plot_file_png)
    }

    return(plot_file_pdf)

  }, error = function(e) {
    warning("Dotplot generation failed: ", e$message)
    return(NULL)
  })
}


#' Export R-M Comprehensive Table to xlsx
#'
#' @description Exports comprehensive table to Excel using xlsx package
#'
#' @param rm_data Comprehensive R-M data
#' @param output_file Output xlsx file path
#' @param verbose Print messages
#'
#' @return Path to exported file
#' @export
export_rm_xlsx <- function(rm_data, output_file, verbose = TRUE) {

  if (is.null(rm_data) || nrow(rm_data) == 0) {
    if (verbose) message("  No data to export")
    return(NULL)
  }

  # Try xlsx package first, then writexl as fallback
  if (requireNamespace("xlsx", quietly = TRUE)) {
    tryCatch({
      # Create workbook
      wb <- xlsx::createWorkbook()

      # Sheet 1: All data
      sheet1 <- xlsx::createSheet(wb, sheetName = "RM_Comprehensive")
      xlsx::addDataFrame(as.data.frame(rm_data), sheet1, row.names = FALSE)

      # Sheet 2: High identity hits (>= 50%)
      high_identity <- rm_data %>%
        dplyr::filter(!is.na(blast_identity) & blast_identity >= 0.5)

      if (nrow(high_identity) > 0) {
        sheet2 <- xlsx::createSheet(wb, sheetName = "High_Identity_50pct")
        xlsx::addDataFrame(as.data.frame(high_identity), sheet2, row.names = FALSE)
      }

      # Sheet 3: Summary by rm_type
      if ("rm_type" %in% names(rm_data)) {
        type_summary <- rm_data %>%
          dplyr::filter(!is.na(blast_identity) & blast_identity >= 0.5) %>%
          dplyr::mutate(
            subunit_clean = dplyr::case_when(
              is.na(subunit) | subunit == "" | tolower(subunit) == "unknown" ~ "R",
              TRUE ~ subunit
            )
          ) %>%
          dplyr::count(rm_type, subunit_clean, name = "count") %>%
          dplyr::arrange(rm_type, subunit_clean)

        if (nrow(type_summary) > 0) {
          sheet3 <- xlsx::createSheet(wb, sheetName = "Type_Summary")
          xlsx::addDataFrame(as.data.frame(type_summary), sheet3, row.names = FALSE)
        }
      }

      # Save workbook
      xlsx::saveWorkbook(wb, output_file)

      if (verbose) message("  Exported: ", output_file, " (xlsx package)")
      return(output_file)

    }, error = function(e) {
      warning("xlsx export failed: ", e$message, ". Trying writexl...")
    })
  }

  # Fallback to writexl
  if (requireNamespace("writexl", quietly = TRUE)) {
    tryCatch({
      sheets <- list(
        RM_Comprehensive = rm_data
      )

      # Add high identity sheet
      high_identity <- rm_data %>%
        dplyr::filter(!is.na(blast_identity) & blast_identity >= 0.5)
      if (nrow(high_identity) > 0) {
        sheets$High_Identity_50pct <- high_identity
      }

      # Add type summary
      if ("rm_type" %in% names(rm_data)) {
        type_summary <- rm_data %>%
          dplyr::filter(!is.na(blast_identity) & blast_identity >= 0.5) %>%
          dplyr::mutate(
            subunit_clean = dplyr::case_when(
              is.na(subunit) | subunit == "" | tolower(subunit) == "unknown" ~ "R",
              TRUE ~ subunit
            )
          ) %>%
          dplyr::count(rm_type, subunit_clean, name = "count")

        if (nrow(type_summary) > 0) {
          sheets$Type_Summary <- type_summary
        }
      }

      writexl::write_xlsx(sheets, output_file)

      if (verbose) message("  Exported: ", output_file, " (writexl package)")
      return(output_file)

    }, error = function(e) {
      warning("writexl export also failed: ", e$message)
      return(NULL)
    })
  } else {
    warning("Neither 'xlsx' nor 'writexl' package available for Excel export")
    return(NULL)
  }
}

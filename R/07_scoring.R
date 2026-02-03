#' @title R-M System Scoring Module
#' @description Functions for comprehensive scoring of R-M system predictions
#' @name rm_scoring
NULL


# ============================================================
# Catalytic Motif Detection
# ============================================================

#' Get R-M Catalytic Motifs
#'
#' @description Returns known catalytic motifs for R-M system components.
#'
#' Key R-M motifs:
#' - PDxDxK: Nuclease catalytic site (REase)
#' - DPPY/DPPF/DPPW: Methyltransferase catalytic motif (MTase)
#' - PxGxG: SAM-binding motif (MTase)
#'
#' For MTase confirmation, BOTH DPPY variants AND PxGxG should be present.
#'
#' @param type Motif type: "all", "mtase", "rease"
#'
#' @return Tibble with motif names, patterns, and descriptions
#' @export
get_rm_catalytic_motifs <- function(type = "all") {

  motifs <- dplyr::tibble(
    motif_name = c(
      # ===== REase (Nuclease) motifs - Type II =====
      "PD-(D/E)xK",
      "PD-ExK",
      "HNH",
      "GIY-YIG",

      # ===== REase (Type I HsdR) - Helicase motifs =====
      # Type I REase has DEAD-box helicase domain (SF2 superfamily)
      "GTGKT",      # Walker A (P-loop) - ATP binding
      "DEAD",       # Walker B variant (DEAD-box)
      "DEAH",       # Walker B variant (DEAH-box)
      "DExx",       # Walker B - ATP hydrolysis (general)
      "SAT/SAH",    # Motif III - DNA translocation coupling (SAT, SAH, SAA, etc.)
      "TAN",        # Motif III variant - common in some helicases
      "QxxR",       # Motif VI - ATP binding/hydrolysis sensor
      "ARGID",      # Motif V - found in SF2 helicases

      # ===== MTase catalytic motifs (Motif IV) =====
      # Pattern: [D/N/S]PP[Y/F/W/H]
      "DPPY",   # D-PP-Y (classic m6A/m4C)
      "DPPF",   # D-PP-F
      "DPPW",   # D-PP-W
      "DPPH",   # D-PP-H
      "NPPY",   # N-PP-Y (m5C type)
      "NPPF",   # N-PP-F
      "NPPW",   # N-PP-W
      "NPPH",   # N-PP-H
      "SPPY",   # S-PP-Y (rare)
      "SPPF",   # S-PP-F
      "DPY_single",  # D-P-Y (single P)
      "NPY_single",  # N-P-Y (single P)

      # ===== MTase SAM-binding motifs (Motif I) =====
      # Pattern: [F/G/P/L/V/I]xGxG
      "FxGxG",
      "GxGxG",
      "LxGxG",
      "VxGxG",
      "IxGxG"
    ),
    # Regex patterns
    pattern = c(
      # Type II REase patterns - PD-(D/E)xK superfamily
      # PD and (D/E)xK are typically separated by 5-30 residues (Steczkiewicz et al. 2012)
      "PD.{5,30}[DE].K",         # PD--(5-30aa)--(D/E)-x-K (canonical spacing)
      "H.{0,4}N.{0,4}H",         # HNH endonuclease
      "G[IL]Y.{5,50}Y[IL]G",     # GIY-YIG nuclease

      # Type I REase (HsdR) helicase motifs - SF2 helicase conserved motifs
      # Walker A (P-loop): GxGKT/S where x is typically T, S, A, or V
      "G[TSAV]GK[TS]",           # Walker A: GTGKT, GSGKS, etc.
      # Walker B: hhhhDE (classic) or DExxH/D for SF2 helicases
      "DE.{1,2}[HD]",            # Walker B: DExxH or DExxD (SF2 helicase core)
      "DEAD",                    # DEAD-box (exact) - SF2 helicase subgroup
      "DEAH",                    # DEAH-box (exact) - SF2 helicase subgroup
      "S[ATV][THPSAQK]",         # Motif III: SAT, SAH, SAA, SAQ, etc.
      "T[AT][NT]",               # Motif III TAN variant: TAN, TAT, TTN, etc.
      "Q.{2}R",                  # Motif VI: QxxR
      "[AR]G[IL]D",              # Motif V: ARGID or RGID

      # MTase Motif IV: [D/N/S]PP[Y/F/W/H]
      "DPPY",
      "DPPF",
      "DPPW",
      "DPPH",
      "NPPY",
      "NPPF",
      "NPPW",
      "NPPH",
      "SPPY",
      "SPPF",
      "DP[YFW]",                  # Single P variant
      "NP[YFW]",

      # MTase Motif I (SAM-binding)
      "F[AGVLIST]G[AGVLIST]G",
      "G[AGVLIST]G[AGVLIST]G",
      "L[AGVLIST]G[AGVLIST]G",
      "V[AGVLIST]G[AGVLIST]G",
      "I[AGVLIST]G[AGVLIST]G"
    ),
    motif_type = c(
      # Type II REase (3 motifs: PD-(D/E)xK, HNH, GIY-YIG)
      "rease_typeII", "rease_typeII", "rease_typeII",
      # Type I REase helicase (8 motifs: Walker A, Walker B, DEAD, DEAH, SAT, TAN, QxxR, ARGID)
      rep("rease_typeI_helicase", 8),
      # MTase catalytic (12 motifs)
      rep("mtase_catalytic", 12),
      # MTase SAM-binding (5 motifs)
      rep("mtase_sam", 5)
    ),
    component = c(
      rep("REase", 11),   # 3 Type II + 8 helicase
      rep("MTase", 17)    # 12 catalytic + 5 SAM
    ),
    description = c(
      # Type II REase (3 motifs)
      "Type II REase PD-(D/E)xK catalytic site (5-30aa spacing)",
      "HNH endonuclease domain",
      "GIY-YIG nuclease domain",

      # Type I REase helicase (8 motifs) - order matches pattern array
      "Type I REase Walker A (GxGKT) - ATP binding P-loop",
      "Type I REase Walker B (DExxH/D) - ATP hydrolysis",
      "Type I REase DEAD-box helicase",
      "Type I REase DEAH-box helicase",
      "Type I REase Motif III (SAT/SAH) - DNA translocation",
      "Type I REase Motif III (TAN) variant",
      "Type I REase Motif VI (QxxR) - RNA/DNA binding",
      "Type I REase Motif V (ARGID)",

      # MTase Motif IV
      "MTase Motif IV DPPY (m6A/m4C)",
      "MTase Motif IV DPPF",
      "MTase Motif IV DPPW",
      "MTase Motif IV DPPH",
      "MTase Motif IV NPPY (m5C)",
      "MTase Motif IV NPPF",
      "MTase Motif IV NPPW",
      "MTase Motif IV NPPH",
      "MTase Motif IV SPPY",
      "MTase Motif IV SPPF",
      "MTase Motif IV DPY (single P)",
      "MTase Motif IV NPY (single P)",

      # MTase Motif I
      "MTase Motif I FxGxG (SAM)",
      "MTase Motif I GxGxG (SAM)",
      "MTase Motif I LxGxG (SAM)",
      "MTase Motif I VxGxG (SAM)",
      "MTase Motif I IxGxG (SAM)"
    ),
    essential = c(
      rep(TRUE, 12),  # REase motifs (4 Type II + 8 helicase)
      rep(TRUE, 12),  # MTase Motif IV (catalytic)
      rep(TRUE, 5)    # MTase Motif I (SAM-binding)
    )
  )

  if (type == "mtase") {
    motifs <- dplyr::filter(motifs, component == "MTase")
  } else if (type == "rease") {
    motifs <- dplyr::filter(motifs, component == "REase")
  }

  return(motifs)
}


#' Detect Catalytic Motifs in Protein Sequences
#'
#' @description Searches for R-M catalytic motifs in protein sequences
#' and returns positions (start, end) for each match.
#'
#' For MTase confirmation:
#' - Must have at least one DPPY-family motif (catalytic site)
#' - Must have at least one GxGxG-family motif (SAM-binding)
#' - Both are required for high confidence MTase classification
#'
#' @param data Data frame with protein sequences
#' @param seq_col Column containing protein sequences
#' @param id_col ID column (auto-detect if NULL)
#' @param component_type "MTase", "REase", or "all"
#'
#' @return Tibble with motif detection results including positions
#' @export
detect_catalytic_motifs <- function(data,
                                     seq_col = "translation",
                                     id_col = NULL,
                                     component_type = "all") {

  if (is.null(data) || nrow(data) == 0) {
    return(dplyr::tibble())
  }

  # Get ID column
  if (is.null(id_col)) {
    id_col <- intersect(c("locus_tag", "protein_id"), names(data))[1]
  }

  if (!seq_col %in% names(data)) {
    warning("Sequence column not found: ", seq_col)
    return(dplyr::tibble())
  }

  # Get motifs for the specified component type
  motifs <- get_rm_catalytic_motifs(type = tolower(component_type))

  # Process each sequence
  results <- list()

  for (i in seq_len(nrow(data))) {
    gene_id <- data[[id_col]][i]
    sequence <- data[[seq_col]][i]

    if (is.na(sequence) || nchar(sequence) < 20) next

    sequence <- toupper(sequence)

    # Search for each motif
    gene_motifs <- list()

    for (j in seq_len(nrow(motifs))) {
      motif_name <- motifs$motif_name[j]
      pattern <- motifs$pattern[j]
      motif_type <- motifs$motif_type[j]

      # Find all matches with positions
      matches <- gregexpr(pattern, sequence, perl = TRUE)[[1]]

      if (matches[1] != -1) {
        match_lengths <- attr(matches, "match.length")

        for (k in seq_along(matches)) {
          start_pos <- matches[k]
          end_pos <- start_pos + match_lengths[k] - 1
          matched_seq <- substr(sequence, start_pos, end_pos)

          gene_motifs[[length(gene_motifs) + 1]] <- dplyr::tibble(
            gene_id = gene_id,
            motif_name = motif_name,
            motif_type = motif_type,
            component = motifs$component[j],
            matched_sequence = matched_seq,
            start_pos = start_pos,
            end_pos = end_pos,
            essential = motifs$essential[j]
          )
        }
      }
    }

    if (length(gene_motifs) > 0) {
      results[[length(results) + 1]] <- dplyr::bind_rows(gene_motifs)
    }
  }

  if (length(results) == 0) {
    message("No catalytic motifs found")
    return(dplyr::tibble())
  }

  result <- dplyr::bind_rows(results)

  # Rename gene_id to match input id_col
  if (id_col != "gene_id") {
    result <- result %>%
      dplyr::rename(!!id_col := gene_id)
  }

  message(sprintf("Found %d catalytic motif matches in %d sequences",
                  nrow(result), dplyr::n_distinct(result[[id_col]])))

  return(result)
}


#' Summarize Catalytic Motifs per Gene
#'
#' @description Creates a summary table with motif info per gene.
#' Useful for adding motif columns to comprehensive table.
#'
#' @param motif_data Results from detect_catalytic_motifs()
#' @param id_col ID column name
#'
#' @return Tibble with one row per gene, summarizing motif findings
#' @export
summarize_catalytic_motifs <- function(motif_data, id_col = "locus_tag") {

  if (is.null(motif_data) || nrow(motif_data) == 0) {
    return(dplyr::tibble())
  }

  # Get the actual id column used
  if (!id_col %in% names(motif_data)) {
    id_col <- intersect(c("locus_tag", "protein_id", "gene_id"), names(motif_data))[1]
  }

  summary_data <- motif_data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      # ===== INTERNAL: Count motifs by type (for evaluation) =====
      .n_nuclease = sum(motif_type == "rease_typeII", na.rm = TRUE),
      .n_helicase = sum(motif_type == "rease_typeI_helicase", na.rm = TRUE),
      .n_mtase_cat = sum(motif_type == "mtase_catalytic", na.rm = TRUE),
      .n_mtase_sam = sum(motif_type == "mtase_sam", na.rm = TRUE),

      # INTERNAL: Specific helicase motif detection (for evaluation)
      .has_walker = any(motif_name == "GTGKT", na.rm = TRUE),
      .has_dead = any(motif_name %in% c("DEAD", "DEAH", "DExx"), na.rm = TRUE),
      .has_sat = any(motif_name %in% c("SAT/SAH", "TAN"), na.rm = TRUE),
      .has_qxxr = any(motif_name == "QxxR", na.rm = TRUE),
      .has_argid = any(motif_name == "ARGID", na.rm = TRUE),

      # ===== REase OUTPUT COLUMNS =====
      # Nuclease motifs (PD-D/E-xK, HNH, GIY-YIG)
      rease_nuclease_motifs = paste(unique(motif_name[motif_type == "rease_typeII"]), collapse = ";"),

      # Walker/DEAD helicase motifs (GTGKT, DEAD, DExx, SAT, QxxR)
      rease_helicase_motifs = paste(unique(motif_name[motif_type == "rease_typeI_helicase"]), collapse = ";"),

      # Combined REase motifs (all)
      rease_motifs = paste(unique(motif_name[component == "REase"]), collapse = ";"),

      # ★ REase motif positions (REQUIRED - all positions combined)
      rease_motif_positions = paste(
        paste0(matched_sequence[component == "REase"], ":",
               start_pos[component == "REase"], "-",
               end_pos[component == "REase"]),
        collapse = ";"
      ),

      # ===== MTase OUTPUT COLUMNS =====
      # Catalytic motifs (DPPY family - Motif IV)
      mtase_catalytic_motifs = paste(unique(motif_name[motif_type == "mtase_catalytic"]), collapse = ";"),

      # SAM-binding motifs (Motif I - xGxG)
      mtase_sam_motifs = paste(unique(motif_name[motif_type == "mtase_sam"]), collapse = ";"),

      # Combined MTase motifs
      mtase_motifs = paste(unique(motif_name[component == "MTase"]), collapse = ";"),

      # ★ MTase motif positions (REQUIRED - all positions combined)
      mtase_motif_positions = paste(
        paste0(matched_sequence[component == "MTase"], ":",
               start_pos[component == "MTase"], "-",
               end_pos[component == "MTase"]),
        collapse = ";"
      ),

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # ===== EVALUATION FLAGS =====
      # Core helicase = Walker + DEAD/DEAH (both required)
      .core_helicase = .has_walker & .has_dead,

      # MTase complete: BOTH catalytic (DPPY) AND SAM-binding (GxGxG)
      mtase_complete = .n_mtase_cat > 0 & .n_mtase_sam > 0,

      # REase complete: has nuclease motif
      rease_nuclease_present = .n_nuclease > 0,

      # REase helicase: Walker + DEAD/DEAH core
      rease_helicase_present = .core_helicase,

      # REase helicase complete: Walker + DEAD + (SAT or QxxR or ARGID)
      rease_helicase_complete = .core_helicase & (.has_sat | .has_qxxr | .has_argid),

      # ===== MOTIF CONFIDENCE =====
      motif_confidence = dplyr::case_when(
        # HIGH: Complete sets
        .n_mtase_cat > 0 & .n_mtase_sam > 0 ~ "high",          # MTase complete
        .core_helicase & .n_nuclease > 0 ~ "high",             # REase helicase + nuclease
        .n_nuclease > 0 & !.has_walker & !.has_dead ~ "high",  # Pure nuclease (Type II)

        # MEDIUM: Partial
        .core_helicase ~ "medium",                              # Helicase core only
        .n_mtase_cat > 0 | .n_mtase_sam > 0 ~ "medium",        # Partial MTase
        .n_nuclease > 0 ~ "medium",                             # Nuclease with partial helicase

        # LOW: Single motifs
        .has_walker | .has_dead ~ "low",                        # Single helicase motif
        .has_qxxr | .has_sat | .has_argid ~ "low",              # Accessory only
        TRUE ~ "none"
      ),

      # ===== CLEAN UP =====
      # Remove empty strings -> NA
      rease_nuclease_motifs = dplyr::if_else(rease_nuclease_motifs == "", NA_character_, rease_nuclease_motifs),
      rease_helicase_motifs = dplyr::if_else(rease_helicase_motifs == "", NA_character_, rease_helicase_motifs),
      rease_motifs = dplyr::if_else(rease_motifs == "", NA_character_, rease_motifs),
      rease_motif_positions = dplyr::if_else(rease_motif_positions == ":" | rease_motif_positions == "", NA_character_, rease_motif_positions),
      mtase_catalytic_motifs = dplyr::if_else(mtase_catalytic_motifs == "", NA_character_, mtase_catalytic_motifs),
      mtase_sam_motifs = dplyr::if_else(mtase_sam_motifs == "", NA_character_, mtase_sam_motifs),
      mtase_motifs = dplyr::if_else(mtase_motifs == "", NA_character_, mtase_motifs),
      mtase_motif_positions = dplyr::if_else(mtase_motif_positions == ":" | mtase_motif_positions == "", NA_character_, mtase_motif_positions)
    ) %>%
    # Remove internal columns (starting with .)
    dplyr::select(-dplyr::starts_with("."))

  return(summary_data)
}


# ============================================================
# R-M System Type Classification (Type I/II/III/IV)
# ============================================================

#' Get R-M System Type Keywords
#'
#' @description Returns keyword patterns for classifying R-M system types
#' @return List of keyword patterns by type
#' @export
get_rm_type_keywords <- function() {
  list(
    # ===== Type I =====
    type_I = list(
      strong = c(
        "type I restriction", "type I R-M", "type I RM",
        "HsdR", "HsdM", "HsdS",
        "type I restriction enzyme",
        "type I modification methylase",
        "type I specificity"
      ),
      subunits = c(
        "restriction enzyme R subunit", "R subunit",
        "modification enzyme M subunit", "M subunit",
        "specificity subunit", "S subunit",
        "DNA specificity domain"
      )
    ),

    # ===== Type II =====
    type_II = list(
      strong = c(
        "type II restriction", "type II R-M", "type II RM",
        "type II restriction enzyme",
        "type II methyltransferase",
        "type IIP", "type IIS", "type IIG", "type IIB"
      ),
      # Common Type II enzyme name patterns
      enzyme_names = c(
        "EcoRI", "BamHI", "HindIII", "PstI", "SalI", "XbaI", "XhoI",
        "NotI", "SmaI", "KpnI", "SacI", "SpeI", "NheI", "AvrII",
        "BglII", "NcoI", "NdeI", "EcoRV", "HpaI", "ClaI", "MluI"
      )
    ),

    # ===== Type III =====
    type_III = list(
      strong = c(
        "type III restriction", "type III R-M", "type III RM",
        "type III restriction enzyme",
        "type III modification"
      ),
      subunits = c(
        "Mod subunit", "modification subunit",
        "Res subunit", "restriction subunit",
        "type III mod", "type III res"
      )
    ),

    # ===== Type IV (Methylation-dependent) =====
    type_IV = list(
      strong = c(
        "type IV restriction", "type IV R-M",
        "methyl-dependent restriction",
        "methylation-dependent restriction",
        "modified DNA restriction"
      ),
      enzymes = c(
        "McrA", "McrB", "McrC", "McrBC",
        "Mrr", "MrrA",
        "GmrSD", "GmrS", "GmrD"
      ),
      keywords = c(
        "methylated cytosine restriction",
        "hydroxymethylcytosine restriction",
        "modified base restriction"
      )
    )
  )
}


#' Classify R-M System Type
#'
#' @description Classifies R-M system components into Type I, II, III, or IV
#' based on multiple evidence sources:
#' - Product annotation keywords
#' - Catalytic motif presence
#' - Protein size
#' - Operon context
#'
#' @param data Data frame with gene annotations (must have 'product' column)
#' @param motif_summary Catalytic motif summary from summarize_catalytic_motifs()
#' @param operon_data Operon analysis results (optional)
#' @param id_col ID column name
#'
#' @return Data frame with rm_system_type and rm_type_evidence columns added
#' @export
classify_rm_system_type <- function(data,
                                     motif_summary = NULL,
                                     operon_data = NULL,
                                     id_col = "locus_tag") {

  if (is.null(data) || nrow(data) == 0) {
    return(data)
  }

  # Get ID column
  if (!id_col %in% names(data)) {
    id_col <- intersect(c("locus_tag", "protein_id", "gene_id"), names(data))[1]
  }

  keywords <- get_rm_type_keywords()

  # Initialize type scores for each gene
  data <- data %>%
    dplyr::mutate(
      # Score columns
      .type_I_score = 0L,
      .type_II_score = 0L,
      .type_III_score = 0L,
      .type_IV_score = 0L,
      .type_evidence = ""
    )

  # Get product column (case-insensitive)
  product_col <- intersect(c("product", "Product", "description"), names(data))[1]
  if (is.na(product_col)) product_col <- "product"

  # ===== 1. KEYWORD-BASED SCORING =====
  for (i in seq_len(nrow(data))) {
    product <- tolower(data[[product_col]][i])
    if (is.na(product)) product <- ""

    evidence <- c()

    # ----- Type I -----
    # Strong keywords (+30)
    if (any(sapply(tolower(keywords$type_I$strong), function(k) grepl(k, product, fixed = TRUE)))) {
      data$.type_I_score[i] <- data$.type_I_score[i] + 30L
      evidence <- c(evidence, "TypeI_keyword")
    }
    # Subunit keywords (+20)
    if (any(sapply(tolower(keywords$type_I$subunits), function(k) grepl(k, product, fixed = TRUE)))) {
      data$.type_I_score[i] <- data$.type_I_score[i] + 20L
      evidence <- c(evidence, "TypeI_subunit")
    }
    # HsdR/M/S specific (+25)
    if (grepl("hsdr|hsdm|hsds", product)) {
      data$.type_I_score[i] <- data$.type_I_score[i] + 25L
      evidence <- c(evidence, "Hsd_subunit")
    }

    # ----- Type II -----
    # Strong keywords (+30)
    if (any(sapply(tolower(keywords$type_II$strong), function(k) grepl(k, product, fixed = TRUE)))) {
      data$.type_II_score[i] <- data$.type_II_score[i] + 30L
      evidence <- c(evidence, "TypeII_keyword")
    }
    # Known enzyme names (+25)
    if (any(sapply(keywords$type_II$enzyme_names, function(k) grepl(k, product, ignore.case = TRUE)))) {
      data$.type_II_score[i] <- data$.type_II_score[i] + 25L
      evidence <- c(evidence, "TypeII_enzyme_name")
    }

    # ----- Type III -----
    # Strong keywords (+30)
    if (any(sapply(tolower(keywords$type_III$strong), function(k) grepl(k, product, fixed = TRUE)))) {
      data$.type_III_score[i] <- data$.type_III_score[i] + 30L
      evidence <- c(evidence, "TypeIII_keyword")
    }
    # Mod/Res subunit (+25)
    if (any(sapply(tolower(keywords$type_III$subunits), function(k) grepl(k, product, fixed = TRUE)))) {
      data$.type_III_score[i] <- data$.type_III_score[i] + 25L
      evidence <- c(evidence, "TypeIII_subunit")
    }
    # "mod" or "res" at word boundary (+15)
    if (grepl("\\bmod\\b|\\bres\\b", product)) {
      data$.type_III_score[i] <- data$.type_III_score[i] + 15L
      evidence <- c(evidence, "Mod_Res_keyword")
    }

    # ----- Type IV -----
    # Strong keywords (+35)
    if (any(sapply(tolower(keywords$type_IV$strong), function(k) grepl(k, product, fixed = TRUE)))) {
      data$.type_IV_score[i] <- data$.type_IV_score[i] + 35L
      evidence <- c(evidence, "TypeIV_keyword")
    }
    # Specific enzymes (+40) - very strong indicator
    if (any(sapply(tolower(keywords$type_IV$enzymes), function(k) grepl(k, product, fixed = TRUE)))) {
      data$.type_IV_score[i] <- data$.type_IV_score[i] + 40L
      evidence <- c(evidence, "TypeIV_enzyme")
    }
    # Methylation-dependent keywords (+25)
    if (any(sapply(tolower(keywords$type_IV$keywords), function(k) grepl(k, product, fixed = TRUE)))) {
      data$.type_IV_score[i] <- data$.type_IV_score[i] + 25L
      evidence <- c(evidence, "TypeIV_methyldep")
    }

    data$.type_evidence[i] <- paste(evidence, collapse = ";")
  }

  # ===== 2. MOTIF-BASED SCORING =====
  if (!is.null(motif_summary) && nrow(motif_summary) > 0) {
    # Get matching ID column in motif_summary
    motif_id_col <- intersect(c(id_col, "locus_tag", "protein_id", "gene_id"), names(motif_summary))[1]

    if (!is.na(motif_id_col)) {
      for (i in seq_len(nrow(data))) {
        gene_id <- data[[id_col]][i]
        motif_row <- motif_summary[motif_summary[[motif_id_col]] == gene_id, ]

        if (nrow(motif_row) > 0) {
          evidence <- strsplit(data$.type_evidence[i], ";")[[1]]

          # ===== MOTIF COMPLETENESS SCORING =====
          # Use NEW column names from simplified summarize_catalytic_motifs()
          # Type I & III REase: BOTH helicase AND nuclease required
          # Type II REase: nuclease ONLY (no helicase)

          has_helicase <- isTRUE(motif_row$rease_helicase_present[1])
          has_nuclease <- isTRUE(motif_row$rease_nuclease_present[1])
          helicase_complete <- isTRUE(motif_row$rease_helicase_complete[1])
          product <- tolower(data[[product_col]][i])
          if (is.na(product)) product <- ""

          # Check for Type I/III keywords
          is_type_III_keyword <- grepl("\\bres\\b|type.?iii|type.?3|\\bmod\\b", product)
          is_type_I_keyword <- grepl("hsdr|hsdm|hsds|type.?i[^iv]|type.?1\\b", product)

          if (has_helicase && has_nuclease) {
            # ★ BOTH helicase + nuclease = Type I or Type III
            evidence <- c(evidence, "helicase+nuclease")

            if (is_type_III_keyword) {
              # Type III Res: helicase + nuclease + "Res" keyword
              data$.type_III_score[i] <- data$.type_III_score[i] + 40L
              evidence <- c(evidence, "TypeIII_Res_complete")
            } else if (is_type_I_keyword) {
              # Type I HsdR: helicase + nuclease + "Hsd" keyword
              data$.type_I_score[i] <- data$.type_I_score[i] + 40L
              evidence <- c(evidence, "TypeI_HsdR_complete")
            } else {
              # Has both but no clear keyword → slight favor to Type I (more common)
              data$.type_I_score[i] <- data$.type_I_score[i] + 25L
              data$.type_III_score[i] <- data$.type_III_score[i] + 20L
              evidence <- c(evidence, "helicase+nuclease_ambiguous")
            }

            # Additional points for helicase completeness (Walker + DEAD + SAT/QxxR)
            if (helicase_complete) {
              data$.type_I_score[i] <- data$.type_I_score[i] + 10L
              data$.type_III_score[i] <- data$.type_III_score[i] + 10L
              evidence <- c(evidence, "helicase_complete")
            }

          } else if (has_helicase && !has_nuclease) {
            # Helicase only (partial - unusual for R-M)
            evidence <- c(evidence, "helicase_only")
            data$.type_I_score[i] <- data$.type_I_score[i] + 15L
            evidence <- c(evidence, "helicase_partial")

          } else if (!has_helicase && has_nuclease) {
            # ★ Nuclease ONLY (no helicase) = Type II
            data$.type_II_score[i] <- data$.type_II_score[i] + 30L
            evidence <- c(evidence, "TypeII_nuclease_only")
          }

          # MTase motif completeness (NEW: mtase_complete instead of mtase_motif_complete)
          if (isTRUE(motif_row$mtase_complete[1])) {
            evidence <- c(evidence, "mtase_complete")
          }

          data$.type_evidence[i] <- paste(unique(evidence[evidence != ""]), collapse = ";")
        }
      }
    }
  }

  # ===== 3. COMPONENT-BASED SCORING =====
  if ("rm_component" %in% names(data)) {
    for (i in seq_len(nrow(data))) {
      comp <- tolower(data$rm_component[i])
      if (is.na(comp)) next

      evidence <- strsplit(data$.type_evidence[i], ";")[[1]]

      # S subunit (specificity) → strong Type I indicator
      if (grepl("specificity|s.subunit|hsds|trd", comp)) {
        data$.type_I_score[i] <- data$.type_I_score[i] + 30L
        evidence <- c(evidence, "S_subunit_component")
      }

      data$.type_evidence[i] <- paste(unique(evidence[evidence != ""]), collapse = ";")
    }
  }

  # ===== 4. OPERON CONTEXT SCORING =====
  if (!is.null(operon_data) && nrow(operon_data) > 0 && "operon_id" %in% names(data)) {
    # Type I typically has 3 genes (R, M, S)
    # Type II typically has 2 genes (R, M) or fused
    # Type III typically has 2 genes (Mod, Res)

    for (i in seq_len(nrow(data))) {
      operon_id <- data$operon_id[i]
      if (is.na(operon_id)) next

      evidence <- strsplit(data$.type_evidence[i], ";")[[1]]

      # Count genes in same operon
      same_operon <- data[!is.na(data$operon_id) & data$operon_id == operon_id, ]
      n_genes <- nrow(same_operon)

      # Check for S subunit in operon → Type I
      products_in_operon <- tolower(same_operon[[product_col]])
      products_in_operon <- products_in_operon[!is.na(products_in_operon)]
      if (length(products_in_operon) > 0 &&
          any(grepl("specificity|hsds|s.subunit", products_in_operon))) {
        data$.type_I_score[i] <- data$.type_I_score[i] + 20L
        evidence <- c(evidence, "S_in_operon")
      }

      # 3 genes with R+M+S pattern → Type I
      if (n_genes == 3) {
        data$.type_I_score[i] <- data$.type_I_score[i] + 10L
        evidence <- c(evidence, "3gene_operon")
      }

      data$.type_evidence[i] <- paste(unique(evidence[evidence != ""]), collapse = ";")
    }
  }

  # ===== 5. DETERMINE FINAL TYPE =====
  data <- data %>%
    dplyr::mutate(
      # Get max score and determine type
      .max_score = pmax(.type_I_score, .type_II_score, .type_III_score, .type_IV_score),

      rm_system_type = dplyr::case_when(
        .max_score == 0 ~ "Unknown",
        .type_IV_score == .max_score & .type_IV_score >= 25 ~ "Type_IV",
        .type_I_score == .max_score & .type_I_score >= 20 ~ "Type_I",
        .type_III_score == .max_score & .type_III_score >= 20 ~ "Type_III",
        .type_II_score == .max_score & .type_II_score >= 15 ~ "Type_II",
        .type_I_score >= 15 ~ "Type_I",  # Default for helicase+ with lower threshold
        .type_II_score >= 10 ~ "Type_II", # Default for most R-M
        TRUE ~ "Unknown"
      ),

      rm_type_evidence = .type_evidence,

      # Also store individual scores for debugging
      rm_type_scores = paste0("I:", .type_I_score, ";II:", .type_II_score,
                               ";III:", .type_III_score, ";IV:", .type_IV_score)
    ) %>%
    dplyr::select(-starts_with(".type_"), -.max_score)

  # Summary message
  type_counts <- table(data$rm_system_type)
  message(sprintf("R-M System Type Classification: %s",
                  paste(names(type_counts), type_counts, sep = "=", collapse = ", ")))

  return(data)
}


#' Calculate R-M System Confidence Score
#'
#' @description Integrates all evidence sources to calculate confidence scores
#'
#' @param mtase_data MTase detection results
#' @param rease_data REase detection results
#' @param operon_data Operon analysis results
#' @param trd_data TRD extraction results
#' @param rebase_results REBASE comparison results
#'
#' @return Tibble with integrated scores for each R-M candidate
#' @export
calculate_rm_scores <- function(mtase_data = NULL,
                                 rease_data = NULL,
                                 operon_data = NULL,
                                 trd_data = NULL,
                                 rebase_results = NULL) {
  
  # Define scoring weights
  weights <- list(
    mtase_pfam = 20,
    mtase_keyword = 10,
    rease_pfam = 20,
    rease_keyword = 10,
    rease_motif = 15,
    operon_complete = 25,
    operon_partial = 10,
    trd_pfam = 15,
    trd_motif = 10,
    rebase_high = 30,
    rebase_medium = 20,
    rebase_low = 10
  )
  
  # Get all unique gene IDs
  all_ids <- unique(c(
    if (!is.null(mtase_data) && nrow(mtase_data) > 0) get_id_col_values(mtase_data) else character(),
    if (!is.null(rease_data) && nrow(rease_data) > 0) get_id_col_values(rease_data) else character(),
    if (!is.null(operon_data) && nrow(operon_data) > 0) get_id_col_values(operon_data) else character(),
    if (!is.null(trd_data) && nrow(trd_data) > 0) get_id_col_values(trd_data) else character()
  ))
  
  if (length(all_ids) == 0) {
    message("No R-M candidates to score")
    return(dplyr::tibble())
  }
  
  # Initialize score table
  scores <- dplyr::tibble(
    gene_id = all_ids,
    mtase_score = 0,
    rease_score = 0,
    operon_score = 0,
    trd_score = 0,
    rebase_score = 0
  )
  
  # MTase scores
  if (!is.null(mtase_data) && nrow(mtase_data) > 0) {
    id_col <- get_id_col(mtase_data)
    mtase_ids <- mtase_data[[id_col]]
    
    scores <- scores %>%
      dplyr::mutate(
        mtase_score = dplyr::case_when(
          gene_id %in% mtase_ids[mtase_data$detection_method == "pfam+keyword"] ~ 
            weights$mtase_pfam + weights$mtase_keyword,
          gene_id %in% mtase_ids[mtase_data$detection_method == "pfam_only"] ~ 
            weights$mtase_pfam,
          gene_id %in% mtase_ids[mtase_data$detection_method == "keyword_only"] ~ 
            weights$mtase_keyword,
          TRUE ~ 0
        )
      )
  }
  
  # REase scores
  if (!is.null(rease_data) && nrow(rease_data) > 0) {
    id_col <- get_id_col(rease_data)
    rease_ids <- rease_data[[id_col]]
    
    # Check for motif evidence
    has_motif <- if ("catalytic_motifs" %in% names(rease_data)) {
      !is.na(rease_data$catalytic_motifs)
    } else {
      rep(FALSE, nrow(rease_data))
    }
    
    scores <- scores %>%
      dplyr::mutate(
        rease_score = dplyr::case_when(
          gene_id %in% rease_ids[has_motif & !is.na(rease_data$re_type_pfam)] ~ 
            weights$rease_pfam + weights$rease_motif,
          gene_id %in% rease_ids[!is.na(rease_data$re_type_pfam)] ~ 
            weights$rease_pfam,
          gene_id %in% rease_ids[has_motif] ~ 
            weights$rease_motif,
          gene_id %in% rease_ids ~ 
            weights$rease_keyword,
          TRUE ~ 0
        )
      )
  }
  
  # Operon scores
  if (!is.null(operon_data) && nrow(operon_data) > 0) {
    id_col <- get_id_col(operon_data)
    
    complete_rm_ids <- operon_data %>%
      dplyr::filter(operon_type == "complete_RM") %>%
      dplyr::pull(!!rlang::sym(id_col))
    
    partial_rm_ids <- operon_data %>%
      dplyr::filter(operon_type %in% c("solitary_MTase", "solitary_REase")) %>%
      dplyr::pull(!!rlang::sym(id_col))
    
    scores <- scores %>%
      dplyr::mutate(
        operon_score = dplyr::case_when(
          gene_id %in% complete_rm_ids ~ weights$operon_complete,
          gene_id %in% partial_rm_ids ~ weights$operon_partial,
          TRUE ~ 0
        )
      )
  }
  
  # TRD scores
  if (!is.null(trd_data) && nrow(trd_data) > 0) {
    id_col <- get_id_col(trd_data)
    trd_ids <- trd_data[[id_col]]
    
    # Check detection method
    has_pfam <- if ("trd_pfam_str" %in% names(trd_data)) {
      !is.na(trd_data$trd_pfam_str)
    } else {
      rep(FALSE, nrow(trd_data))
    }
    
    has_motif <- if ("trd_motifs_found" %in% names(trd_data)) {
      !is.na(trd_data$trd_motifs_found) & trd_data$trd_motifs_found != ""
    } else {
      rep(FALSE, nrow(trd_data))
    }
    
    scores <- scores %>%
      dplyr::mutate(
        trd_score = dplyr::case_when(
          gene_id %in% trd_ids[has_pfam & has_motif] ~ 
            weights$trd_pfam + weights$trd_motif,
          gene_id %in% trd_ids[has_pfam] ~ 
            weights$trd_pfam,
          gene_id %in% trd_ids[has_motif] ~ 
            weights$trd_motif,
          gene_id %in% trd_ids ~ 
            5,  # minimal score for keyword-only
          TRUE ~ 0
        )
      )
  }
  
  # REBASE comparison scores
  if (!is.null(rebase_results) && nrow(rebase_results) > 0) {
    # Handle different column names from REBASE results
    # Could be: query_id + identity (raw BLAST) or query_id + best_match_identity (processed)

    if ("best_match_identity" %in% names(rebase_results)) {
      # Processed REBASE results (from get_best_rebase_match)
      best_matches <- rebase_results %>%
        dplyr::select(query_id, rebase_identity = best_match_identity) %>%
        dplyr::distinct()
    } else if ("pct_identity" %in% names(rebase_results)) {
      # Raw BLAST results with pct_identity column - get best per query
      best_matches <- rebase_results %>%
        dplyr::group_by(query_id) %>%
        dplyr::slice_max(order_by = pct_identity, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(query_id, rebase_identity = pct_identity)
    } else {
      best_matches <- dplyr::tibble(query_id = character(), rebase_identity = numeric())
    }

    if (nrow(best_matches) > 0) {
      scores <- scores %>%
        dplyr::left_join(best_matches, by = c("gene_id" = "query_id")) %>%
        dplyr::mutate(
          rebase_score = dplyr::case_when(
            !is.na(rebase_identity) & rebase_identity >= 0.7 ~ weights$rebase_high,
            !is.na(rebase_identity) & rebase_identity >= 0.5 ~ weights$rebase_medium,
            !is.na(rebase_identity) & rebase_identity >= 0.3 ~ weights$rebase_low,
            TRUE ~ 0
          )
        ) %>%
        dplyr::select(-rebase_identity)
    }
  }
  
  # Calculate total score and confidence
  scores <- scores %>%
    dplyr::mutate(
      total_score = mtase_score + rease_score + operon_score + trd_score + rebase_score,
      max_possible = sum(unlist(weights)),
      normalized_score = total_score / max_possible * 100,
      confidence_level = dplyr::case_when(
        normalized_score >= 60 ~ "very_high",
        normalized_score >= 40 ~ "high",
        normalized_score >= 25 ~ "medium",
        normalized_score >= 10 ~ "low",
        TRUE ~ "very_low"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(total_score))
  
  return(scores)
}


#' Helper: Get ID Column Values
#' @keywords internal
get_id_col_values <- function(data) {
  id_col <- intersect(c("locus_tag", "protein_id", "gene_id"), names(data))[1]
  if (is.na(id_col)) return(character())
  data[[id_col]]
}


#' Helper: Get ID Column Name
#' @keywords internal
get_id_col <- function(data) {
  intersect(c("locus_tag", "protein_id", "gene_id"), names(data))[1]
}


#' Classify R-M System Type
#'
#' @description Assigns final R-M system classification based on all evidence
#'
#' @param scored_data Scored data from calculate_rm_scores()
#' @param mtase_data MTase data with type information
#' @param rease_data REase data with type information
#' @param operon_data Operon data with predicted type
#'
#' @return Tibble with final classifications
#' @export
classify_rm_systems <- function(scored_data,
                                 mtase_data = NULL,
                                 rease_data = NULL,
                                 operon_data = NULL) {

  # Handle empty or invalid scored_data
 if (is.null(scored_data) || nrow(scored_data) == 0) {
    message("No scored data to classify")
    return(dplyr::tibble(
      gene_id = character(),
      mtase_score = numeric(),
      rease_score = numeric(),
      operon_score = numeric(),
      trd_score = numeric(),
      rebase_score = numeric(),
      total_score = numeric(),
      confidence_level = character(),
      is_mtase = logical(),
      is_rease = logical(),
      is_trd = logical(),
      rm_classification = character(),
      rm_type = character()
    ))
  }

  # Check required columns exist
  required_cols <- c("mtase_score", "rease_score", "trd_score")
  missing_cols <- setdiff(required_cols, names(scored_data))
  if (length(missing_cols) > 0) {
    # Add missing columns with 0
    for (col in missing_cols) {
      scored_data[[col]] <- 0
    }
  }

  result <- scored_data

  # Initialize classification columns
  result <- result %>%
    dplyr::mutate(
      is_mtase = mtase_score > 0,
      is_rease = rease_score > 0,
      is_trd = trd_score > 0,
      rm_classification = "unknown",
      rm_type = "unknown"
    )
  
  # Add MTase type
  if (!is.null(mtase_data) && nrow(mtase_data) > 0) {
    id_col <- get_id_col(mtase_data)
    mtase_types <- mtase_data %>%
      dplyr::select(dplyr::all_of(id_col), mtase_type) %>%
      dplyr::rename(gene_id = !!id_col, mtase_type_detail = mtase_type)
    
    result <- result %>%
      dplyr::left_join(mtase_types, by = "gene_id")
  }
  
  # Add REase type
  if (!is.null(rease_data) && nrow(rease_data) > 0) {
    id_col <- get_id_col(rease_data)
    rease_types <- rease_data %>%
      dplyr::select(dplyr::all_of(id_col), re_type) %>%
      dplyr::rename(gene_id = !!id_col, rease_type_detail = re_type)
    
    result <- result %>%
      dplyr::left_join(rease_types, by = "gene_id")
  }
  
  # Add operon predicted type
  if (!is.null(operon_data) && nrow(operon_data) > 0 && "predicted_type" %in% names(operon_data)) {
    id_col <- get_id_col(operon_data)
    operon_types <- operon_data %>%
      dplyr::select(dplyr::all_of(id_col), predicted_type) %>%
      dplyr::distinct() %>%
      dplyr::rename(gene_id = !!id_col, operon_predicted_type = predicted_type)
    
    result <- result %>%
      dplyr::left_join(operon_types, by = "gene_id")
  }
  
  # Final classification
  result <- result %>%
    dplyr::mutate(
      rm_classification = dplyr::case_when(
        is_mtase & is_rease ~ "RM_complete",
        is_mtase & !is_rease ~ "MTase_only",
        !is_mtase & is_rease ~ "REase_only",
        is_trd ~ "TRD_component",
        TRUE ~ "unknown"
      ),
      rm_type = dplyr::coalesce(
        if ("operon_predicted_type" %in% names(.)) operon_predicted_type else NA_character_,
        if ("rease_type_detail" %in% names(.)) rease_type_detail else NA_character_,
        if ("mtase_type_detail" %in% names(.)) mtase_type_detail else NA_character_,
        "unknown"
      )
    )
  
  return(result)
}


#' Generate R-M System Summary Report
#'
#' @description Creates comprehensive summary of R-M system analysis
#'
#' @param classified_data Classified R-M data
#' @param dnmb_data Original DNMB data for context
#'
#' @return List with summary statistics and visualizations
#' @export
generate_rm_summary <- function(classified_data, dnmb_data = NULL) {

  # Handle empty data
  if (is.null(classified_data) || nrow(classified_data) == 0) {
    message("No data to summarize")
    return(list(
      classification_counts = dplyr::tibble(rm_classification = character(), n = integer()),
      type_counts = dplyr::tibble(rm_type = character(), n = integer()),
      confidence_distribution = dplyr::tibble(confidence_level = character(), n = integer()),
      score_statistics = dplyr::tibble(n_total = 0L, mean_score = NA_real_,
                                       median_score = NA_real_, max_score = NA_real_,
                                       n_high_confidence = 0L),
      top_candidates = dplyr::tibble()
    ))
  }

  # Count by classification (handle missing column)
  if ("rm_classification" %in% names(classified_data)) {
    classification_counts <- classified_data %>%
      dplyr::count(rm_classification, name = "n")
  } else {
    classification_counts <- dplyr::tibble(rm_classification = "unknown", n = nrow(classified_data))
  }

  # Count by type (handle missing column)
  if ("rm_type" %in% names(classified_data)) {
    type_counts <- classified_data %>%
      dplyr::filter(rm_type != "unknown") %>%
      dplyr::count(rm_type, name = "n")
  } else {
    type_counts <- dplyr::tibble(rm_type = character(), n = integer())
  }

  # Confidence distribution (handle missing column)
  if ("confidence_level" %in% names(classified_data)) {
    confidence_dist <- classified_data %>%
      dplyr::count(confidence_level, name = "n") %>%
      dplyr::mutate(
        confidence_level = factor(confidence_level,
                                   levels = c("very_high", "high", "medium", "low", "very_low"))
      ) %>%
      dplyr::arrange(confidence_level)
  } else {
    confidence_dist <- dplyr::tibble(confidence_level = character(), n = integer())
  }

  # Score statistics (handle missing columns)
  if ("total_score" %in% names(classified_data)) {
    score_stats <- classified_data %>%
      dplyr::summarise(
        n_total = dplyr::n(),
        mean_score = mean(total_score, na.rm = TRUE),
        median_score = median(total_score, na.rm = TRUE),
        max_score = max(total_score, na.rm = TRUE),
        n_high_confidence = if ("confidence_level" %in% names(.)) {
          sum(confidence_level %in% c("very_high", "high"), na.rm = TRUE)
        } else 0L
      )
  } else {
    score_stats <- dplyr::tibble(n_total = nrow(classified_data), mean_score = NA_real_,
                                 median_score = NA_real_, max_score = NA_real_,
                                 n_high_confidence = 0L)
  }

  # Top candidates
  if ("confidence_level" %in% names(classified_data) && "total_score" %in% names(classified_data)) {
    top_candidates <- classified_data %>%
      dplyr::filter(confidence_level %in% c("very_high", "high")) %>%
      dplyr::arrange(dplyr::desc(total_score)) %>%
      dplyr::slice_head(n = 20)
  } else {
    top_candidates <- classified_data %>%
      dplyr::slice_head(n = 20)
  }
  
  # Print summary
  message("\n", paste(rep("=", 50), collapse = ""))
  message("R-M SYSTEM ANALYSIS SUMMARY")
  message(paste(rep("=", 50), collapse = ""))
  
  message("\nClassification Summary:")
  for (i in seq_len(nrow(classification_counts))) {
    message(sprintf("  %s: %d", 
                    classification_counts$rm_classification[i],
                    classification_counts$n[i]))
  }
  
  message("\nType Distribution:")
  for (i in seq_len(nrow(type_counts))) {
    message(sprintf("  %s: %d", type_counts$rm_type[i], type_counts$n[i]))
  }
  
  message("\nConfidence Distribution:")
  for (i in seq_len(nrow(confidence_dist))) {
    message(sprintf("  %s: %d", confidence_dist$confidence_level[i], confidence_dist$n[i]))
  }
  
  message("\nScore Statistics:")
  message(sprintf("  Total candidates: %d", score_stats$n_total))
  message(sprintf("  High confidence: %d (%.1f%%)", 
                  score_stats$n_high_confidence,
                  100 * score_stats$n_high_confidence / score_stats$n_total))
  message(sprintf("  Mean score: %.1f", score_stats$mean_score))
  message(sprintf("  Max score: %.1f", score_stats$max_score))
  
  return(list(
    classification_counts = classification_counts,
    type_counts = type_counts,
    confidence_distribution = confidence_dist,
    score_statistics = score_stats,
    top_candidates = top_candidates
  ))
}


#' Export R-M Analysis Results
#'
#' @description Exports all R-M analysis results to files
#'
#' @param classified_data Final classified data
#' @param output_dir Output directory
#' @param prefix File name prefix
#' @param formats Export formats: "csv", "xlsx", "rds"
#'
#' @return Paths to exported files
#' @export
export_rm_results <- function(classified_data,
                               output_dir = ".",
                               prefix = "rm_analysis",
                               formats = c("csv", "rds")) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  exported_files <- character()
  
  # CSV export
  if ("csv" %in% formats) {
    csv_file <- file.path(output_dir, paste0(prefix, "_results.csv"))
    utils::write.csv(classified_data, csv_file, row.names = FALSE)
    exported_files <- c(exported_files, csv_file)
    message("Exported: ", csv_file)
  }
  
  # RDS export
  if ("rds" %in% formats) {
    rds_file <- file.path(output_dir, paste0(prefix, "_results.rds"))
    saveRDS(classified_data, rds_file)
    exported_files <- c(exported_files, rds_file)
    message("Exported: ", rds_file)
  }
  
  # Excel export
  if ("xlsx" %in% formats) {
    if (requireNamespace("writexl", quietly = TRUE)) {
      xlsx_file <- file.path(output_dir, paste0(prefix, "_results.xlsx"))
      writexl::write_xlsx(classified_data, xlsx_file)
      exported_files <- c(exported_files, xlsx_file)
      message("Exported: ", xlsx_file)
    } else {
      warning("Package 'writexl' needed for Excel export. Skipping.")
    }
  }
  
  return(exported_files)
}


#' Visualize R-M Scores
#'
#' @description Creates visualization of R-M system scores
#'
#' @param scored_data Scored data
#'
#' @return ggplot object
#' @export
plot_rm_scores <- function(scored_data) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting")
  }
  
  # Prepare data for stacked bar chart
  plot_data <- scored_data %>%
    dplyr::select(gene_id, mtase_score, rease_score, operon_score, trd_score, rebase_score) %>%
    tidyr::pivot_longer(
      cols = -gene_id,
      names_to = "score_type",
      values_to = "score"
    ) %>%
    dplyr::mutate(
      score_type = factor(score_type, 
                          levels = c("rebase_score", "trd_score", "operon_score", 
                                     "rease_score", "mtase_score"),
                          labels = c("REBASE", "TRD", "Operon", "REase", "MTase"))
    )
  
  # Order genes by total score
  gene_order <- scored_data %>%
    dplyr::arrange(dplyr::desc(total_score)) %>%
    dplyr::pull(gene_id)
  
  plot_data <- plot_data %>%
    dplyr::mutate(gene_id = factor(gene_id, levels = gene_order))
  
  # Limit to top N for readability
  if (length(gene_order) > 30) {
    plot_data <- plot_data %>%
      dplyr::filter(gene_id %in% gene_order[1:30])
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = gene_id, y = score, fill = score_type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Evidence") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      x = "Gene ID",
      y = "Score",
      title = "R-M System Prediction Scores",
      subtitle = "Breakdown by evidence type"
    )
  
  return(p)
}

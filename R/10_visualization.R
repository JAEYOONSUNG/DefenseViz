#' @title Visualization Module
#' @description Functions for generating heatmaps and visualizations of R-M system analysis
#' @name visualization
NULL

#' Generate R-M System Heatmap
#'
#' @description Creates a heatmap visualization of R-M system candidates showing
#' detection evidence, type classification, and scores
#'
#' @param classified_data Classified/scored data from pipeline
#' @param mtase_data MTase detection results
#' @param rease_data REase detection results
#' @param output_file Output file path (PNG)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param max_genes Maximum genes to show
#'
#' @return ggplot object (also saves to file if output_file provided)
#' @export
generate_rm_heatmap <- function(classified_data,
                                 mtase_data = NULL,
                                 rease_data = NULL,
                                 output_file = NULL,
                                 width = 12,
                                 height = 10,
                                 max_genes = 50) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  if (nrow(classified_data) == 0) {
    message("No data to plot")
    return(NULL)
  }

  # Prepare heatmap matrix data
  heatmap_data <- prepare_heatmap_matrix(
    classified_data = classified_data,
    mtase_data = mtase_data,
    rease_data = rease_data,
    max_genes = max_genes
  )

  if (is.null(heatmap_data) || nrow(heatmap_data) == 0) {
    message("No data to plot after preparation")
    return(NULL)
  }

  # Create the heatmap
  p <- create_evidence_heatmap(heatmap_data)

  # Save to file if specified
  if (!is.null(output_file)) {
    ggplot2::ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )
    message("Heatmap saved to: ", output_file)
  }

  return(p)
}


#' Prepare Heatmap Matrix Data
#'
#' @description Prepares data in long format for heatmap visualization
#'
#' @keywords internal
prepare_heatmap_matrix <- function(classified_data, mtase_data, rease_data, max_genes) {

  # Get top genes by score
  top_genes <- classified_data %>%
    dplyr::arrange(dplyr::desc(total_score)) %>%
    dplyr::slice_head(n = max_genes)

  if (nrow(top_genes) == 0) return(NULL)

  # Build evidence matrix
  evidence_cols <- c(
    "mtase_score", "rease_score", "operon_score", "trd_score", "rebase_score"
  )

  # Normalize scores for visualization
  matrix_data <- top_genes %>%
    dplyr::select(gene_id, dplyr::any_of(evidence_cols)) %>%
    tidyr::pivot_longer(
      cols = -gene_id,
      names_to = "evidence_type",
      values_to = "score"
    ) %>%
    dplyr::mutate(
      evidence_type = factor(
        evidence_type,
        levels = c("mtase_score", "rease_score", "operon_score", "trd_score", "rebase_score"),
        labels = c("MTase", "REase", "Operon", "TRD", "REBASE")
      ),
      # Normalize to 0-1 for color scale
      normalized = dplyr::case_when(
        score == 0 ~ 0,
        TRUE ~ pmin(score / 30, 1)  # Cap at 30 for visualization
      )
    )

  # Add classification info
  class_info <- top_genes %>%
    dplyr::select(gene_id, rm_classification, rm_type, confidence_level, total_score)

  matrix_data <- matrix_data %>%
    dplyr::left_join(class_info, by = "gene_id")

  # Order genes by total score
  gene_order <- top_genes %>%
    dplyr::arrange(dplyr::desc(total_score)) %>%
    dplyr::pull(gene_id)

  matrix_data$gene_id <- factor(matrix_data$gene_id, levels = rev(gene_order))

  return(matrix_data)
}


#' Create Evidence Heatmap
#'
#' @description Creates the actual heatmap plot
#'
#' @keywords internal
create_evidence_heatmap <- function(heatmap_data) {

  # Main heatmap
  p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = evidence_type, y = gene_id)) +
    ggplot2::geom_tile(ggplot2::aes(fill = normalized), color = "white", size = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "white",
      mid = "#74ADD1",
      high = "#313695",
      midpoint = 0.5,
      limits = c(0, 1),
      name = "Evidence\nStrength"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 7),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "Evidence Type",
      y = "Gene ID",
      title = "R-M System Detection Evidence Heatmap",
      subtitle = "Evidence strength across detection methods"
    )

  # Add confidence level as side annotation if available
  if ("confidence_level" %in% names(heatmap_data)) {
    # Create side bar data
    conf_data <- heatmap_data %>%
      dplyr::select(gene_id, confidence_level) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        x = "Confidence",
        conf_score = dplyr::case_when(
          confidence_level == "very_high" ~ 1.0,
          confidence_level == "high" ~ 0.75,
          confidence_level == "medium" ~ 0.5,
          confidence_level == "low" ~ 0.25,
          TRUE ~ 0.1
        )
      )

    p <- p +
      ggplot2::geom_tile(
        data = conf_data,
        ggplot2::aes(x = x, y = gene_id, fill = conf_score),
        color = "white", size = 0.3
      )
  }

  return(p)
}


#' Generate R-M Type Distribution Plot
#'
#' @description Creates a bar plot showing distribution of R-M system types
#'
#' @param classified_data Classified data
#' @param output_file Output file path (optional)
#'
#' @return ggplot object
#' @export
plot_rm_type_distribution <- function(classified_data, output_file = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  # Count by type
  type_counts <- classified_data %>%
    dplyr::filter(rm_type != "unknown") %>%
    dplyr::count(rm_type, name = "count") %>%
    dplyr::arrange(dplyr::desc(count))

  if (nrow(type_counts) == 0) {
    message("No typed R-M systems to plot")
    return(NULL)
  }

  p <- ggplot2::ggplot(type_counts, ggplot2::aes(x = reorder(rm_type, count), y = count)) +
    ggplot2::geom_bar(stat = "identity", fill = "#3288BD", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = count), hjust = -0.2, size = 4) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = "R-M System Type",
      y = "Count",
      title = "Distribution of R-M System Types"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 8, height = 6, dpi = 300)
    message("Type distribution plot saved to: ", output_file)
  }

  return(p)
}


#' Generate Operon Structure Visualization
#'
#' @description Creates a visualization of R-M operon structures
#'
#' @param operon_data Operon data from identify_rm_operons()
#' @param output_file Output file path (optional)
#' @param max_operons Maximum operons to show
#'
#' @return ggplot object
#' @export
plot_operon_structure <- function(operon_data, output_file = NULL, max_operons = 20) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  if (nrow(operon_data) == 0) {
    message("No operon data to plot")
    return(NULL)
  }

  # Limit to top operons
  top_operons <- operon_data %>%
    dplyr::filter(operon_id %in% unique(operon_id)[1:min(max_operons, dplyr::n_distinct(operon_id))])

  # Prepare gene arrow data
  plot_data <- top_operons %>%
    dplyr::mutate(
      gene_color = dplyr::case_when(
        is_mtase ~ "#E41A1C",   # Red for MTase
        is_rease ~ "#377EB8",   # Blue for REase
        TRUE ~ "#999999"        # Grey for other
      ),
      operon_label = paste("Operon", operon_id, "-", system_type)
    )

  # Create gene arrow plot
  p <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = start, xend = end,
        y = factor(operon_id), yend = factor(operon_id),
        color = gene_color
      ),
      size = 8,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"))
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::facet_wrap(~ operon_label, scales = "free_x", ncol = 1) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Genome Position (bp)",
      y = "",
      title = "R-M Operon Structures",
      subtitle = "Red = MTase, Blue = REase"
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 12, height = max_operons * 0.5 + 2, dpi = 300)
    message("Operon structure plot saved to: ", output_file)
  }

  return(p)
}


#' Generate Score Breakdown Plot
#'
#' @description Creates a stacked bar chart showing score breakdown for each gene
#'
#' @param scored_data Scored data from calculate_rm_scores()
#' @param output_file Output file path (optional)
#' @param top_n Number of top genes to show
#'
#' @return ggplot object
#' @export
plot_score_breakdown <- function(scored_data, output_file = NULL, top_n = 30) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }

  if (nrow(scored_data) == 0) {
    message("No scored data to plot")
    return(NULL)
  }

  # Get top genes
  top_genes <- scored_data %>%
    dplyr::arrange(dplyr::desc(total_score)) %>%
    dplyr::slice_head(n = top_n)

  # Prepare long format for stacking
  plot_data <- top_genes %>%
    dplyr::select(gene_id, mtase_score, rease_score, operon_score, trd_score, rebase_score) %>%
    tidyr::pivot_longer(
      cols = -gene_id,
      names_to = "score_type",
      values_to = "score"
    ) %>%
    dplyr::mutate(
      score_type = factor(
        score_type,
        levels = c("rebase_score", "trd_score", "operon_score", "rease_score", "mtase_score"),
        labels = c("REBASE", "TRD", "Operon", "REase", "MTase")
      )
    )

  # Order by total score
  gene_order <- top_genes$gene_id
  plot_data$gene_id <- factor(plot_data$gene_id, levels = gene_order)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = gene_id, y = score, fill = score_type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_brewer(palette = "Set2", name = "Evidence") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      x = "Gene ID",
      y = "Score",
      title = "R-M System Prediction Scores",
      subtitle = "Breakdown by evidence type"
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 12, height = 8, dpi = 300)
    message("Score breakdown plot saved to: ", output_file)
  }

  return(p)
}


#' Generate Recognition Sequence Logo
#'
#' @description Creates a sequence logo for predicted recognition sequences
#' (Requires ggseqlogo package)
#'
#' @param recognition_data Recognition sequence predictions
#' @param rm_type Filter by R-M type (optional)
#' @param output_file Output file path (optional)
#'
#' @return ggplot object or NULL if ggseqlogo not available
#' @export
plot_recognition_logo <- function(recognition_data, rm_type = NULL, output_file = NULL) {

  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    message("Package 'ggseqlogo' is required for sequence logo. Install with: install.packages('ggseqlogo')")
    return(NULL)
  }

  if (nrow(recognition_data) == 0 || !"predicted_recognition_seq" %in% names(recognition_data)) {
    message("No recognition sequences to plot")
    return(NULL)
  }

  # Filter by type if specified
  seqs <- recognition_data
  if (!is.null(rm_type)) {
    seqs <- seqs %>% dplyr::filter(rebase_type == rm_type)
  }

  # Get unique sequences
  unique_seqs <- unique(seqs$predicted_recognition_seq)
  unique_seqs <- unique_seqs[!is.na(unique_seqs) & unique_seqs != ""]

  if (length(unique_seqs) == 0) {
    message("No valid recognition sequences")
    return(NULL)
  }

  # Create logo
  p <- ggseqlogo::ggseqlogo(unique_seqs) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Recognition Sequence Logo", if (!is.null(rm_type)) paste("-", rm_type) else ""),
      x = "Position",
      y = "Bits"
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 8, height = 4, dpi = 300)
    message("Recognition logo saved to: ", output_file)
  }

  return(p)
}


#' Generate Complete Visualization Report
#'
#' @description Generates all visualization plots and saves them
#'
#' @param results Results list from rmscan_pipeline()
#' @param output_dir Output directory
#'
#' @return List of generated plot objects
#' @export
generate_visualization_report <- function(results, output_dir = NULL) {

  if (is.null(output_dir)) {
    output_dir <- results$output_dir
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  plots <- list()

  # 1. Main heatmap
  message("Generating evidence heatmap...")
  plots$heatmap <- tryCatch({
    generate_rm_heatmap(
      classified_data = results$classified_data,
      mtase_data = results$mtase_data,
      rease_data = results$rease_data,
      output_file = file.path(output_dir, "rm_evidence_heatmap.png")
    )
  }, error = function(e) {
    warning("Heatmap generation failed: ", e$message)
    NULL
  })

  # 2. Type distribution
  message("Generating type distribution plot...")
  plots$type_dist <- tryCatch({
    plot_rm_type_distribution(
      classified_data = results$classified_data,
      output_file = file.path(output_dir, "rm_type_distribution.png")
    )
  }, error = function(e) {
    warning("Type distribution plot failed: ", e$message)
    NULL
  })

  # 3. Score breakdown
  message("Generating score breakdown plot...")
  plots$score_breakdown <- tryCatch({
    plot_score_breakdown(
      scored_data = results$scored_data,
      output_file = file.path(output_dir, "rm_score_breakdown.png")
    )
  }, error = function(e) {
    warning("Score breakdown plot failed: ", e$message)
    NULL
  })

  # 4. Operon structure
  if (!is.null(results$operon_data) && nrow(results$operon_data) > 0) {
    message("Generating operon structure plot...")
    plots$operon_structure <- tryCatch({
      plot_operon_structure(
        operon_data = results$operon_data,
        output_file = file.path(output_dir, "rm_operon_structure.png")
      )
    }, error = function(e) {
      warning("Operon structure plot failed: ", e$message)
      NULL
    })
  }

  message("Visualization report complete. Files saved to: ", output_dir)
  return(plots)
}

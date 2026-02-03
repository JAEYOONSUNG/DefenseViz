#' @title Data Loading and Preprocessing Module
#' @description Functions for loading DNMB output and extracting R-M system candidates
#' @name data_loading
NULL

#' Load DNMB Output File
#'
#' @description Reads DNMB annotation output (Excel or CSV format) and standardizes column names
#'
#' @param file_path Path to DNMB output file (.xlsx, .xls, or .csv)
#' @param sheet Sheet name or number for Excel files (default: 1)
#'
#' @return A tibble with standardized DNMB annotation data
#' @export
#'
#' @examples
#' \dontrun{
#' dnmb_data <- load_dnmb("path/to/annotation.xlsx")
#' }
load_dnmb <- function(file_path, sheet = 1) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Determine file type and read
  ext <- tolower(tools::file_ext(file_path))
  
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required. Install with: install.packages('readxl')")
    }
    data <- readxl::read_excel(file_path, sheet = sheet)
  } else if (ext == "csv") {
    data <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  } else if (ext == "tsv" || ext == "txt") {
    data <- utils::read.delim(file_path, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file format: ", ext, ". Use .xlsx, .xls, .csv, or .tsv")
  }
  
  # Convert to tibble
  data <- dplyr::as_tibble(data)
  
  # Validate required columns
  required_cols <- c("locus_tag", "start", "end", "direction", "translation", "product")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    warning("Missing recommended columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Add metadata
  attr(data, "source_file") <- basename(file_path)
  attr(data, "load_time") <- Sys.time()
  attr(data, "n_genes") <- nrow(data)
  
  message(sprintf("Loaded %d genes from %s", nrow(data), basename(file_path)))
  
  return(data)
}


#' Extract R-M System Candidates by Keywords
#'
#' @description Identifies potential R-M system genes based on product annotation keywords
#'
#' @param dnmb_data DNMB data frame from load_dnmb()
#' @param keywords Character vector of keywords to search in product column
#' @param product_col Column name containing product description (default: "product")
#' @param case_sensitive Logical, whether to use case-sensitive matching (default: FALSE)
#'
#' @return Filtered tibble containing potential R-M system genes
#' @export
#'
#' @examples
#' \dontrun{
#' rm_candidates <- extract_rm_candidates(dnmb_data)
#' }
extract_rm_candidates <- function(dnmb_data, 
                                   keywords = NULL,
                                   product_col = "product",
                                   case_sensitive = FALSE) {
  
  # Default R-M system keywords
  if (is.null(keywords)) {
    keywords <- c(
      # Methyltransferases
      "methyltransferase", "methylase", "MTase",
      "N-6 adenine", "N6-adenine", "N4-cytosine", "C-5 cytosine",
      "dam", "dcm", "hsdM",
      # Restriction enzymes
      "restriction", "endonuclease", "REase",
      "hsdR", "hsdS",
      # Type-specific
      "type I", "type II", "type III", "type IV",
      "specificity", "modification"
    )
  }
  
  if (!product_col %in% names(dnmb_data)) {
    stop("Product column '", product_col, "' not found in data")
  }
  
  # Build regex pattern
  if (case_sensitive) {
    pattern <- paste(keywords, collapse = "|")
  } else {
    pattern <- paste0("(?i)", paste(keywords, collapse = "|"))
  }
  
  # Filter by keywords
  candidates <- dnmb_data %>%
    dplyr::filter(stringr::str_detect(.data[[product_col]], pattern))
  
  message(sprintf("Found %d R-M candidate genes (%.1f%% of total)", 
                  nrow(candidates), 
                  100 * nrow(candidates) / nrow(dnmb_data)))
  
  return(candidates)
}


#' Detect Available Annotation Columns
#'
#' @description Checks which annotation sources are available in DNMB data
#'
#' @param dnmb_data DNMB data frame
#'
#' @return A named list indicating available annotation types
#' @export
detect_annotation_sources <- function(dnmb_data) {
  
  cols <- names(dnmb_data)
  
  sources <- list(
    # PFAM annotations
    pfam = any(stringr::str_detect(cols, "(?i)pfam")),
    pfam_cols = cols[stringr::str_detect(cols, "(?i)pfam")],
    
    # InterProScan signatures
    interpro = any(stringr::str_detect(cols, "(?i)signature")),
    
    # CDD (Conserved Domain Database)
    cdd = any(stringr::str_detect(cols, "(?i)cdd")),
    cdd_cols = cols[stringr::str_detect(cols, "(?i)cdd")],
    
    # TIGRFAM
    tigrfam = any(stringr::str_detect(cols, "(?i)tigrfam")),
    
    # PANTHER
    panther = any(stringr::str_detect(cols, "(?i)panther")),
    
    # COG/eggNOG
    cog = any(stringr::str_detect(cols, "(?i)cog|eggnog")),
    
    # Gene3D
    gene3d = any(stringr::str_detect(cols, "(?i)gene3d")),
    
    # Hamap
    hamap = any(stringr::str_detect(cols, "(?i)hamap")),
    
    # Sequence data
    has_protein_seq = "translation" %in% cols,
    has_nt_seq = any(c("nt_seq", "rearranged_nt_seq") %in% cols),
    
    # Position data
    has_coordinates = all(c("start", "end") %in% cols),
    has_direction = "direction" %in% cols
  )
  
  # Summary message
  available <- names(sources)[sapply(sources, function(x) {
    if (is.logical(x)) x else FALSE
  })]
  
  message("Available annotation sources: ", 
          paste(available[!grepl("_cols|has_", available)], collapse = ", "))
  
  return(sources)
}


#' Get DNMB Column Mapping
#'
#' @description Returns standard column name mappings for DNMB data
#'
#' @return Named list of column mappings
#' @export
get_column_mapping <- function() {
  list(
    # Core identifiers
    id = c("locus_tag", "protein_id", "gene"),
    
    # Position
    position = c("start", "end", "direction", "contig"),
    
    # Sequences
    protein_seq = c("translation"),
    nucleotide_seq = c("nt_seq", "rearranged_nt_seq"),
    
    # PFAM (from eggNOG-mapper)
    pfam_eggnog = c("PFAMs"),
    
    # PFAM (from InterProScan)
    pfam_interpro = c("Signature accession_Pfam", "Signature description_Pfam"),
    
    # CDD
    cdd = c("Signature accession_CDD", "Signature description_CDD"),
    
    # Product annotation
    product = c("product", "Description", "Preferred_name")
  )
}

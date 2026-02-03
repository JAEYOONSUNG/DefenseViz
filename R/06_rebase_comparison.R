#' @title REBASE Comparison Module
#' @description Functions for comparing sequences with REBASE database
#' @name rebase_comparison
NULL

#' REBASE FTP URLs
#'
#' @description Returns URLs for REBASE data files
#'
#' @return Named list of REBASE URLs
#' @export
get_rebase_urls <- function() {
  list(
    # FTP files (actual data)
    protein_gold = "ftp://ftp.neb.com/pub/rebase/protein_seqs.txt",
    allenz = "ftp://ftp.neb.com/pub/rebase/allenz.txt",
    bairoch = "ftp://ftp.neb.com/pub/rebase/bairoch.txt",

    # Type-specific (Gold standards from REBASE website)
    type_i = "http://rebase.neb.com/rebase/rebase.typeI.html",
    type_ii = "http://rebase.neb.com/rebase/rebase.typeII.html",
    type_iii = "http://rebase.neb.com/rebase/rebase.typeIII.html"
  )
}


#' Get REBASE Cache Directory
#'
#' @description Returns the directory where REBASE data is cached
#'
#' @return Path to cache directory
#' @export
get_rebase_cache_dir <- function() {
  cache_dir <- tools::R_user_dir("DefenseViz", which = "cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  return(cache_dir)
}


#' Get REBASE Data as R Table
#'
#' @description Downloads (if needed) and returns REBASE gold standard data as a tibble.
#' Contains enzyme names, recognition sequences, protein sequences, and type classifications.
#'
#' REBASE updates on the 1st of each month. This function:
#' - Downloads fresh data if no local cache exists
#' - Re-downloads if local cache is from a previous month
#' - Uses cached data if it's from the current month
#'
#' @param force_download Force re-download even if cache exists
#' @param verbose Print status messages
#'
#' @return Tibble with columns: enzyme_name, enz_type, rec_seq (recognition sequence),
#'         org_name, sequence (protein), seq_length, rm_type, subunit
#' @export
#'
#' @examples
#' \dontrun{
#' # Get REBASE data table
#' rebase <- get_rebase_data()
#'
#' # Filter for Type I enzymes with recognition sequences
#' type_i <- rebase %>%
#'   filter(rm_type == "Type_I", rec_seq != "")
#'
#' # Get unique recognition sequences
#' rec_seqs <- rebase %>%
#'   filter(rec_seq != "") %>%
#'   distinct(enzyme_name, rec_seq, rm_type)
#' }
get_rebase_data <- function(force_download = FALSE, verbose = TRUE) {

  cache_dir <- get_rebase_cache_dir()
  rds_file <- file.path(cache_dir, "rebase_data.rds")
  txt_file <- file.path(cache_dir, "REBASE_protein_seqs.txt")

  if (verbose) {
    message("=== REBASE Data Location ===")
    message("Cache directory: ", cache_dir)
    message("RDS cache file: ", rds_file)
    message("Raw text file:  ", txt_file)
  }

  # Check if we need to update
  need_download <- force_download

  if (!need_download && file.exists(rds_file)) {
    file_info <- file.info(rds_file)
    file_date <- as.Date(file_info$mtime)
    current_date <- Sys.Date()
    current_month_start <- as.Date(format(current_date, "%Y-%m-01"))

    if (file_date >= current_month_start) {
      if (verbose) {
        message("Using cached REBASE data (current month: ", format(current_date, "%Y-%m"), ")")
        message("Last updated: ", file_date)
      }
      rebase_data <- readRDS(rds_file)
      if (verbose) message("Loaded ", nrow(rebase_data), " REBASE sequences")
      return(rebase_data)
    } else {
      if (verbose) {
        message("Cache is from previous month (", file_date, "). Checking for update...")
      }
      need_download <- TRUE
    }
  } else if (!file.exists(rds_file)) {
    need_download <- TRUE
  }

  # Download and parse
  if (need_download) {
    if (verbose) message("\nDownloading REBASE data...")

    txt_file <- download_rebase_gold(cache_dir, force_download = TRUE)

    if (is.null(txt_file)) {
      # Try to use existing RDS if download failed
      if (file.exists(rds_file)) {
        warning("Download failed, using existing (older) cached data")
        return(readRDS(rds_file))
      }
      stop("Failed to download REBASE data and no cache available")
    }

    if (verbose) message("Parsing REBASE data...")
    rebase_data <- parse_rebase_sequences(txt_file)

    # Save RDS cache
    if (nrow(rebase_data) > 0) {
      if (verbose) message("Saving parsed data to cache: ", rds_file)
      saveRDS(rebase_data, rds_file)
      # Update timestamp
      Sys.setFileTime(rds_file, Sys.time())
    }
  } else {
    # Load from RDS
    rebase_data <- readRDS(rds_file)
  }

  if (verbose) message("Total REBASE sequences: ", nrow(rebase_data))
  return(rebase_data)
}


#' Get REBASE Recognition Sequences Table
#'
#' @description Returns a simplified table of REBASE enzymes with their recognition sequences.
#' Filters to only include entries with known recognition sequences.
#'
#' @param force_download Force re-download
#'
#' @return Tibble with enzyme_name, rec_seq, rm_type, subunit, org_name
#' @export
#'
#' @examples
#' \dontrun{
#' rec_seqs <- get_rebase_recognition_sequences()
#' print(rec_seqs)
#' }
get_rebase_recognition_sequences <- function(force_download = FALSE) {

  rebase_data <- get_rebase_data(force_download = force_download, verbose = FALSE)

  result <- rebase_data %>%
    dplyr::filter(!is.na(rec_seq) & rec_seq != "" & rec_seq != "?") %>%
    dplyr::select(enzyme_name, rec_seq, rm_type, subunit, enz_type, org_name) %>%
    dplyr::distinct() %>%
    dplyr::arrange(rm_type, enzyme_name)

  message(sprintf("Found %d enzymes with recognition sequences", nrow(result)))
  return(result)
}


#' Show REBASE Cache Status
#'
#' @description Displays information about the local REBASE cache
#'
#' @return List with cache status information
#' @export
rebase_cache_status <- function() {

  cache_dir <- get_rebase_cache_dir()
  rds_file <- file.path(cache_dir, "rebase_data.rds")
  txt_file <- file.path(cache_dir, "REBASE_protein_seqs.txt")
  fasta_file <- file.path(cache_dir, "rebase_db.fasta")

  status <- list(
    cache_directory = cache_dir,
    rds_exists = file.exists(rds_file),
    txt_exists = file.exists(txt_file),
    fasta_exists = file.exists(fasta_file)
  )

  message("=== REBASE Cache Status ===")
  message("Cache directory: ", cache_dir)
  message("")

  if (status$rds_exists) {
    info <- file.info(rds_file)
    status$rds_size_mb <- round(info$size / 1e6, 1)
    status$rds_date <- as.Date(info$mtime)
    message("RDS cache:  EXISTS")
    message("  Size: ", status$rds_size_mb, " MB")
    message("  Date: ", status$rds_date)

    # Check if current month
    current_month_start <- as.Date(format(Sys.Date(), "%Y-%m-01"))
    status$rds_is_current <- status$rds_date >= current_month_start
    message("  Current: ", ifelse(status$rds_is_current, "YES (up to date)", "NO (needs update)"))
  } else {
    message("RDS cache:  NOT FOUND")
  }

  message("")

  if (status$txt_exists) {
    info <- file.info(txt_file)
    status$txt_size_mb <- round(info$size / 1e6, 1)
    status$txt_date <- as.Date(info$mtime)
    message("Raw text:   EXISTS")
    message("  Size: ", status$txt_size_mb, " MB")
    message("  Date: ", status$txt_date)
  } else {
    message("Raw text:   NOT FOUND")
  }

  message("")

  if (status$fasta_exists) {
    info <- file.info(fasta_file)
    status$fasta_size_mb <- round(info$size / 1e6, 1)
    message("BLAST FASTA: EXISTS (", status$fasta_size_mb, " MB)")
  } else {
    message("BLAST FASTA: NOT FOUND (will be created on first BLAST run)")
  }

  invisible(status)
}


#' Clear REBASE Cache
#'
#' @description Removes all cached REBASE files to force fresh download.
#' Use confirm = FALSE to clear without interactive prompt (for scripts).
#'
#' @param confirm Require confirmation before deleting. Set FALSE for non-interactive use.
#'
#' @return Invisible TRUE if cleared
#' @export
#'
#' @examples
#' \dontrun{
#' # Interactive - will ask for confirmation
#' clear_rebase_cache()
#'
#' # Non-interactive - clears immediately
#' clear_rebase_cache(confirm = FALSE)
#' }
clear_rebase_cache <- function(confirm = TRUE) {

  cache_dir <- get_rebase_cache_dir()

  files_to_delete <- c(
    file.path(cache_dir, "rebase_data.rds"),
    file.path(cache_dir, "REBASE_protein_seqs.txt"),
    file.path(cache_dir, "rebase_db.fasta"),
    file.path(cache_dir, "rebase_db.fasta.phr"),
    file.path(cache_dir, "rebase_db.fasta.pin"),
    file.path(cache_dir, "rebase_db.fasta.psq")
  )

  existing <- files_to_delete[file.exists(files_to_delete)]

  if (length(existing) == 0) {
    message("No REBASE cache files found")
    return(invisible(FALSE))
  }

  message("Files to delete:")
  for (f in existing) {
    message("  ", f)
  }

  if (confirm) {
    response <- readline("Delete these files? (y/n): ")
    if (tolower(response) != "y") {
      message("Cancelled")
      return(invisible(FALSE))
    }
  }

  for (f in existing) {
    unlink(f)
  }
  message("Cache cleared successfully")
  message("Next call to get_rebase_data() will download fresh data")

  invisible(TRUE)
}


#' Force Refresh REBASE Data
#'
#' @description Clears cache and re-downloads REBASE data in one step.
#' Useful after updating parsing logic or when data is corrupted.
#'
#' @param verbose Print progress messages
#'
#' @return Refreshed REBASE tibble
#' @export
#'
#' @examples
#' \dontrun{
#' # Clear cache and get fresh data
#' rebase <- refresh_rebase_data()
#' }
refresh_rebase_data <- function(verbose = TRUE) {

  if (verbose) message("=== Refreshing REBASE Data ===\n")

  # Step 1: Clear cache
  if (verbose) message("Step 1: Clearing cache...")
  clear_rebase_cache(confirm = FALSE)

  # Step 2: Download and parse fresh data
  if (verbose) message("\nStep 2: Downloading fresh data...")
  rebase_data <- get_rebase_data(force_download = TRUE, verbose = verbose)

  if (verbose) {
    message("\n=== Refresh Complete ===")
    message("Total sequences: ", nrow(rebase_data))

    # Show rec_seq quality check
    n_with_rec <- sum(!is.na(rebase_data$rec_seq) &
                      rebase_data$rec_seq != "" &
                      rebase_data$rec_seq != "?", na.rm = TRUE)
    message("With recognition sequences: ", n_with_rec)

    # Check for parsing errors (rec_seq should NOT contain "GenBank")
    n_bad <- sum(grepl("GenBank", rebase_data$rec_seq, ignore.case = TRUE), na.rm = TRUE)
    if (n_bad > 0) {
      warning("Parsing issue: ", n_bad, " rec_seq entries contain 'GenBank' - check parse_rebase_sequences()")
    } else {
      message("Quality check: rec_seq parsing looks correct (no GenBank contamination)")
    }
  }

  return(rebase_data)
}


#' Download REBASE Gold Standard Sequences
#'
#' @description Downloads REBASE protein sequences from FTP.
#' Uses monthly caching - REBASE updates on the 1st of each month.
#'
#' @param output_dir Directory to save files (default: package cache directory)
#' @param force_download Re-download even if current month's cache exists
#'
#' @return Path to downloaded file
#' @export
download_rebase_gold <- function(output_dir = NULL, force_download = FALSE) {

  url <- "ftp://ftp.neb.com/pub/rebase/protein_seqs.txt"

  # Use persistent package cache if no directory provided
  if (is.null(output_dir)) {
    output_dir <- get_rebase_cache_dir()
  }

  output_file <- file.path(output_dir, "REBASE_protein_seqs.txt")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  do_download <- force_download

  # Smart monthly caching logic
  if (file.exists(output_file) && !force_download) {
    file_info <- file.info(output_file)
    file_date <- as.Date(file_info$mtime)

    # REBASE updates on the 1st of each month
    current_date <- Sys.Date()
    current_month_start <- as.Date(format(current_date, "%Y-%m-01"))

    # Check if file is valid (>1MB)
    if (file_info$size > 1000000) {

       # If file is from current month, use it
       if (file_date >= current_month_start) {
         message("Using cached REBASE file (current month)")
         message("  Location: ", output_file)
         message("  Date: ", file_date)
         message("  Size: ", round(file_info$size / 1e6, 1), " MB")
         return(output_file)
       } else {
         message(sprintf("Local REBASE file (dated %s) is from previous month.", file_date))
         message(sprintf("REBASE updates monthly on the 1st. Current month: %s",
                         format(current_date, "%Y-%m")))
         message("Downloading new version...")
         do_download <- TRUE
       }
    } else {
       message("Local REBASE file is incomplete/empty. Re-downloading...")
       do_download <- TRUE
    }
  } else if (!file.exists(output_file)) {
    message("No local REBASE file found. Downloading...")
    do_download <- TRUE
  }

  if (!do_download) return(output_file)

  message("Downloading REBASE gold standard sequences from FTP...")
  message("URL: ", url)

  # Download to temporary file first to avoid corrupting existing good file
  temp_dl <- paste0(output_file, ".tmp")

  # Increase timeout for large file
  old_timeout <- getOption("timeout")
  options(timeout = 3600)
  on.exit(options(timeout = old_timeout))

  tryCatch({
    utils::download.file(
      url = url,
      destfile = temp_dl,
      mode = "wb",
      quiet = FALSE,
      method = "auto"
    )

    if (file.exists(temp_dl) && file.info(temp_dl)$size > 1000000) {
      if (file.exists(output_file)) unlink(output_file)
      file.rename(temp_dl, output_file)

      # IMPORTANT: Update the timestamp to NOW (not server's timestamp)
      Sys.setFileTime(output_file, Sys.time())

      message("Success! REBASE file saved to: ", output_file)
      message("File size: ", round(file.info(output_file)$size / 1e6, 1), " MB")

      # Invalidate RDS and FASTA caches when raw file is updated
      rds_file <- file.path(output_dir, "rebase_data.rds")
      fasta_file <- file.path(output_dir, "rebase_db.fasta")

      if (file.exists(rds_file)) {
        message("Invalidating old RDS cache...")
        unlink(rds_file)
      }
      if (file.exists(fasta_file)) {
        message("Invalidating old BLAST database...")
        unlink(fasta_file)
        unlink(paste0(fasta_file, ".phr"))
        unlink(paste0(fasta_file, ".pin"))
        unlink(paste0(fasta_file, ".psq"))
      }

      return(output_file)
    } else {
      warning("Download completed but file seems too small (<1MB)")
      if (file.exists(temp_dl)) unlink(temp_dl)
      # Fallback to existing file if available
      if (file.exists(output_file)) {
        warning("Reverting to existing (older) local file.")
        return(output_file)
      }
      return(NULL)
    }
  }, error = function(e) {
    warning("REBASE download failed: ", e$message)
    if (file.exists(temp_dl)) unlink(temp_dl)
    # Fallback
    if (file.exists(output_file)) {
      message("Using existing local file despite download failure.")
      return(output_file)
    }
    return(NULL)
  })
}


#' Parse REBASE Protein Sequences File
#'
#' @description Parses REBASE protein_seqs.txt format with Gold Standard priority.
#' Header format: >REBASE:EnzName EnzType:TYPE RecSeq:SEQUENCE GenBank:... OrgName:...
#'
#' Gold Standard entries (experimentally verified) are marked and prioritized.
#' These typically have more complete information (recognition sequences, etc.)
#'
#' @param rebase_file Path to REBASE file
#'
#' @return Tibble with parsed sequences and metadata including gold_standard flag
#' @export
parse_rebase_sequences <- function(rebase_file) {

  if (!file.exists(rebase_file)) {
    stop("File not found: ", rebase_file)
  }

  message("Parsing REBASE file: ", rebase_file)

  # Read all headers first for consistent parsing
  all_lines <- readLines(rebase_file, warn = FALSE)
  header_lines <- all_lines[grep("^>", all_lines)]
  clean_headers <- sub("^>", "", header_lines)

  message("Found ", length(clean_headers), " sequence entries")

  # ===== ROBUST FIELD EXTRACTION =====
  # Header format: >REBASE:EnzName EnzType:TYPE RecSeq:SEQUENCE GenBank:... OrgName:...
  # Fields are SPACE-separated, values end at next FieldName: or end of line

  # Helper function to extract field value safely
  extract_field <- function(headers, field_name) {
    values <- rep("", length(headers))
    has_field <- grepl(paste0(field_name, ":"), headers)

    if (any(has_field)) {
      # Extract value after FieldName: until next field (word followed by :) or end
      # Pattern: FieldName:VALUE where VALUE ends at whitespace+Word: or end
      pattern <- paste0(field_name, ":([^:]+?)(?=\\s+[A-Za-z]+:|$)")
      matches <- regmatches(headers[has_field], regexpr(pattern, headers[has_field], perl = TRUE))

      # Extract just the value part
      extracted <- sub(paste0("^", field_name, ":"), "", matches)
      extracted <- trimws(extracted)
      values[has_field] <- extracted
    }
    return(values)
  }

  # ===== EXTRACT ENZYME NAME =====
  # Format: REBASE:EnzymeName EnzType:... or REBASE:EnzymeName (space before EnzType)
  # Enzyme name is between "REBASE:" and the first space/field
  enzyme_names <- rep("", length(clean_headers))

  for (i in seq_along(clean_headers)) {
    header <- clean_headers[i]

    # Remove REBASE: prefix
    if (grepl("^REBASE:", header)) {
      header <- sub("^REBASE:", "", header)
    }

    # Extract enzyme name: everything up to first space or first FieldName:
    # EnzType:, RecSeq:, GenBank:, OrgName: are known field markers
    enz_match <- regmatches(header, regexpr("^[^\\s]+", header, perl = TRUE))
    if (length(enz_match) > 0) {
      enzyme_names[i] <- enz_match
    }
  }
  enzyme_names <- trimws(enzyme_names)

  # ===== EXTRACT OTHER FIELDS =====
  enz_types <- extract_field(clean_headers, "EnzType")
  org_names <- extract_field(clean_headers, "OrgName")

  # Debug: show sample parsed values
  message("  Sample enzyme names: ", paste(head(enzyme_names[enzyme_names != ""], 3), collapse = ", "))
  message("  Sample EnzTypes: ", paste(head(enz_types[enz_types != ""], 3), collapse = ", "))

  # ===== CAREFUL REC_SEQ EXTRACTION =====
  # rec_seq should ONLY contain sequence characters (A-Z, ^, /, parentheses)
  # Must stop BEFORE GenBank: or any other field
  rec_seqs <- rep("", length(clean_headers))
  has_recseq <- grepl("RecSeq:", clean_headers)

  if (any(has_recseq)) {
    # Extract everything after RecSeq: until next field or tab
    raw_recseq <- sub(".*RecSeq:([^\\t]+).*", "\\1", clean_headers[has_recseq], perl = TRUE)

    # Clean: keep only valid sequence characters
    # Valid: A-Z (IUPAC), ^(methylation), /(alternative), ()
    # Stop at: space followed by uppercase word (next field), digits, or GenBank
    cleaned <- gsub("\\s+GenBank.*$", "", raw_recseq, perl = TRUE)  # Remove GenBank and after
    cleaned <- gsub("\\s+[A-Z][a-z]+:.*$", "", cleaned, perl = TRUE)  # Remove next Field:
    cleaned <- gsub("\\s+\\d+.*$", "", cleaned, perl = TRUE)  # Remove numbers
    cleaned <- trimws(cleaned)
    cleaned <- toupper(cleaned)

    # Final validation: only keep if it looks like a valid recognition sequence
    # Valid patterns: letters, ^, /, (), numbers for position
    valid_pattern <- "^[A-Z^/()+0-9]+$"
    cleaned[!grepl(valid_pattern, cleaned)] <- ""

    rec_seqs[has_recseq] <- cleaned
  }

  # ===== GOLD STANDARD DETECTION =====
  # Gold standard entries are experimentally verified
  # They typically have: complete rec_seq, specific EnzType (not putative)
  # Check for indicators in header or enzyme naming
  is_gold <- rep(FALSE, length(clean_headers))

  # Gold standard indicators:
  # 1. Has valid recognition sequence (experimentally determined)
  # 2. EnzType is specific (Type I, II, III, IV) not "putative" or "probable"
  # 3. No "putative" or "probable" or "predicted" in the header
  has_valid_recseq <- rec_seqs != "" & rec_seqs != "?" & nchar(rec_seqs) >= 2
  is_putative <- grepl("(?i)putative|probable|predicted|hypothetical", clean_headers)
  has_specific_type <- grepl("(?i)Type\\s*(I|II|III|IV)|^(I|II|III|IV)$", enz_types)

  is_gold <- has_valid_recseq & !is_putative & has_specific_type

  message(sprintf("  Gold standard entries: %d (%.1f%%)",
                  sum(is_gold), 100 * sum(is_gold) / length(is_gold)))

  # ===== SEQUENCE EXTRACTION =====
  if (requireNamespace("Biostrings", quietly = TRUE)) {
    message("Using Biostrings for sequence parsing...")
    tryCatch({
      aa <- Biostrings::readAAStringSet(rebase_file)
      sequences <- as.character(aa)
      seq_lengths <- Biostrings::width(aa)
    }, error = function(e) {
      message("Biostrings failed, using base R: ", e$message)
      sequences <- NULL
    })
  }

  if (!exists("sequences") || is.null(sequences)) {
    # Base R sequence extraction
    message("Using base R for sequence parsing...")
    header_idx <- grep("^>", all_lines)
    sequences <- character(length(header_idx))
    seq_start <- header_idx + 1
    seq_end <- c(header_idx[-1] - 1, length(all_lines))

    for (i in seq_along(header_idx)) {
      if (seq_start[i] <= seq_end[i]) {
        sequences[i] <- paste(all_lines[seq_start[i]:seq_end[i]], collapse = "")
      }
    }
    sequences <- gsub("\\s+", "", sequences)
    seq_lengths <- nchar(sequences)
  }

  # ===== BUILD RESULT TABLE =====
  result <- dplyr::tibble(
    enzyme_name = enzyme_names,
    enz_type = enz_types,
    rec_seq = rec_seqs,
    org_name = org_names,
    sequence = sequences,
    seq_length = seq_lengths,
    is_gold_standard = is_gold
  )

  # Clean up
  result$enz_type[is.na(result$enz_type) | result$enz_type == ""] <- "unknown"
  result$rec_seq[is.na(result$rec_seq)] <- ""

  # Filter valid sequences (>50 aa)
  result <- result %>% dplyr::filter(!is.na(sequence) & seq_length > 50)

  # ===== R-M TYPE CLASSIFICATION =====
  result <- result %>%
    dplyr::mutate(
      rm_type = dplyr::case_when(
        grepl("(?i)^I$|Type\\s*I[^IV]|Type\\s*I$", enz_type) ~ "Type_I",
        grepl("(?i)^II|Type\\s*II", enz_type) ~ "Type_II",
        grepl("(?i)^III|Type\\s*III", enz_type) ~ "Type_III",
        grepl("(?i)^IV|Type\\s*IV|Mcr|Mrr", enz_type) ~ "Type_IV",
        grepl("(?i)^M\\.", enzyme_name) ~ "Type_II",  # M. prefix usually Type II
        grepl("(?i)^R\\.", enzyme_name) ~ "Type_II",  # R. prefix usually Type II
        grepl("(?i)^Hsd[MRS]", enzyme_name) ~ "Type_I",
        grepl("(?i)^Mod\\.|^Res\\.", enzyme_name) ~ "Type_III",
        TRUE ~ "unknown"
      ),
      subunit = dplyr::case_when(
        grepl("(?i)^M\\.|Mtase|HsdM|^Mod\\.", enzyme_name) ~ "M",
        grepl("(?i)^R\\.|Rease|HsdR|^Res\\.", enzyme_name) ~ "R",
        grepl("(?i)^S\\.|HsdS", enzyme_name) ~ "S",
        TRUE ~ "unknown"
      )
    )

  # ===== PRIORITIZE GOLD STANDARD =====
  # For duplicate enzyme names, keep gold standard version if available
  # This ensures we get the most complete information
  result <- result %>%
    dplyr::arrange(enzyme_name, dplyr::desc(is_gold_standard), dplyr::desc(nchar(rec_seq)))

  message(sprintf("Parsed %d valid REBASE sequences", nrow(result)))
  message(sprintf("  With recognition sequences: %d",
                  sum(result$rec_seq != "" & result$rec_seq != "?")))
  message(sprintf("  Gold standard: %d", sum(result$is_gold_standard)))

  # Type distribution
  type_dist <- table(result$rm_type)
  message("  Type distribution: ", paste(names(type_dist), type_dist, sep = "=", collapse = ", "))

  return(result)
}


#' Debug REBASE Parsing
#'
#' @description Shows sample parsed entries for debugging enzyme_name and enz_type extraction
#'
#' @param rebase_data Parsed REBASE data from parse_rebase_sequences()
#' @param n Number of samples to show
#'
#' @return Invisible NULL (prints debug info)
#' @export
debug_rebase_parsing <- function(rebase_data, n = 10) {
  if (is.null(rebase_data) || nrow(rebase_data) == 0) {
    message("No REBASE data to debug")
    return(invisible(NULL))
  }

  message("\n=== REBASE Parsing Debug ===\n")

  # Show column names
  message("Columns: ", paste(names(rebase_data), collapse = ", "))

  # Show sample entries
  message("\nSample entries:")
  sample_data <- head(rebase_data, n)

  for (i in seq_len(nrow(sample_data))) {
    message(sprintf("\n[%d] enzyme_name: %s", i, sample_data$enzyme_name[i]))
    message(sprintf("    enz_type: %s", sample_data$enz_type[i]))
    message(sprintf("    rm_type: %s", sample_data$rm_type[i]))
    message(sprintf("    rec_seq: %s", substr(sample_data$rec_seq[i], 1, 30)))
    message(sprintf("    subunit: %s", sample_data$subunit[i]))
    message(sprintf("    gold_standard: %s", sample_data$is_gold_standard[i]))
  }

  # Check for parsing issues
  message("\n=== Potential Issues ===")

  # Check if enzyme_name contains EnzType (parsing bug)
  has_enztype_in_name <- grepl("EnzType:", rebase_data$enzyme_name)
  if (any(has_enztype_in_name)) {
    message(sprintf("⚠ %d enzyme names contain 'EnzType:' (parsing bug!)",
                    sum(has_enztype_in_name)))
    message("  Examples: ", paste(head(rebase_data$enzyme_name[has_enztype_in_name], 3), collapse = ", "))
  } else {
    message("✓ No enzyme names contain 'EnzType:' (parsing OK)")
  }

  # Check enz_type values
  enz_type_dist <- table(rebase_data$enz_type, useNA = "ifany")
  message("\nEnzType distribution:")
  for (type in names(enz_type_dist)) {
    message(sprintf("  %s: %d", ifelse(is.na(type), "NA", type), enz_type_dist[type]))
  }

  # Check rm_type values
  rm_type_dist <- table(rebase_data$rm_type, useNA = "ifany")
  message("\nrm_type distribution:")
  for (type in names(rm_type_dist)) {
    message(sprintf("  %s: %d", ifelse(is.na(type), "NA", type), rm_type_dist[type]))
  }

  return(invisible(NULL))
}


#' Run Local BLASTP Against REBASE
#'
#' @description Runs blastp against REBASE database
#'
#' @param query_fasta Path to query FASTA file
#' @param rebase_db Path to REBASE BLAST database (or FASTA to create db)
#' @param output_file Output file path
#' @param evalue E-value threshold
#' @param num_threads Number of threads
#' @param max_target_seqs Maximum hits per query
#'
#' @return Path to BLAST results file
#' @export
run_blastp_rebase <- function(query_fasta,
                               rebase_db,
                               output_file = NULL,
                               evalue = 1e-5,
                               num_threads = 4,
                               max_target_seqs = 5) {

  # Check BLAST availability
  blastp_check <- suppressWarnings(system("which blastp", intern = TRUE))
  if (length(blastp_check) == 0) {
    stop("blastp not found. Please install NCBI BLAST+")
  }

  if (is.null(output_file)) {
    output_file <- tempfile(fileext = ".blast.txt")
  }

  # Check if we need to create BLAST database
  db_files <- paste0(rebase_db, c(".phr", ".pin", ".psq"))
  if (!all(file.exists(db_files))) {
    message("  [1/2] Creating BLAST database from REBASE FASTA...")
    message("        This takes ~2-3 minutes for 850K sequences...")
    start_db <- Sys.time()

    makeblastdb_cmd <- sprintf(
      "makeblastdb -in %s -dbtype prot -out %s -title 'REBASE' 2>&1",
      shQuote(rebase_db), shQuote(rebase_db)
    )
    db_result <- system(makeblastdb_cmd, intern = TRUE)

    elapsed_db <- round(difftime(Sys.time(), start_db, units = "secs"), 1)
    message("        Database created in ", elapsed_db, " seconds")
  } else {
    message("  [1/2] Using existing BLAST database")
  }

  # Count query sequences
  query_count <- length(grep("^>", readLines(query_fasta, warn = FALSE)))
  message("  [2/2] Running BLASTP: ", query_count, " queries vs REBASE database...")
  message("        E-value threshold: ", evalue)
  message("        Threads: ", num_threads)
  message("        This may take 5-15 minutes depending on query count...")

  start_blast <- Sys.time()

  blast_cmd <- sprintf(
    "blastp -query %s -db %s -out %s -evalue %s -num_threads %d -max_target_seqs %d -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'",
    shQuote(query_fasta),
    shQuote(rebase_db),
    shQuote(output_file),
    evalue,
    num_threads,
    max_target_seqs
  )

  system(blast_cmd)

  elapsed_blast <- round(difftime(Sys.time(), start_blast, units = "secs"), 1)

  if (file.exists(output_file) && file.info(output_file)$size > 0) {
    n_hits <- length(readLines(output_file, warn = FALSE))
    message("        BLAST completed in ", elapsed_blast, " seconds")
    message("        Total hits: ", n_hits)
    return(output_file)
  } else {
    warning("BLAST produced no results (", elapsed_blast, " seconds)")
    return(NULL)
  }
}


#' Parse BLAST Tabular Output
#'
#' @description Parses BLAST format 6 output with improved enzyme name extraction.
#'
#' @param blast_file BLAST output file
#' @param verbose Print diagnostic info
#'
#' @return Tibble with parsed results
#' @export
parse_blast_results <- function(blast_file, verbose = TRUE) {

  if (!file.exists(blast_file) || file.info(blast_file)$size == 0) {
    message("No BLAST results to parse")
    return(dplyr::tibble())
  }

  col_names <- c("query_id", "subject_id", "pident", "length",
                 "mismatch", "gapopen", "qstart", "qend",
                 "sstart", "send", "evalue", "bitscore")

  result <- utils::read.delim(
    blast_file,
    header = FALSE,
    col.names = col_names,
    stringsAsFactors = FALSE
  )

  result <- dplyr::as_tibble(result) %>%
    dplyr::mutate(
      # Use pct_identity to avoid conflict with base::identity function
      pct_identity = pident / 100,
      # Extract enzyme name from subject - take the full subject_id
      # (should already be clean from sanitized FASTA writing)
      rebase_enzyme = subject_id
    )

  if (verbose) {
    message(sprintf("Parsed %d BLAST hits", nrow(result)))

    # ===== DIAGNOSTIC CHECKS =====
    n_queries <- dplyr::n_distinct(result$query_id)
    n_subjects <- dplyr::n_distinct(result$subject_id)
    n_enzymes <- dplyr::n_distinct(result$rebase_enzyme)

    message(sprintf("  Unique queries: %d", n_queries))
    message(sprintf("  Unique subjects (REBASE enzymes): %d", n_subjects))

    # Check for suspicious patterns
    if (n_subjects == 1 && n_queries > 1) {
      message("  ⚠ WARNING: All queries matched to the SAME enzyme!")
      message("    This is unusual - check BLAST database and query files")
      message("    Single subject: ", unique(result$subject_id)[1])
    }

    # Show e-value distribution
    evalue_ranges <- cut(result$evalue,
                         breaks = c(0, 1e-50, 1e-20, 1e-10, 1e-5, 1, Inf),
                         labels = c("<1e-50", "1e-50 to 1e-20", "1e-20 to 1e-10",
                                    "1e-10 to 1e-5", "1e-5 to 1", ">1"))
    message("  E-value distribution:")
    ev_table <- table(evalue_ranges)
    for (i in seq_along(ev_table)) {
      if (ev_table[i] > 0) {
        message(sprintf("    %s: %d hits", names(ev_table)[i], ev_table[i]))
      }
    }

    # Sample of matches
    if (nrow(result) > 0) {
      message("  Sample matches (first 5):")
      sample_matches <- result %>%
        dplyr::select(query_id, rebase_enzyme, pident, evalue) %>%
        dplyr::slice_head(n = 5)
      for (i in seq_len(nrow(sample_matches))) {
        message(sprintf("    %s -> %s (%.1f%% identity, e=%.2e)",
                        sample_matches$query_id[i],
                        sample_matches$rebase_enzyme[i],
                        sample_matches$pident[i],
                        sample_matches$evalue[i]))
      }
    }
  }

  return(result)
}


#' Debug BLAST Matching Results
#'
#' @description Diagnostic function to check BLAST matching quality.
#' Use this when results seem suspicious (e.g., all queries matching same enzyme).
#'
#' @param blast_results BLAST results from parse_blast_results()
#' @param rebase_data REBASE data for additional context
#' @param query_fasta Path to query FASTA (optional, for verification)
#' @param rebase_fasta Path to REBASE FASTA (optional, for verification)
#'
#' @return List with diagnostic information
#' @export
debug_blast_matching <- function(blast_results,
                                  rebase_data = NULL,
                                  query_fasta = NULL,
                                  rebase_fasta = NULL) {

  message("=== BLAST Matching Diagnostics ===\n")

  diagnostics <- list()

  if (is.null(blast_results) || nrow(blast_results) == 0) {
    message("No BLAST results to diagnose")
    return(diagnostics)
  }

  # 1. Basic statistics
  diagnostics$n_total_hits <- nrow(blast_results)
  diagnostics$n_unique_queries <- dplyr::n_distinct(blast_results$query_id)
  diagnostics$n_unique_subjects <- dplyr::n_distinct(blast_results$subject_id)

  message("1. Basic Statistics:")
  message(sprintf("   Total hits: %d", diagnostics$n_total_hits))
  message(sprintf("   Unique queries: %d", diagnostics$n_unique_queries))
  message(sprintf("   Unique REBASE matches: %d", diagnostics$n_unique_subjects))

  # 2. Check for over-concentration
  subject_counts <- blast_results %>%
    dplyr::count(subject_id, sort = TRUE)

  diagnostics$top_subjects <- head(subject_counts, 10)
  diagnostics$subject_concentration <- subject_counts$n[1] / nrow(blast_results)

  message("\n2. Top 5 matched REBASE enzymes:")
  for (i in seq_len(min(5, nrow(subject_counts)))) {
    message(sprintf("   %s: %d hits (%.1f%%)",
                    subject_counts$subject_id[i],
                    subject_counts$n[i],
                    100 * subject_counts$n[i] / nrow(blast_results)))
  }

  if (diagnostics$subject_concentration > 0.5 && diagnostics$n_unique_queries > 3) {
    message("   ⚠ HIGH CONCENTRATION: Top enzyme has >50% of all hits")
    message("     This may indicate a database or query issue")
  }

  # 3. Check identity distribution
  diagnostics$identity_stats <- list(
    min = min(blast_results$pct_identity, na.rm = TRUE),
    max = max(blast_results$pct_identity, na.rm = TRUE),
    mean = mean(blast_results$pct_identity, na.rm = TRUE),
    median = median(blast_results$pct_identity, na.rm = TRUE)
  )

  message("\n3. Identity distribution:")
  message(sprintf("   Min: %.1f%%", diagnostics$identity_stats$min * 100))
  message(sprintf("   Max: %.1f%%", diagnostics$identity_stats$max * 100))
  message(sprintf("   Mean: %.1f%%", diagnostics$identity_stats$mean * 100))
  message(sprintf("   Median: %.1f%%", diagnostics$identity_stats$median * 100))

  # 4. Check e-value distribution
  diagnostics$evalue_stats <- list(
    min = min(blast_results$evalue, na.rm = TRUE),
    max = max(blast_results$evalue, na.rm = TRUE),
    n_significant = sum(blast_results$evalue < 1e-5, na.rm = TRUE)
  )

  message("\n4. E-value distribution:")
  message(sprintf("   Best (lowest): %.2e", diagnostics$evalue_stats$min))
  message(sprintf("   Worst (highest): %.2e", diagnostics$evalue_stats$max))
  message(sprintf("   Significant hits (< 1e-5): %d (%.1f%%)",
                  diagnostics$evalue_stats$n_significant,
                  100 * diagnostics$evalue_stats$n_significant / nrow(blast_results)))

  # 5. Check FASTA files if provided
  if (!is.null(query_fasta) && file.exists(query_fasta)) {
    query_headers <- grep("^>", readLines(query_fasta, warn = FALSE), value = TRUE)
    diagnostics$n_query_seqs <- length(query_headers)
    message("\n5. Query FASTA check:")
    message(sprintf("   Sequences in file: %d", diagnostics$n_query_seqs))
    message("   Sample headers:")
    for (h in head(query_headers, 3)) {
      message(sprintf("     %s", h))
    }
  }

  if (!is.null(rebase_fasta) && file.exists(rebase_fasta)) {
    rebase_headers <- grep("^>", readLines(rebase_fasta, n = 100, warn = FALSE), value = TRUE)
    diagnostics$n_rebase_sample <- length(rebase_headers)
    message("\n6. REBASE FASTA check (first 100 lines):")
    message("   Sample headers:")
    for (h in head(rebase_headers, 5)) {
      message(sprintf("     %s", h))
    }
  }

  # 6. Cross-check with REBASE data
  if (!is.null(rebase_data) && nrow(rebase_data) > 0) {
    matched_enzymes <- unique(blast_results$subject_id)
    in_rebase <- matched_enzymes %in% rebase_data$enzyme_name

    diagnostics$n_matched_in_rebase <- sum(in_rebase)
    diagnostics$n_not_in_rebase <- sum(!in_rebase)

    message("\n7. Cross-check with REBASE database:")
    message(sprintf("   Matched enzymes found in REBASE: %d", diagnostics$n_matched_in_rebase))
    message(sprintf("   Matched enzymes NOT in REBASE: %d", diagnostics$n_not_in_rebase))

    if (diagnostics$n_not_in_rebase > 0) {
      message("   Enzymes not found:")
      for (e in head(matched_enzymes[!in_rebase], 5)) {
        message(sprintf("     - %s", e))
      }
    }
  }

  message("\n=== End Diagnostics ===")

  invisible(diagnostics)
}


#' Filter BLAST Results
#'
#' @description Filters BLAST results to remove low-quality matches before scoring.
#' Removes matches with:
#' - identity < min_identity (default 10%)
#' - alignment length < min_length (default 50 aa)
#' - e-value > max_evalue (default 1e-3)
#'
#' @param blast_results BLAST results from parse_blast_results()
#' @param min_identity Minimum identity threshold (0-1). Default 0.10 (10%)
#' @param min_length Minimum alignment length. Default 50
#' @param max_evalue Maximum e-value. Default 1e-3
#' @param verbose Print filtering statistics
#'
#' @return Filtered tibble
#' @export
filter_blast_results <- function(blast_results,
                                  min_identity = 0.10,
                                  min_length = 50,
                                  max_evalue = 1e-3,
                                  verbose = TRUE) {

  if (is.null(blast_results) || nrow(blast_results) == 0) {
    if (verbose) message("No BLAST results to filter")
    return(dplyr::tibble())
  }

  n_before <- nrow(blast_results)
  n_queries_before <- dplyr::n_distinct(blast_results$query_id)

  # Apply filters
  filtered <- blast_results %>%
    dplyr::filter(
      pct_identity >= min_identity,
      length >= min_length,
      evalue <= max_evalue
    )

  n_after <- nrow(filtered)
  n_queries_after <- dplyr::n_distinct(filtered$query_id)

  if (verbose) {
    message(sprintf("BLAST filtering (identity >= %.0f%%, length >= %d, evalue <= %.0e):",
                    min_identity * 100, min_length, max_evalue))
    message(sprintf("  Hits: %d -> %d (removed %d)",
                    n_before, n_after, n_before - n_after))
    message(sprintf("  Queries with hits: %d -> %d",
                    n_queries_before, n_queries_after))
  }

  return(filtered)
}


#' Write Sequences to FASTA for BLAST
#'
#' @description Writes sequences to FASTA format for BLAST input.
#' Sanitizes IDs to be BLAST-friendly (removes spaces, special chars).
#'
#' @param data Data frame with sequences
#' @param output_file Output FASTA path
#' @param id_col ID column name
#' @param seq_col Sequence column name
#'
#' @return Path to output file
#' @export
write_fasta_for_blast <- function(data, output_file, id_col = "locus_tag", seq_col = "translation") {

  # ===== VALIDATE INPUT DATA =====
  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
    warning("No data to write to FASTA")
    return(NULL)
  }

  message("  write_fasta_for_blast: Input data has ", nrow(data), " rows, ", ncol(data), " cols")
  message("    Available columns: ", paste(names(data), collapse = ", "))

  if (!id_col %in% names(data)) {
    id_col <- intersect(c("locus_tag", "protein_id", "gene_id", "id"), names(data))[1]
    if (!is.na(id_col)) {
      message("    Using detected ID column: ", id_col)
    }
  }

  if (is.na(id_col) || !id_col %in% names(data)) {
    stop("ID column not found. Available columns: ", paste(names(data), collapse = ", "))
  }

  if (!seq_col %in% names(data)) {
    stop("Sequence column not found: ", seq_col, ". Available: ", paste(names(data), collapse = ", "))
  }

  # ===== CHECK ID COLUMN CONTENT =====
  id_values <- data[[id_col]]
  n_empty_ids <- sum(is.na(id_values) | id_values == "")
  if (n_empty_ids > 0) {
    message("    WARNING: ", n_empty_ids, " rows have empty/NA ID values")
  }

  # Filter valid sequences (using base R for reliability)
  seq_values <- data[[seq_col]]
  valid_mask <- !is.na(seq_values) & nchar(seq_values) > 50
  valid_data <- data[valid_mask, , drop = FALSE]

  if (nrow(valid_data) == 0) {
    warning("No valid sequences to write")
    return(NULL)
  }

  # ===== SANITIZE IDS FOR BLAST =====
  # BLAST uses first word (up to space) as ID
  # Enzyme names like "M.AaaI" should work, but ensure no problematic chars
  sanitize_id <- function(x) {
    # Replace spaces with underscores
    x <- gsub("\\s+", "_", x)
    # Remove problematic characters that could break BLAST parsing
    x <- gsub("[|;,()\\[\\]\\{\\}]", "_", x)
    # Keep only alphanumeric, underscore, dot, hyphen
    x <- gsub("[^A-Za-z0-9_.-]", "", x)
    # Ensure unique - add index if duplicated
    x
  }

  # Use base R to avoid dplyr .data pronoun issues
  valid_data$.blast_id <- sanitize_id(valid_data[[id_col]])
  valid_data$.original_id <- valid_data[[id_col]]

  # Check for duplicates after sanitization and make unique
  dup_mask <- duplicated(valid_data$.blast_id) | duplicated(valid_data$.blast_id, fromLast = TRUE)
  if (any(dup_mask)) {
    # Make unique by appending row index
    valid_data$.blast_id[dup_mask] <- paste0(valid_data$.blast_id[dup_mask], "_", seq_len(nrow(valid_data))[dup_mask])
  }

  n_seqs <- nrow(valid_data)

  # For large files, write in chunks to show progress
  if (n_seqs > 10000) {
    message(sprintf("  Writing %d sequences to FASTA (chunked)...", n_seqs))

    # Open file connection
    con <- file(output_file, "w")
    on.exit(close(con))

    chunk_size <- 50000
    n_chunks <- ceiling(n_seqs / chunk_size)

    for (chunk in seq_len(n_chunks)) {
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, n_seqs)

      chunk_data <- valid_data[start_idx:end_idx, ]

      # Build FASTA for this chunk - use sanitized ID
      fasta_lines <- character(nrow(chunk_data) * 2)
      for (i in seq_len(nrow(chunk_data))) {
        fasta_lines[(i-1)*2 + 1] <- paste0(">", chunk_data$.blast_id[i])
        fasta_lines[(i-1)*2 + 2] <- chunk_data[[seq_col]][i]
      }

      writeLines(fasta_lines, con)

      # Progress
      pct <- round(100 * end_idx / n_seqs)
      message(sprintf("  Progress: %d/%d sequences (%d%%)", end_idx, n_seqs, pct))
    }

    message(sprintf("  Wrote %d sequences to %s", n_seqs, output_file))

  } else {
    # Small file - write all at once
    fasta_lines <- character(n_seqs * 2)
    for (i in seq_len(n_seqs)) {
      fasta_lines[(i-1)*2 + 1] <- paste0(">", valid_data$.blast_id[i])
      fasta_lines[(i-1)*2 + 2] <- valid_data[[seq_col]][i]
    }

    writeLines(fasta_lines, output_file)
    message(sprintf("Wrote %d sequences to %s", n_seqs, output_file))
  }

  return(output_file)
}


#' Extract R-M Type from REBASE Enzyme Name
#'
#' @description Determines R-M type from enzyme naming convention
#'
#' @param enzyme_name REBASE enzyme name
#'
#' @return Character: Type_I, Type_II, etc.
#' @export
extract_rm_type_from_rebase <- function(enzyme_name) {
  if (is.na(enzyme_name) || enzyme_name == "") return("unknown")

  dplyr::case_when(
    # Type I: HsdM, HsdR, HsdS
    stringr::str_detect(enzyme_name, "(?i)Hsd[MRS]") ~ "Type_I",
    stringr::str_detect(enzyme_name, "(?i)EcoK|EcoB|EcoA|StySB") ~ "Type_I",

    # Type III: Mod, Res
    stringr::str_detect(enzyme_name, "(?i)^Mod\\.") ~ "Type_III",
    stringr::str_detect(enzyme_name, "(?i)^Res\\.") ~ "Type_III",
    stringr::str_detect(enzyme_name, "(?i)EcoP") ~ "Type_III",

    # Type IV: Mcr, Mrr
    stringr::str_detect(enzyme_name, "(?i)^Mcr|^Mrr") ~ "Type_IV",

    # Type II: M. or R. prefix (most common)
    stringr::str_detect(enzyme_name, "(?i)^M\\.") ~ "Type_II",
    stringr::str_detect(enzyme_name, "(?i)^R\\.") ~ "Type_II",

    TRUE ~ "unknown"
  )
}


#' Get Best REBASE Match with Type
#'
#' @description Returns best match per query including type info.
#' Prioritizes Gold Standard entries when multiple matches exist.
#'
#' @param blast_results BLAST results from parse_blast_results()
#' @param rebase_data Parsed REBASE data (with is_gold_standard column)
#' @param min_identity Minimum identity threshold
#'
#' @return Tibble with best matches including rec_seq
#' @export
get_best_rebase_match <- function(blast_results, rebase_data = NULL, min_identity = 0.3) {

  if (nrow(blast_results) == 0) {
    return(dplyr::tibble(
      query_id = character(),
      best_rebase_match = character(),
      best_match_identity = numeric(),
      rm_type = character(),
      subunit = character(),
      rec_seq = character(),
      is_gold_standard = logical()
    ))
  }

  # Filter by minimum identity
  result <- blast_results %>%
    dplyr::filter(pct_identity >= min_identity)

  if (nrow(result) == 0) {
    message("No BLAST hits passing identity threshold")
    return(dplyr::tibble())
  }

  # Join REBASE info BEFORE selecting best match (to use gold standard for tiebreaking)
  if (!is.null(rebase_data) && nrow(rebase_data) > 0) {
    # Check if is_gold_standard column exists
    has_gold_standard <- "is_gold_standard" %in% names(rebase_data)

    # ★ Strip _\d+ suffix from rebase_enzyme (added for duplicate handling in FASTA)
    # e.g., "M.Gst45ORF7625P_123" -> "M.Gst45ORF7625P"
    result <- result %>%
      dplyr::mutate(
        .rebase_enzyme_clean = sub("_\\d+$", "", rebase_enzyme)
      )

    # Create lookup table prioritizing Gold Standard
    # ★ enz_type에서 rm_type 추출, enzyme_name에서 subunit 추출
    if (has_gold_standard) {
      rebase_lookup <- rebase_data %>%
        dplyr::select(enzyme_name, enz_type, rec_seq, is_gold_standard) %>%
        dplyr::group_by(enzyme_name) %>%
        dplyr::arrange(
          dplyr::desc(dplyr::if_else(is.na(is_gold_standard), FALSE, is_gold_standard)),
          dplyr::desc(nchar(rec_seq))
        ) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()
    } else {
      # No gold standard column - just use rec_seq length
      rebase_lookup <- rebase_data %>%
        dplyr::select(enzyme_name, enz_type, rec_seq) %>%
        dplyr::mutate(is_gold_standard = FALSE) %>%
        dplyr::group_by(enzyme_name) %>%
        dplyr::arrange(dplyr::desc(nchar(rec_seq))) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()
    }

    # ★ enz_type에서 rm_type 추출, enzyme_name에서 subunit 추출
    rebase_lookup <- rebase_lookup %>%
      dplyr::mutate(
        # enz_type에서 Type I~IV만 추출 (Type IIG 등 서브타입 제외)
        rm_type = dplyr::case_when(
          grepl("Type I[^IV]|Type I$", enz_type) ~ "Type I",
          grepl("Type II[^I]|Type II$", enz_type) ~ "Type II",
          grepl("Type III", enz_type) ~ "Type III",
          grepl("Type IV", enz_type) ~ "Type IV",
          TRUE ~ NA_character_
        ),
        # enzyme_name 앞의 접두사에서 subunit 추출 (M., R., S. -> M, R, S, 없으면 R)
        subunit = dplyr::case_when(
          grepl("^M\\.", enzyme_name) ~ "M",
          grepl("^R\\.", enzyme_name) ~ "R",
          grepl("^S\\.", enzyme_name) ~ "S",
          TRUE ~ "R"  # 접두사 없으면 R
        )
      )

    # ★ Join using cleaned enzyme name (without _\d+ suffix)
    result <- result %>%
      dplyr::left_join(
        rebase_lookup,
        by = c(".rebase_enzyme_clean" = "enzyme_name")
      ) %>%
      dplyr::select(-.rebase_enzyme_clean)

    # Now select best match per query
    # Priority: 1) Higher identity, 2) Gold standard, 3) Has rec_seq
    result <- result %>%
      dplyr::group_by(query_id) %>%
      dplyr::arrange(
        dplyr::desc(pct_identity),
        dplyr::desc(dplyr::if_else(is.na(is_gold_standard), FALSE, is_gold_standard)),
        dplyr::desc(nchar(dplyr::if_else(is.na(rec_seq), "", rec_seq)))
      ) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
  } else {
    # No REBASE data - just select by identity
    result <- result %>%
      dplyr::group_by(query_id) %>%
      dplyr::slice_max(order_by = pct_identity, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        rm_type = NA_character_,
        subunit = NA_character_,
        rec_seq = NA_character_,
        is_gold_standard = NA
      )
  }

  # Final output - enz_type에서 추출한 rm_type, enzyme_name에서 추출한 subunit
  result %>%
    dplyr::select(
      query_id,
      best_rebase_match = rebase_enzyme,
      best_match_identity = pct_identity,
      rm_type,           # ★ enz_type에서 추출 (Type I~IV만)
      subunit,           # ★ enzyme_name 접두사에서 추출 (M., R., S. -> M/R/S, 없으면 R)
      rec_seq,
      dplyr::any_of(c("is_gold_standard", "evalue", "bitscore", "length"))
    )
}


#' Annotate R-M Candidates with REBASE Results
#'
#' @description Adds REBASE match info to MTase/REase data including recognition sequences
#'
#' @param rm_data MTase or REase candidate data
#' @param rebase_results Results from get_best_rebase_match()
#' @param id_col ID column name
#'
#' @return Data with REBASE columns added (best_rebase_match, best_match_identity, rm_type, subunit, rec_seq, etc.)
#' @export
annotate_with_rebase <- function(rm_data, rebase_results, id_col = "locus_tag") {

  # Get actual id column
  if (!id_col %in% names(rm_data)) {
    id_col <- intersect(c("locus_tag", "protein_id"), names(rm_data))[1]
  }

  if (is.null(rebase_results) || nrow(rebase_results) == 0) {
    # Add empty columns when no REBASE results
    rm_data$best_rebase_match <- NA_character_
    rm_data$best_match_identity <- NA_real_
    rm_data$rm_type <- NA_character_
    rm_data$subunit <- NA_character_
    rm_data$rec_seq <- NA_character_
    return(rm_data)
  }

  # Select columns to join (rm_type and subunit directly)
  rebase_cols <- rebase_results %>%
    dplyr::select(
      query_id,
      dplyr::any_of(c("best_rebase_match", "best_match_identity",
                      "rm_type", "subunit",
                      "rec_seq", "evalue", "bitscore"))
    ) %>%
    dplyr::rename(!!id_col := query_id)

  # Ensure required columns exist
  if (!"rec_seq" %in% names(rebase_cols)) rebase_cols$rec_seq <- NA_character_
  if (!"rm_type" %in% names(rebase_cols)) rebase_cols$rm_type <- NA_character_
  if (!"subunit" %in% names(rebase_cols)) rebase_cols$subunit <- NA_character_

  result <- rm_data %>%
    dplyr::left_join(rebase_cols, by = id_col)

  # Fill NA for unmatched entries
  if (!"best_rebase_match" %in% names(result)) result$best_rebase_match <- NA_character_
  if (!"best_match_identity" %in% names(result)) result$best_match_identity <- NA_real_
  if (!"rm_type" %in% names(result)) result$rm_type <- NA_character_
  if (!"subunit" %in% names(result)) result$subunit <- NA_character_
  if (!"rec_seq" %in% names(result)) result$rec_seq <- NA_character_

  return(result)
}


#' Full REBASE Comparison Pipeline
#'
#' @description Complete pipeline: download REBASE, run BLAST, annotate results.
#' REBASE data is cached in a persistent directory that survives R sessions.
#'
#' Cache location: tools::R_user_dir("DefenseViz", which = "cache")
#' - REBASE_protein_seqs.txt: Raw downloaded file
#' - rebase_data.rds: Parsed R object
#' - rebase_db.fasta: FASTA for BLAST
#' - rebase_db.fasta.phr/pin/psq: BLAST database files
#'
#' @param query_data Data frame with sequences (needs locus_tag and translation)
#' @param output_dir Directory for query-specific intermediate files
#' @param evalue BLAST e-value threshold
#' @param min_identity Minimum identity for reporting
#' @param use_local_blast Use local BLAST (TRUE) or simple comparison (FALSE)
#' @param verbose Print detailed progress messages
#'
#' @return Tibble with REBASE annotations
#' @export
#'
#' @examples
#' \dontrun{
#' # Check cache status first
#' rebase_cache_status()
#'
#' # Run comparison
#' results <- run_rebase_comparison(my_proteins)
#'
#' # Get REBASE data as table
#' rebase <- get_rebase_data()
#' }
run_rebase_comparison <- function(query_data,
                                   output_dir = tempdir(),
                                   evalue = 1e-5,
                                   min_identity = 0.3,
                                   use_local_blast = TRUE,
                                   verbose = TRUE) {

  # Find ID column
  id_col <- intersect(c("locus_tag", "protein_id", "gene_id", "id"), names(query_data))[1]
  if (is.na(id_col)) {
    stop("No ID column found in query_data. Expected: locus_tag, protein_id, gene_id, or id. ",
         "Available: ", paste(names(query_data), collapse = ", "))
  }

  # Get persistent cache directory for REBASE files
  cache_dir <- get_rebase_cache_dir()

  # Step 1: Get REBASE data (uses monthly caching)
  if (verbose) {
    message("\n=== Step 1: Get REBASE Data ===")
    message("Cache directory: ", cache_dir)
  }

  rebase_data <- get_rebase_data(verbose = verbose)

  if (nrow(rebase_data) == 0) {
    warning("No REBASE sequences available")
    return(dplyr::tibble())
  }

  # ===== VALIDATE REBASE DATA STRUCTURE =====
  required_cols <- c("enzyme_name", "sequence")
  missing_cols <- setdiff(required_cols, names(rebase_data))
  if (length(missing_cols) > 0) {
    stop("REBASE data missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nAvailable columns: ", paste(names(rebase_data), collapse = ", "),
         "\nTry running: clear_rebase_cache() then re-run the pipeline")
  }

  # Check for empty enzyme_name (parsing bug indicator)
  n_empty_enzyme <- sum(is.na(rebase_data$enzyme_name) | rebase_data$enzyme_name == "")
  if (n_empty_enzyme > nrow(rebase_data) * 0.5) {
    stop("More than 50% of REBASE entries have empty enzyme_name. ",
         "This indicates a parsing issue. ",
         "Try running: clear_rebase_cache() then re-run the pipeline")
  }

  if (verbose) {
    message("  REBASE data validated: ", nrow(rebase_data), " sequences")
    message("    enzyme_name non-empty: ", sum(rebase_data$enzyme_name != ""), "/", nrow(rebase_data))
  }

  # Step 2: Write query sequences (to output_dir - these are project-specific)
  if (verbose) message("\n=== Step 2: Prepare query sequences ===")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  query_fasta <- file.path(output_dir, "query_sequences.fasta")
  write_fasta_for_blast(query_data, query_fasta, id_col)

  # Step 3: Create/Use REBASE FASTA for BLAST (in cache_dir - persistent)
  rebase_fasta <- file.path(cache_dir, "rebase_db.fasta")

  # Check if REBASE FASTA already exists and is valid (>100MB for ~850K sequences)
  if (file.exists(rebase_fasta) && file.info(rebase_fasta)$size > 100000000) {
    if (verbose) {
      message("\n=== Step 3: Using cached REBASE FASTA ===")
      message("Location: ", rebase_fasta)
      message("Size: ", round(file.info(rebase_fasta)$size / 1e6, 1), " MB")
    }
  } else {
    if (verbose) {
      message("\n=== Step 3: Writing REBASE FASTA (", nrow(rebase_data), " sequences) ===")
      message("Destination: ", rebase_fasta)
      message("This may take 1-2 minutes for ~850K sequences...")
    }

    # ===== SAFE RENAME: Check columns exist before renaming =====
    rebase_for_fasta <- rebase_data
    if (!"locus_tag" %in% names(rebase_for_fasta) && "enzyme_name" %in% names(rebase_for_fasta)) {
      rebase_for_fasta <- rebase_for_fasta %>% dplyr::rename(locus_tag = enzyme_name)
    }
    if (!"translation" %in% names(rebase_for_fasta) && "sequence" %in% names(rebase_for_fasta)) {
      rebase_for_fasta <- rebase_for_fasta %>% dplyr::rename(translation = sequence)
    }

    write_fasta_for_blast(
      rebase_for_fasta,
      rebase_fasta,
      id_col = "locus_tag",
      seq_col = "translation"
    )

    if (verbose) message("REBASE FASTA saved to cache: ", rebase_fasta)
  }

  # Step 4: Run BLAST or comparison
  if (verbose) message("\n=== Step 4: Run BLASTP ===")

  if (use_local_blast) {
    # BLAST results go to output_dir (project-specific)
    blast_output <- file.path(output_dir, "blast_results.txt")

    # BLAST database is in cache_dir (persistent)
    blast_file <- run_blastp_rebase(
      query_fasta, rebase_fasta, blast_output,
      evalue = evalue
    )

    if (!is.null(blast_file)) {
      blast_results <- parse_blast_results(blast_file)
    } else {
      blast_results <- dplyr::tibble()
    }
  } else {
    # Fallback: simple k-mer comparison (slower but no BLAST needed)
    if (verbose) message("Using k-mer based comparison (BLAST not available)")
    blast_results <- compare_by_kmer(query_data, rebase_data, id_col)
  }

  # Step 5: Get best matches
  if (verbose) message("\n=== Step 5: Get best matches ===")
  best_matches <- get_best_rebase_match(blast_results, rebase_data, min_identity)

  if (verbose) {
    message(sprintf("Found %d queries with REBASE hits (>= %.0f%% identity)",
                    nrow(best_matches), min_identity * 100))
    message("\n=== REBASE Comparison Complete ===")
    message("Cache location: ", cache_dir)
  }

  return(best_matches)
}


#' Search REBASE by Recognition Sequence
#'
#' @description Find enzymes with matching or similar recognition sequences
#'
#' @param pattern Recognition sequence pattern (supports IUPAC ambiguity codes)
#' @param exact If TRUE, exact match only. If FALSE, pattern matching
#'
#' @return Tibble with matching enzymes
#' @export
#'
#' @examples
#' \dontrun{
#' # Find enzymes recognizing GAATTC (EcoRI site)
#' search_rebase_by_recognition("GAATTC")
#'
#' # Find enzymes with recognition sequences containing GATC
#' search_rebase_by_recognition("GATC", exact = FALSE)
#' }
search_rebase_by_recognition <- function(pattern, exact = FALSE) {

  rebase_data <- get_rebase_data(verbose = FALSE)

  if (exact) {
    matches <- rebase_data %>%
      dplyr::filter(rec_seq == toupper(pattern))
  } else {
    matches <- rebase_data %>%
      dplyr::filter(grepl(toupper(pattern), rec_seq, ignore.case = TRUE))
  }

  if (nrow(matches) == 0) {
    message("No enzymes found with recognition sequence matching: ", pattern)
    return(dplyr::tibble())
  }

  result <- matches %>%
    dplyr::select(enzyme_name, rec_seq, rm_type, subunit, enz_type, org_name) %>%
    dplyr::distinct() %>%
    dplyr::arrange(rm_type, enzyme_name)

  message(sprintf("Found %d enzymes with recognition sequence matching '%s'",
                  nrow(result), pattern))

  return(result)
}


#' Get REBASE Summary Statistics
#'
#' @description Returns summary statistics about the REBASE database
#'
#' @return List with summary statistics
#' @export
get_rebase_summary <- function() {

  rebase_data <- get_rebase_data(verbose = FALSE)

  summary_stats <- list(
    total_sequences = nrow(rebase_data),
    unique_enzymes = dplyr::n_distinct(rebase_data$enzyme_name),
    with_recognition_seq = sum(!is.na(rebase_data$rec_seq) & rebase_data$rec_seq != "" & rebase_data$rec_seq != "?"),
    by_rm_type = table(rebase_data$rm_type),
    by_subunit = table(rebase_data$subunit),
    unique_organisms = dplyr::n_distinct(rebase_data$org_name)
  )

  message("=== REBASE Database Summary ===")
  message("Total sequences:       ", summary_stats$total_sequences)
  message("Unique enzyme names:   ", summary_stats$unique_enzymes)
  message("With recognition seq:  ", summary_stats$with_recognition_seq)
  message("Unique organisms:      ", summary_stats$unique_organisms)
  message("")
  message("By R-M Type:")
  print(summary_stats$by_rm_type)
  message("")
  message("By Subunit:")
  print(summary_stats$by_subunit)

  invisible(summary_stats)
}


# ============================================================
# HMMER vs BLAST for R-M System Detection
# ============================================================
#
# HMMER (Hidden Markov Models) vs BLAST comparison:
#
# BLAST Advantages:
# - Faster for sequence-to-sequence comparison
# - Good for finding close homologs (>30% identity)
# - REBASE is already a sequence database (BLAST-ready)
# - Simple to interpret results (% identity, e-value)
#
# HMMER Advantages:
# - Higher sensitivity for remote homologs (<30% identity)
# - Better for detecting divergent R-M systems
# - PFAM domains are HMM-based (already used in MTase/REase detection)
# - Can detect partial/truncated domains
#
# Recommendation for R-M Analysis:
# 1. Use PFAM HMM search (already implemented) for INITIAL detection
#    - This catches divergent R-M genes that BLAST might miss
# 2. Use BLAST against REBASE for TYPE ASSIGNMENT and rec_seq
#    - REBASE has specific type/recognition sequence annotations
# 3. Consider HMMER against REBASE for difficult cases
#
# The current pipeline already uses a hybrid approach:
# - Step 2-3: PFAM domain detection (HMM-based via InterProScan/HMMER)
# - Step 5-6: BLAST against REBASE (sequence similarity)
#
# Future Enhancement: Could add hmmsearch against REBASE HMM profiles
# for more sensitive detection. Would require building HMM profiles
# from REBASE sequence families.
# ============================================================


#' Check HMMER Availability
#'
#' @description Checks if HMMER is installed and available
#'
#' @return TRUE if hmmsearch is available, FALSE otherwise
#' @export
check_hmmer_available <- function() {
  hmmer_check <- suppressWarnings(system("which hmmsearch", intern = TRUE))
  return(length(hmmer_check) > 0)
}


#' Run HMMER Search Against PFAM (Alternative to BLAST)
#'
#' @description Uses hmmscan against PFAM for more sensitive R-M detection.
#' This is an alternative/complement to BLAST when detecting divergent R-M systems.
#'
#' Note: This requires PFAM database to be installed locally.
#' Download from: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/
#'
#' @param query_fasta Path to query FASTA file
#' @param pfam_db Path to PFAM-A.hmm database
#' @param output_file Output file path
#' @param evalue E-value threshold
#' @param num_threads Number of threads
#'
#' @return Path to HMMER results file or NULL if not available
#' @export
run_hmmscan_pfam <- function(query_fasta,
                              pfam_db,
                              output_file = NULL,
                              evalue = 1e-5,
                              num_threads = 4) {

  if (!check_hmmer_available()) {
    message("HMMER not found. Install with: conda install -c bioconda hmmer")
    message("Falling back to BLAST-based detection")
    return(NULL)
  }

  if (!file.exists(pfam_db)) {
    message("PFAM database not found: ", pfam_db)
    message("Download from: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz")
    return(NULL)
  }

  if (is.null(output_file)) {
    output_file <- tempfile(fileext = ".hmmer.txt")
  }

  message("Running HMMER scan against PFAM...")
  message("  Query: ", query_fasta)
  message("  Database: ", pfam_db)

  hmmer_cmd <- sprintf(
    "hmmscan --cpu %d -E %s --tblout %s %s %s",
    num_threads,
    evalue,
    output_file,
    shQuote(pfam_db),
    shQuote(query_fasta)
  )

  system(hmmer_cmd)

  if (file.exists(output_file) && file.info(output_file)$size > 0) {
    message("HMMER scan complete: ", output_file)
    return(output_file)
  } else {
    message("HMMER produced no results")
    return(NULL)
  }
}


#' Simple K-mer Based Comparison (Fallback)
#'
#' @keywords internal
compare_by_kmer <- function(query_data, rebase_data, id_col, k = 5, top_n = 3) {

  message("Running k-mer comparison...")

  get_kmers <- function(seq, k) {
    if (is.na(seq) || nchar(seq) < k) return(character())
    substring(seq, 1:(nchar(seq) - k + 1), k:nchar(seq))
  }

  kmer_similarity <- function(seq1, seq2, k) {
    kmers1 <- unique(get_kmers(toupper(seq1), k))
    kmers2 <- unique(get_kmers(toupper(seq2), k))
    if (length(kmers1) == 0 || length(kmers2) == 0) return(0)
    length(intersect(kmers1, kmers2)) / length(union(kmers1, kmers2))
  }

  results <- list()

  for (i in seq_len(nrow(query_data))) {
    query_id <- query_data[[id_col]][i]
    query_seq <- query_data$translation[i]

    if (is.na(query_seq)) next

    # Compare to each REBASE entry
    similarities <- sapply(seq_len(nrow(rebase_data)), function(j) {
      kmer_similarity(query_seq, rebase_data$sequence[j], k)
    })

    # Get top hits
    top_idx <- order(similarities, decreasing = TRUE)[1:min(top_n, length(similarities))]
    top_idx <- top_idx[similarities[top_idx] > 0.1]

    if (length(top_idx) > 0) {
      for (idx in top_idx) {
        results[[length(results) + 1]] <- dplyr::tibble(
          query_id = query_id,
          rebase_enzyme = rebase_data$enzyme_name[idx],
          identity = similarities[idx],
          pident = similarities[idx] * 100
        )
      }
    }

    if (i %% 10 == 0) message("  Processed ", i, "/", nrow(query_data))
  }

  if (length(results) == 0) return(dplyr::tibble())
  dplyr::bind_rows(results)
}

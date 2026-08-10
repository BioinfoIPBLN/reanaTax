##############################################################################
# exploreMetaTax — Kraken2 / KrakenUniq / Bracken / MetaPhlAn Report Viewer
#
# Interactive Shiny application for visualizing metagenomics taxonomy
# classification results from Kraken2, KrakenUniq, Bracken2, and MetaPhlAn.
#
# Features:
#   - Upload Kraken2 reports, KrakenUniq reports, Bracken per-sample files, combined Bracken, or MetaPhlAn profiles
#   - Optional metadata upload for group-based analyses
#   - In-app sample renaming
#   - Visualization tabs: Composition, Heatmap, Krona, Distribution,
#     Alpha Diversity, Rarefaction, PCA/Beta, LifemapR Tree, Data Table
#
# Dependencies:
#   shiny, ggplot2, plotly, dplyr, tidyr, DT, RColorBrewer,
#   vegan, LifemapR, scales, ggdendro, shinyjs, data.table, Rcpp
##############################################################################

# Local dev convenience: uncomment to make ?data= URL loads work outside
# ShinyProxy. In production, ShinyProxy sets SHINYPROXY_USERNAME for the
# authenticated B-Fabric user; an empty value means the public instance
# (no auth), where gstore loads are intentionally denied.
## Sys.setenv(SHINYPROXY_USERNAME = Sys.getenv("USER"))

# Increase upload limit to 500 MB (Shiny default is 5 MB)
options(shiny.maxRequestSize = 500 * 1024^2)

# Use temp dir for sass cache (avoids /root/.cache permission errors in Docker)
Sys.setenv(R_USER_CACHE_DIR = tempdir())

library(shiny)
library(shinydashboard)
library(fresh)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(DT)
library(RColorBrewer)
library(scales)
library(ggdendro)
library(shinyjs)
library(data.table)
library(Rcpp)


# ── FAST KRAKEN2 LINEAGE PARSER ──
# Two implementations of the same function. The C++ one is compiled at startup
# and used wherever a toolchain exists; the R one is the fallback for webR /
# shinylive builds, which run R in WebAssembly and have no C++ compiler.
#
# Both walk a Kraken2 report top-down carrying the current ancestor at each of
# the nine main ranks: a row at rank X sets the ancestor for X and clears every
# rank below it, and every row is then stamped with the ancestor state as it
# stands. Sub-ranks ("S1", "D1", ...) match no main rank, so they neither set
# nor clear anything and simply inherit the state - which is what puts a strain
# under its species.
#
# The R version reaches that result without a row loop, which matters in webR:
# for each rank it finds the most recent row that SET that rank and the most
# recent row that CLEARED it, and keeps the value only when the set is the more
# recent of the two.
KRAKEN_MAIN_RANKS <- c("R", "D", "K", "P", "C", "O", "F", "G", "S")
KRAKEN_RANK_LABELS <- c("Root", "Domain", "Kingdom", "Phylum", "Class",
                        "Order", "Family", "Genus", "Species")

build_kraken_lineage_r <- function(rank, name) {
  rank <- as.character(rank)
  name <- as.character(name)
  n <- length(rank)
  if (n == 0L) {
    empty <- rep(list(character(0)), length(KRAKEN_RANK_LABELS))
    names(empty) <- KRAKEN_RANK_LABELS
    return(empty)
  }
  idx <- match(rank, KRAKEN_MAIN_RANKS)   # NA for sub-ranks such as "S1"
  pos <- seq_len(n)
  out <- lapply(seq_along(KRAKEN_MAIN_RANKS), function(j) {
    last_set   <- cummax(ifelse(!is.na(idx) & idx == j, pos, 0L))
    last_clear <- cummax(ifelse(!is.na(idx) & idx <  j, pos, 0L))
    keep <- last_set > 0L & last_set > last_clear
    ifelse(keep, name[pmax(last_set, 1L)], NA_character_)
  })
  names(out) <- KRAKEN_RANK_LABELS
  out
}

# cppFunction() shells out to a compiler, so it throws under webR. Failing over
# to the R implementation keeps a shinylive build working instead of dying at
# startup; a normal deployment still gets the compiled version.
cpp_lineage_ok <- tryCatch({
  cppFunction('
List build_kraken_lineage_cpp(CharacterVector rank, CharacterVector name) {
    int n = rank.size();
    std::vector<std::string> main_ranks = {"R", "D", "K", "P", "C", "O", "F", "G", "S"};
    int num_ranks = main_ranks.size();
    
    std::vector<CharacterVector> out(num_ranks);
    for(int i = 0; i < num_ranks; i++) {
        out[i] = CharacterVector(n, NA_STRING);
    }
    
    std::vector<String> current_ancestors(num_ranks, NA_STRING);
    
    for (int i = 0; i < n; i++) {
        String r = rank[i];
        String nm = name[i];
        
        int r_idx = -1;
        for(int j = 0; j < num_ranks; j++) {
            if (r == main_ranks[j]) { r_idx = j; break; }
        }
        
        if (r_idx != -1) {
            current_ancestors[r_idx] = nm;
            for (int j = r_idx + 1; j < num_ranks; j++) {
                current_ancestors[j] = NA_STRING;
            }
        }
        
        for(int j = 0; j < num_ranks; j++) {
            out[j][i] = current_ancestors[j];
        }
    }
    
    return List::create(
        _["Root"] = out[0], _["Domain"] = out[1], _["Kingdom"] = out[2],
        _["Phylum"] = out[3], _["Class"] = out[4], _["Order"] = out[5],
        _["Family"] = out[6], _["Genus"] = out[7], _["Species"] = out[8]
    );
}
')
  TRUE
}, error = function(e) {
  message("[exploreMetaTax] No C++ toolchain (", conditionMessage(e),
          "); using the pure-R Kraken2 lineage parser.")
  FALSE
})

if (!isTRUE(cpp_lineage_ok)) {
  build_kraken_lineage_cpp <- build_kraken_lineage_r
}

# Optional: graceful fallback if not installed
HAS_VEGAN <- requireNamespace("vegan", quietly = TRUE)
if (HAS_VEGAN) library(vegan)

# sortable powers the drag-and-drop Visible/Hidden sample buckets on the Upload
# tab. If it's missing from the image we fall back to the plain checkbox list.
HAS_SORTABLE <- requireNamespace("sortable", quietly = TRUE)
if (HAS_SORTABLE) library(sortable)

# Turn on SortableJS "MultiDrag": click several sample chips to select them
# (each gets the default `sortable-selected` class — see the CSS below), then
# drag any one to move the whole selection at once, including across the
# Visible/Hidden buckets. The `sortable` package ships the plugin and mounts it
# itself, so setting the flag is all that's needed (same as exploreDE's buckets).
.multidrag_opts <- if (HAS_SORTABLE)
  sortable_options(multiDrag = TRUE, animation = 150) else NULL

HAS_LIFEMAPR <- requireNamespace("LifemapR", quietly = TRUE)
if (HAS_LIFEMAPR) library(LifemapR)

# taxplore (the Krona/sunburst tab) is GitHub-only and has no WebAssembly build
# in the webR CRAN mirror, so it cannot ship in a shinylive bundle. Hard-coding
# FALSE - rather than probing for it - is also what keeps shinylive's dependency
# scanner from trying to fetch it. The tab already handles the missing package
# and shows its install hint; the pipeline writes standalone Krona charts to
# `krona/` anyway.
HAS_TAXPLORE <- FALSE

# ═══════════════════════════════════════════════════════════════════════════
# PARSING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

#' Detect file format from content
#' @return One of "kraken2", "krakenuniq", "bracken", "combined_bracken", "metaphlan", "merged_metaphlan", or "unknown"
detect_format <- function(filepath) {
  first_lines <- readLines(filepath, n = 10, warn = FALSE)
  if (length(first_lines) == 0) return("unknown")

  # Skip comment lines (MetaPhlAn headers start with #)
  data_lines <- first_lines[!grepl("^#", first_lines)]
  comment_lines <- first_lines[grepl("^#", first_lines)]

  # KrakenUniq detection: comment lines contain "KrakenUniq"
  if (any(grepl("KrakenUniq", comment_lines, ignore.case = FALSE))) {
    return("krakenuniq")
  }

  # Also detect by KrakenUniq header: first data line starts with %<tab>
  if (length(data_lines) > 0 && grepl("^%\\t", data_lines[1])) {
    return("krakenuniq")
  }

  # MetaPhlAn detection: comment lines with #mpa or #clade_name, or
  # header/data containing pipe-separated taxonomy (k__|p__)
  if (any(grepl("#mpa_v|#clade_name|#SampleID", comment_lines, ignore.case = TRUE))) {
    # Check if merged (multiple sample columns) or single
    header_line <- data_lines[1]
    if (!is.null(header_line)) {
      header_fields <- strsplit(header_line, "\t")[[1]]
      # Merged MetaPhlAn: clade_name + NCBI_tax_id + multiple sample cols
      # Single MetaPhlAn: clade_name + NCBI_tax_id + relative_abundance [+ additional_species]
      if (length(header_fields) > 4 &&
          !any(grepl("relative_abundance|additional_species", header_fields, ignore.case = TRUE))) {
        return("merged_metaphlan")
      }
    }
    return("metaphlan")
  }

  # Also detect by data pattern: MetaPhlAn pipe-separated lineage, e.g.
  # "k__Bacteria|p__Firmicutes|c__Clostridia". We must require the PIPE here,
  # not just a bare "p__"/"s__" prefix: GTDB-based Kraken/Bracken databases
  # (e.g. HRGM) carry those same rank prefixes in single, non-piped taxon
  # names like "p__Firmicutes_A", and a prefix-only match misclassified those
  # Kraken2 reports as merged_metaphlan — which then parsed to nothing.
  if (length(data_lines) >= 2) {
    if (any(grepl("\\|[a-z]__", data_lines[1:min(3, length(data_lines))]))) {
      header_fields <- strsplit(data_lines[1], "\t")[[1]]
      if (length(header_fields) > 4 &&
          !any(grepl("relative_abundance|additional_species", header_fields, ignore.case = TRUE))) {
        return("merged_metaphlan")
      }
      return("metaphlan")
    }
  }

  first_line <- if (length(data_lines) > 0) data_lines[1] else first_lines[1]
  fields <- strsplit(first_line, "\t")[[1]]

  # Combined Bracken: header has "name" AND multiple sample columns with _num or _frac
  if (grepl("^name\t", first_line, ignore.case = TRUE) &&
      (sum(grepl("_num$|_frac$", fields)) >= 2)) {
    return("combined_bracken")
  }

  # Bracken per-sample: header starts with "name" and has ~7 columns
  if (grepl("^name\t", first_line, ignore.case = TRUE) &&
      length(fields) >= 6 && length(fields) <= 8) {
    return("bracken")
  }

  # Kraken2 report: no header, first field is numeric (percentage)
  first_field <- trimws(fields[1])
  if (!is.na(suppressWarnings(as.numeric(first_field)))) {
    if (length(fields) == 6 || length(fields) == 8) {
      return("kraken2")
    }
  }

  return("unknown")
}


#' Parse a Kraken2 report file (6 or 8 columns, no header)
parse_kraken2_report <- function(filepath, sample_name = NULL) {
  if (is.null(sample_name)) {
    sample_name <- sub("_report\\.txt$|_kraken2\\.txt$|\\.kreport2?$|\\.txt$", "",
                       basename(filepath), ignore.case = TRUE)
  }

  df <- tryCatch(
    as.data.frame(data.table::fread(filepath, sep = "\t", header = FALSE, quote = "",
                                    fill = TRUE, showProgress = FALSE, nThread = 4)),
    error = function(e) return(NULL)
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  ncols <- ncol(df)
  if (ncols == 6) {
    colnames(df) <- c("percent", "reads_clade", "reads_taxon", "rank", "taxid", "name")
  } else if (ncols == 8) {
    colnames(df) <- c("percent", "reads_clade", "reads_taxon",
                       "minimizers_clade", "minimizers_taxon", "rank", "taxid", "name")
  } else {
    # Try to read as 6-column anyway
    colnames(df)[1:min(ncols, 6)] <- c("percent", "reads_clade", "reads_taxon",
                                         "rank", "taxid", "name")[1:min(ncols, 6)]
  }

  # Compute indent-based depth from the raw name column (before trimws)
  # Kraken2 uses 2-space indentation per level
  raw_names <- df$name
  df$indent_depth <- (nchar(raw_names) - nchar(trimws(raw_names, which = "left"))) / 2

  df$name <- trimws(df$name)
  df$percent <- as.numeric(df$percent)
  df$reads_clade <- as.numeric(df$reads_clade)
  df$reads_taxon <- as.numeric(df$reads_taxon)
  df$taxid <- as.integer(df$taxid)
  df$Sample <- sample_name

  # Standardise rank codes
  df$rank <- trimws(df$rank)

  # ── Build lineage columns fast using C++ ──
  lineage_list <- build_kraken_lineage_cpp(df$rank, df$name)
  df <- cbind(df, as.data.frame(lineage_list, stringsAsFactors = FALSE))

  # Keep essential columns (include minimizers, indent_depth, and lineage if present)
  rank_labels <- c("Root", "Domain", "Kingdom", "Phylum", "Class",
                   "Order", "Family", "Genus", "Species")
  cols_keep <- intersect(
    c("percent", "reads_clade", "reads_taxon",
      "minimizers_clade", "minimizers_taxon",
      "rank", "taxid", "name", "Sample", "indent_depth",
      rank_labels),
    colnames(df)
  )
  df[, cols_keep, drop = FALSE]
}



#' Parse a KrakenUniq report file (9 columns, with comment header)
#' KrakenUniq reports have: %, reads, taxReads, kmers, dup, cov, taxID, rank, taxName
#' Comment lines start with # and the header line starts with %
parse_krakenuniq_report <- function(filepath, sample_name = NULL) {
  if (is.null(sample_name)) {
    sample_name <- sub("_report\\.txt$|_krakenuniq\\.txt$|\\.kreport2?$|\\.txt$", "",
                       basename(filepath), ignore.case = TRUE)
  }

  # Read all lines
  all_lines <- readLines(filepath, warn = FALSE)

  # Skip comment lines (start with #) and the header line (starts with %)
  comment_mask <- grepl("^#", all_lines)
  header_mask <- grepl("^%\\t", all_lines)
  data_mask <- !comment_mask & !header_mask & nchar(trimws(all_lines)) > 0
  data_lines <- all_lines[data_mask]

  if (length(data_lines) == 0) return(NULL)

  # Parse tab-separated data
  df <- tryCatch(
    as.data.frame(data.table::fread(text = paste(data_lines, collapse = "\n"),
                                    sep = "\t", header = FALSE, quote = "",
                                    fill = TRUE, showProgress = FALSE, nThread = 4)),
    error = function(e) return(NULL)
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # KrakenUniq reports have 9 columns:
  # %  reads  taxReads  kmers  dup  cov  taxID  rank  taxName
  if (ncol(df) < 9) return(NULL)
  colnames(df)[1:9] <- c("percent", "reads_clade", "reads_taxon",
                          "kmers", "dup", "cov", "taxid", "rank_full", "name")

  # Compute indent-based depth from the raw name column (before trimws)
  # KrakenUniq uses 2-space indentation per level like Kraken2
  raw_names <- df$name
  df$indent_depth <- (nchar(raw_names) - nchar(trimws(raw_names, which = "left"))) / 2

  df$name    <- trimws(df$name)
  df$percent <- as.numeric(df$percent)
  df$reads_clade <- as.numeric(df$reads_clade)
  df$reads_taxon <- as.numeric(df$reads_taxon)
  df$kmers   <- as.numeric(df$kmers)
  df$dup     <- as.numeric(df$dup)
  df$cov     <- as.numeric(df$cov)
  df$taxid   <- as.integer(df$taxid)
  df$Sample  <- sample_name

  # Map KrakenUniq full rank names to standard single-letter codes
  krakenuniq_rank_map <- c(
    "superkingdom" = "D", "kingdom" = "K", "phylum" = "P",
    "subphylum" = "P1", "superclass" = "C1", "class" = "C",
    "subclass" = "C2", "superorder" = "O1", "order" = "O",
    "suborder" = "O2", "infraorder" = "O3", "parvorder" = "O4",
    "superfamily" = "F1", "family" = "F", "subfamily" = "F2",
    "genus" = "G", "species" = "S", "subspecies" = "S1",
    "species group" = "G1", "species subgroup" = "G2",
    "tribe" = "F3", "subtribe" = "F4",
    "forma" = "S2", "varietas" = "S3"
  )

  # Convert rank names to standard codes
  rank_lower <- tolower(trimws(df$rank_full))
  df$rank <- ifelse(
    rank_lower %in% names(krakenuniq_rank_map),
    krakenuniq_rank_map[rank_lower],
    ifelse(df$name == "unclassified", "U",
           ifelse(df$name == "root", "R", "U1"))
  )

  # ── Build lineage columns fast using C++ ──
  lineage_list <- build_kraken_lineage_cpp(df$rank, df$name)
  df <- cbind(df, as.data.frame(lineage_list, stringsAsFactors = FALSE))

  # Keep essential columns (include KrakenUniq-specific cols, indent_depth, and lineage)
  rank_labels <- c("Root", "Domain", "Kingdom", "Phylum", "Class",
                   "Order", "Family", "Genus", "Species")
  cols_keep <- intersect(
    c("percent", "reads_clade", "reads_taxon",
      "kmers", "dup", "cov",
      "rank", "taxid", "name", "Sample", "indent_depth",
      rank_labels),
    colnames(df)
  )
  df[, cols_keep, drop = FALSE]
}


#' Parse a Bracken per-sample output file (7 columns, with header)
parse_bracken_output <- function(filepath, sample_name = NULL) {
  if (is.null(sample_name)) {
    sample_name <- sub("_bracken.*$|\\.bracken$|\\.txt$", "",
                       basename(filepath), ignore.case = TRUE)
  }

  df <- tryCatch(
    as.data.frame(data.table::fread(filepath, sep = "\t", header = TRUE, quote = "",
                                    fill = TRUE, showProgress = FALSE, nThread = 4,
                                    check.names = FALSE)),
    error = function(e) return(NULL)
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Standardise column names
  col_lower <- tolower(colnames(df))
  colnames(df) <- col_lower

  # Map to our standard schema
  name_col <- intersect(c("name"), col_lower)[1]
  taxid_col <- intersect(c("taxonomy_id", "taxid", "tax_id"), col_lower)[1]
  rank_col <- intersect(c("taxonomy_lvl", "level", "rank"), col_lower)[1]
  reads_col <- intersect(c("new_est_reads", "est_reads", "kraken_assigned_reads"), col_lower)[1]
  frac_col <- intersect(c("fraction_total_reads", "fraction"), col_lower)[1]

  result <- data.frame(
    name = if (!is.na(name_col)) df[[name_col]] else "unknown",
    taxid = if (!is.na(taxid_col)) as.integer(df[[taxid_col]]) else NA_integer_,
    rank = if (!is.na(rank_col)) df[[rank_col]] else "S",
    reads_clade = if (!is.na(reads_col)) as.numeric(df[[reads_col]]) else 0,
    reads_taxon = if (!is.na(reads_col)) as.numeric(df[[reads_col]]) else 0,
    percent = if (!is.na(frac_col)) as.numeric(df[[frac_col]]) * 100 else 0,
    Sample = sample_name,
    stringsAsFactors = FALSE
  )

  result
}


#' Parse a combined Bracken table (from combine_bracken_outputs.py)
parse_combined_bracken <- function(filepath) {
  df <- tryCatch(
    as.data.frame(data.table::fread(filepath, sep = "\t", header = TRUE, quote = "",
                                    fill = TRUE, showProgress = FALSE, nThread = 4,
                                    check.names = FALSE)),
    error = function(e) return(NULL)
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  col_lower <- tolower(colnames(df))

  # Identify taxonomy columns
  name_col <- which(col_lower == "name")[1]
  taxid_col <- which(col_lower %in% c("taxonomy_id", "taxid", "tax_id"))[1]
  rank_col <- which(col_lower %in% c("taxonomy_lvl", "level", "rank"))[1]

  # Identify sample columns: columns ending in _num or _frac
  num_cols <- grep("_num$", colnames(df), value = TRUE)
  frac_cols <- grep("_frac$", colnames(df), value = TRUE)

  # If no _num/_frac pattern, treat remaining numeric columns as sample columns
  meta_cols_idx <- na.omit(c(name_col, taxid_col, rank_col))
  if (length(num_cols) == 0 && length(frac_cols) == 0) {
    remaining <- setdiff(seq_len(ncol(df)), meta_cols_idx)
    num_cols <- colnames(df)[remaining[sapply(df[remaining], is.numeric)]]
  }

  results <- list()

  if (length(num_cols) > 0) {
    for (col in num_cols) {
      sample_name <- sub("_num$", "", col)
      result <- data.frame(
        name = if (!is.na(name_col)) df[[name_col]] else "unknown",
        taxid = if (!is.na(taxid_col)) as.integer(df[[taxid_col]]) else NA_integer_,
        rank = if (!is.na(rank_col)) df[[rank_col]] else "S",
        reads_clade = as.numeric(df[[col]]),
        reads_taxon = as.numeric(df[[col]]),
        Sample = sample_name,
        stringsAsFactors = FALSE
      )
      # Add fraction if available
      frac_name <- paste0(sample_name, "_frac")
      if (frac_name %in% colnames(df)) {
        result$percent <- as.numeric(df[[frac_name]]) * 100
      } else {
        total <- sum(as.numeric(df[[col]]), na.rm = TRUE)
        result$percent <- if (total > 0) (result$reads_clade / total) * 100 else 0
      }
      results[[length(results) + 1]] <- result
    }
  }

  if (length(results) > 0) {
    return(bind_rows(results))
  }

  return(NULL)
}


#' Parse a MetaPhlAn profile (single-sample or merged multi-sample)
#' Handles MetaPhlAn 3/4 output with clade_name|NCBI_tax_id|relative_abundance columns
#' Also handles merged tables from merge_metaphlan_tables.py
parse_metaphlan_output <- function(filepath, sample_name = NULL, merged = FALSE) {
  if (is.null(sample_name)) {
    sample_name <- sub("_profile.*$|_metaphlan.*$|\\.txt$|\\.tsv$", "",
                       basename(filepath), ignore.case = TRUE)
  }

  # Read all lines, identify comments and header
  all_lines <- readLines(filepath, warn = FALSE)
  comment_mask <- grepl("^#", all_lines)

  # Find the header: prefer #clade_name over #SampleID (MetaPhlAn files may
  # have #SampleID before #clade_name, but #clade_name is the real data header)
  header_line <- NULL
  header_idx <- NULL

  # Pass 1: look for #clade_name (the actual data header)
  for (i in seq_along(all_lines)) {
    if (grepl("^#?clade_name", all_lines[i], ignore.case = TRUE)) {
      header_line <- sub("^#", "", all_lines[i])
      header_idx <- i
      break
    }
  }

  # Pass 2: if no #clade_name, try #SampleID (merged MetaPhlAn tables)
  if (is.null(header_line)) {
    for (i in seq_along(all_lines)) {
      if (grepl("^#SampleID", all_lines[i], ignore.case = TRUE)) {
        header_line <- sub("^#", "", all_lines[i])
        header_idx <- i
        break
      }
    }
  }

  # Pass 3: fallback — use first non-comment line if it contains known column names
  if (is.null(header_line)) {
    first_data <- which(!comment_mask)[1]
    if (!is.null(first_data)) {
      candidate <- all_lines[first_data]
      if (grepl("clade_name|NCBI_tax_id|relative_abundance", candidate, ignore.case = TRUE)) {
        header_line <- candidate
        header_idx <- first_data
      }
    }
  }

  # Read data
  if (!is.null(header_idx)) {
    data_lines <- all_lines[(header_idx + 1):length(all_lines)]
    data_lines <- data_lines[data_lines != "" & !grepl("^#", data_lines)]
    if (length(data_lines) == 0) return(NULL)

    header_fields <- strsplit(header_line, "\t")[[1]]
    data_split <- strsplit(data_lines, "\t")

    # Build data.frame
    df <- do.call(rbind, lapply(data_split, function(x) {
      length(x) <- length(header_fields)  # pad with NA if shorter
      x
    }))
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    colnames(df) <- header_fields
  } else {
    # Fallback: try read.table skipping comments
    df <- tryCatch(
      as.data.frame(data.table::fread(filepath, sep = "\t", header = TRUE, quote = "",
                                      fill = TRUE, showProgress = FALSE, nThread = 4,
                                      check.names = FALSE, skip = "#")),
      error = function(e) return(NULL)
    )
  }

  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Identify columns
  clade_col <- grep("clade_name", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  taxid_col <- grep("NCBI_tax_id|tax_id|taxid", colnames(df), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(clade_col)) {
    # Maybe first column is the clade
    clade_col <- colnames(df)[1]
  }

  # MetaPhlAn rank mapping from prefix to standard code
  mpa_rank_map <- c("k__" = "D", "p__" = "P", "c__" = "C", "o__" = "O",
                    "f__" = "F", "g__" = "G", "s__" = "S", "t__" = "T")

  # Extract the deepest taxon name and rank from clade_name
  extract_taxon_info <- function(clade_str) {
    parts <- strsplit(clade_str, "\\|")[[1]]
    deepest <- parts[length(parts)]

    # Get rank prefix
    prefix <- sub("^([a-z]__)(.*)", "\\1", deepest)
    taxon_name <- sub("^[a-z]__", "", deepest)
    taxon_name <- gsub("_", " ", taxon_name)

    # Normalise UNCLASSIFIED → unclassified for consistency with filter logic
    if (toupper(taxon_name) == "UNCLASSIFIED") taxon_name <- "unclassified"

    rank_code <- if (prefix %in% names(mpa_rank_map)) mpa_rank_map[prefix] else "U"
    list(name = taxon_name, rank = unname(rank_code))
  }

  # Extract last taxid from pipe-separated NCBI_tax_id column
  extract_last_taxid <- function(taxid_str) {
    if (is.na(taxid_str) || taxid_str == "") return(NA_integer_)
    parts <- strsplit(as.character(taxid_str), "\\|")[[1]]
    as.integer(parts[length(parts)])
  }

  # Determine if this is a merged table (multiple sample columns)
  known_meta_cols <- c(clade_col, taxid_col, "additional_species")
  known_meta_cols <- known_meta_cols[!is.na(known_meta_cols)]
  abundance_col <- grep("relative_abundance", colnames(df), ignore.case = TRUE, value = TRUE)[1]

  if (!merged && !is.na(abundance_col)) {
    # ── Single-sample MetaPhlAn profile ──
    taxon_info <- lapply(df[[clade_col]], extract_taxon_info)
    result <- data.frame(
      name = sapply(taxon_info, `[[`, "name"),
      taxid = if (!is.na(taxid_col)) sapply(df[[taxid_col]], extract_last_taxid) else NA_integer_,
      rank = sapply(taxon_info, `[[`, "rank"),
      reads_clade = 0,
      reads_taxon = 0,
      percent = as.numeric(df[[abundance_col]]),
      Sample = sample_name,
      stringsAsFactors = FALSE
    )
    return(result)

  } else {
    # ── Merged MetaPhlAn table (multiple samples as columns) ──
    sample_cols <- setdiff(colnames(df), known_meta_cols)
    # Remove non-numeric columns
    sample_cols <- sample_cols[sapply(df[sample_cols], function(x) {
      !all(is.na(suppressWarnings(as.numeric(x))))
    })]

    if (length(sample_cols) == 0) return(NULL)

    results <- list()
    taxon_info <- lapply(df[[clade_col]], extract_taxon_info)

    for (scol in sample_cols) {
      result <- data.frame(
        name = sapply(taxon_info, `[[`, "name"),
        taxid = if (!is.na(taxid_col)) sapply(df[[taxid_col]], extract_last_taxid) else NA_integer_,
        rank = sapply(taxon_info, `[[`, "rank"),
        reads_clade = 0,
        reads_taxon = 0,
        percent = as.numeric(df[[scol]]),
        Sample = scol,
        stringsAsFactors = FALSE
      )
      results[[length(results) + 1]] <- result
    }
    return(bind_rows(results))
  }
}


# ═══════════════════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

#' Get standard rank choices
rank_choices <- c("Domain" = "D", "Kingdom" = "K", "Phylum" = "P",
                  "Class" = "C", "Order" = "O", "Family" = "F",
                  "Genus" = "G", "Species" = "S")

#' Palette choices
palette_choices <- c("Set3", "Paired", "Set1", "Set2", "Dark2",
                     "Accent", "Pastel1", "Pastel2", "Spectral", "RdYlBu")

#' Generate N distinct colours from a palette
get_palette_colors <- function(n, palette_name = "Set3") {
  max_pal <- brewer.pal.info[palette_name, "maxcolors"]
  if (n <= max_pal) {
    brewer.pal(max(3, n), palette_name)[1:n]
  } else {
    colorRampPalette(brewer.pal(max_pal, palette_name))(n)
  }
}

#' Compute diversity indices manually (fallback if vegan not available)
compute_diversity <- function(x) {
  x <- x[x > 0]
  total <- sum(x)
  if (total == 0) return(data.frame(Shannon = 0, Simpson = 0, InvSimpson = 0))
  p <- x / total
  shannon <- -sum(p * log(p))
  simpson <- 1 - sum(p^2)
  inv_simpson <- 1 / sum(p^2)
  data.frame(Shannon = shannon, Simpson = simpson, InvSimpson = inv_simpson)
}

# ─── Group-comparison statistics (non-parametric) ──────────────────────────
# Abundance data is typically non-normal and heavily skewed; a rank-based
# test is a safer default than a t-test/ANOVA. Two groups → Wilcoxon
# rank-sum; ≥3 → Kruskal–Wallis. Any group with n < 2 is skipped (the
# strict minimum wilcox.test / kruskal.test accept). When 2 ≤ min(n) < 3,
# the test still runs but the returned note flags the low-n caveat: the
# smallest achievable p is capped (e.g. p ≥ ~0.33 for n=2 vs n=2), so a
# non-significant result on very small samples means very little.
# Returns a one-row data.frame (or NULL when the test can't be run).
compute_group_stat <- function(values, groups) {
  ok <- !is.na(values) & !is.na(groups) & nzchar(as.character(groups))
  v <- values[ok]; g <- as.character(groups[ok])
  if (length(v) == 0) return(NULL)
  ns <- table(g)
  levs <- names(ns)
  if (length(levs) < 2 || any(ns < 2)) {
    return(data.frame(
      test = NA_character_, statistic = NA_real_,
      p = NA_real_, note = "n < 2 in ≥1 group",
      stringsAsFactors = FALSE
    ))
  }
  low_n_note <- if (any(ns < 3)) {
    sprintf("low n (min group size = %d): min achievable p is limited",
            min(ns))
  } else NA_character_
  if (length(levs) == 2) {
    t <- tryCatch(wilcox.test(v ~ factor(g), exact = FALSE),
                  error = function(e) NULL)
    if (is.null(t)) return(NULL)
    data.frame(test = "Wilcoxon rank-sum",
               statistic = unname(t$statistic),
               p = t$p.value, note = low_n_note,
               stringsAsFactors = FALSE)
  } else {
    t <- tryCatch(kruskal.test(v ~ factor(g)), error = function(e) NULL)
    if (is.null(t)) return(NULL)
    data.frame(test = "Kruskal-Wallis",
               statistic = unname(t$statistic),
               p = t$p.value, note = low_n_note,
               stringsAsFactors = FALSE)
  }
}

# Pairwise Wilcoxon with BH FDR correction, returned as a long data.frame.
# Requires ≥3 groups (with n ≥ 2 each after filtering); NULL otherwise.
pairwise_wilcox_bh <- function(values, groups) {
  ok <- !is.na(values) & !is.na(groups) & nzchar(as.character(groups))
  v <- values[ok]; g <- as.character(groups[ok])
  ns <- table(g)
  g_ok <- names(ns)[ns >= 2]
  if (length(g_ok) < 3) return(NULL)
  keep <- g %in% g_ok
  v <- v[keep]; g <- g[keep]
  res <- tryCatch(
    suppressWarnings(pairwise.wilcox.test(v, g, p.adjust.method = "BH",
                                          exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(res) || is.null(res$p.value)) return(NULL)
  m <- res$p.value
  out <- expand.grid(group1 = rownames(m), group2 = colnames(m),
                     stringsAsFactors = FALSE)
  out$p_adj <- as.vector(m)
  out <- out[!is.na(out$p_adj), , drop = FALSE]
  out[order(out$p_adj), , drop = FALSE]
}

# Per-taxon group-comparison table over a set of taxa, with BH-adjusted
# p-values across the taxa panel. `df` must be a long data.frame with
# columns Sample, name, percent (or the value column requested), plus the
# grouping column. Returns a data.frame ordered by adjusted p-value (or
# NULL if no comparable groups).
per_taxon_group_table <- function(df, taxa, group_col, value_col = "percent") {
  if (is.null(df) || nrow(df) == 0 || is.null(group_col) ||
      !(group_col %in% colnames(df))) return(NULL)
  rows <- lapply(taxa, function(tx) {
    sub <- df[df$name == tx, , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)
    res <- compute_group_stat(sub[[value_col]], sub[[group_col]])
    if (is.null(res)) return(NULL)
    medians <- tapply(sub[[value_col]], sub[[group_col]],
                      function(x) median(x, na.rm = TRUE))
    med_str <- paste(sprintf("%s: %.3g", names(medians), medians),
                     collapse = "; ")
    cbind(data.frame(taxon = tx, stringsAsFactors = FALSE),
          res,
          data.frame(medians = med_str, stringsAsFactors = FALSE))
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)
  out$p_adj <- p.adjust(out$p, method = "BH")
  out <- out[, c("taxon", "test", "statistic", "p", "p_adj", "medians", "note"),
             drop = FALSE]
  out[order(out$p_adj, out$p), , drop = FALSE]
}

#' Regex for filenames typical of a given metagenomics format
format_filename_pattern <- function(format) {
  switch(format,
    "kraken2"          = "\\.report\\.txt$|_report\\.txt$|\\.kreport2?$|_kraken2\\.txt$",
    "krakenuniq"       = "_krakenuniq.*\\.txt$|\\.krakenuniq$|_krakenuniq_report\\.txt$",
    "bracken"          = "\\.bracken$|_bracken[^/]*\\.(txt|tsv)$",
    "combined_bracken" = "combined.*bracken.*\\.(txt|tsv)$",
    "metaphlan"        = "_profile\\.txt$|_metaphlan.*\\.(txt|tsv)$|_metaphlan_bugs_list\\.tsv$",
    "merged_metaphlan" = "merged_abundance.*\\.(txt|tsv)$|merged_metaphlan.*\\.(txt|tsv)$",
    "humann"           = "pathabund(ance)?.*\\.tsv$|reactions?.*\\.tsv$|gene[_-]?famil(ies|y).*\\.tsv$|(^|[_-])ko(_cpm)?\\.tsv$|kegg[_-]?kos?.*\\.tsv$",
    NULL
  )
}

#' Classify a HUMAnN per-sample TSV filename into a table slot
#' (pathways / reactions / kegg_kos / gene_families) or NA if not recognised.
#' The SUSHI HUMAnN app writes *_4_pathabundance.tsv / *_4_reactions.tsv /
#' *_4_ko_cpm.tsv / *_4_genefamilies_cpm.tsv, while stock HUMAnN drops
#' the "_4_" prefix.
detect_humann_slot <- function(name) {
  low <- tolower(name)
  if (grepl("pathabund", low))                     return("pathways")
  if (grepl("gene[_-]?famil", low))                return("gene_families")
  if (grepl("reaction", low))                      return("reactions")
  if (grepl("(^|[_-])ko(_cpm)?\\.tsv$|kegg", low)) return("kegg_kos")
  NA_character_
}

#' Strip the HUMAnN table suffix from a filename to recover the sample name.
strip_humann_suffix <- function(name, slot) {
  sfx <- switch(slot,
    pathways      = "(_4)?[_-]?pathabund(ance)?(_[a-z]+)?\\.tsv$",
    reactions     = "(_4)?[_-]?reactions?(_[a-z]+)?\\.tsv$",
    kegg_kos      = "(_4)?[_-]?(ko(_cpm)?|kegg[_-]?kos?)\\.tsv$",
    gene_families = "(_4)?[_-]?gene[_-]?famil(ies|y)(_[a-z]+)?\\.tsv$"
  )
  s <- sub(sfx, "", name, ignore.case = TRUE, perl = TRUE)
  sub("[_-]+$", "", s)
}

#' Parse a user-uploaded set of HUMAnN per-sample TSVs into the same
#' shape produced by load_humann_tables() for the SUSHI dataset.tsv path.
#' Returns a named list of joined wide matrices keyed by slot
#' (pathways/reactions/kegg_kos/gene_families) or NULL if nothing recognised.
parse_humann_uploaded_files <- function(paths, names) {
  if (length(paths) == 0) return(NULL)
  slots   <- vapply(names, detect_humann_slot, character(1))
  keep    <- !is.na(slots)
  if (!any(keep)) return(NULL)
  paths   <- paths[keep]; names <- names[keep]; slots <- slots[keep]
  samples <- mapply(strip_humann_suffix, names, slots,
                    USE.NAMES = FALSE, SIMPLIFY = TRUE)

  out <- list()
  for (slot in unique(slots)) {
    idx <- which(slots == slot)
    per_sample <- mapply(read_humann_tsv, paths[idx], samples[idx],
                         SIMPLIFY = FALSE)
    per_sample <- per_sample[!vapply(per_sample, is.null, logical(1))]
    if (length(per_sample) == 0) next
    joined <- Reduce(function(a, b) merge(a, b, by = "Feature", all = TRUE),
                     per_sample)
    for (j in setdiff(colnames(joined), "Feature")) {
      v <- joined[[j]]; v[is.na(v)] <- 0
      joined[[j]] <- v
    }
    out[[slot]] <- joined
  }
  canon <- c("pathways", "reactions", "kegg_kos", "gene_families")
  out[intersect(canon, names(out))]
}

#' Logical vector — which `names` match a given format's filename patterns
filter_names_for_format <- function(names, format) {
  pattern <- format_filename_pattern(format)
  if (is.null(pattern) || length(names) == 0) return(rep(FALSE, length(names)))
  grepl(pattern, names, ignore.case = TRUE)
}

#' Server-side: absolute paths in `folder` matching the format's patterns
list_files_for_format <- function(folder, format) {
  if (!isTRUE(nzchar(folder)) || !dir.exists(folder)) return(character(0))
  pattern <- format_filename_pattern(format)
  if (is.null(pattern)) return(character(0))
  list.files(folder, pattern = pattern, ignore.case = TRUE,
             full.names = TRUE, recursive = FALSE)
}

#' Server-side: list any plausible report files in `folder` (any common
#' metagenomics extension) for downstream per-file format auto-detection.
list_candidate_report_files <- function(folder) {
  if (!isTRUE(nzchar(folder)) || !dir.exists(folder)) return(character(0))
  list.files(folder,
             pattern = "\\.(txt|tsv|kreport|kreport2|bracken|profile)$",
             ignore.case = TRUE, full.names = TRUE, recursive = FALSE)
}

#' Resolve a URL-supplied path against known gstore project roots.
#' Mirrors the exploreDE convention: a relative path like
#' "p1234/Kraken_2025-.../" is resolved under /srv/gstore/projects first,
#' then under the course gstore root. An absolute path is honored as-is.
#' Returns the first matching existing directory, or NULL if none.
GSTORE_ROOTS <- c(
  "/srv/gstore/projects",
  "/srv/GT/analysis/course_sushi/public/gstore/projects"
)
resolve_gstore_folder <- function(path) {
  if (!isTRUE(nzchar(path))) return(NULL)
  candidates <- unique(c(path, file.path(GSTORE_ROOTS, path)))
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) NULL else hit[1]
}

#' Resolve a path written in a SUSHI dataset.tsv to an absolute file path.
#' Such paths are typically relative to a gstore root (e.g.
#' "p41135/o41225_.../foo.report.txt"). Absolute paths are honored as-is.
#' Returns NA_character_ if no candidate exists on disk.
resolve_gstore_file <- function(rel_or_abs) {
  if (is.na(rel_or_abs) || !nzchar(rel_or_abs)) return(NA_character_)
  candidates <- unique(c(rel_or_abs, file.path(GSTORE_ROOTS, rel_or_abs)))
  hit <- candidates[file.exists(candidates) & !file.info(candidates)$isdir]
  if (length(hit) == 0) NA_character_ else hit[1]
}

#' Read a SUSHI dataset.tsv and return the report-file manifest + metadata.
#'
#' SUSHI dataset.tsv conventions (verified against the kraken2 SUSHI app):
#'   - Each row is one sample.
#'   - `Name` column → human-readable sample name.
#'   - One column ends in `Report [File]` (Kraken/Bracken/KrakenUniq, e.g.
#'     "KrakenReport [File]") OR `Profile [File]` (MetaPhlAn, e.g.
#'     "MetaPhlAnProfile [File]" — "profile" is the canonical upstream term)
#'     and holds the path to the per-sample report, relative to a gstore root.
#'   - `[File]` / `[Link]` columns are inputs/outputs; everything else
#'     (Condition, Sample Id, Order Id, …) is sample metadata that the
#'     downstream UI joins on `Sample`.
#'
#' Returns list(report_paths, sample_names, metadata, missing) where
#'   - report_paths: resolved absolute paths (NA where the file is missing)
#'   - sample_names: aligned `Name` values
#'   - metadata:     data.frame keyed by `Sample` (= Name) with annotation cols
#'   - missing:      character vector of report values that didn't resolve
#'
#' Stops returning NULL only when the file is unreadable; downstream code is
#' responsible for hard-erroring on missing files / missing report column.
read_sushi_dataset_tsv <- function(dir_path) {
  tsv_path <- file.path(dir_path, "dataset.tsv")
  if (!file.exists(tsv_path)) return(NULL)
  df <- tryCatch(
    as.data.frame(data.table::fread(tsv_path, sep = "\t", header = TRUE,
                                    quote = "", fill = TRUE,
                                    showProgress = FALSE, nThread = 4,
                                    check.names = FALSE)),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  cn <- colnames(df)
  name_col <- cn[match("Name", cn)]
  # Pick the first column whose header matches "<something>Report [File]"
  # (KrakenReport, KrakenUniqReport, BrackenReport) OR
  # "<something>Profile [File]" (MetaPhlAnProfile — canonical MetaPhlAn term).
  report_col <- grep("(Report|Profile)\\s*\\[File\\]\\s*$", cn,
                     ignore.case = TRUE, value = TRUE)[1]

  # HUMAnN-mode detection: if the dataset.tsv carries per-sample functional
  # tables (any of PathAbundance / ReactionsCPM / KEGGKO_CPM / GeneFamiliesCPM
  # as a [File] column), classify as humann mode. PathAbundance is the
  # canonical anchor — required if any humann column is present.
  HUMANN_TABLES <- c(
    pathways      = "PathAbundance",
    reactions     = "ReactionsCPM",
    kegg_kos      = "KEGGKO_CPM",
    gene_families = "GeneFamiliesCPM"
  )
  humann_cols <- vapply(HUMANN_TABLES, function(stem) {
    hit <- grep(paste0("^", stem, "\\s*\\[File\\]\\s*$"), cn,
                ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) hit[1] else NA_character_
  }, FUN.VALUE = NA_character_)
  any_humann <- !all(is.na(humann_cols))
  # Mode: prefer the [File] column (per-sample text file with content
  # "full" or "translated"). Fall back to the [Characteristic] scheme
  # if an older upstream job wrote it inline.
  mode_file_col <- grep("^Mode\\s*\\[File\\]\\s*$", cn,
                        ignore.case = TRUE, value = TRUE)[1]
  mode_col      <- grep("^Mode\\s*\\[Characteristic\\]\\s*$", cn,
                        ignore.case = TRUE, value = TRUE)[1]
  app_mode <- if (any_humann) "humann" else "taxonomy"

  list(
    name_col      = name_col,
    report_col    = report_col,
    humann_cols   = humann_cols,
    mode_file_col = mode_file_col,
    mode_col      = mode_col,
    app_mode      = app_mode,
    raw           = df
  )
}

# ─── HUMAnN per-sample-TSV joiner ───────────────────────────────────────────
# HUMAnN writes one TSV per sample per table. The first column is the feature
# ID ("# Gene Family", "# Reaction", "# Pathway") followed by a single sample
# abundance column. We join all per-sample files for a given table by feature
# ID, using full_join semantics (NA for samples that don't carry the feature).
#
# Stratified rows (feature|species) are kept in the joined matrix so the
# full DT table can display them; the plotting layer filters them out for
# top-N / heatmap / PCoA.
read_humann_tsv <- function(path, sample_name, quote = "\"", colClasses = NULL,
                            nThread = 4L) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  # quote = "\"" (fread's default) is critical: `humann_rename_table -n
  # kegg-orthology` CSV-escapes KO names that contain a literal `"` (e.g.
  # K00984's "streptomycin 3\"-adenylyltransferase") inside a TSV — the
  # field arrives on disk as `"K00984: ...3""-..."` and must be dequoted
  # here or the surrounding `"` and doubled internal `""` propagate into
  # every display, filter, and CSV/TSV export downstream. Other HUMAnN
  # slot files never emit `"` so this doesn't affect them — the caller passes
  # quote = "" for the (large, quote-free) gene-family table, which skips
  # fread's quote handling, and colClasses to skip column-type detection.
  df <- tryCatch(
    as.data.frame(data.table::fread(path, sep = "\t", header = TRUE,
                                    quote = quote, check.names = FALSE,
                                    colClasses = colClasses,
                                    nThread = nThread,
                                    showProgress = FALSE)),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0 || ncol(df) < 2) return(NULL)
  # First column = feature id; second column = abundance for this sample.
  # The sample col name written by HUMAnN includes "_Abundance"; we replace
  # it with the SUSHI Name.
  colnames(df)[1] <- "Feature"
  colnames(df)[2] <- sample_name
  df[, 1:2, drop = FALSE]
}

#' Fast wide join of per-sample HUMAnN tables.
#'
#' Each element of `per_sample` is a 2-column data frame (Feature, <sample>)
#' as returned by read_humann_tsv(). Produces one wide frame: rows = the union
#' of all feature ids (first-seen order), cols = Feature + one numeric column
#' per sample, absent features filled with 0.
#'
#' Replaces a 43-way `Reduce(merge(..., all = TRUE))` fold. That fold re-sorted
#' a multi-million-row string key on every step (super-linear) and was the
#' dominant cost when loading gene-family tables (~1.9M rows x 44 samples).
#' Here we build the feature union once and fill a preallocated numeric matrix
#' by hashed match() — a single O(total rows) pass, no repeated sorting.
join_humann_wide <- function(per_sample) {
  if (length(per_sample) == 0) return(NULL)
  # Sample name is the 2nd column name of each per-sample frame (read_humann_tsv
  # set it to the SUSHI Name); mapply keys `per_sample` by path, not sample.
  sample_names <- vapply(per_sample, function(d) colnames(d)[2], character(1))
  master <- unique(unlist(lapply(per_sample, `[[`, "Feature"),
                          use.names = FALSE))
  mat <- matrix(0, nrow = length(master), ncol = length(per_sample),
                dimnames = list(NULL, sample_names))
  for (j in seq_along(per_sample)) {
    d    <- per_sample[[j]]
    vals <- as.numeric(d[[2]])
    vals[is.na(vals)] <- 0
    mat[match(d$Feature, master), j] <- vals
  }
  data.frame(Feature = master, mat, check.names = FALSE,
             stringsAsFactors = FALSE)
}

#' Load + join HUMAnN per-sample tables referenced by a parsed dataset.tsv.
#' Returns a named list of data frames keyed by humann_cols slot name
#' (pathways / reactions / kegg_kos / gene_families); each frame is wide
#' (rows = features, cols = Feature + one column per sample), NA-filled as 0.
#' Tables whose column is missing in the dataset.tsv get a NULL slot.
#'
#' `slots_wanted` restricts loading to a subset of slot names. The large
#' gene-family table is loaded lazily (only when its tab is opened) via this
#' argument, so the initial load reads only the small pathways/reactions/KO
#' tables. NULL loads every slot.
load_humann_tables <- function(parsed_ds, slots_wanted = NULL) {
  ds <- parsed_ds$raw
  if (is.null(ds) || nrow(ds) == 0) return(list())
  sample_names <- as.character(ds[[parsed_ds$name_col]])

  slot_names <- names(parsed_ds$humann_cols)
  if (!is.null(slots_wanted)) slot_names <- intersect(slot_names, slots_wanted)

  out <- list()
  for (slot in slot_names) {
    col <- parsed_ds$humann_cols[[slot]]
    if (is.na(col)) { out[[slot]] <- NULL; next }
    rels <- as.character(ds[[col]])
    paths <- vapply(rels, resolve_gstore_file,
                    FUN.VALUE = NA_character_, USE.NAMES = FALSE)
    # gene_families is huge and never contains quotes; disable quote handling,
    # fix column types, and give fread several threads so the parse (the
    # first-open cost) is parallelised with the most threads. The other (small)
    # slots keep quote handling — kegg_kos genuinely needs it (KO names with ")
    # — and read with 4 threads.
    is_gf    <- identical(slot, "gene_families")
    rd_quote <- if (is_gf) "" else "\""
    rd_cc    <- if (is_gf) c("character", "numeric") else NULL
    rd_thr   <- if (is_gf) HUMANN_GF_READ_THREADS else 4L
    per_sample <- mapply(read_humann_tsv, paths, sample_names,
                         MoreArgs = list(quote = rd_quote, colClasses = rd_cc,
                                         nThread = rd_thr),
                         SIMPLIFY = FALSE)
    per_sample <- per_sample[!vapply(per_sample, is.null, logical(1))]
    if (length(per_sample) == 0) { out[[slot]] <- NULL; next }
    out[[slot]] <- join_humann_wide(per_sample)
  }
  out
}

# ─── Async helper: run a function in a background R process, return a promise ─
# Backs the Shiny ExtendedTask that loads the large gene-family table off the
# main session, so the UI stays responsive while it reads. callr (already a
# dependency) runs `func` in a separate R process; we poll it on the `later`
# loop (ships with Shiny) and resolve/reject the promise when it exits. No new
# package is needed. `func` must be self-contained (a fresh R session): pass
# everything it needs via `args`; it may use `pkg::fn` for installed packages.
# The child returns its value through callr's own RDS hand-off — that is process
# IPC, not an on-disk cache.
callr_promise <- function(func, args = list()) {
  promises::promise(function(resolve, reject) {
    proc <- tryCatch(callr::r_bg(func = func, args = args),
                     error = function(e) { reject(e); NULL })
    if (is.null(proc)) return(invisible(NULL))
    poll <- function() {
      if (isTRUE(proc$is_alive())) {
        later::later(poll, delay = 0.3)
      } else {
        res <- tryCatch(proc$get_result(), error = function(e) e)
        if (inherits(res, "condition")) reject(res) else resolve(res)
      }
    }
    poll()
  })
}

# Background worker for the gene-family read (runs in a fresh callr R session).
# Kept top-level so callr serialises it with a clean (global) environment. The
# child has none of the app's globals, so the reader/joiner are passed in as
# arguments (read_fn = read_humann_tsv, join_fn = join_humann_wide — both are
# self-contained, using only base R + data.table::fread).
gf_read_join_bg <- function(paths, sample_names, read_fn, join_fn,
                            quote, colClasses, nThread) {
  per_sample <- mapply(read_fn, paths, sample_names,
                       MoreArgs = list(quote = quote, colClasses = colClasses,
                                       nThread = nThread),
                       SIMPLIFY = FALSE)
  per_sample <- per_sample[!vapply(per_sample, is.null, logical(1))]
  if (length(per_sample) == 0) return(NULL)
  join_fn(per_sample)
}

# ─── HUMAnN plotting helpers (inlined — exploreMetaTax has no ezRun dep) ────

# Max stratified rows fed to the Stratified-tab plot/table when the query is
# empty (or very broad). Rows are ranked by summed CPM and a banner warns when
# the cap is hit; a specific organism query normally narrows well below it.
# Without this cap an unfiltered gene-family view reshapes ~1M rows x N samples.
HUMANN_STRAT_ROW_CAP <- 2000L

# fread threads for the one large read (the gene-family table: ~1.9M rows x N
# samples, read lazily on first Gene-families/Stratified tab open). Parsing the
# gene-family table IS the first-open cost, so it gets the most threads; the
# small taxonomy / pathway / reaction / KO reads use 4 (see read_humann_tsv and
# the taxonomy readers). fread caps the request to the cores actually available,
# so these are safe even if a container is allocated fewer. Tune here if needed.
HUMANN_GF_READ_THREADS <- 24L

# Community-level rows only (drop "feature|species" stratified rows). HUMAnN
# uses "|" as the SGB separator inside the Feature ID.
humann_community_only <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df[!grepl("|", df$Feature, fixed = TRUE), , drop = FALSE]
}

# Drop HUMAnN "special" features (UNMAPPED, READS_UNMAPPED, UNINTEGRATED,
# UNGROUPED) — they aren't biological signal and dominate top-N otherwise.
humann_drop_specials <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  specials <- c("UNMAPPED", "READS_UNMAPPED", "UNINTEGRATED", "UNGROUPED")
  feat_root <- sub("\\:.*$", "", df$Feature)  # strip "PWY-NNN: name" suffix
  df[!(df$Feature %in% specials | feat_root %in% specials), , drop = FALSE]
}

# Top-N features by mean abundance across samples
humann_top_n <- function(df, n = 20) {
  if (is.null(df) || nrow(df) == 0) return(df)
  sample_cols <- setdiff(colnames(df), "Feature")
  means <- rowMeans(df[, sample_cols, drop = FALSE], na.rm = TRUE)
  ord <- order(means, decreasing = TRUE)
  df[ord[seq_len(min(n, length(ord)))], , drop = FALSE]
}

# MetaCyc reaction id → description lookup. Bundled from HUMAnN 4's
# utility_mapping (map_metacyc-rxn_name.txt.gz). Descriptions are prefixed
# with the source in the raw file (`(expasy) …` or `(metacyc) …`) — those
# tags are dropped here so the DT column shows just the human text.
# Missing on disk → empty vector; downstream code treats that as "no
# description available" and leaves the column blank rather than erroring.
metacyc_rxn_desc <- local({
  path <- "map_metacyc-rxn_name.txt.gz"
  if (!file.exists(path)) return(character(0))
  m <- tryCatch(
    read.table(gzfile(path), header = FALSE, sep = "\t", quote = "",
               comment.char = "", stringsAsFactors = FALSE,
               col.names = c("id", "name")),
    error = function(e) NULL
  )
  if (is.null(m) || nrow(m) == 0) return(character(0))
  m$name <- sub("^\\((expasy|metacyc)\\)\\s+", "", m$name)
  setNames(m$name, m$id)
})

# Per-slot database-link + description enrichment for HUMAnN tables. Feature
# ids are the leftmost token before ":" after stripping any "|s__..." tail;
# UniRef90_* accessions link to UniProtKB, UniClust90_* rows get "-" since
# HUMAnN's UniClust cluster numbers don't resolve on any public web page.
humann_link_specs <- list(
  pathways      = list(col   = "MetaCyc",
                       url   = "https://metacyc.org/pathway?orgid=META&id=",
                       id_of = function(base) trimws(sub(":.*$", "", base)),
                       blank = ""),
  kegg_kos      = list(col   = "KEGG",
                       url   = "https://www.kegg.jp/entry/",
                       id_of = function(base) trimws(sub(":.*$", "", base)),
                       blank = ""),
  gene_families = list(col   = "UniProt",
                       url   = "https://www.uniprot.org/uniprotkb/",
                       id_of = function(base) ifelse(grepl("^UniRef90_", base),
                                                     sub("^UniRef90_", "", base), ""),
                       blank = "-")
)

# Enrich a HUMAnN slot data frame with a database-link column and, for the
# reactions slot, a description column. Shared by the on-screen DT and the
# R-side CSV/TSV download handlers so both stay in sync.
humann_slot_enrich <- function(df, slot_name) {
  if (is.null(df) || nrow(df) == 0) return(df)
  spec <- humann_link_specs[[slot_name]]
  if (!is.null(spec)) {
    base_feat <- sub("\\|.*$", "", df$Feature)
    id        <- spec$id_of(base_feat)
    specials  <- c("UNMAPPED", "UNINTEGRATED", "UNGROUPED",
                   "READS_UNMAPPED")
    linkable  <- nzchar(id) & !id %in% specials
    df[[spec$col]] <- ifelse(linkable, id, spec$blank)
    df <- df[, c("Feature", spec$col,
                 setdiff(colnames(df), c("Feature", spec$col)))]
  }
  if (identical(slot_name, "reactions") && length(metacyc_rxn_desc) > 0) {
    base_feat <- sub("\\|.*$", "", df$Feature)
    desc <- unname(metacyc_rxn_desc[base_feat])
    desc[is.na(desc)] <- ""
    df$Description <- desc
    df <- df[, c("Feature", "Description",
                 setdiff(colnames(df), c("Feature", "Description")))]
  }
  df
}

# Shared exportOptions for every DT copy/CSV/Excel button. Enforces two
# things DataTables Buttons doesn't do by default:
#   modifier.page = "all"  — export every row across pages, not just the
#                            currently visible page.
#   format.body            — strip any residual HTML that leaks into a cell
#                            so createdCell-generated <a> nodes don't corrupt
#                            CSV quoting on rows whose plain data itself
#                            contains a `"` (e.g. KEGG K00984's name).
# orthogonal = "export" bypasses the render/createdCell pipeline for exports.
dt_export_options <- list(
  modifier   = list(page = "all", search = "applied"),
  orthogonal = "export",
  format     = list(
    body = DT::JS(
      "function(data){",
      "  if (data == null) return '';",
      "  var s = String(data);",
      "  if (s.indexOf('<') === -1) return s;",
      "  var tmp = document.createElement('div');",
      "  tmp.innerHTML = s;",
      "  return tmp.textContent || tmp.innerText || '';",
      "}"
    )
  )
)

# Bray-Curtis dissimilarity (pure base R; sample x sample matrix)
humann_bray <- function(mat) {
  # mat: rows = features, cols = samples
  n <- ncol(mat)
  d <- matrix(0, n, n, dimnames = list(colnames(mat), colnames(mat)))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      x <- mat[, i]; y <- mat[, j]
      num <- 2 * sum(pmin(x, y), na.rm = TRUE)
      den <- sum(x, na.rm = TRUE) + sum(y, na.rm = TRUE)
      d[i, j] <- d[j, i] <- if (den > 0) 1 - num / den else 0
    }
  }
  as.dist(d)
}

#' LDAP-based access check for a gstore URL load. Mirrors exploreDE.
#' Returns a list(allowed, role, reason, projects).
#'   - Empty SHINYPROXY_USERNAME → denied. Empty username means we're on the
#'     public ShinyProxy instance (no auth) or running outside ShinyProxy
#'     entirely; in neither case should a gstore ?data= URL be honoured.
#'     For local dev, uncomment a line at the top of this app.R:
#'         Sys.setenv(SHINYPROXY_USERNAME = Sys.getenv("USER"))
#'   - course_sushi paths → allowed without LDAP (course data is public).
#'   - FGCZ Employee (LDAP role R_2) → allowed to any project.
#'   - FGCZ User (R_3) → allowed only if extracted p<NNNN> ∈ memberOf cn=P_<NNNN>.
#'   - Other / no role → denied.
ldap_access_check <- function(raw_path, resolved) {
  username <- Sys.getenv("SHINYPROXY_USERNAME")
  # Defensive: usernames at FGCZ are simple ASCII; strip anything else
  # before splicing into a shell command.
  username <- gsub("[^a-zA-Z0-9._-]", "", username)

  if (!nzchar(username)) {
    return(list(allowed = FALSE, role = "anonymous",
                projects = character(0),
                reason = "no SHINYPROXY_USERNAME — gstore URL loads are not available on the public instance"))
  }
  if (grepl("/srv/GT/analysis/course_sushi/public/gstore/projects",
            resolved, fixed = TRUE)) {
    return(list(allowed = TRUE, role = "course",
                projects = character(0),
                reason = "course_sushi path — no LDAP check"))
  }

  m <- regmatches(raw_path, regexec("p[0-9]{4,}", raw_path))[[1]]
  projectFromUrl <- if (length(m) > 0) m[1] else NA_character_

  roles <- tryCatch(
    system(paste0(
      "ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' ",
      "'(cn=", username, ")' memberof | grep Roles | sed 's/,ou=.*//g;s,.*cn=,,g'"
    ), intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character(0)
  )
  allowedProjects <- tryCatch(
    system(paste0(
      "ldapsearch -x -H ldaps://fgcz-bfabric-ldap:636 -b 'dc=bfabric,dc=org' ",
      "'(cn=", username, ")' memberof | grep Projec | sed 's/,ou=.*//g;s,.*cn=P_,p,g' | sort | uniq"
    ), intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character(0)
  )

  if ("R_2" %in% roles) {
    list(allowed = TRUE, role = "employee", projects = allowedProjects,
         reason = "FGCZ employee — full access")
  } else if ("R_3" %in% roles) {
    granted <- !is.na(projectFromUrl) && projectFromUrl %in% allowedProjects
    list(allowed = granted, role = "user", projects = allowedProjects,
         reason = if (granted) paste0("member of ", projectFromUrl)
                  else paste0("not a member of ", projectFromUrl %||% "?"))
  } else {
    list(allowed = FALSE, role = "unauthorized", projects = allowedProjects,
         reason = "no recognised FGCZ role")
  }
}


# ═══════════════════════════════════════════════════════════════════════════
# UI — shinydashboard layout (matching exploreDE aesthetics)
# ═══════════════════════════════════════════════════════════════════════════

my_theme <- create_theme(
  adminlte_color(
    light_blue = "#C4692A"
  )
)

ui <- dashboardPage(
  dashboardHeader(
    title = "exploreMetaTax"
  ),
  dashboardSidebar(
    shinyjs::useShinyjs(),
    # Sidebar contents are rendered server-side so we can swap between the
    # taxonomy menu (default) and the HUMAnN-only menu once a HUMAnN dataset
    # has been loaded via ?data=... URL.
    sidebarMenuOutput("dynamicSidebar")
  ),
  dashboardBody(
    use_theme(my_theme),
    tags$head(
      # Custom message handler used by the Load & Export Settings tab to
      # restore sortable::bucket_list state. The R side can't drive these
      # widgets via update_*() (none exists) or sendInputMessage (their JS
      # binding ignores it), so instead we ship the saved order to the
      # client and reorder DOM items across the two buckets ourselves,
      # then push the new value back into Shiny via setInputValue.
      tags$script(HTML(paste(sep = "\n",
        "Shiny.addCustomMessageHandler('restore_sortable_buckets', function(msg) {",
        "  var collect = function(root) {",
        "    // rank_list items have class 'rank-list-item'; fall back to",
        "    // direct children if the class ever changes.",
        "    var nodes = root.querySelectorAll('.rank-list-item');",
        "    if (nodes.length === 0) nodes = root.children;",
        "    var map = {};",
        "    Array.prototype.slice.call(nodes).forEach(function(el) {",
        "      var t = (el.textContent || '').trim();",
        "      if (t) map[t] = el;",
        "    });",
        "    return map;",
        "  };",
        "  var applyPair = function() {",
        "    var v = document.getElementById(msg.visible_id);",
        "    var h = document.getElementById(msg.hidden_id);",
        "    if (!v || !h) return false;",
        "    var itemsV = collect(v), itemsH = collect(h);",
        "    var items = Object.assign({}, itemsV, itemsH);",
        "    var vOrd = (msg.visible_order || []).filter(function(s) { return items[s]; });",
        "    var hOrd = (msg.hidden_order  || []).filter(function(s) { return items[s]; });",
        "    if (vOrd.length + hOrd.length === 0) return false;",
        "    // Move items to their target bucket in the target order. Since",
        "    // appendChild moves DOM nodes across parents, this handles",
        "    // both cross-bucket transfers and within-bucket reordering.",
        "    vOrd.forEach(function(s) { v.appendChild(items[s]); });",
        "    hOrd.forEach(function(s) { h.appendChild(items[s]); });",
        "    Shiny.setInputValue(msg.visible_id, vOrd, { priority: 'event' });",
        "    Shiny.setInputValue(msg.hidden_id,  hOrd, { priority: 'event' });",
        "    return true;",
        "  };",
        "  var tries = 0;",
        "  var attempt = function() {",
        "    if (applyPair() || tries >= 60) return;",
        "    tries += 1;",
        "    setTimeout(attempt, 150);",
        "  };",
        "  attempt();",
        "})"
      ))),
      tags$style(HTML("
        .shiny-split-layout > div {
          overflow: visible;
        }
        /* SortableJS MultiDrag: highlight click-selected sample chips so it's
           clear which ones will move together on the next drag. */
        .rank-list-item.sortable-selected,
        .sortable-selected {
          background: #F6D9BE !important;
          outline: 2px solid #C4692A;
          outline-offset: -2px;
        }
        .box {
          box-shadow: 0 4px 8px rgba(0,0,0,0.1);
          transition: all 0.3s ease;
        }
        .box:hover {
          box-shadow: 0 8px 16px rgba(0,0,0,0.2);
        }
        .box.box-solid.box-primary>.box-header {
          color:#fff;
          background: linear-gradient(135deg, #C4692A 0%, #E6A76A 100%);
        }
        .box.box-solid.box-primary{
          border-bottom-color:#C4692A;
          border-left-color:#C4692A;
          border-right-color:#C4692A;
          border-top-color:#C4692A;
        }
        .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
          background-color: #D68C48;
          transition: background-color 0.3s ease;
        }
        .skin-blue .sidebar-menu > li:hover > a {
          border-left-color: #9A4F17;
        }
        /* Main Header Gradient */
        .skin-blue .main-header .logo {
          background: linear-gradient(135deg, #C4692A 0%, #E6A76A 100%) !important;
        }
        .skin-blue .main-header .navbar {
          background: linear-gradient(135deg, #C4692A 0%, #E6A76A 100%) !important;
        }
        /* body */
        .content-wrapper, .right-side {
          background-color: #FFFFFF;
        }
        .sidebar-menu > li.header {
          font-size: 16px;
          font-weight: bold;
          color: #b8c7ce;
          padding: 10px 25px 10px 15px;
          text-shadow: 1px 1px 2px rgba(0,0,0,0.1);
        }
        .status-text {
          background: #FDF3EA;
          border: 1px solid #E6A76A;
          border-radius: 6px;
          padding: 8px 12px;
          font-size: 0.85em;
          color: #C4692A;
          margin-top: 8px;
        }
        .status-warning {
          background: #fef9e7;
          border-color: #f4d03f;
          color: #7d6608;
        }
        .btn-load {
          background: linear-gradient(135deg, #C4692A 0%, #E6A76A 100%);
          color: #fff;
          border: none;
          border-radius: 4px;
          font-weight: 500;
          padding: 8px 20px;
          transition: all 0.2s ease;
        }
        .btn-load:hover {
          background: linear-gradient(135deg, #9A4F17 0%, #D68C48 100%);
          box-shadow: 0 4px 12px rgba(196, 105, 42, 0.3);
          color: #fff;
        }
        /* Loading spinner for the Krona widget. Shiny auto-toggles the
           `recalculating` class on outputs while their render runs; we
           paint a centered spinner + label overlay while it's there.
           Scoped to #krona_plot so other outputs aren't affected. */
        #krona_plot {
          position: relative;
        }
        #krona_plot.recalculating::before {
          content: '';
          position: absolute;
          z-index: 10;
          top: 50%;
          left: 50%;
          width: 56px;
          height: 56px;
          margin: -28px 0 0 -28px;
          border: 5px solid #E6A76A;
          border-top-color: #C4692A;
          border-radius: 50%;
          animation: krona-spin 0.9s linear infinite;
        }
        #krona_plot.recalculating::after {
          content: 'Rendering Krona chart…';
          position: absolute;
          z-index: 10;
          top: calc(50% + 36px);
          left: 0;
          right: 0;
          text-align: center;
          color: #C4692A;
          font-size: 14px;
          font-weight: 500;
        }
        @keyframes krona-spin {
          to { transform: rotate(360deg); }
        }
        /* Explicit loading banner shown for a 5s floor on each Krona
           (re)generation — covers the client-side SVG build, which lags after
           the server value arrives and the `recalculating` spinner clears. */
        .krona-loading-banner {
          display: flex;
          align-items: center;
          gap: 12px;
          background: linear-gradient(135deg, #C4692A 0%, #E6A76A 100%);
          color: #fff;
          padding: 12px 18px;
          border-radius: 6px;
          margin-bottom: 12px;
          font-weight: 500;
          box-shadow: 0 2px 6px rgba(0,0,0,0.15);
        }
        .krona-loading-spinner {
          flex: 0 0 auto;
          width: 22px;
          height: 22px;
          border: 3px solid rgba(255,255,255,0.4);
          border-top-color: #fff;
          border-radius: 50%;
          animation: krona-spin 0.9s linear infinite;
        }
        /* Full-viewport modal overlay used during blocking loads
           (URL-triggered HUMAnN + file-upload paths). Kept simple — a
           centered spinner + label on a translucent backdrop. Toggled by
           the `appLoading` custom message from the server. */
        #app_loading_overlay {
          position: fixed;
          inset: 0;
          background: rgba(0,0,0,0.55);
          backdrop-filter: blur(2px);
          -webkit-backdrop-filter: blur(2px);
          z-index: 9999;
          display: none;
          align-items: center;
          justify-content: center;
        }
        #app_loading_overlay.show { display: flex; }
        #app_loading_overlay .app-loading-box {
          color: #fff;
          text-align: center;
          padding: 24px 32px;
        }
        #app_loading_overlay .app-loading-spinner {
          margin: 0 auto 20px;
          width: 72px;
          height: 72px;
          border: 6px solid rgba(255,255,255,0.35);
          border-top-color: #fff;
          border-radius: 50%;
          animation: krona-spin 0.9s linear infinite;
        }
        #app_loading_overlay .app-loading-text {
          font-size: 20px;
          font-weight: 500;
          letter-spacing: 0.4px;
        }
      ")),
      tags$script(HTML(
        "Shiny.addCustomMessageHandler('appLoading', function(msg){",
        "  var el = document.getElementById('app_loading_overlay');",
        "  if (!el) return;",
        "  var state = (typeof msg === 'object') ? msg.state : msg;",
        "  if (state === 'start') {",
        "    if (typeof msg === 'object' && msg.text) {",
        "      var t = el.querySelector('.app-loading-text');",
        "      if (t) t.textContent = msg.text;",
        "    }",
        "    el.classList.add('show');",
        "  } else {",
        "    el.classList.remove('show');",
        "  }",
        "});"
      ))
    ),
    tags$div(id = "app_loading_overlay",
      tags$div(class = "app-loading-box",
        tags$div(class = "app-loading-spinner"),
        tags$div(class = "app-loading-text", "Loading…")
      )
    ),
    tabItems(

      # ════════════════════════════════════════════════════════════════════
      # TAB: Upload Data
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "uploadTab",
        fluidPage(
          fluidRow(
            column(
              width = 6,
              box(
                title = "Upload Report Files",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                radioButtons("input_format", "Input Format:",
                             choices = c("Kraken2 Reports" = "kraken2",
                                         "KrakenUniq Reports" = "krakenuniq",
                                         "Bracken Per-Sample" = "bracken",
                                         "Combined Bracken Table" = "combined_bracken",
                                         "MetaPhlAn Profile" = "metaphlan",
                                         "Merged MetaPhlAn Table" = "merged_metaphlan",
                                         "HUMAnN 4 Per-Sample TSVs" = "humann"),
                             selected = "kraken2", inline = FALSE),
                helpText("Format is auto-detected when possible. Select manually if detection fails."),
                conditionalPanel(
                  condition = "input.input_format == 'humann'",
                  helpText(
                    tags$b("HUMAnN uploads:"),
                    " select (or folder-browse) the per-sample TSVs produced by HUMAnN 4. ",
                    "The app groups files by table using the filename suffix — recognised ",
                    "patterns are ",
                    tags$code("*_pathabundance.tsv"), ", ",
                    tags$code("*_reactions.tsv"), ", ",
                    tags$code("*_ko_cpm.tsv"), " (or ", tags$code("*_kegg_kos.tsv"), "), and ",
                    tags$code("*_genefamilies*.tsv"), ". ",
                    "The sample name is inferred from the filename prefix (the ",
                    tags$code("_4_"), " infix used by the SUSHI HUMAnN app is stripped)."
                  )
                ),
                checkboxInput("upload_folder",
                              "Browse a folder (uploads all matching reports inside)",
                              value = FALSE),
                fileInput("data_files", "Upload Report Files or Folder:",
                          multiple = TRUE,
                          accept = c(".txt", ".tsv", ".kreport", ".kreport2", ".bracken",
                                     ".profile", ".biom")),
                actionButton("load_data", "Load Data", icon = icon("play"),
                             class = "btn-load", width = "100%"),
                uiOutput("status_text")
              ),
              box(
                title = "Metadata (Optional)",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                helpText("Upload a TSV or CSV with a 'Sample' column matching your sample names.",
                         "Additional columns (Group, Treatment, etc.) enable group-based analyses."),
                downloadButton("download_example_metadata", "See Example",
                               icon = icon("eye")),
                hr(),
                fileInput("metadata_file", "Upload Metadata (TSV/CSV):",
                          multiple = FALSE,
                          accept = c(".tsv", ".csv", ".txt")),
                actionButton("load_metadata", "Load Metadata", icon = icon("plus-circle"),
                             class = "btn-load", width = "100%"),
                uiOutput("metadata_status")
              )
            ),
            column(
              width = 6,
              box(
                title = "Sample Names",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                helpText("Rename your samples here. Changes propagate to all tabs."),
                uiOutput("rename_ui")
              ),
              box(
                title = "Filters",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                fluidRow(
                  column(4,
                    selectInput("rank_select", "Taxonomic Rank:",
                                choices = rank_choices, selected = "S")
                  ),
                  column(4,
                    sliderInput("min_abundance", "Min Abundance (%):",
                                min = 0, max = 50, value = 0, step = 0.1)
                  ),
                  column(4,
                    numericInput("min_counts", "Min Counts:",
                                value = 0, min = 0, step = 1)
                  )
                ),
                uiOutput("minimizer_filters_ui"),
                uiOutput("sample_checkboxes"),
                hr(),
                tags$b("Entries to discard"),
                helpText(
                  "Remove any taxon matching these queries from every tab ",
                  tags$b("except the Data Table"), ". Cascading match per query: ",
                  tags$b("exact name"), " → ",
                  tags$b("name substring"), " → ",
                  tags$b("Genus lineage substring"), ". ",
                  "Applied globally across all ranks — click ",
                  tags$b("Apply discards"), " to commit."
                ),
                textAreaInput(
                  "discard_query_text",
                  label = NULL,
                  value = "",
                  width = "100%", height = "90px",
                  placeholder = "e.g.\nHomo\nunclassified\nBacteroides fragilis"
                ),
                actionButton("discard_apply_btn", "Apply discards",
                             icon = icon("filter"), class = "btn-default"),
                actionButton("discard_clear_btn", "Clear"),
                uiOutput("discard_summary")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Composition Barplot
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "compositionTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                fluidRow(
                  column(3, numericInput("top_n", "Top N Taxa:", value = 10, min = 1, max = 100)),
                  column(3, selectInput("display_mode", "Display Mode:",
                                        choices = c("Percentages" = "percent", "Counts" = "counts"),
                                        selected = "percent")),
                  column(3, selectInput("color_palette", "Color Palette:",
                                        choices = palette_choices, selected = "Set3")),
                  column(3, checkboxInput("show_labels", "Show Labels on Bars", value = TRUE))
                ),
                fluidRow(
                  column(6, uiOutput("group_by_selector"))
                ),
                hr(),
                fluidRow(
                  column(3, numericInput("plot_width", "Export Width (in):", value = 12, min = 4, max = 24)),
                  column(3, numericInput("plot_height", "Export Height (in):", value = 7, min = 3, max = 20)),
                  column(6,
                    downloadButton("download_pdf", "PDF"),
                    downloadButton("download_png", "PNG"),
                    downloadButton("download_svg", "SVG"),
                    downloadButton("download_csv", "CSV")
                  )
                )
              ),
              box(
                title = "Taxonomic Composition",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                plotlyOutput("composition_plot", height = "650px")
              ),
              box(
                title = "Group comparison (top-N taxa)",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                helpText("Per-taxon Wilcoxon rank-sum (2 groups) or ",
                         "Kruskal–Wallis (≥3 groups) on the top-N taxa, ",
                         "with Benjamini–Hochberg FDR correction across taxa. ",
                         "Requires a metadata Group By with ≥2 samples per group (n<3 flagged in the note column — min achievable p is limited). ",
                         "Values tested are those shown on the y-axis (percent or counts)."),
                uiOutput("composition_stats_status"),
                DTOutput("composition_stats_table"),
                br(),
                downloadButton("download_composition_stats_csv", "Download stats (CSV)")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Organisms of Interest
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "organismsTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Search",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                fluidRow(
                  column(
                    width = 7,
                    textAreaInput(
                      "organism_query_text",
                      label = "Taxon names (one per line, or comma-separated):",
                      value = "",
                      width = "100%", height = "120px",
                      placeholder = "e.g.\nBacteroides fragilis\nEscherichia coli\nFaecalibacterium"
                    )
                  ),
                  column(
                    width = 5,
                    div(
                      style = "padding: 5px 0;",
                      tags$p(tags$strong("Match strategy (cascading; stops at first non-empty tier):")),
                      tags$ol(
                        tags$li(tags$b("exact"), " — case-insensitive equality on the current-rank ", tags$code("name"), "."),
                        tags$li(tags$b("substring"), " — ", tags$code("name"), " contains the query (e.g. ", tags$code("Bacteroides"), " finds ", tags$em("Bacteroides fragilis"), ")."),
                        tags$li(tags$b("genus fallback"), " — first word of the query, substring-matched against the ", tags$code("Genus"), " column.")
                      ),
                      tags$p(tags$em("All other tab filters (rank, samples, min abundance/counts) still apply."))
                    )
                  )
                ),
                fluidRow(
                  column(12,
                    actionButton("organism_search_btn", "Search", icon = icon("search"),
                                 class = "btn-primary"),
                    actionButton("organism_clear_btn", "Clear"),
                    uiOutput("organism_external_links", inline = TRUE)
                  )
                )
              ),
              box(
                title = "Matches",
                width = NULL, solidHeader = TRUE, status = "info",
                collapsible = TRUE,
                uiOutput("organism_match_summary")
              ),
              box(
                title = "Composition (matched taxa only)",
                width = NULL, solidHeader = TRUE, status = "primary",
                fluidRow(
                  column(3, selectInput("organism_display_mode", "Display Mode:",
                                        choices = c("Percentages" = "percent", "Counts" = "counts"),
                                        selected = "percent")),
                  column(3, selectInput("organism_color_palette", "Color Palette:",
                                        choices = palette_choices, selected = "Set3")),
                  column(3, checkboxInput("organism_show_labels", "Show Labels on Bars", value = TRUE)),
                  column(3, uiOutput("organism_group_selector"))
                ),
                hr(),
                fluidRow(
                  column(3, numericInput("organism_plot_width", "Export Width (in):",
                                         value = 12, min = 4, max = 24)),
                  column(3, numericInput("organism_plot_height", "Export Height (in):",
                                         value = 7, min = 3, max = 20)),
                  column(6,
                    downloadButton("download_organism_pdf", "PDF"),
                    downloadButton("download_organism_png", "PNG"),
                    downloadButton("download_organism_svg", "SVG"),
                    downloadButton("organism_download_csv", "Matched rows (CSV)")
                  )
                ),
                plotlyOutput("organism_plot", height = "550px")
              ),
              box(
                title = "Group comparison (matched taxa)",
                width = NULL, solidHeader = TRUE, status = "primary",
                collapsible = TRUE, collapsed = FALSE,
                helpText("Per-taxon Wilcoxon rank-sum (2 groups) or Kruskal–Wallis ",
                         "(≥3 groups) on the matched taxa, with Benjamini–Hochberg ",
                         "FDR across taxa. Requires a metadata Group By with ≥3 ",
                         "samples per group."),
                uiOutput("organism_stats_status"),
                DTOutput("organism_stats_table"),
                br(),
                downloadButton("download_organism_stats_csv", "Download stats (CSV)")
              ),
              box(
                title = "Distribution (matched taxa only)",
                width = NULL, solidHeader = TRUE, status = "primary",
                collapsible = TRUE, collapsed = FALSE,
                fluidRow(
                  column(3, selectInput("organism_view_mode", "View Mode:",
                                        choices = c("Single Taxon" = "single",
                                                    "Faceted (All matched)" = "faceted"),
                                        selected = "single")),
                  column(3, conditionalPanel(
                    condition = "input.organism_view_mode == 'single'",
                    selectInput("organism_violin_taxon", "Select Taxon:", choices = NULL)
                  )),
                  column(3, selectInput("organism_plot_type", "Plot Type:",
                                        choices = c("Violin" = "violin",
                                                    "Box + Strip" = "box",
                                                    "Strip Plot" = "dot"))),
                  column(3, uiOutput("organism_distrib_group_selector"))
                ),
                conditionalPanel(
                  condition = "input.organism_view_mode == 'single'",
                  uiOutput("organism_taxon_search_links")
                ),
                plotlyOutput("organism_violin_plot", height = "600px"),
                br(),
                downloadButton("download_organism_distrib_pdf", "PDF"),
                downloadButton("download_organism_distrib_png", "PNG"),
                downloadButton("download_organism_distrib_svg", "SVG")
              ),
              box(
                title = "Statistics (Distribution — matched taxa)",
                width = NULL, solidHeader = TRUE, status = "primary",
                collapsible = TRUE, collapsed = FALSE,
                helpText("Wilcoxon rank-sum (2 groups) or Kruskal–Wallis (≥3 groups) ",
                         "on abundance across the selected metadata Group By. Faceted ",
                         "mode compares each matched taxon with BH FDR across taxa."),
                uiOutput("organism_distrib_stats_status"),
                DTOutput("organism_distrib_stats_table"),
                br(),
                uiOutput("organism_distrib_pairwise_header"),
                DTOutput("organism_distrib_stats_pairwise"),
                br(),
                downloadButton("download_organism_distrib_stats_csv",
                               "Download stats (CSV)")
              ),
              box(
                title = "Matched rows",
                width = NULL, solidHeader = TRUE, status = "primary",
                collapsible = TRUE,
                DTOutput("organism_table")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Heatmap
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "heatmapTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Abundance Heatmap",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                fluidRow(
                  column(6, uiOutput("heatmap_group_selector")),
                  column(6,
                    div(style = "margin-top: 25px;",
                      downloadButton("download_heatmap_pdf", "PDF"),
                      downloadButton("download_heatmap_png", "PNG"),
                      downloadButton("download_heatmap_svg", "SVG")
                    )
                  )
                ),
                br(),
                plotlyOutput("heatmap_plot", height = "700px")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Distribution (Violin/Box)
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "distributionTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                fluidRow(
                  column(3, selectInput("violin_view_mode", "View Mode:",
                                        choices = c("Single Taxon" = "single",
                                                    "Faceted (All Top N)" = "faceted"),
                                        selected = "single")),
                  column(3, conditionalPanel(
                    condition = "input.violin_view_mode == 'single'",
                    selectInput("violin_taxon", "Select Taxon:", choices = NULL)
                  )),
                  column(3, selectInput("violin_plot_type", "Plot Type:",
                                        choices = c("Violin" = "violin",
                                                    "Box + Strip" = "box",
                                                    "Strip Plot" = "dot"))),
                  column(3, uiOutput("violin_group_selector"))
                ),
                # Search links (shown when single taxon is selected)
                conditionalPanel(
                  condition = "input.violin_view_mode == 'single'",
                  uiOutput("taxon_search_links")
                )
              ),
              box(
                title = "Distribution Plot",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                fluidRow(
                  column(6,
                    downloadButton("download_distrib_pdf", "PDF"),
                    downloadButton("download_distrib_png", "PNG"),
                    downloadButton("download_distrib_svg", "SVG")
                  )
                ),
                br(),
                plotlyOutput("violin_plot", height = "650px")
              ),
              box(
                title = "Statistics",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                helpText("Wilcoxon rank-sum (2 groups) or Kruskal–Wallis (≥3 groups) ",
                         "on abundance across the selected metadata Group By. ",
                         "For ≥3 groups an additional pairwise Wilcoxon table with ",
                         "Benjamini–Hochberg FDR is shown. Faceted mode compares each ",
                         "top-N taxon separately with BH correction across taxa."),
                uiOutput("distrib_stats_status"),
                DTOutput("distrib_stats_table"),
                br(),
                uiOutput("distrib_pairwise_header"),
                DTOutput("distrib_stats_pairwise"),
                br(),
                downloadButton("download_distrib_stats_csv", "Download stats (CSV)")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Alpha Diversity
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "diversityTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                fluidRow(
                  column(4, selectInput("diversity_index", "Diversity Index:",
                                        choices = c("Shannon", "Simpson", "Inverse Simpson"),
                                        selected = "Shannon")),
                  column(4, uiOutput("diversity_color_selector"))
                )
              ),
              box(
                title = "Alpha Diversity",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                fluidRow(
                  column(6,
                    downloadButton("download_alpha_pdf", "PDF"),
                    downloadButton("download_alpha_png", "PNG"),
                    downloadButton("download_alpha_svg", "SVG")
                  )
                ),
                br(),
                plotlyOutput("diversity_plot", height = "550px")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Rarefaction Curves
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "rarefactionTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                helpText("Rarefaction curves plot observed taxon richness against",
                         "sequencing depth, repeatedly subsampling each sample down",
                         "in fixed steps. A curve that plateaus means the sample was",
                         "sequenced deeply enough to capture most of its taxa.",
                         "Requires read counts — not available for MetaPhlAn",
                         "relative-abundance profiles. Uses the rank, samples and",
                         "filters selected in the sidebar."),
                fluidRow(
                  column(3, numericInput("rarefy_points", "Curve resolution (points/sample):",
                                         value = 150, min = 20, max = 1000, step = 10)),
                  column(3, uiOutput("rarefaction_color_selector")),
                  column(3, selectInput("rarefy_palette", "Color Palette:",
                                        choices = palette_choices, selected = "Set3")),
                  column(3,
                    tags$label("Common depth", class = "control-label"),
                    checkboxInput("rarefy_normalize",
                                  "Cap curves at shallowest sample",
                                  value = FALSE))
                )
              ),
              box(
                title = "Rarefaction Curves",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                fluidRow(
                  column(6,
                    downloadButton("download_rarefaction_pdf", "PDF"),
                    downloadButton("download_rarefaction_png", "PNG"),
                    downloadButton("download_rarefaction_svg", "SVG")
                  )
                ),
                br(),
                plotlyOutput("rarefaction_plot", height = "600px")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: PCA / Beta Diversity
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "pcaTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                fluidRow(
                  column(4, uiOutput("pca_color_selector")),
                  column(4, selectInput("pca_dist_method", "Distance:",
                                        choices = c("Euclidean" = "euclidean",
                                                    "Bray-Curtis" = "bray",
                                                    "Jaccard" = "jaccard"),
                                        selected = "euclidean"))
                )
              )
            )
          ),
          fluidRow(
            column(
              width = 6,
              box(
                title = "PCA Ordination",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                fluidRow(
                  column(6,
                    downloadButton("download_pca_pdf", "PDF"),
                    downloadButton("download_pca_png", "PNG"),
                    downloadButton("download_pca_svg", "SVG")
                  )
                ),
                br(),
                plotlyOutput("pca_plot", height = "500px")
              )
            ),
            column(
              width = 6,
              box(
                title = "Hierarchical Clustering",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                fluidRow(
                  column(6,
                    downloadButton("download_dendro_pdf", "PDF"),
                    downloadButton("download_dendro_png", "PNG"),
                    downloadButton("download_dendro_svg", "SVG")
                  )
                ),
                br(),
                plotlyOutput("dendro_plot", height = "500px")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: LifemapR Tree
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "lifemapTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "LifemapR Taxonomy Tree",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                if (HAS_LIFEMAPR) {
                  tagList(
                    uiOutput("lifemap_warning"),
                    helpText("Project detected taxa onto the NCBI Tree of Life.",
                             "Requires internet access and NCBI TaxIDs.",
                             "Uses the samples and filters selected in the sidebar."),
                    actionButton("run_lifemap", "Generate Tree", icon = icon("seedling"),
                                 class = "btn-load"),
                    br(), br(),
                    uiOutput("lifemap_status")
                  )
                } else {
                  div(class = "status-text status-warning",
                      icon("exclamation-triangle"),
                      " LifemapR package not installed. Install with: ",
                      code('remotes::install_github("damiendevienne/LifemapR")'))
                }
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Krona
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "kronaTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Settings",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                fluidRow(
                  column(4, uiOutput("krona_sample_selector")),
                  column(4, checkboxInput("krona_use_abundance", "Scale by abundance", value = TRUE)),
                  column(4, selectInput("krona_value", "Abundance from:",
                                        choices = c("Percentage" = "percent",
                                                    "Counts" = "reads_clade"),
                                        selected = "percent"))
                ),
                uiOutput("krona_group_ui")
              ),
              box(
                title = "Krona Chart",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                if (HAS_TAXPLORE) {
                  tagList(
                    shinyjs::hidden(
                      div(id = "krona_loading", class = "krona-loading-banner",
                          div(class = "krona-loading-spinner"),
                          div(class = "krona-loading-text",
                              "Generating Krona chart — please wait while it renders…"))
                    ),
                    KronaChartOutput("krona_plot", height = "700px")
                  )
                } else {
                  div(class = "status-text status-warning",
                      icon("exclamation-triangle"),
                      " The ", code("taxplore"), " package is required for Krona charts. ",
                      "Install with: ", code('remotes::install_github("markschl/taxplore")'))
                }
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Data Table
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "tableTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = "Full Data Table",
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                downloadButton("download_table_csv", "Export Full Table (CSV)"),
                br(), br(),
                DTOutput("full_data_table")
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: HUMAnN (only shown when mode == "humann")
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "humannTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = tagList(icon("sliders-h"),
                                "Filters & renaming (plots only)"),
                width = NULL, status = "info", solidHeader = FALSE,
                collapsible = TRUE, collapsed = TRUE,
                helpText("Applies to HUMAnN plot bars only — the 'Full table' ",
                         "DTs always show the raw, unfiltered data."),
                fluidRow(
                  column(6,
                    uiOutput("humann_sample_filter_ui"),
                    tags$br(),
                    uiOutput("humann_group_by_ui"),
                    uiOutput("humann_meta_filter_ui")
                  ),
                  column(6,
                    tags$h5("Rename samples for plots"),
                    helpText("Double-click a cell in the 'Display' column ",
                             "to edit; clear it to reset to the original."),
                    DTOutput("humann_rename_dt", height = "300px")
                  )
                )
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              box(
                title = "HUMAnN functional profiles",
                width = NULL,
                solidHeader = TRUE, status = "primary",
                tabsetPanel(
                  id = "humann_subtab",
                  tabPanel("Info",              uiOutput("humann_info_ui")),
                  tabPanel("Pathways",          uiOutput("humann_pathways_ui")),
                  tabPanel("Reactions",         uiOutput("humann_reactions_ui")),
                  tabPanel("KEGG KOs",          uiOutput("humann_kegg_kos_ui")),
                  tabPanel("Gene families",     uiOutput("humann_gene_families_ui")),
                  tabPanel("Stratified by organisms",
                                                uiOutput("humann_stratified_ui")),
                  tabPanel("Method notes",      uiOutput("humann_method_ui"))
                )
              )
            )
          )
        )
      ),

      # ════════════════════════════════════════════════════════════════════
      # TAB: Load & Export Settings
      # ════════════════════════════════════════════════════════════════════
      tabItem(
        tabName = "settingsTab",
        fluidPage(
          fluidRow(
            column(
              width = 12,
              box(
                title = tagList(icon("save"),
                                "Save current app state"),
                width = NULL, status = "primary", solidHeader = TRUE,
                helpText("Download a single .rds bundle that captures the ",
                         "loaded data, sample renames, discard filters, and ",
                         "every widget value across all tabs. The file stays ",
                         "on your machine — nothing is uploaded to the server."),
                downloadButton("download_settings",
                               "Save settings (.rds)",
                               class = "btn-primary"),
                br(), br(),
                verbatimTextOutput("settings_save_status")
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              box(
                title = tagList(icon("upload"),
                                "Load a previously saved state"),
                width = NULL, status = "primary", solidHeader = TRUE,
                helpText("Upload an .rds bundle exported from a previous ",
                         "session to restore the same data, filters and ",
                         "widget values across every tab."),
                fileInput("load_settings_file",
                          "Choose settings file (.rds)",
                          accept = c(".rds", "application/octet-stream")),
                verbatimTextOutput("load_settings_status")
              )
            )
          )
        )
      )
    )
  )
)



# ═══════════════════════════════════════════════════════════════════════════
# SERVER
# ═══════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Reactive values ──
  rv <- reactiveValues(
    raw_data = NULL,       # Parsed data (all samples, long format)
    metadata = NULL,       # Optional metadata
    sample_renames = NULL, # Named vector: old_name => new_name
    loaded = FALSE,
    load_time = NULL,      # Timer for parsing execution
    lifemap_obj = NULL,    # LifemapR build_Lifemap result
    lifemap_port = NULL,   # Port for lifemap sub-app
    lifemap_bg = NULL,     # Background process for lifemap sub-app
    # ── HUMAnN mode ──
    app_mode = "taxonomy", # "taxonomy" | "humann"
    humann = NULL,         # list of joined data frames per HUMAnN table
    humann_parsed_ds = NULL, # parsed dataset.tsv, kept so gene_families can be
                             # loaded lazily on first tab open (see below)
    humann_modes = NULL,   # named char vector: sample -> "full"/"translated"
    humann_source_dir = NULL, # absolute path of the dataset.tsv parent
    discard_committed = character(0),  # taxa-discard queries (committed on Apply)
    # ── Restore-from-file state (used by the Load & Export Settings tab) ──
    # pending_ui_state holds sortable-widget snapshots (samples_visible,
    # samples_hidden, humann_visible_samples, humann_hidden_samples,
    # sample_filter) that the relevant renderUIs consult on next render.
    # pending_inputs holds the flat input$ snapshot; a retry observer
    # sendInputMessage()s each entry until the corresponding widget exists.
    pending_ui_state = list(),
    pending_inputs = list(),
    # restore_inputs mirrors the input snapshot but PERSISTS until a new dataset
    # is loaded. pending_inputs drives a ~6 s sendInputMessage retry loop that
    # only reaches widgets already in the DOM; every "Group by"/"Color by"
    # dropdown, the per-column HUMAnN filters and the minimizer filters live
    # inside a renderUI on a tab that is suspended (hidden) at load time, so they
    # miss that window. Those widgets instead seed their saved value from this
    # store via restored_input() whenever they first render — even on a tab
    # opened long after the file was loaded.
    restore_inputs = list()
  )

  # Seed value for a dynamically-rendered (renderUI / update*Input) widget from
  # a loaded settings file. Returns the saved value for `id`, or `default` when
  # nothing was saved. `valid` (optional) restricts to values still present in
  # the widget's current choices — for a single-select a stale value falls back
  # to `default`; for a multi-select the surviving subset is kept. isolate() so
  # reading the store never makes the calling renderUI re-run when the store
  # changes (it re-runs off its data dependency, e.g. rv$metadata, instead).
  restored_input <- function(id, default = NULL, valid = NULL) {
    v <- isolate(rv$restore_inputs[[id]])
    if (is.null(v)) return(default)
    if (!is.null(valid)) {
      v <- v[v %in% valid]
      if (length(v) == 0) return(default)
    }
    v
  }

  # ── Dynamic sidebar — taxonomy menu by default, HUMAnN-only when a
  # functional dataset is loaded. The kraken/metaphlan tabs and the HUMAnN
  # tab are mutually exclusive: showing both would imply the same `rv`
  # carries both kinds of data, which it does not.
  output$dynamicSidebar <- renderMenu({
    if (isTRUE(rv$app_mode == "humann")) {
      sidebarMenu(
        id = "tabs",
        tags$li(class = "header", "HUMAnN"),
        menuItem("Functional profiles", tabName = "humannTab",
                 icon = icon("dna"), selected = TRUE),
        tags$li(class = "header", "Settings"),
        menuItem("Load & Export Settings", tabName = "settingsTab",
                 icon = icon("save"))
      )
    } else {
      sidebarMenu(
        id = "tabs",
        tags$li(class = "header", "Data Input"),
        menuItem("Upload and Filter Data", tabName = "uploadTab", icon = icon("upload")),
        tags$li(class = "header", "Visualisations"),
        menuItem("Composition", tabName = "compositionTab", icon = icon("chart-bar")),
        menuItem("Heatmap", tabName = "heatmapTab", icon = icon("th")),
        menuItem("Krona", tabName = "kronaTab", icon = icon("circle-notch")),
        menuItem("Distribution", tabName = "distributionTab", icon = icon("chart-area")),
        menuItem("Organisms of Interest", tabName = "organismsTab", icon = icon("search")),
        menuItem("Alpha Diversity", tabName = "diversityTab", icon = icon("calculator")),
        menuItem("Rarefaction", tabName = "rarefactionTab", icon = icon("chart-line")),
        menuItem("PCA / Beta", tabName = "pcaTab", icon = icon("project-diagram")),
        menuItem("LifemapR Tree", tabName = "lifemapTab", icon = icon("tree")),
        tags$li(class = "header", "Tables"),
        menuItem("Data Table", tabName = "tableTab", icon = icon("database")),
        tags$li(class = "header", "Settings"),
        menuItem("Load & Export Settings", tabName = "settingsTab",
                 icon = icon("save"))
      )
    }
  })

  # ═══════════════════════════════════════════════════════════════════════
  # DATA LOADING
  # ═══════════════════════════════════════════════════════════════════════

  # Shared loader used by the file upload button, the "Load Directory" button,
  # and the ?folder=... URL auto-load. When `skip_unknown` is TRUE, files for
  # which detect_format() returns "unknown" are silently skipped (used by the
  # URL handler where no user-selected format is available).
  parse_files_and_update <- function(file_paths, file_names,
                                     fallback_format, skip_unknown = FALSE,
                                     sample_names = NULL) {
    if (length(file_paths) == 0) {
      showNotification("✗ No matching files found.",
                       type = "error", duration = 8)
      return(invisible(FALSE))
    }

    start_time <- Sys.time()
    all_data <- list()

    withProgress(message = "Loading data...", value = 0, {
      for (i in seq_along(file_paths)) {
        incProgress(1 / length(file_paths), detail = basename(file_names[i]))
        fpath <- file_paths[i]
        # sample_names (when provided by the URL/dataset.tsv loader) wins over
        # the filename-derived fallback so the Krona dropdown / metadata join
        # use the SUSHI `Name` column.
        sname <- if (!is.null(sample_names) && i <= length(sample_names) &&
                     !is.na(sample_names[i]) && nzchar(sample_names[i])) {
          sample_names[i]
        } else {
          sub("\\.[^.]+$", "", file_names[i])
        }

        detected <- detect_format(fpath)
        use_format <- if (detected != "unknown") {
          detected
        } else if (skip_unknown) {
          NULL
        } else {
          fallback_format
        }
        if (is.null(use_format)) next

        parsed <- switch(use_format,
          "kraken2"          = parse_kraken2_report(fpath, sname),
          "krakenuniq"       = parse_krakenuniq_report(fpath, sname),
          "bracken"          = parse_bracken_output(fpath, sname),
          "combined_bracken" = parse_combined_bracken(fpath),
          "metaphlan"        = parse_metaphlan_output(fpath, sname, merged = FALSE),
          "merged_metaphlan" = parse_metaphlan_output(fpath, sname, merged = TRUE),
          NULL
        )

        if (!is.null(parsed) && nrow(parsed) > 0) {
          all_data[[length(all_data) + 1]] <- parsed
        }
      }
    })

    if (length(all_data) == 0) {
      showNotification("✗ No data could be parsed. Check file format.",
                       type = "error", duration = 8)
      return(invisible(FALSE))
    }

    # Normalise columns across all parsed data frames before binding.
    # Different formats (or Kraken2 with/without lineage cols) may have
    # different column sets. Ensure all share the same columns.
    all_cols <- unique(unlist(lapply(all_data, colnames)))
    all_data <- lapply(all_data, function(d) {
      missing <- setdiff(all_cols, colnames(d))
      for (m in missing) d[[m]] <- NA
      d[, all_cols, drop = FALSE]
    })
    combined <- tryCatch(
      as.data.frame(data.table::rbindlist(all_data, use.names = TRUE, fill = TRUE)),
      error = function(e) {
        showNotification(paste("Error combining data:", e$message),
                         type = "error", duration = 10)
        NULL
      }
    )
    if (is.null(combined)) return(invisible(FALSE))

    rv$raw_data <- combined
    rv$loaded <- TRUE

    samples <- sort(unique(combined$Sample))
    rv$sample_renames <- setNames(samples, samples)

    # If the current rank filter has zero rows in the newly loaded data,
    # fall back to the deepest rank that IS present. This is the common case
    # for genus-level Bracken reports (no "S" rows) auto-loaded from a
    # ?data=... Live Report link — otherwise every plot renders blank.
    present_ranks <- unique(combined$rank)
    current_rank  <- isolate(input$rank_select)
    if (!is.null(current_rank) && !(current_rank %in% present_ranks)) {
      priority <- c("S", "G", "F", "O", "C", "P", "K", "D")
      fallback <- priority[priority %in% present_ranks][1]
      if (!is.na(fallback)) {
        updateSelectInput(session, "rank_select", selected = fallback)
      }
    }

    end_time <- Sys.time()
    load_time_val <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
    rv$load_time <- load_time_val

    showNotification(
      paste0("✓ Loaded ", length(unique(combined$Sample)), " sample(s), ",
             nrow(combined), " rows in ", load_time_val, "s"),
      type = "message", duration = 5
    )
    invisible(TRUE)
  }

  # Toggle the underlying <input type="file"> between file-picker and
  # folder-picker mode via the webkitdirectory attribute. The same Browse
  # button is reused — only the OS dialog changes.
  observeEvent(input$upload_folder, {
    js <- if (isTRUE(input$upload_folder)) {
      "var el = document.getElementById('data_files'); if (el) { el.setAttribute('webkitdirectory',''); el.setAttribute('directory',''); el.setAttribute('mozdirectory',''); }"
    } else {
      "var el = document.getElementById('data_files'); if (el) { el.removeAttribute('webkitdirectory'); el.removeAttribute('directory'); el.removeAttribute('mozdirectory'); }"
    }
    shinyjs::runjs(js)
  }, ignoreInit = FALSE)

  observeEvent(input$load_data, {
    req(input$data_files)
    files <- input$data_files
    # A fresh manual load supersedes any selections queued from a settings file.
    rv$restore_inputs <- list()
    rv$pending_inputs <- list()

    # Folder-mode: a webkitdirectory browse dumps every file in the picked
    # folder into input$data_files. Pre-filter to filenames matching the
    # selected radio format's patterns (e.g. *.report.txt for Kraken2).
    if (isTRUE(input$upload_folder)) {
      keep <- filter_names_for_format(files$name, input$input_format)
      n_total <- nrow(files)
      n_kept  <- sum(keep)
      if (n_kept == 0) {
        showNotification(
          paste0("No files matching '", input$input_format,
                 "' patterns in the selected folder (", n_total, " file(s) inspected)."),
          type = "warning", duration = 8
        )
        return()
      }
      files <- files[keep, , drop = FALSE]
      if (n_kept < n_total) {
        showNotification(
          paste0("Loading ", n_kept, " of ", n_total,
                 " files (filtered by '", input$input_format, "' patterns)."),
          type = "message", duration = 5
        )
      }
    }

    # ── HUMAnN branch: build rv$humann from per-sample TSVs directly ──
    if (identical(input$input_format, "humann")) {
      session$sendCustomMessage("appLoading",
        list(state = "start", text = "Loading HUMAnN reports…"))
      on.exit(session$sendCustomMessage("appLoading", "stop"), add = TRUE)
      humann <- withProgress(
        message = "Loading HUMAnN reports…",
        detail  = "Grouping per-sample TSVs by table",
        value   = 0.1,
        {
          res <- parse_humann_uploaded_files(files$datapath, files$name)
          setProgress(value = 0.9, detail = "Assembling reactive state")
          res
        }
      )
      if (is.null(humann) || length(humann) == 0) {
        showNotification(
          paste0("✗ No HUMAnN per-sample TSVs recognised. Expected filenames ",
                 "like *_pathabundance.tsv, *_reactions.tsv, *_ko_cpm.tsv, ",
                 "*_genefamilies*.tsv."),
          type = "error", duration = 10
        )
        return()
      }
      all_samples <- unique(unlist(
        lapply(humann, function(df) setdiff(colnames(df), "Feature")),
        use.names = FALSE
      ))
      rv$humann            <- humann
      rv$humann_parsed_ds  <- NULL   # uploads are fully parsed in memory already
      rv$humann_modes      <- setNames(rep(NA_character_, length(all_samples)),
                                       all_samples)
      rv$humann_source_dir <- NA_character_
      rv$sample_renames    <- setNames(all_samples, all_samples)
      rv$loaded            <- TRUE
      rv$app_mode          <- "humann"
      showNotification(
        sprintf("HUMAnN dataset loaded: %d sample(s), %d table(s)",
                length(all_samples), length(humann)),
        type = "message", duration = 5
      )
      return()
    }

    parse_files_and_update(
      file_paths       = files$datapath,
      file_names       = files$name,
      fallback_format  = input$input_format,
      skip_unknown     = FALSE
    )
  })

  # ── URL parameter ?data=<path> (or ?folder=<path>) auto-load. Runs once at
  # session start. Server-side path only — matches the exploreDE convention:
  # a relative value like "p1234/Kraken_2025-.../" is resolved under
  # /srv/gstore/projects first, then under the course gstore root. Absolute
  # paths are honored as-is. Format is auto-detected per file; unrecognised
  # files are skipped (no radio-button fallback in URL mode).
  session$onFlushed(function() {
    query <- parseQueryString(isolate(session$clientData$url_search))
    raw_path <- query$data %||% query$folder
    if (is.null(raw_path) || !nzchar(raw_path)) return()

    # Immediate feedback while LDAP + dataset.tsv parse run silently (they
    # precede withProgress). Cleared on error return or once loading starts.
    showNotification(
      "Loading dataset from URL — resolving path and manifest…",
      id = "url_boot", type = "default", duration = NULL, closeButton = FALSE
    )
    on.exit(removeNotification("url_boot"), add = TRUE)

    resolved <- resolve_gstore_folder(raw_path)
    if (is.null(resolved)) {
      showNotification(
        paste0("URL folder not found under known gstore roots: ", raw_path),
        type = "error", duration = 10
      )
      return()
    }

    # ── LDAP access control (mirrors exploreDE) ──
    ldap <- ldap_access_check(raw_path, resolved)
    message(sprintf("[exploreMetaTax] user=%s role=%s allowed=%s reason=%s",
                    Sys.getenv("SHINYPROXY_USERNAME"),
                    ldap$role, ldap$allowed, ldap$reason))
    if (!isTRUE(ldap$allowed)) {
      showModal(modalDialog(
        title = "Access denied",
        paste0("You are not authorised to view this dataset.\n",
               "User: ", Sys.getenv("SHINYPROXY_USERNAME"), "\n",
               "Role: ", ldap$role, "\n",
               "Reason: ", ldap$reason, "\n\n",
               "If you believe this is a mistake, contact FGCZ."),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return()
    }

    # ── dataset.tsv-driven load (SUSHI convention) ──
    # We deliberately do NOT glob for *.report.txt anymore. The canonical
    # source of truth for "what reports belong to this analysis" is the
    # dataset.tsv at the root of the project folder. We hard-error if it
    # is missing or if any KrakenReport path doesn't resolve — silently
    # falling back to a glob hid metadata mismatches and double-counted
    # stray .report.txt files left behind by previous runs.
    parsed_ds <- read_sushi_dataset_tsv(resolved)
    if (is.null(parsed_ds)) {
      showModal(modalDialog(
        title = "dataset.tsv not found",
        sprintf(paste0("No readable dataset.tsv at:\n  %s\n\n",
                       "URL-based loading requires a SUSHI dataset.tsv ",
                       "manifest. Re-run the upstream Kraken SUSHI app, ",
                       "or upload reports manually via the Upload tab."),
                file.path(resolved, "dataset.tsv")),
        easyClose = TRUE, footer = modalButton("Close")
      ))
      return()
    }

    # ── HUMAnN branch: dataset.tsv carries per-sample functional tables ──
    # We short-circuit here so the rest of this handler (which is wired up
    # around the kraken/metaphlan Report/Profile column) doesn't run.
    if (identical(parsed_ds$app_mode, "humann")) {
      ds <- parsed_ds$raw
      sample_names <- as.character(ds[[parsed_ds$name_col]])

      # Wrap the on-disk join in withProgress so the browser shows a bar
      # while the per-sample TSVs are read (can take several seconds on a
      # large gene_families table).
      session$sendCustomMessage("appLoading",
        list(state = "start", text = "Loading HUMAnN reports…"))
      on.exit(session$sendCustomMessage("appLoading", "stop"), add = TRUE)
      humann <- withProgress(
        message = "Loading HUMAnN reports…",
        detail  = "Reading per-sample TSVs and joining tables",
        value   = 0.1,
        {
          # Gene families (the large ~1.9M-row table) is loaded lazily and
          # asynchronously the first time its tab is opened — see the gf_task
          # ExtendedTask below. Reading only the small tables here keeps the
          # initial open fast even for 40+ sample runs.
          res <- load_humann_tables(
            parsed_ds,
            slots_wanted = c("pathways", "reactions", "kegg_kos"))
          setProgress(value = 0.9, detail = "Assembling reactive state")
          res
        }
      )
      if (length(humann) == 0) {
        showModal(modalDialog(
          title = "No HUMAnN tables could be loaded",
          paste0("dataset.tsv references HUMAnN [File] columns but none ",
                 "of the per-sample TSVs resolved on disk. Check that the ",
                 "result_dir paths are still valid."),
          easyClose = TRUE, footer = modalButton("Close")
        ))
        return()
      }

      # Per-sample mode. Prefer Mode [File] (canonical, written by
      # HUMAnNApp) → read the file for each sample. Fall back to the
      # older Mode [Characteristic] inline column. Otherwise NA → the
      # Info tab shows "unknown".
      modes_vec <- rep(NA_character_, length(sample_names))
      names(modes_vec) <- sample_names
      if (!is.na(parsed_ds$mode_file_col) && nzchar(parsed_ds$mode_file_col)) {
        mode_rels  <- as.character(ds[[parsed_ds$mode_file_col]])
        mode_paths <- vapply(mode_rels, resolve_gstore_file,
                             FUN.VALUE = NA_character_, USE.NAMES = FALSE)
        for (i in seq_along(mode_paths)) {
          p <- mode_paths[i]
          if (!is.na(p) && file.exists(p)) {
            v <- tryCatch(trimws(readLines(p, n = 1, warn = FALSE)[1]),
                          error = function(e) NA_character_)
            if (!is.na(v) && nzchar(v)) modes_vec[i] <- v
          }
        }
      } else if (!is.na(parsed_ds$mode_col) && nzchar(parsed_ds$mode_col)) {
        m <- as.character(ds[[parsed_ds$mode_col]])
        modes_vec[seq_along(m)] <- m
      }

      # Metadata cleanup (same logic as the kraken/metaphlan branch below)
      cn <- colnames(ds)
      is_path_col <- grepl("\\[(File|Link)\\]\\s*$", cn)
      is_name_col <- cn == parsed_ds$name_col
      is_mode_col <- (!is.na(parsed_ds$mode_col)      & cn == parsed_ds$mode_col) |
                    (!is.na(parsed_ds$mode_file_col) & cn == parsed_ds$mode_file_col)
      meta_cols   <- cn[!is_path_col & !is_name_col & !is_mode_col]
      if (length(meta_cols) > 0) {
        meta_df <- ds[, meta_cols, drop = FALSE]
        colnames(meta_df) <- trimws(sub("\\s*\\[[^]]+\\]\\s*$", "",
                                        colnames(meta_df)))
        meta_df <- cbind(Sample = sample_names, meta_df, stringsAsFactors = FALSE)
        keep <- vapply(meta_df, function(x) !all(is.na(x) | !nzchar(as.character(x))),
                       logical(1))
        rv$metadata <- meta_df[, keep, drop = FALSE]
      }

      rv$humann            <- humann
      rv$humann_parsed_ds  <- parsed_ds
      rv$humann_modes      <- modes_vec
      rv$humann_source_dir <- resolved
      rv$sample_renames    <- setNames(sample_names, sample_names)
      rv$loaded            <- TRUE
      rv$app_mode          <- "humann"
      showNotification(
        sprintf("HUMAnN dataset loaded: %d sample(s), %d table(s)",
                length(sample_names),
                sum(!vapply(humann, is.null, logical(1)))),
        type = "message", duration = 5
      )
      return()
    }

    if (is.na(parsed_ds$report_col) || !nzchar(parsed_ds$report_col)) {
      showModal(modalDialog(
        title = "No report column in dataset.tsv",
        sprintf(paste0("dataset.tsv at %s has no column matching ",
                       "'<Tool>Report [File]' or '<Tool>Profile [File]' ",
                       "(e.g. 'KrakenReport [File]', 'MetaPhlAnProfile [File]'). ",
                       "Columns present: %s"),
                file.path(resolved, "dataset.tsv"),
                paste(colnames(parsed_ds$raw), collapse = ", ")),
        easyClose = TRUE, footer = modalButton("Close")
      ))
      return()
    }
    if (is.na(parsed_ds$name_col)) {
      showModal(modalDialog(
        title = "No Name column in dataset.tsv",
        sprintf("dataset.tsv at %s is missing the required 'Name' column.",
                file.path(resolved, "dataset.tsv")),
        easyClose = TRUE, footer = modalButton("Close")
      ))
      return()
    }

    ds <- parsed_ds$raw
    report_rel    <- as.character(ds[[parsed_ds$report_col]])
    sample_names  <- as.character(ds[[parsed_ds$name_col]])
    report_paths  <- vapply(report_rel, resolve_gstore_file,
                            FUN.VALUE = NA_character_, USE.NAMES = FALSE)

    # Hard-error if any listed report file is missing on disk — better to
    # tell the user exactly which rows broke than to silently load a subset.
    missing_idx <- which(is.na(report_paths))
    if (length(missing_idx) > 0) {
      detail <- paste(sprintf("  - %s  (row %d, sample '%s')",
                              report_rel[missing_idx],
                              missing_idx,
                              sample_names[missing_idx]),
                      collapse = "\n")
      showModal(modalDialog(
        title = "Missing report files",
        tags$div(
          tags$p(sprintf("%d of %d report file(s) listed in dataset.tsv ",
                         length(missing_idx), nrow(ds)),
                 "could not be found on disk under any known gstore root:"),
          tags$pre(detail)
        ),
        easyClose = TRUE, footer = modalButton("Close"),
        size = "l"
      ))
      return()
    }

    # ── Pull non-file/non-link columns as sample metadata ──
    # SUSHI annotates types in the header, e.g. "Condition [Factor]" /
    # "Sample Id [B-Fabric]". Strip those suffixes for cleaner column names
    # in the diversity / PCA Color By dropdowns. Exclude Name (it becomes
    # the keying `Sample` column) and any [File]/[Link] columns (those are
    # data paths, not phenotype).
    cn <- colnames(ds)
    is_path_col <- grepl("\\[(File|Link)\\]\\s*$", cn)
    is_name_col <- cn == parsed_ds$name_col
    meta_cols   <- cn[!is_path_col & !is_name_col]
    if (length(meta_cols) > 0) {
      meta_df <- ds[, meta_cols, drop = FALSE]
      colnames(meta_df) <- trimws(sub("\\s*\\[[^]]+\\]\\s*$", "",
                                      colnames(meta_df)))
      meta_df <- cbind(Sample = sample_names, meta_df, stringsAsFactors = FALSE)
      # Drop columns that are entirely empty after the strip
      keep <- vapply(meta_df, function(x) !all(is.na(x) | !nzchar(as.character(x))),
                     logical(1))
      rv$metadata <- meta_df[, keep, drop = FALSE]
    }

    parse_files_and_update(
      file_paths      = report_paths,
      file_names      = basename(report_paths),
      fallback_format = NULL,
      skip_unknown    = TRUE,
      sample_names    = sample_names
    )
  }, once = TRUE)

  # ── Status text ──
  output$status_text <- renderUI({
    if (rv$loaded) {
      n_samples <- length(unique(rv$raw_data$Sample))
      n_rows <- nrow(rv$raw_data)
      div(class = "status-text",
          icon("check-circle"),
          sprintf(" %d sample(s) loaded | %s rows | in %s sec",
                  n_samples, format(n_rows, big.mark = ","), rv$load_time))
    }
  })

  # ═══════════════════════════════════════════════════════════════════════
  # METADATA LOADING
  # ═══════════════════════════════════════════════════════════════════════

  observeEvent(input$load_metadata, {
    req(input$metadata_file)
    # Fresh metadata supersedes any group/colour selections queued from a
    # settings file (their columns may no longer exist).
    rv$restore_inputs <- list()
    rv$pending_inputs <- list()
    fpath <- input$metadata_file$datapath
    fname <- input$metadata_file$name

    sep <- if (grepl("\\.csv$", fname, ignore.case = TRUE)) "," else "\t"

    meta <- tryCatch(
      read.table(fpath, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                 fill = TRUE, check.names = FALSE),
      error = function(e) NULL
    )

    if (!is.null(meta) && "Sample" %in% colnames(meta)) {
      rv$metadata <- meta
      showNotification(
        paste0("✓ Metadata loaded: ", nrow(meta), " rows, ",
               ncol(meta), " columns"),
        type = "message", duration = 5
      )
    } else {
      showNotification(
        "✗ Metadata must have a 'Sample' column. Check format.",
        type = "error", duration = 8)
    }
  })

  output$metadata_status <- renderUI({
    if (!is.null(rv$metadata)) {
      div(class = "status-text",
          icon("check-circle"),
          sprintf(" Metadata: %d samples × %d columns",
                  nrow(rv$metadata), ncol(rv$metadata) - 1))
    }
  })

  # ── Download example metadata template ──
  output$download_example_metadata <- downloadHandler(
    filename = function() { "example_metadata.tsv" },
    content = function(file) {
      example_path <- file.path(dirname(getwd()), "example_metadata.tsv")
      # Try app directory first, then fall back to inline content
      if (file.exists("example_metadata.tsv")) {
        file.copy("example_metadata.tsv", file)
      } else if (file.exists(example_path)) {
        file.copy(example_path, file)
      } else {
        # Generate inline
        writeLines(
          c("Sample\tGroup\tTreatment\tReplicate",
            "SampleA\tControl\tNone\t1",
            "SampleB\tControl\tNone\t2",
            "SampleC\tTreatment1\tDrugX\t1",
            "SampleD\tTreatment1\tDrugX\t2",
            "SampleE\tTreatment2\tDrugY\t1",
            "SampleF\tTreatment2\tDrugY\t2"),
          file
        )
      }
    }
  )

  # ═══════════════════════════════════════════════════════════════════════
  # SAMPLE RENAMING
  # ═══════════════════════════════════════════════════════════════════════

  output$rename_ui <- renderUI({
    req(rv$loaded, rv$sample_renames)
    renames <- rv$sample_renames
    n <- length(renames)

    if (n == 0) return(div("No samples loaded."))

    inputs <- lapply(seq_along(renames), function(i) {
      old <- names(renames)[i]
      new_val <- renames[i]
      fluidRow(
        column(6, tags$small(tags$b(old))),
        column(6, textInput(
          inputId = paste0("rename_", i),
          label = NULL,
          value = new_val,
          width = "100%"
        ))
      )
    })

    tagList(
      tags$div(style = "max-height: 300px; overflow-y: auto; padding-right: 6px;",
        inputs
      ),
      actionButton("apply_renames", "Apply Names", icon = icon("check"),
                   class = "btn-default", width = "100%",
                   style = "margin-top: 8px;")
    )
  })

  observeEvent(input$apply_renames, {
    req(rv$sample_renames)
    renames <- rv$sample_renames
    n <- length(renames)

    new_names <- sapply(seq_len(n), function(i) {
      val <- input[[paste0("rename_", i)]]
      if (is.null(val) || val == "") names(renames)[i] else val
    })

    rv$sample_renames <- setNames(new_names, names(renames))
    showNotification("✓ Sample names updated", type = "message", duration = 3)
  })

  # ═══════════════════════════════════════════════════════════════════════
  # REACTIVE: Renamed + Filtered Data
  # ═══════════════════════════════════════════════════════════════════════

  renamed_data <- reactive({
    req(rv$raw_data, rv$sample_renames)
    df <- rv$raw_data
    renames <- rv$sample_renames
    df$Sample <- renames[df$Sample]
    df
  })

  # Merged with metadata
  merged_data <- reactive({
    df <- renamed_data()
    if (!is.null(rv$metadata)) {
      meta <- rv$metadata
      df <- merge(df, meta, by = "Sample", all.x = TRUE)
    }
    df
  })

  # Filtered data (rank, abundance, sample selection)
  filtered_data <- reactive({
    df <- merged_data()
    req(nrow(df) > 0)

    # Pre-allocate a single logical vector for highly memory-efficient filtering
    keep <- rep(TRUE, nrow(df))

    # Filter by rank
    rank_val <- input$rank_select
    if (!is.null(rank_val) && rank_val != "" && "rank" %in% colnames(df)) {
      keep <- keep & (df$rank == rank_val | df$name == "unclassified")
    }

    # Sample selection + order. With the drag buckets the "Visible" list gives
    # both which samples to keep and the order to show them in; otherwise the
    # checkbox group's selection is used (order = original sort). Empty/NULL
    # means "no explicit selection" → keep all, matching the old behaviour.
    selected_samples <- if (HAS_SORTABLE) input$samples_visible else input$sample_filter
    if (!is.null(selected_samples) && length(selected_samples) > 0) {
      keep <- keep & (df$Sample %in% selected_samples)
    }

    # Filter by minimum abundance (%)
    if (!is.null(input$min_abundance) && input$min_abundance > 0) {
      keep <- keep & (df$percent >= input$min_abundance | df$name == "unclassified")
    }

    # Filter by minimum absolute counts
    if (!is.null(input$min_counts) && input$min_counts > 0 &&
        "reads_clade" %in% colnames(df)) {
      keep <- keep & (df$reads_clade >= input$min_counts | df$name == "unclassified")
    }

    # Filter by minimizer counts (only present in Kraken2 --report-minimizer-data)
    if (!is.null(input$min_minimizers_clade) && input$min_minimizers_clade > 0 &&
        "minimizers_clade" %in% colnames(df)) {
      keep <- keep & (df$minimizers_clade >= input$min_minimizers_clade |
                df$name == "unclassified")
    }
    if (!is.null(input$min_minimizers_taxon) && input$min_minimizers_taxon > 0 &&
        "minimizers_taxon" %in% colnames(df)) {
      keep <- keep & (df$minimizers_taxon >= input$min_minimizers_taxon |
                df$name == "unclassified")
    }

    # Filter by KrakenUniq unique k-mers (only present in KrakenUniq reports)
    if (!is.null(input$min_kmers) && input$min_kmers > 0 &&
        "kmers" %in% colnames(df)) {
      keep <- keep & (df$kmers >= input$min_kmers | df$name == "unclassified")
    }

    # Filter by KrakenUniq genome coverage
    if (!is.null(input$min_cov) && input$min_cov > 0 &&
        "cov" %in% colnames(df)) {
      keep <- keep & (df$cov >= input$min_cov | df$name == "unclassified")
    }

    # ── Discard filter (Upload/Filters → "Entries to discard") ──
    # Removes any row whose taxon name is in the committed discard set.
    # Computed globally (across all ranks) from rv$raw_data so a query like
    # "Homo" hits both Genus rows and Species rows underneath it. Data Table
    # bypasses filtered_data() and stays unaffected.
    disc <- discard_names()
    if (length(disc) > 0) {
      keep <- keep & !(df$name %in% disc)
    }

    # Subset the huge dataframe exactly once
    out <- df[keep, , drop = FALSE]

    # Propagate the chosen sample order: make Sample an ordered factor whose
    # levels follow the Visible bucket. Every downstream plot/table that puts
    # Sample on an axis or in a legend then respects the drag order (tabs that
    # deliberately re-sort, e.g. Alpha Diversity by value or the clustered
    # Heatmap, still override this on purpose).
    if (!is.null(selected_samples) && length(selected_samples) > 0) {
      lev <- selected_samples[selected_samples %in% out$Sample]
      if (length(lev) > 0) out$Sample <- factor(out$Sample, levels = lev)
    }
    out
  })

  # ── "Entries to discard" (Upload/Filters box) ──
  # Same cascading matcher as the Organisms-of-interest search, but applied
  # against rv$raw_data (all ranks + full lineage) so a single query removes
  # its matches globally. Only committed on the Apply button so the user isn't
  # penalised while typing. The Data Table renderer uses merged_data() and
  # is intentionally NOT filtered by this.
  discard_queries <- reactive({
    txt <- input$discard_query_text
    if (is.null(txt) || !nzchar(txt)) return(character(0))
    qs <- unlist(strsplit(txt, "[,\n]", perl = TRUE), use.names = FALSE)
    qs <- trimws(qs); qs[nzchar(qs)]
  })

  observeEvent(input$discard_apply_btn, {
    rv$discard_committed <- discard_queries()
    showNotification(sprintf("✓ Discard filter applied (%d quer%s).",
                             length(rv$discard_committed),
                             if (length(rv$discard_committed) == 1) "y" else "ies"),
                     type = "message", duration = 3)
  })

  observeEvent(input$discard_clear_btn, {
    updateTextAreaInput(session, "discard_query_text", value = "")
    rv$discard_committed <- character(0)
  })

  # Resolve committed queries into a concrete taxon-name blacklist.
  # Cascade (per query): exact name → name substring → Genus-lineage substring.
  discard_names <- reactive({
    qs <- rv$discard_committed
    if (length(qs) == 0) return(character(0))
    if (is.null(rv$raw_data)) return(character(0))
    df <- rv$raw_data
    cols_want <- intersect(c("name", "Genus", "Species", "Phylum", "Family", "rank"),
                           colnames(df))
    pool <- unique(df[, cols_want, drop = FALSE])
    pool$.name_nrm  <- .organism_nrm(pool$name)
    pool$.genus_nrm <- if ("Genus" %in% colnames(pool)) .organism_nrm(pool$Genus) else ""
    hit <- character(0)
    for (q in qs) {
      qn <- .organism_nrm(q)
      if (!nzchar(qn)) next
      m <- pool[pool$.name_nrm == qn, , drop = FALSE]
      if (nrow(m) == 0)
        m <- pool[grepl(qn, pool$.name_nrm, fixed = TRUE), , drop = FALSE]
      if (nrow(m) == 0 && "Genus" %in% colnames(pool)) {
        first <- strsplit(qn, " ", fixed = TRUE)[[1]][1]
        if (length(first) && nzchar(first)) {
          m <- pool[grepl(first, pool$.genus_nrm, fixed = TRUE), , drop = FALSE]
        }
      }
      if (nrow(m) > 0) hit <- c(hit, m$name)
    }
    unique(hit)
  })

  # Compact status line under the discard input: how many queries committed
  # and how many taxa they resolve to. Uses the same nzchar/status-text style
  # as the other Upload-tab widgets.
  output$discard_summary <- renderUI({
    qs <- rv$discard_committed
    if (length(qs) == 0) return(NULL)
    nms <- discard_names()
    unmatched <- length(qs) - sum(vapply(qs, function(q) {
      any(grepl(.organism_nrm(q),
                .organism_nrm(nms), fixed = TRUE))
    }, logical(1)))
    class_extra <- if (length(nms) == 0) " status-warning" else ""
    header_line <- sprintf(" Discarding %d taxon%s across %d quer%s%s.",
                           length(nms), if (length(nms) == 1) "" else "s",
                           length(qs), if (length(qs) == 1) "y" else "ies",
                           if (unmatched > 0)
                             sprintf(" — %d quer%s had no matches",
                                     unmatched, if (unmatched == 1) "y" else "ies")
                           else "")
    div(class = paste0("status-text", class_extra),
        icon("filter"), header_line,
        # Show the concrete taxa removed below the header. Sorted for stable
        # display; scrollable when the list is long so the Filters box doesn't
        # grow unbounded.
        if (length(nms) > 0) tagList(
          tags$div(style = "margin-top: 6px; font-weight: 600;", "Removed taxa:"),
          tags$div(
            style = paste("max-height: 140px; overflow-y: auto;",
                          "padding: 4px 8px; margin-top: 2px;",
                          "background: rgba(255,255,255,0.6); border-radius: 4px;",
                          "font-family: monospace; font-size: 12px;"),
            paste(sort(nms), collapse = ", ")
          )
        )
    )
  })

  # ── Dynamic UI: minimizer / KrakenUniq filter inputs (shown only if data has those columns) ──
  output$minimizer_filters_ui <- renderUI({
    req(rv$loaded, rv$raw_data)
    df <- rv$raw_data
    has_min_clade <- "minimizers_clade" %in% colnames(df)
    has_min_taxon <- "minimizers_taxon" %in% colnames(df)
    has_kmers <- "kmers" %in% colnames(df)
    has_cov   <- "cov" %in% colnames(df)

    ui_elements <- list()

    # Kraken2 minimizer filters
    if (has_min_clade || has_min_taxon) {
      ui_elements <- c(ui_elements, list(
        hr(),
        helpText(icon("info-circle"),
                 "Minimizer columns detected (--report-minimizer-data)."),
        fluidRow(
          if (has_min_clade) column(6,
            numericInput("min_minimizers_clade", "Min Minimizers (clade):",
                         value = restored_input("min_minimizers_clade", 0),
                         min = 0, step = 1)
          ),
          if (has_min_taxon) column(6,
            numericInput("min_minimizers_taxon", "Min Unique Minimizers (taxon):",
                         value = restored_input("min_minimizers_taxon", 0),
                         min = 0, step = 1),
            actionLink("kmer_help", label = NULL,
                       icon = icon("question-circle"),
                       style = "color: #C4692A; font-size: 14px; margin-top: -10px;")
          )
        )
      ))
    }

    # KrakenUniq-specific filters
    if (has_kmers || has_cov) {
      ui_elements <- c(ui_elements, list(
        hr(),
        helpText(icon("info-circle"),
                 "KrakenUniq columns detected (unique k-mers, genome coverage)."),
        fluidRow(
          if (has_kmers) column(6,
            numericInput("min_kmers", "Min Unique K-mers:",
                         value = restored_input("min_kmers", 0),
                         min = 0, step = 1),
            actionLink("krakenuniq_kmer_help", label = NULL,
                       icon = icon("question-circle"),
                       style = "color: #C4692A; font-size: 14px; margin-top: -10px;")
          ),
          if (has_cov) column(6,
            numericInput("min_cov", "Min Genome Coverage:",
                         value = restored_input("min_cov", 0),
                         min = 0, step = 0.0001)
          )
        )
      ))
    }

    if (length(ui_elements) > 0) do.call(tagList, ui_elements)
  })

  # Modal for k-mer help
  observeEvent(input$kmer_help, {
    showModal(modalDialog(
      title = "Why filter by unique k-mers?",
      HTML(paste0(
        "<p>Each read of length <b>L</b> can produce up to <b>L &minus; k + 1</b> ",
        "distinct k-mers (e.g. ~70 for 100 bp reads with k=31).</p>",
        "<p>If a taxon assignment has far fewer unique k-mers than expected, it likely reflects ",
        "<b>contamination</b> or <b>low-complexity sequence</b> rather than a true match.</p>",
        "<h4>Rule of thumb</h4>",
        "<p>The unique k-mer count should be at least <b>~5&times; the read count</b> ",
        "for a confident classification.</p>",
        "<p>Please also consider an absolute threshold.</p>",
        "<h4>Exception</h4>",
        "<p>For small genomes (e.g. viruses), the unique k-mer count cannot exceed the genome size, ",
        "so lower ratios may still be valid.</p>"
      )),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # Modal for KrakenUniq k-mer help
  observeEvent(input$krakenuniq_kmer_help, {
    showModal(modalDialog(
      title = "KrakenUniq: Why filter by unique k-mers?",
      HTML(paste0(
        "<p>KrakenUniq reports the number of <b>unique k-mers</b> mapping to each taxon, ",
        "providing a more reliable signal than read counts alone.</p>",
        "<p>A taxon with many reads but very few unique k-mers likely represents ",
        "<b>contamination</b>, <b>low-complexity sequence</b>, or <b>misclassification</b>.</p>",
        "<h4>Rule of thumb</h4>",
        "<p>Set a minimum unique k-mer threshold to remove spurious classifications. ",
        "Higher thresholds increase specificity at the cost of sensitivity.</p>",
        "<h4>Genome coverage</h4>",
        "<p>The <b>coverage</b> column estimates the fraction of the reference genome ",
        "covered by unique k-mers. Very low coverage (<0.001) for a reported taxon ",
        "suggests a false positive.</p>"
      )),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # ── Dynamic UI: sample checkboxes ──
  output$sample_checkboxes <- renderUI({
    req(rv$loaded, rv$sample_renames)
    samples <- sort(unique(unname(rv$sample_renames)))

    # Preserve saved split across re-renders. Priority: a .rds restore
    # (rv$pending_ui_state, isolated so it doesn't re-invalidate when
    # cleared post-flush), then whatever the widget last reported
    # (input$… via isolate), then defaults (all visible).
    pending <- isolate(rv$pending_ui_state)
    visible_saved <- pending$samples_visible
    hidden_saved  <- pending$samples_hidden
    if (is.null(visible_saved) && is.null(hidden_saved)) {
      visible_saved <- isolate(input$samples_visible)
      hidden_saved  <- isolate(input$samples_hidden)
    }
    if (length(intersect(c(visible_saved, hidden_saved), samples)) > 0) {
      visible <- intersect(visible_saved, samples)
      hidden  <- intersect(hidden_saved,  samples)
      # New samples default to visible, appended at the end.
      visible <- c(visible, setdiff(samples, c(visible, hidden)))
    } else {
      visible <- samples
      hidden  <- character(0)
    }

    if (HAS_SORTABLE) {
      # Two drag-and-drop buckets. Samples in "Visible" are shown in every
      # plot/table, in the order they appear here; dragging a sample to
      # "Hidden" drops it. Both lists post their current contents to shiny
      # inputs (samples_visible / samples_hidden); filtered_data() reads the
      # visible list for BOTH selection and ordering.
      sortable::bucket_list(
        header = "Drag to hide/show samples, or reorder within Visible:",
        group_name = "sample_buckets",
        orientation = "horizontal",
        sortable::add_rank_list(
          text = "Visible (order applies to plots)",
          labels = visible,
          input_id = "samples_visible",
          options = .multidrag_opts
        ),
        sortable::add_rank_list(
          text = "Hidden",
          labels = hidden,
          input_id = "samples_hidden",
          options = .multidrag_opts
        )
      )
    } else {
      sel <- pending$sample_filter
      if (is.null(sel)) sel <- visible
      sel <- intersect(sel, samples)
      if (length(sel) == 0) sel <- samples
      checkboxGroupInput("sample_filter", "Select Samples:",
                         choices = samples, selected = sel)
    }
  })

  # ── Dynamic UI: group by selector ──
  output$group_by_selector <- renderUI({
    if (!is.null(rv$metadata)) {
      group_cols <- setdiff(colnames(rv$metadata), "Sample")
      selectInput("group_by", "Group By (metadata):",
                  choices = c("None", group_cols),
                  selected = restored_input("group_by", "None",
                                            valid = c("None", group_cols)))
    }
  })

  # ── Helper: get grouping column ──
  get_group_col <- reactive({
    if (!is.null(input$group_by) && input$group_by != "None" &&
        input$group_by %in% colnames(merged_data())) {
      input$group_by
    } else {
      NULL
    }
  })

  # ── Helper: metadata factor columns available for grouping ──
  # Shared by the Composition / Heatmap / Krona "Group by" selectors. Returns
  # the metadata column names (excluding the Sample key), or character(0) when
  # no metadata was uploaded.
  metadata_group_cols <- reactive({
    if (is.null(rv$metadata)) return(character(0))
    setdiff(colnames(rv$metadata), "Sample")
  })

  # Resolve a per-tab "Group by" selectInput value to a usable metadata column,
  # or NULL when set to "None" / absent / not a real metadata column.
  resolve_group_col <- function(val) {
    if (!is.null(val) && val != "None" && val %in% metadata_group_cols()) val
    else NULL
  }

  # Uniform message shown by every tab when filtered_data() is empty after
  # the current rank / abundance / sample filters. Kept in one place so all
  # tabs read the same and users know what to change.
  empty_filter_message <- function() {
    rname <- names(rank_choices)[match(input$rank_select, rank_choices)] %||%
             input$rank_select
    sprintf(paste0("No taxa at rank '%s' passed the current filters. ",
                   "Try a different rank or relax the abundance/count filters."),
            rname)
  }

  # ═══════════════════════════════════════════════════════════════════════
  # TAB 1: COMPOSITION BARPLOT
  # ═══════════════════════════════════════════════════════════════════════

  composition_reactive <- reactive({
    df <- filtered_data()
    req(nrow(df) > 0)

    top_n <- input$top_n

    # Find global top N taxa (by mean abundance)
    clean <- df[df$name != "unclassified", , drop = FALSE]
    if (nrow(clean) > 0) {
      taxa_means <- aggregate(percent ~ name, data = clean, FUN = mean)
      taxa_means <- taxa_means[order(taxa_means$percent, decreasing = TRUE), ]
      top_taxa <- head(taxa_means$name, top_n)
    } else {
      top_taxa <- character(0)
    }

    # Group non-top taxa into "Other"
    df$name_grouped <- ifelse(
      df$name %in% c("unclassified", top_taxa),
      df$name, "Other"
    )

    # Aggregate — keep both percent and reads_clade
    group_col <- get_group_col()

    agg_cols <- "percent"
    has_reads <- "reads_clade" %in% colnames(df)
    if (has_reads) agg_cols <- c(agg_cols, "reads_clade")

    group_vars <- c("Sample", "name_grouped")
    if (!is.null(group_col)) group_vars <- c(group_vars, group_col)

    # Use aggregate with cbind to aggregate multiple value columns
    if (has_reads) {
      agg_formula <- as.formula(paste("cbind(percent, reads_clade) ~",
                                       paste(group_vars, collapse = " + ")))
    } else {
      agg_formula <- as.formula(paste("percent ~",
                                       paste(group_vars, collapse = " + ")))
    }
    plot_df <- aggregate(agg_formula, data = df, sum)

    # Order taxa by total abundance (descending) for meaningful stacking
    taxa_totals <- aggregate(percent ~ name_grouped, data = plot_df, sum)
    taxa_totals <- taxa_totals[order(taxa_totals$percent), ]  # ascending = bottom of stack gets largest
    # Keep unclassified and Other at the beginning
    special <- c("unclassified", "Other")
    ordered_taxa <- taxa_totals$name_grouped[!taxa_totals$name_grouped %in% special]
    level_order <- c(special, as.character(ordered_taxa))
    plot_df$name_grouped <- factor(plot_df$name_grouped, levels = level_order)

    # Colors
    n_colors <- length(top_taxa)
    if (n_colors > 0) {
      dyn_colors <- get_palette_colors(n_colors, input$color_palette)
      names(dyn_colors) <- top_taxa
      all_colors <- c("unclassified" = "#95a5a6", "Other" = "#d5dbdb", dyn_colors)
    } else {
      all_colors <- c("unclassified" = "#95a5a6", "Other" = "#d5dbdb")
    }

    list(plot_df = plot_df, all_colors = all_colors, group_col = group_col)
  })

  output$composition_plot <- renderPlotly({
    shiny::validate(shiny::need(nrow(filtered_data()) > 0, empty_filter_message()))
    comp <- composition_reactive()
    plot_df <- comp$plot_df
    all_colors <- comp$all_colors
    group_col <- comp$group_col
    req(nrow(plot_df) > 0)

    colors_to_use <- all_colors[names(all_colors) %in% levels(plot_df$name_grouped)]

    if (input$display_mode == "counts" && "reads_clade" %in% colnames(plot_df)) {
      y_val <- "reads_clade"
      y_lab <- "Counts"
    } else {
      y_val <- "percent"
      y_lab <- "% Reads"
    }

    # Create a standard y column to avoid .data[[]] issues
    plot_df$y_value <- plot_df[[y_val]]

    p <- ggplot(plot_df, aes(x = Sample, y = y_value, fill = name_grouped,
                              text = paste0("Sample: ", Sample,
                                            "\nTaxon: ", name_grouped,
                                            "\nValue: ", sprintf("%.2f", y_value)))) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_manual(values = colors_to_use, drop = FALSE) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 16),
        panel.grid.major.x = element_blank()
      ) +
      labs(
        title = paste("Taxonomic Composition (Rank:", names(rank_choices[rank_choices == input$rank_select]), ")"),
        x = NULL, y = y_lab, fill = "Taxon"
      )

    # Labels
    if (input$show_labels) {
      label_df <- plot_df[plot_df$y_value > 1, , drop = FALSE]
      if (nrow(label_df) > 0) {
        p <- p + geom_text(
          data = label_df,
          aes(label = sprintf("%.1f", y_value)),
          position = position_stack(vjust = 0.5),
          size = 2.8, color = "black", fontface = "bold"
        )
      }
    }

    # Faceting
    if (!is.null(group_col)) {
      p <- p + facet_grid(cols = vars(!!sym(group_col)), scales = "free_x", space = "free_x")
    }

    ggplotly(p, tooltip = "text") %>%
      layout(
        legend = list(font = list(size = 10)),
        margin = list(b = 120)
      )
  })

  # ── Composition Group comparison ──
  # Non-parametric per-taxon test across the metadata Group By, over the
  # same top-N taxa shown on the composition plot. Uses the underlying
  # long-form filtered_data() (unaggregated) so we can compare per-sample
  # abundances properly.
  composition_stats_df <- reactive({
    df <- filtered_data()
    shiny::validate(shiny::need(nrow(df) > 0, empty_filter_message()))
    gcol <- get_group_col()
    if (is.null(gcol)) return(NULL)
    clean <- df[df$name != "unclassified", , drop = FALSE]
    if (nrow(clean) == 0) return(NULL)
    taxa_means <- aggregate(percent ~ name, data = clean, FUN = mean)
    taxa_means <- taxa_means[order(taxa_means$percent, decreasing = TRUE), ]
    top_taxa <- head(taxa_means$name, input$top_n)
    val_col <- if (identical(input$display_mode, "counts") &&
                    "reads_clade" %in% colnames(clean)) "reads_clade" else "percent"
    per_taxon_group_table(clean, top_taxa, gcol, value_col = val_col)
  })

  output$composition_stats_status <- renderUI({
    gcol <- get_group_col()
    if (is.null(gcol)) {
      return(div(class = "status-text status-warning",
                 icon("info-circle"),
                 " Choose a Group By value above to enable per-taxon comparisons."))
    }
    res <- tryCatch(composition_stats_df(), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) {
      return(div(class = "status-text status-warning",
                 icon("exclamation-triangle"),
                 " No comparable groups (need ≥2 samples per group)."))
    }
    NULL
  })

  output$composition_stats_table <- renderDT({
    res <- composition_stats_df()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Nothing to report — set a Group By with ≥2 samples per group."))
    show <- res
    for (col in c("statistic", "p", "p_adj")) {
      if (col %in% colnames(show))
        show[[col]] <- signif(as.numeric(show[[col]]), 3)
    }
    datatable(
      show, rownames = FALSE, filter = "top",
      extensions = "Buttons",
      options = list(pageLength = 10, dom = "Blfrtip",
                     buttons = list(
                       list(extend = "copyHtml5",  exportOptions = dt_export_options),
                       list(extend = "csvHtml5",   exportOptions = dt_export_options),
                       list(extend = "excelHtml5", exportOptions = dt_export_options)
                     ))
    )
  })

  output$download_composition_stats_csv <- downloadHandler(
    filename = function() paste0("composition_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- composition_stats_df()
      if (is.null(res)) res <- data.frame()
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ═══════════════════════════════════════════════════════════════════════
  # TAB: ORGANISMS OF INTEREST
  # Cascading match (exact -> substring -> genus-fallback) on user-supplied
  # taxon names, producing a composition-style barplot and a DT of the
  # matched rows. Mirrors the "Organisms of interest" section in the HTML
  # differential-abundance reports.
  # ═══════════════════════════════════════════════════════════════════════

  # Normalise a string for case- and punctuation-insensitive matching.
  .organism_nrm <- function(s) {
    s <- as.character(s); s[is.na(s)] <- ""
    s <- tolower(s)
    s <- gsub("[^a-z0-9 ]", " ", s, perl = TRUE)
    s <- gsub("\\s+", " ", s, perl = TRUE)
    trimws(s)
  }

  # Split the textarea content into clean queries (newline + comma separated).
  organism_queries <- reactive({
    txt <- input$organism_query_text
    if (is.null(txt) || !nzchar(txt)) return(character(0))
    qs <- unlist(strsplit(txt, "[,\n]", perl = TRUE), use.names = FALSE)
    qs <- trimws(qs); qs[nzchar(qs)]
  })

  # Clear button: wipe the input and reset the output.
  observeEvent(input$organism_clear_btn, {
    updateTextAreaInput(session, "organism_query_text", value = "")
  })

  # External search buttons (NCBI, Scholar, PubMed, ...). With one query the
  # button is a direct link; with multiple it becomes a dropdown listing one
  # anchor per organism so the user can open each in its own tab. Direct anchor
  # clicks bypass the pop-up blocker that nukes batched window.open() calls.
  output$organism_external_links <- renderUI({
    qs <- organism_queries()
    if (length(qs) == 0) return(NULL)

    clean <- function(s) {
      s <- tolower(s)
      s <- gsub("[^a-z0-9 _]", "", s)
      s <- gsub("[_ ]+", "+", trimws(s))
      s
    }
    cleaned <- clean(qs)
    keep <- nzchar(cleaned)
    cleaned <- cleaned[keep]
    qs <- qs[keep]
    if (length(cleaned) == 0) return(NULL)

    ncbi_urls <- paste0(
      "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?searchTerm=",
      cleaned,
      "&searchMode=complete+name&lock=1&unlock=1&command=search",
      "&curr_id=1427378&lvl=3&filter="
    )
    scholar_urls <- paste0(
      "https://scholar.google.es/scholar?hl=en&as_sdt=0%2C5&q=", cleaned, "&btnG="
    )
    pubmed_urls <- paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", cleaned)
    google_ai_urls <- paste0(
      "https://www.google.com/search?aep=11&udm=50&q=", cleaned
    )

    make_btn <- function(urls, queries, icon_name, label) {
      # Single query — direct anchor, no dropdown.
      if (length(urls) == 1) {
        return(tags$a(
          href = urls[[1]], target = "_blank", rel = "noopener",
          class = "btn btn-default", style = "margin-left: 4px;",
          icon(icon_name), " ", label
        ))
      }
      # Multi-query — dropdown listing one direct anchor per organism.
      # No "Open all" shortcut: browser pop-up blockers reliably nuke any
      # batched window.open() from a single click, so it's misleading UX.
      items <- mapply(function(url, q) {
        tags$li(tags$a(href = url, target = "_blank", rel = "noopener", q))
      }, urls, queries, SIMPLIFY = FALSE, USE.NAMES = FALSE)
      tags$div(
        class = "btn-group", style = "margin-left: 4px;",
        tags$button(
          type = "button",
          class = "btn btn-default dropdown-toggle",
          `data-toggle` = "dropdown",
          `aria-haspopup` = "true", `aria-expanded` = "false",
          icon(icon_name), " ", label, " ",
          tags$span(class = "badge", length(urls)), " ",
          tags$span(class = "caret")
        ),
        tags$ul(class = "dropdown-menu", items)
      )
    }

    tagList(
      make_btn(ncbi_urls,      qs, "search",         "NCBI Taxonomy"),
      make_btn(scholar_urls,   qs, "graduation-cap", "Google Scholar"),
      make_btn(pubmed_urls,    qs, "book-medical",   "PubMed"),
      make_btn(google_ai_urls, qs, "robot",          "Google AI Search"),
      tags$a(
        href   = "https://scholar.google.com/scholar_labs/search",
        target = "_blank",
        class  = "btn btn-default",
        style  = "margin-left: 4px;",
        icon("flask"), " Google Scholar Labs"
      )
    )
  })

  # Match each query against the current filtered_data() taxa universe using
  # the three-tier strategy. Recomputed only when the Search button is hit so
  # the output doesn't redraw on every keystroke.
  organism_matches <- eventReactive(input$organism_search_btn, {
    qs <- organism_queries()
    if (length(qs) == 0) return(NULL)
    df <- filtered_data()
    shiny::validate(shiny::need(nrow(df) > 0, "No data after current filters."))
    cols_want <- intersect(c("name", "Genus", "Species", "Phylum", "Family", "rank"),
                           colnames(df))
    pool <- unique(df[, cols_want, drop = FALSE])
    pool$.name_nrm  <- .organism_nrm(pool$name)
    pool$.genus_nrm <- if ("Genus" %in% colnames(pool)) .organism_nrm(pool$Genus) else ""
    lapply(qs, function(q) {
      qn <- .organism_nrm(q); strat <- "none"; m <- pool[0L, , drop = FALSE]
      if (nzchar(qn)) {
        m <- pool[pool$.name_nrm == qn, , drop = FALSE]; strat <- "exact"
        if (nrow(m) == 0) {
          m <- pool[grepl(qn, pool$.name_nrm, fixed = TRUE), , drop = FALSE]
          strat <- "substring"
        }
        if (nrow(m) == 0 && "Genus" %in% colnames(pool)) {
          first <- strsplit(qn, " ", fixed = TRUE)[[1]][1]
          if (length(first) && nzchar(first)) {
            m <- pool[grepl(first, pool$.genus_nrm, fixed = TRUE), , drop = FALSE]
            strat <- "genus"
          }
        }
        if (nrow(m) == 0) strat <- "none"
      }
      list(query = q, strategy = strat,
           taxa  = m[, setdiff(colnames(m), c(".name_nrm", ".genus_nrm")), drop = FALSE])
    })
  })

  # Union of matched taxon names across all queries.
  organism_matched_names <- reactive({
    ms <- organism_matches()
    if (is.null(ms) || length(ms) == 0) return(character(0))
    unique(unlist(lapply(ms, function(x) as.character(x$taxa$name)), use.names = FALSE))
  })

  # Rows of filtered_data() restricted to the matched taxa — the basis for
  # both the composition plot and the DT below.
  organism_subset <- reactive({
    nms <- organism_matched_names()
    if (length(nms) == 0) return(NULL)
    df <- filtered_data()
    df[df$name %in% nms, , drop = FALSE]
  })

  # Per-query match summary with a colored strategy tag.
  output$organism_match_summary <- renderUI({
    ms <- organism_matches()
    if (is.null(ms))
      return(tags$p(tags$em("Enter taxon names above and click Search.")))
    badge_style <- function(strat) {
      bg <- switch(strat,
                   exact     = "#c6efce", substring = "#fff2cc",
                   genus     = "#ffd9c4", none      = "#eee", "#eef")
      fg <- switch(strat,
                   exact     = "#2c5934", substring = "#7f5f00",
                   genus     = "#803b00", none      = "#666", "#335")
      sprintf("display:inline-block; padding:1px 8px; border-radius:3px; font-size:.85em; background:%s; color:%s; margin-left:.5em;", bg, fg)
    }
    chunks <- lapply(ms, function(x) {
      header <- tags$h4(
        "Query: ", tags$code(x$query),
        tags$span(style = badge_style(x$strategy), x$strategy)
      )
      if (x$strategy == "none" || nrow(x$taxa) == 0)
        return(tagList(header,
          tags$p(tags$em("No matches. Try a shorter substring or just the Genus."))))
      tagList(
        header,
        tags$p(sprintf("%d taxon match%s — using all of these for the plot/table below.",
                       nrow(x$taxa), if (nrow(x$taxa) > 1) "es" else "")),
        HTML(renderTable(x$taxa, striped = TRUE, hover = TRUE, bordered = TRUE,
                         width = "100%")())
      )
    })
    do.call(tagList, chunks)
  })

  output$organism_group_selector <- renderUI({
    cols <- metadata_group_cols()
    if (length(cols) == 0) return(NULL)
    selectInput("organism_group_by", "Group By (metadata):",
                choices = c("None", cols),
                selected = restored_input("organism_group_by", "None",
                                          valid = c("None", cols)))
  })

  # Composition-style barplot restricted to the matched taxa (no Top-N rollup;
  # every matched taxon gets its own stack segment). Shared by the on-screen
  # plotly view and the PDF/PNG/SVG download handlers.
  organism_gg <- reactive({
    df <- organism_subset()
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0,
                                "No matched rows for the current filters/queries."))
    group_col <- resolve_group_col(input$organism_group_by)
    if (!is.null(group_col) && !group_col %in% colnames(df)) group_col <- NULL

    # Aggregate per sample x name (+ group_col when faceting).
    has_reads <- "reads_clade" %in% colnames(df)
    group_vars <- c("Sample", "name")
    if (!is.null(group_col)) group_vars <- c(group_vars, group_col)
    rhs <- paste(group_vars, collapse = " + ")
    if (has_reads) {
      f <- as.formula(paste("cbind(percent, reads_clade) ~", rhs))
    } else {
      f <- as.formula(paste("percent ~", rhs))
    }
    plot_df <- aggregate(f, data = df, sum)

    # Order taxa by total (ascending so largest sits at bottom of stack).
    totals <- aggregate(percent ~ name, data = plot_df, sum)
    totals <- totals[order(totals$percent), , drop = FALSE]
    plot_df$name <- factor(plot_df$name, levels = as.character(totals$name))

    if (input$organism_display_mode == "counts" && has_reads) {
      plot_df$y_value <- plot_df$reads_clade; y_lab <- "Counts"
    } else {
      plot_df$y_value <- plot_df$percent;     y_lab <- "% Reads"
    }
    n_taxa <- nlevels(plot_df$name)
    cols <- get_palette_colors(max(n_taxa, 1), input$organism_color_palette)
    names(cols) <- levels(plot_df$name)

    p <- ggplot(plot_df, aes(x = Sample, y = y_value, fill = name,
                              text = paste0("Sample: ", Sample,
                                            "\nTaxon: ", name,
                                            "\nValue: ", sprintf("%.3f", y_value)))) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_manual(values = cols, drop = FALSE) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.position = "right",
            panel.grid.major.x = element_blank()) +
      labs(title = "Organisms of interest — composition",
           x = NULL, y = y_lab, fill = "Taxon")

    if (isTRUE(input$organism_show_labels)) {
      label_df <- plot_df[plot_df$y_value > 0.5, , drop = FALSE]
      if (nrow(label_df) > 0)
        p <- p + geom_text(data = label_df,
                           aes(label = sprintf("%.1f", y_value)),
                           position = position_stack(vjust = 0.5),
                           size = 2.8, color = "black", fontface = "bold")
    }
    if (!is.null(group_col)) {
      p <- p + facet_grid(cols = vars(!!sym(group_col)),
                          scales = "free_x", space = "free_x")
    }
    p
  })

  output$organism_plot <- renderPlotly({
    ggplotly(organism_gg(), tooltip = "text") %>%
      layout(legend = list(font = list(size = 10)), margin = list(b = 120))
  })

  # DT of the matched rows. Columns ordered for readability; everything from
  # the underlying filtered_data() carries through so users can sort/search.
  output$organism_table <- renderDT({
    df <- organism_subset()
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0,
                                "Run a search to populate the table."))
    pref <- c("Sample", "name", "rank", "percent", "reads_clade",
              "Genus", "Species", "Family", "Phylum")
    cols <- c(intersect(pref, colnames(df)),
              setdiff(colnames(df), pref))
    datatable(
      df[, cols, drop = FALSE],
      rownames = FALSE, filter = "top",
      extensions = "Buttons",
      options = list(scrollX = TRUE, dom = "lftBip",
                     buttons = list(
                       list(extend = "copyHtml5",  exportOptions = dt_export_options),
                       list(extend = "csvHtml5",   exportOptions = dt_export_options),
                       list(extend = "excelHtml5", exportOptions = dt_export_options)
                     ),
                     pageLength = 20))
  })

  output$organism_download_csv <- downloadHandler(
    filename = function()
      paste0("organisms_of_interest_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      df <- organism_subset()
      if (is.null(df)) df <- data.frame()
      write.csv(df, file, row.names = FALSE)
    }
  )

  # User-controllable export dimensions (mirrors main Composition tab).
  organism_export_dims <- reactive({
    list(
      width  = if (is.null(input$organism_plot_width))  12 else input$organism_plot_width,
      height = if (is.null(input$organism_plot_height))  7 else input$organism_plot_height
    )
  })

  output$download_organism_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_organisms_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) {
      d <- organism_export_dims()
      save_gg(file, organism_gg(), "pdf", width = d$width, height = d$height)
    }
  )
  output$download_organism_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_organisms_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      d <- organism_export_dims()
      save_gg(file, organism_gg(), "png", width = d$width, height = d$height)
    }
  )
  output$download_organism_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_organisms_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) {
      d <- organism_export_dims()
      save_gg(file, organism_gg(), "svg", width = d$width, height = d$height)
    }
  )

  # ── Organisms: Composition Group comparison ──
  # Same non-parametric per-taxon logic as the main Composition tab, but
  # scoped to the matched taxa only. Uses the shared organism_group_by
  # selector for the grouping variable.
  organism_stats_df <- reactive({
    df <- organism_subset()
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0,
                                "No matched rows for the current filters/queries."))
    gcol <- resolve_group_col(input$organism_group_by)
    if (is.null(gcol) || !(gcol %in% colnames(df))) return(NULL)
    val_col <- if (identical(input$organism_display_mode, "counts") &&
                    "reads_clade" %in% colnames(df)) "reads_clade" else "percent"
    per_taxon_group_table(df, unique(df$name), gcol, value_col = val_col)
  })

  output$organism_stats_status <- renderUI({
    gcol <- resolve_group_col(input$organism_group_by)
    if (is.null(gcol)) {
      return(div(class = "status-text status-warning",
                 icon("info-circle"),
                 " Set a Group By value above to enable per-taxon comparisons."))
    }
    res <- tryCatch(organism_stats_df(), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) {
      return(div(class = "status-text status-warning",
                 icon("exclamation-triangle"),
                 " No comparable groups (need ≥2 samples per group)."))
    }
    NULL
  })

  output$organism_stats_table <- renderDT({
    res <- organism_stats_df()
    shiny::validate(shiny::need(!is.null(res) && nrow(res) > 0,
                                "Nothing to report — set a Group By with ≥2 samples per group."))
    show <- res
    for (col in c("statistic", "p", "p_adj")) {
      if (col %in% colnames(show))
        show[[col]] <- signif(as.numeric(show[[col]]), 3)
    }
    datatable(show, rownames = FALSE, filter = "top",
              extensions = "Buttons",
              options = list(pageLength = 10, dom = "Blfrtip",
                             buttons = list(
                               list(extend = "copyHtml5",  exportOptions = dt_export_options),
                               list(extend = "csvHtml5",   exportOptions = dt_export_options),
                               list(extend = "excelHtml5", exportOptions = dt_export_options)
                             )))
  })

  output$download_organism_stats_csv <- downloadHandler(
    filename = function() paste0("organism_composition_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- organism_stats_df()
      if (is.null(res)) res <- data.frame()
      write.csv(res, file, row.names = FALSE)
    }
  )

  # ── Organisms: Distribution panel (matched taxa only) ──
  # Own Group By selector (independent from the Composition panel).
  output$organism_distrib_group_selector <- renderUI({
    cols <- metadata_group_cols()
    if (length(cols) == 0) return(NULL)
    selectInput("organism_distrib_group", "Group By (metadata):",
                choices = c("None", cols),
                selected = restored_input("organism_distrib_group", "None",
                                          valid = c("None", cols)))
  })

  # Keep the taxon selector in sync with the matched taxa list.
  observe({
    df <- organism_subset()
    taxa <- if (is.null(df) || nrow(df) == 0) character(0)
            else sort(unique(df$name[df$name != "unclassified"]))
    # Prefer a taxon restored from a settings file (consumed once); else first.
    sel <- restored_input("organism_violin_taxon", NULL, valid = taxa)
    if (!is.null(sel)) {
      isolate({ ri <- rv$restore_inputs; ri[["organism_violin_taxon"]] <- NULL
                rv$restore_inputs <- ri })
    } else {
      sel <- if (length(taxa) > 0) taxa[1] else NULL
    }
    updateSelectizeInput(session, "organism_violin_taxon", choices = taxa,
                        selected = sel, server = TRUE)
  })

  # External search links for the currently selected matched taxon.
  output$organism_taxon_search_links <- renderUI({
    tx <- input$organism_violin_taxon
    if (is.null(tx) || !nzchar(tx)) return(NULL)
    clean_name <- gsub("[^a-z0-9 _]", "", tolower(tx))
    clean_url  <- gsub("[_ ]+", "+", trimws(clean_name))
    ncbi_url <- paste0(
      "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?",
      "searchTerm=", clean_url,
      "&searchMode=complete+name&lock=1&unlock=1&command=search",
      "&curr_id=1427378&lvl=3&filter="
    )
    scholar_url <- paste0("https://scholar.google.es/scholar?hl=en&as_sdt=0%2C5&q=", clean_url, "&btnG=")
    pubmed_url  <- paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", clean_url)
    google_ai_url <- paste0("https://www.google.com/search?aep=11&udm=50&q=", clean_url)
    tags$div(
      style = "display: flex; flex-wrap: wrap; gap: 6px; margin: 6px 0;",
      tags$a(href = ncbi_url,      target = "_blank", class = "btn btn-default btn-sm",
             icon("search"),         " NCBI Taxonomy"),
      tags$a(href = scholar_url,   target = "_blank", class = "btn btn-default btn-sm",
             icon("graduation-cap"), " Google Scholar"),
      tags$a(href = pubmed_url,    target = "_blank", class = "btn btn-default btn-sm",
             icon("book-medical"),   " PubMed"),
      tags$a(href = google_ai_url, target = "_blank", class = "btn btn-default btn-sm",
             icon("robot"),          " Google AI Search")
    )
  })

  # ggplot for the matched-taxa distribution plot. Mirrors output$violin_plot
  # (single vs faceted, Violin / Box+Strip / Strip) but constrained to the
  # matched taxa. Uses the tab-local organism_distrib_group input.
  organism_violin_gg <- reactive({
    df <- organism_subset()
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0,
                                "No matched rows for the current filters/queries."))
    gcol <- resolve_group_col(input$organism_distrib_group)
    view_mode <- if (!is.null(input$organism_view_mode)) input$organism_view_mode else "single"
    plot_type <- input$organism_plot_type %||% "violin"
    pal <- input$organism_color_palette %||% "Set3"

    if (view_mode == "single") {
      tx <- input$organism_violin_taxon
      shiny::validate(shiny::need(!is.null(tx) && nzchar(tx),
                                  "Select a taxon to plot."))
      sub <- df[df$name == tx, , drop = FALSE]
      shiny::validate(shiny::need(nrow(sub) > 0,
                                  "Selected taxon has no rows in the matched set."))
      if (!is.null(gcol)) sub$.grp <- as.character(sub[[gcol]]) else sub$.grp <- "All"
      base <- ggplot(sub, aes(x = .grp, y = percent, fill = .grp))
      if (plot_type == "violin") {
        p <- base + geom_violin(alpha = 0.6) +
          geom_boxplot(width = 0.15, alpha = 0.4) +
          geom_jitter(width = 0.1, size = 2, alpha = 0.6)
      } else if (plot_type == "box") {
        p <- base + geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          geom_jitter(width = 0.15, size = 2, alpha = 0.6)
      } else {
        p <- base + geom_jitter(width = 0.15, size = 2, alpha = 0.7)
      }
      grp_levels <- unique(sub$.grp)
      p + scale_fill_manual(values = get_palette_colors(length(grp_levels), pal)) +
        theme_minimal(base_size = 14) +
        theme(legend.position = if (is.null(gcol)) "none" else "right") +
        labs(title = paste("Distribution:", tx),
             x = if (!is.null(gcol)) gcol else NULL,
             y = "Abundance (%)")
    } else {
      # Faceted: one panel per matched taxon.
      clean <- df[df$name != "unclassified", , drop = FALSE]
      shiny::validate(shiny::need(nrow(clean) > 0, "No matched taxa to plot."))
      taxa <- sort(unique(clean$name))
      clean$name <- factor(clean$name, levels = taxa)
      if (!is.null(gcol)) clean$.grp <- as.character(clean[[gcol]]) else clean$.grp <- "All"
      base <- ggplot(clean, aes(x = .grp, y = percent, fill = .grp))
      if (plot_type == "violin") {
        p <- base + geom_violin(alpha = 0.6) +
          geom_boxplot(width = 0.15, alpha = 0.4)
      } else if (plot_type == "box") {
        p <- base + geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          geom_jitter(width = 0.15, size = 1.5, alpha = 0.5)
      } else {
        p <- base + geom_jitter(width = 0.15, size = 1.5, alpha = 0.7)
      }
      grp_levels <- unique(clean$.grp)
      p + scale_fill_manual(values = get_palette_colors(length(grp_levels), pal)) +
        facet_wrap(~name, scales = "free_y") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              legend.position = if (is.null(gcol)) "none" else "right") +
        labs(title = "Matched taxa — distribution",
             x = if (!is.null(gcol)) gcol else NULL,
             y = "Abundance (%)")
    }
  })

  output$organism_violin_plot <- renderPlotly({
    ggplotly(organism_violin_gg(), tooltip = c("x", "y", "fill")) %>%
      layout(margin = list(b = 90))
  })

  output$download_organism_distrib_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_organism_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) {
      d <- organism_export_dims()
      save_gg(file, organism_violin_gg(), "pdf", width = d$width, height = d$height)
    }
  )
  output$download_organism_distrib_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_organism_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      d <- organism_export_dims()
      save_gg(file, organism_violin_gg(), "png", width = d$width, height = d$height)
    }
  )
  output$download_organism_distrib_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_organism_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) {
      d <- organism_export_dims()
      save_gg(file, organism_violin_gg(), "svg", width = d$width, height = d$height)
    }
  )

  # ── Organisms: Distribution Statistics ──
  organism_distrib_stats <- reactive({
    df <- organism_subset()
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0,
                                "No matched rows for the current filters/queries."))
    gcol <- resolve_group_col(input$organism_distrib_group)
    if (is.null(gcol) || !(gcol %in% colnames(df)))
      return(list(main = NULL, pairwise = NULL, mode = NULL))
    view_mode <- if (!is.null(input$organism_view_mode)) input$organism_view_mode else "single"
    clean <- df[df$name != "unclassified", , drop = FALSE]
    if (nrow(clean) == 0) return(list(main = NULL, pairwise = NULL, mode = view_mode))
    if (view_mode == "single") {
      tx <- input$organism_violin_taxon
      if (is.null(tx) || !nzchar(tx))
        return(list(main = NULL, pairwise = NULL, mode = view_mode))
      sub <- clean[clean$name == tx, , drop = FALSE]
      main <- per_taxon_group_table(sub, tx, gcol, value_col = "percent")
      pw   <- pairwise_wilcox_bh(sub$percent, sub[[gcol]])
      list(main = main, pairwise = pw, mode = view_mode, taxon = tx)
    } else {
      taxa <- unique(clean$name)
      main <- per_taxon_group_table(clean, taxa, gcol, value_col = "percent")
      list(main = main, pairwise = NULL, mode = view_mode)
    }
  })

  output$organism_distrib_stats_status <- renderUI({
    if (is.null(resolve_group_col(input$organism_distrib_group))) {
      return(div(class = "status-text status-warning",
                 icon("info-circle"),
                 " Choose a Group By value to enable group-comparison tests."))
    }
    res <- tryCatch(organism_distrib_stats(), error = function(e) NULL)
    if (is.null(res) || is.null(res$main) || nrow(res$main) == 0) {
      return(div(class = "status-text status-warning",
                 icon("exclamation-triangle"),
                 " No comparable groups (need ≥2 samples per group)."))
    }
    NULL
  })

  output$organism_distrib_stats_table <- renderDT({
    res <- organism_distrib_stats()
    shiny::validate(shiny::need(!is.null(res$main) && nrow(res$main) > 0,
                                "Nothing to report — set a Group By with ≥2 samples per group."))
    show <- res$main
    for (col in c("statistic", "p", "p_adj")) {
      if (col %in% colnames(show))
        show[[col]] <- signif(as.numeric(show[[col]]), 3)
    }
    datatable(show, rownames = FALSE, filter = "top",
              extensions = "Buttons",
              options = list(pageLength = 10, dom = "Blfrtip",
                             buttons = list(
                               list(extend = "copyHtml5",  exportOptions = dt_export_options),
                               list(extend = "csvHtml5",   exportOptions = dt_export_options),
                               list(extend = "excelHtml5", exportOptions = dt_export_options)
                             )))
  })

  output$organism_distrib_pairwise_header <- renderUI({
    res <- tryCatch(organism_distrib_stats(), error = function(e) NULL)
    if (is.null(res) || is.null(res$pairwise) || nrow(res$pairwise) == 0) return(NULL)
    tags$h4(sprintf("Pairwise Wilcoxon (BH-adjusted) — %s", res$taxon %||% ""))
  })

  output$organism_distrib_stats_pairwise <- renderDT({
    res <- organism_distrib_stats()
    if (is.null(res$pairwise) || nrow(res$pairwise) == 0) return(NULL)
    show <- res$pairwise
    show$p_adj <- signif(show$p_adj, 3)
    datatable(show, rownames = FALSE,
              options = list(pageLength = 10, dom = "tp"))
  })

  output$download_organism_distrib_stats_csv <- downloadHandler(
    filename = function() paste0("organism_distribution_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- organism_distrib_stats()
      main <- res$main; if (is.null(main)) main <- data.frame()
      write.csv(main, file, row.names = FALSE)
    }
  )

  # ═══════════════════════════════════════════════════════════════════════
  # TAB 2: HEATMAP
  # ═══════════════════════════════════════════════════════════════════════

  output$heatmap_group_selector <- renderUI({
    cols <- metadata_group_cols()
    if (length(cols) == 0) return(NULL)
    selectInput("heatmap_group_by", "Facet by (metadata):",
                choices = c("None", cols),
                selected = restored_input("heatmap_group_by", "None",
                                          valid = c("None", cols)))
  })

  # Build the (optionally faceted) heatmap ggplot from a heatmap_prep() list.
  # Used for the faceted on-screen view (via ggplotly) and all PDF/PNG exports.
  build_heatmap_gg <- function(prep) {
    mat  <- prep$mat
    long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    colnames(long) <- c("name", "Sample", "percent")
    long$name   <- factor(long$name,   levels = rownames(mat))
    long$Sample <- factor(long$Sample, levels = colnames(mat))

    faceted <- !is.null(prep$group_col) && !is.null(prep$samp_group)
    if (faceted) long$.grp <- factor(prep$samp_group[as.character(long$Sample)])
    title <- if (faceted) {
      paste("Abundance Heatmap (Top", prep$top_n, ") — faceted by", prep$group_col)
    } else {
      paste("Abundance Heatmap (Top", prep$top_n, ")")
    }

    g <- ggplot(long, aes(x = Sample, y = name, fill = percent)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "#f8f9fa", high = "#1a5276") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = title, x = NULL, y = NULL, fill = "%")

    if (faceted) {
      g <- g + facet_grid(cols = vars(.grp), scales = "free_x", space = "free_x")
    }
    g
  }

  # ── Shared heatmap builder ──
  # Returns a list(mat, top_taxa, taxa_order, samp_order, group_col) so the
  # interactive (plotly) and download (ggplot) renderers stay in sync. `mat`
  # is taxa×samples with taxa rows clustered; samples ordered by group (when a
  # metadata facet is chosen) then clustered, otherwise clustered globally.
  heatmap_prep <- reactive({
    df <- filtered_data()
    req(nrow(df) > 0)
    top_n <- input$top_n
    group_col <- resolve_group_col(input$heatmap_group_by)

    clean <- df[df$name != "unclassified", ]
    taxa_means <- aggregate(percent ~ name, data = clean, FUN = mean)
    taxa_means <- taxa_means[order(taxa_means$percent, decreasing = TRUE), ]
    top_taxa <- head(taxa_means$name, top_n)

    hm_data <- df[df$name %in% top_taxa, ]
    req(nrow(hm_data) > 0)

    wide <- hm_data %>%
      select(Sample, name, percent) %>%
      group_by(Sample, name) %>%
      summarise(percent = sum(percent, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Sample, values_from = percent, values_fill = 0)

    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$name

    # Cluster taxa (rows) for a stable, readable vertical order.
    if (nrow(mat) >= 2) {
      row_order <- tryCatch(hclust(dist(mat))$order, error = function(e) seq_len(nrow(mat)))
      mat <- mat[row_order, , drop = FALSE]
    }

    # Map each sample → its group level (for faceting / column ordering).
    samp_group <- NULL
    if (!is.null(group_col)) {
      lookup <- hm_data[!duplicated(hm_data$Sample), c("Sample", group_col)]
      samp_group <- setNames(as.character(lookup[[group_col]]), lookup$Sample)
      samp_group[is.na(samp_group)] <- "NA"
    }

    # Order columns: by group (when faceting), clustering within the full set.
    if (ncol(mat) >= 2) {
      col_order <- tryCatch(hclust(dist(t(mat)))$order, error = function(e) seq_len(ncol(mat)))
      mat <- mat[, col_order, drop = FALSE]
    }
    if (!is.null(samp_group)) {
      g <- samp_group[colnames(mat)]
      mat <- mat[, order(g), drop = FALSE]
    }

    list(mat = mat, top_n = top_n, group_col = group_col, samp_group = samp_group)
  })

  output$heatmap_plot <- renderPlotly({
    shiny::validate(shiny::need(nrow(filtered_data()) > 0, empty_filter_message()))
    prep <- heatmap_prep()
    mat <- prep$mat
    top_n <- prep$top_n

    if (is.null(prep$group_col)) {
      # No metadata facet → original single native-plotly heatmap.
      plot_ly(
        z = mat,
        x = colnames(mat),
        y = rownames(mat),
        type = "heatmap",
        colorscale = list(c(0, "#f8f9fa"), c(0.5, "#3498db"), c(1, "#1a5276")),
        hovertemplate = "Taxon: %{y}<br>Sample: %{x}<br>Abundance: %{z:.2f}%<extra></extra>"
      ) %>%
        layout(
          title = list(text = paste("Abundance Heatmap (Top", top_n, "taxa)"),
                       font = list(size = 16, family = "Inter")),
          xaxis = list(title = "", tickangle = -45, tickfont = list(size = 10)),
          yaxis = list(title = "", tickfont = list(size = 10)),
          margin = list(l = 200, b = 130)
        )
    } else {
      # Faceted heatmap: one panel per metadata level, sharing the taxa axis.
      gg <- build_heatmap_gg(prep)
      ggplotly(gg) %>%
        layout(margin = list(l = 200, b = 130))
    }
  })

  # ═══════════════════════════════════════════════════════════════════════
  # TAB 3: VIOLIN / BOX DISTRIBUTION
  # ═══════════════════════════════════════════════════════════════════════

  # Update taxon choices safely using server-side processing
  observe({
    df <- filtered_data()
    req(nrow(df) > 0)
    taxa <- sort(unique(df$name[df$name != "unclassified"]))

    # Prefer a taxon restored from a settings file (consumed once so later
    # filter changes don't keep re-forcing it); otherwise default to the first.
    sel <- restored_input("violin_taxon", NULL, valid = taxa)
    if (!is.null(sel)) {
      isolate({ ri <- rv$restore_inputs; ri[["violin_taxon"]] <- NULL
                rv$restore_inputs <- ri })
    } else {
      sel <- if (length(taxa) > 0) taxa[1] else NULL
    }
    # Use updateSelectizeInput with server = TRUE to prevent DOM lockup
    updateSelectizeInput(session, "violin_taxon", choices = taxa,
                         selected = sel, server = TRUE)
  })

  output$violin_group_selector <- renderUI({
    if (!is.null(rv$metadata)) {
      group_cols <- setdiff(colnames(rv$metadata), "Sample")
      selectInput("violin_group", "Color/Group By:",
                  choices = c("None", group_cols),
                  selected = restored_input("violin_group", "None",
                                            valid = c("None", group_cols)))
    }
  })

  # ── Taxon search links (NCBI, Scholar, PubMed, Google AI) ──
  output$taxon_search_links <- renderUI({
    req(input$violin_taxon)

    clean_name <- tolower(input$violin_taxon)
    clean_name <- gsub("[^a-z0-9 _]", "", clean_name)
    clean_url  <- gsub("[_ ]+", "+", trimws(clean_name))

    ncbi_url <- paste0(
      "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?",
      "searchTerm=", clean_url,
      "&searchMode=complete+name&lock=1&unlock=1&command=search",
      "&curr_id=1427378&lvl=3&filter="
    )
    scholar_url <- paste0(
      "https://scholar.google.es/scholar?hl=en&as_sdt=0%2C5&q=", clean_url, "&btnG="
    )
    pubmed_url <- paste0(
      "https://pubmed.ncbi.nlm.nih.gov/?term=", clean_url
    )
    google_ai_url <- paste0(
      "https://www.google.com/search?aep=11&udm=50&q=", clean_url
    )

    tagList(
      tags$div(
        style = "display: flex; flex-wrap: wrap; gap: 6px; margin-top: 4px;",
        tags$a(href = ncbi_url,       target = "_blank", class = "btn btn-default btn-sm",
               icon("search"), " NCBI Taxonomy"),
        tags$a(href = scholar_url,    target = "_blank", class = "btn btn-default btn-sm",
               icon("graduation-cap"), " Google Scholar"),
        tags$a(href = pubmed_url,     target = "_blank", class = "btn btn-default btn-sm",
               icon("book-medical"), " PubMed"),
        tags$a(href = google_ai_url,  target = "_blank", class = "btn btn-default btn-sm",
               icon("robot"), " Google AI Search"),
        tags$a(href = "https://scholar.google.com/scholar_labs/search",
               target = "_blank", class = "btn btn-default btn-sm",
               icon("flask"), " Scholar Labs (experimental)")
      )
    )
  })

  output$violin_plot <- renderPlotly({
    df <- filtered_data()
    shiny::validate(shiny::need(nrow(df) > 0, empty_filter_message()))

    # NB: filtered_data() already joins rv$metadata (via merged_data()). Do NOT
    # merge it again here — a second merge on "Sample" renames the shared
    # columns to Condition.x/.y, so input$violin_group ("Condition") no longer
    # matches colnames(df) and grouping silently does nothing.

    view_mode <- if (!is.null(input$violin_view_mode)) input$violin_view_mode else "single"

    group_var <- if (!is.null(input$violin_group) && input$violin_group != "None" &&
                     input$violin_group %in% colnames(df)) {
      input$violin_group
    } else {
      NULL
    }

    plot_type <- input$violin_plot_type

    # ── Faceted mode: all top N taxa ──
    if (view_mode == "faceted") {
      # Get filtered data at current rank (outer validate above already fired
      # if this was empty; the local rebind here is just for clarity).
      filt <- df

      clean <- filt[filt$name != "unclassified", , drop = FALSE]
      if (nrow(clean) == 0) return(NULL)

      taxa_means <- aggregate(percent ~ name, data = clean, mean)
      taxa_means <- taxa_means[order(taxa_means$percent, decreasing = TRUE), ]
      top_n_val <- if (!is.null(input$top_n)) input$top_n else 10
      top_taxa <- head(taxa_means$name, top_n_val)

      plot_data <- df[df$name %in% top_taxa, , drop = FALSE]
      req(nrow(plot_data) > 0)
      plot_data$name <- factor(plot_data$name, levels = top_taxa)

      if (!is.null(group_var)) {
        plot_data$group_val <- plot_data[[group_var]]
        p <- ggplot(plot_data, aes(x = group_val, y = percent, fill = group_val))
      } else {
        plot_data$All <- "All"
        p <- ggplot(plot_data, aes(x = All, y = percent))
      }

      if (plot_type == "violin") {
        # Native plotly violin — build directly
        if (!is.null(group_var)) {
          p_ly <- plot_ly()
          groups <- unique(plot_data$group_val)
          pal <- get_palette_colors(length(groups), input$color_palette)
          for (i in seq_along(groups)) {
            grp_data <- plot_data[plot_data$group_val == groups[i], , drop = FALSE]
            p_ly <- add_trace(p_ly, data = grp_data, type = "violin",
                              x = ~name, y = ~percent, split = ~group_val,
                              name = groups[i],
                              meanline = list(visible = TRUE),
                              box = list(visible = TRUE),
                              marker = list(color = pal[i]),
                              line = list(color = pal[i]),
                              fillcolor = adjustcolor(pal[i], alpha.f = 0.3))
          }
        } else {
          p_ly <- plot_ly(plot_data, type = "violin",
                          x = ~name, y = ~percent,
                          meanline = list(visible = TRUE),
                          box = list(visible = TRUE),
                          marker = list(color = "#C4692A"),
                          line = list(color = "#C4692A"),
                          fillcolor = "rgba(196,105,42,0.3)")
        }
        n_facets <- length(top_taxa)
        plot_height <- max(600, ceiling(n_facets / 3) * 280)
        return(
          p_ly %>% layout(
            title = list(text = paste("Distribution of Top", top_n_val, "Taxa"),
                         font = list(size = 16)),
            xaxis = list(title = ""),
            yaxis = list(title = "Abundance (%)"),
            violinmode = "group",
            height = plot_height
          )
        )
      }

      if (plot_type == "box") {
        p <- p +
          geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          geom_point(position = position_jitter(width = 0.15, height = 0),
                     size = 2, alpha = 0.7)
      } else {
        p <- p +
          geom_point(position = position_jitter(width = 0.2, height = 0),
                     size = 3, alpha = 0.7)
      }

      p <- p +
        facet_wrap(~name, scales = "free_y", ncol = 3) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = if (is.null(group_var)) "none" else "right"
        ) +
        labs(
          title = paste("Distribution of Top", top_n_val, "Taxa"),
          x = if (!is.null(group_var)) group_var else "",
          y = "Abundance (%)"
        )

      if (!is.null(group_var)) {
        pal_colors <- get_palette_colors(length(unique(plot_data[[group_var]])), input$color_palette)
        p <- p + scale_fill_manual(values = pal_colors)
      } else {
        p <- p + scale_fill_brewer(palette = "Set3")
      }

      n_facets <- length(top_taxa)
      plot_height <- max(600, ceiling(n_facets / 3) * 280)
      ggplotly(p, tooltip = "text", height = plot_height)

    } else {
      # ── Single taxon mode (original behavior) ──
      req(input$violin_taxon)
      taxon <- input$violin_taxon
      taxon_data <- df[df$name == taxon, , drop = FALSE]
      req(nrow(taxon_data) > 0)

      if (!is.null(group_var)) {
        taxon_data$group_val <- taxon_data[[group_var]]
        p <- ggplot(taxon_data, aes(x = group_val, y = percent,
                                     fill = group_val,
                                     text = paste0("Sample: ", Sample,
                                                   "\n", group_var, ": ", group_val,
                                                   "\nAbundance: ", sprintf("%.3f%%", percent))))
      } else {
        taxon_data$All <- "All Samples"
        p <- ggplot(taxon_data, aes(x = All, y = percent,
                                     text = paste0("Sample: ", Sample,
                                                   "\nAbundance: ", sprintf("%.3f%%", percent))))
      }

      if (plot_type == "violin") {
        # Native plotly violin for single taxon
        if (!is.null(group_var)) {
          groups <- unique(taxon_data$group_val)
          pal <- get_palette_colors(length(groups), input$color_palette)
          p_ly <- plot_ly()
          for (i in seq_along(groups)) {
            grp_data <- taxon_data[taxon_data$group_val == groups[i], , drop = FALSE]
            p_ly <- add_trace(p_ly, data = grp_data, type = "violin",
                              x = ~group_val, y = ~percent,
                              name = groups[i],
                              meanline = list(visible = TRUE),
                              box = list(visible = TRUE),
                              points = "all",
                              jitter = 0.3,
                              pointpos = -1.8,
                              marker = list(color = pal[i], size = 5),
                              line = list(color = pal[i]),
                              fillcolor = adjustcolor(pal[i], alpha.f = 0.3),
                              text = ~paste0("Sample: ", Sample,
                                             "\n", group_var, ": ", group_val,
                                             "\nAbundance: ", sprintf("%.3f%%", percent)),
                              hoverinfo = "text")
          }
        } else {
          taxon_data$x_label <- "All Samples"
          p_ly <- plot_ly(taxon_data, type = "violin",
                          x = ~x_label, y = ~percent,
                          meanline = list(visible = TRUE),
                          box = list(visible = TRUE),
                          points = "all",
                          jitter = 0.3,
                          pointpos = -1.8,
                          marker = list(color = "#C4692A", size = 5),
                          line = list(color = "#C4692A"),
                          fillcolor = "rgba(196,105,42,0.3)",
                          text = ~paste0("Sample: ", Sample,
                                         "\nAbundance: ", sprintf("%.3f%%", percent)),
                          hoverinfo = "text")
        }
        return(
          p_ly %>% layout(
            title = list(text = paste("Distribution:", taxon),
                         font = list(size = 16, family = "Inter")),
            yaxis = list(title = "Abundance (%)"),
            xaxis = list(title = if (!is.null(group_var)) group_var else "")
          )
        )
      }

      if (plot_type == "box") {
        p <- p +
          geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          geom_point(position = position_jitter(width = 0.15, height = 0),
                     size = 3, alpha = 0.7, color = "#2c3e50")
      } else {
        p <- p +
          geom_point(position = position_jitter(width = 0.15, height = 0),
                     size = 3, alpha = 0.7, color = "#2c3e50")
      }

      p <- p +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 16),
          axis.text.x = element_text(size = 11),
          legend.position = if (is.null(group_var)) "none" else "right"
        ) +
        labs(
          title = paste("Distribution:", taxon),
          y = "Abundance (%)",
          x = if (!is.null(group_var)) group_var else ""
        )

      if (!is.null(group_var)) {
        pal_colors <- get_palette_colors(length(unique(taxon_data[[group_var]])), input$color_palette)
        p <- p + scale_fill_manual(values = pal_colors)
      }

      ggplotly(p, tooltip = "text")
    }
  })

  # ── Distribution Statistics ──
  # Single-taxon mode: one test on the selected taxon across the grouping.
  # Faceted mode: per-taxon test over the top-N taxa with BH FDR across taxa.
  # Empty when there's no metadata grouping or fewer than 2 comparable groups.
  distrib_group_col <- reactive({
    if (!is.null(input$violin_group) && input$violin_group != "None" &&
        input$violin_group %in% colnames(merged_data())) input$violin_group
    else NULL
  })

  distrib_stats <- reactive({
    df <- filtered_data()
    shiny::validate(shiny::need(nrow(df) > 0, empty_filter_message()))
    gcol <- distrib_group_col()
    if (is.null(gcol)) return(list(main = NULL, pairwise = NULL, mode = NULL))
    view_mode <- if (!is.null(input$violin_view_mode)) input$violin_view_mode else "single"
    clean <- df[df$name != "unclassified", , drop = FALSE]
    if (nrow(clean) == 0) return(list(main = NULL, pairwise = NULL, mode = view_mode))
    if (view_mode == "single") {
      tx <- input$violin_taxon
      if (is.null(tx) || !nzchar(tx)) return(list(main = NULL, pairwise = NULL, mode = view_mode))
      sub <- clean[clean$name == tx, , drop = FALSE]
      main <- per_taxon_group_table(sub, tx, gcol, value_col = "percent")
      pw   <- pairwise_wilcox_bh(sub$percent, sub[[gcol]])
      list(main = main, pairwise = pw, mode = view_mode, taxon = tx)
    } else {
      taxa_means <- aggregate(percent ~ name, data = clean, FUN = mean)
      taxa_means <- taxa_means[order(taxa_means$percent, decreasing = TRUE), ]
      top_n_val <- if (!is.null(input$top_n)) input$top_n else 10
      top_taxa  <- head(taxa_means$name, top_n_val)
      main <- per_taxon_group_table(clean, top_taxa, gcol, value_col = "percent")
      list(main = main, pairwise = NULL, mode = view_mode)
    }
  })

  output$distrib_stats_status <- renderUI({
    if (is.null(distrib_group_col())) {
      return(div(class = "status-text status-warning",
                 icon("info-circle"),
                 " Choose a Group By value to enable group-comparison tests."))
    }
    res <- tryCatch(distrib_stats(), error = function(e) NULL)
    if (is.null(res) || is.null(res$main) || nrow(res$main) == 0) {
      return(div(class = "status-text status-warning",
                 icon("exclamation-triangle"),
                 " No comparable groups (need ≥2 samples per group)."))
    }
    NULL
  })

  output$distrib_stats_table <- renderDT({
    res <- distrib_stats()
    shiny::validate(shiny::need(!is.null(res$main) && nrow(res$main) > 0,
                                "Nothing to report — set a Group By with ≥2 samples per group."))
    show <- res$main
    for (col in c("statistic", "p", "p_adj")) {
      if (col %in% colnames(show))
        show[[col]] <- signif(as.numeric(show[[col]]), 3)
    }
    datatable(show, rownames = FALSE, filter = "top",
              extensions = "Buttons",
              options = list(pageLength = 10, dom = "Blfrtip",
                             buttons = list(
                               list(extend = "copyHtml5",  exportOptions = dt_export_options),
                               list(extend = "csvHtml5",   exportOptions = dt_export_options),
                               list(extend = "excelHtml5", exportOptions = dt_export_options)
                             )))
  })

  output$distrib_pairwise_header <- renderUI({
    res <- tryCatch(distrib_stats(), error = function(e) NULL)
    if (is.null(res) || is.null(res$pairwise) || nrow(res$pairwise) == 0) return(NULL)
    tags$h4(sprintf("Pairwise Wilcoxon (BH-adjusted) — %s",
                    res$taxon %||% ""))
  })

  output$distrib_stats_pairwise <- renderDT({
    res <- distrib_stats()
    if (is.null(res$pairwise) || nrow(res$pairwise) == 0) return(NULL)
    show <- res$pairwise
    show$p_adj <- signif(show$p_adj, 3)
    datatable(show, rownames = FALSE,
              options = list(pageLength = 10, dom = "tp"))
  })

  output$download_distrib_stats_csv <- downloadHandler(
    filename = function() paste0("distribution_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- distrib_stats()
      main <- res$main; if (is.null(main)) main <- data.frame()
      write.csv(main, file, row.names = FALSE)
    }
  )

  # ═══════════════════════════════════════════════════════════════════════
  # TAB: KRONA / SUNBURST (taxplore)
  # ═══════════════════════════════════════════════════════════════════════

  output$krona_sample_selector <- renderUI({
    req(rv$loaded, rv$sample_renames)
    samples <- sort(unique(unname(rv$sample_renames)))
    selectInput("krona_sample", "Select Sample:",
                choices = c("All (sum)" = "__ALL__", samples),
                selected = restored_input("krona_sample", "__ALL__",
                                          valid = c("__ALL__", samples)))
  })

  # Optional metadata grouping for Krona. When a factor is chosen, the chart is
  # built with one dataset per level (each summing its samples) — Krona renders
  # a dataset dropdown to switch between them. Overrides the sample selector.
  output$krona_group_ui <- renderUI({
    cols <- metadata_group_cols()
    if (length(cols) == 0) return(NULL)
    tagList(
      fluidRow(
        column(8,
          selectInput("krona_group_by", "Split into datasets by (metadata):",
                      choices = c("None", cols),
                      selected = restored_input("krona_group_by", "None",
                                                valid = c("None", cols))))
      ),
      conditionalPanel(
        condition = "input.krona_group_by && input.krona_group_by != 'None'",
        div(class = "status-text",
            icon("circle-info"),
            " One Krona dataset per level (sum of its samples) — use the chart's ",
            tags$b("dataset dropdown"), " to switch. The sample selector above is ignored.")
      )
    )
  })

  if (HAS_TAXPLORE) {
    # Show the loading banner whenever the Krona chart (re)generates — on tab
    # entry and on any input that drives the render. shinyjs::delay() holds it
    # for a 5s floor so it stays up while the widget builds its SVG client-side
    # (the freeze the built-in `recalculating` spinner doesn't cover).
    observeEvent(
      list(input$tabs, input$krona_sample, input$krona_value, input$krona_use_abundance,
           input$krona_group_by),
      {
        if (isTRUE(input$tabs == "kronaTab") && isTRUE(rv$loaded)) {
          shinyjs::show("krona_loading")
          shinyjs::delay(5000, shinyjs::hide("krona_loading"))
        }
      },
      ignoreNULL = FALSE
    )

    output$krona_plot <- renderKronaChart({
      req(rv$loaded, rv$raw_data, input$krona_sample)

      df <- rv$raw_data

      # Apply sample renames
      if (!is.null(rv$sample_renames)) {
        df$Sample <- rv$sample_renames[df$Sample]
      }

      # Optional metadata grouping → one dataset per factor level (Krona's
      # dataset dropdown). When active it OVERRIDES the single-sample selector:
      # all samples are kept and split by level. Otherwise honour the selector.
      krona_grp <- resolve_group_col(input$krona_group_by)
      if (!is.null(krona_grp)) {
        meta <- rv$metadata[, c("Sample", krona_grp), drop = FALSE]
        df <- merge(df, meta, by = "Sample", all.x = TRUE)
        df$.level <- as.character(df[[krona_grp]])
        df$.level[is.na(df$.level) | df$.level == ""] <- "NA"
      } else if (input$krona_sample != "__ALL__") {
        df <- df[df$Sample == input$krona_sample, , drop = FALSE]
      }

      req(nrow(df) > 0)

      # Determine value column
      val_col <- input$krona_value
      if (!val_col %in% colnames(df)) val_col <- "percent"

      # ── Build taxonomy matrix using pre-computed lineage columns ──
      rank_labels <- c("Domain", "Kingdom", "Phylum", "Class",
                       "Order", "Family", "Genus", "Species")
      available_ranks <- intersect(rank_labels, colnames(df))

      if (length(available_ranks) > 0) {
        # Use only rows at main ranks (have a lineage), exclude unclassified/root
        rank_map <- c("D" = "Domain", "K" = "Kingdom", "P" = "Phylum",
                      "C" = "Class",  "O" = "Order",   "F" = "Family",
                      "G" = "Genus",  "S" = "Species")
        df_main <- df[df$rank %in% names(rank_map) &
                       df$name != "unclassified", , drop = FALSE]
        req(nrow(df_main) > 0)

        # Pick the per-row magnitude. Krona sums child magnitudes into parents,
        # so we must feed it taxon-only counts — feeding clade-cumulative
        # values (Kraken2's reads_clade / percent) makes every internal node
        # contribute its own count AGAIN on top of its descendants', massively
        # inflating parents in the chart.
        #
        # reads_taxon (Kraken2) / taxReads (KrakenUniq) is taxon-only and is
        # the right input. Bracken sets reads_taxon = reads_clade (only
        # species rows, so equivalent). MetaPhlAn has reads_taxon = 0 and
        # only a `percent` relative abundance — fall back to species-only
        # rows in that case.
        has_reads_taxon <- "reads_taxon" %in% colnames(df_main) &&
          any(suppressWarnings(as.numeric(df_main$reads_taxon)) > 0, na.rm = TRUE)

        if (has_reads_taxon) {
          df_main$.krona_mag <- as.numeric(df_main$reads_taxon)
          # In percent mode, normalize to 100. For a single dataset that is over
          # the whole chart here; for grouped datasets we defer to a per-level
          # normalization after the lineage×level pivot below.
          if (val_col == "percent" && is.null(krona_grp)) {
            tot <- sum(df_main$.krona_mag, na.rm = TRUE)
            if (is.finite(tot) && tot > 0) {
              df_main$.krona_mag <- (df_main$.krona_mag / tot) * 100
            }
          }
        } else {
          # MetaPhlAn-style: keep only deepest-rank rows so the chosen
          # column (relative abundance) doesn't double-count up the tree.
          df_main <- df_main[df_main$rank == "S", , drop = FALSE]
          req(nrow(df_main) > 0)
          df_main$.krona_mag <- as.numeric(df_main[[val_col]])
        }

        # Forward-fill helper: coerce to character and replace each NA with the
        # last known (shallower-rank) ancestor name, e.g.
        # (Bacteria, .., Neisseria, NA) → (Bacteria, .., Neisseria, Neisseria).
        # This avoids NAs entirely (newer taxplore asserts is.character), keeps
        # higher-rank reads inside the correct branch instead of a single
        # "_unclassified_" sibling, and produces single-child chains that
        # Krona's "Collapse" toggle hides automatically. Column-wise loop
        # preserves data.frame structure/names across edge cases (1-row td).
        fill_lineage <- function(td) {
          td[] <- lapply(td, as.character)
          for (i in seq_along(td)) {
            na_mask <- is.na(td[[i]])
            if (i == 1L) td[[i]][na_mask] <- "Unknown"
            else td[[i]][na_mask] <- td[[i - 1L]][na_mask]
          }
          td
        }

        sample_label <- if (input$krona_sample == "__ALL__") "All samples" else input$krona_sample

        if (is.null(krona_grp)) {
          # ── Single dataset (all samples summed, or one selected sample) ──
          # `dplyr::group_by` preserves NA groups, which `aggregate.formula`
          # does NOT (it silently drops rows with any NA grouping var) — that
          # was the cause of "No taxonomy data available" on inputs where every
          # leaf-rank row had at least one NA ancestor.
          agg <- df_main %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(available_ranks))) %>%
            dplyr::summarise(magnitude = sum(.data$.krona_mag, na.rm = TRUE),
                             .groups = "drop") %>%
            as.data.frame()

          tax_df <- agg[, available_ranks, drop = FALSE]
          magnitudes <- agg$magnitude
          non_empty <- sapply(tax_df, function(x) !all(is.na(x)))
          tax_df <- tax_df[, non_empty, drop = FALSE]

          # NOTE: keep `shiny::` explicit. jsonlite's validate() masks shiny's
          # when it lands later in the search path, and treats the `need()`
          # object as a JSON string ("is.character(txt) is not TRUE").
          shiny::validate(
            need(nrow(tax_df) > 0 && ncol(tax_df) > 0,
                 sprintf(paste0("No taxonomy data available for Krona chart ",
                                "(input had %d rows, %d at standard ranks, but ",
                                "no resolvable lineage columns)."),
                         nrow(df), nrow(df_main)))
          )

          tax_df <- fill_lineage(tax_df)
          # Surface the sample name inside the chart: taxplore's `root_label`
          # controls the central-circle text (default "Root").
          krona_opts <- list(root_label = sample_label)

          if (input$krona_use_abundance) {
            mag <- as.numeric(magnitudes)
            mag[!is.finite(mag)] <- 0
            plot_krona(tax_df, mag, dataset_group = sample_label, opts = krona_opts)
          } else {
            plot_krona(tax_df, dataset_group = sample_label, opts = krona_opts)
          }

        } else {
          # ── One dataset per metadata level (Krona dataset dropdown) ──
          # Aggregate per lineage × level, then pivot to a taxa×level matrix so
          # plot_krona() gets one magnitude column per level. dataset_group =
          # the level names makes each column its own switchable dataset.
          agg <- df_main %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c(available_ranks, ".level")))) %>%
            dplyr::summarise(magnitude = sum(.data$.krona_mag, na.rm = TRUE),
                             .groups = "drop") %>%
            as.data.frame()

          wide <- as.data.frame(
            tidyr::pivot_wider(agg, names_from = ".level",
                               values_from = "magnitude", values_fill = 0)
          )

          levels_present <- setdiff(colnames(wide), available_ranks)
          tax_df <- wide[, available_ranks, drop = FALSE]
          non_empty <- sapply(tax_df, function(x) !all(is.na(x)))
          tax_df <- tax_df[, non_empty, drop = FALSE]

          shiny::validate(
            need(nrow(tax_df) > 0 && ncol(tax_df) > 0 && length(levels_present) > 0,
                 "No taxonomy data available for the grouped Krona chart.")
          )

          tax_df <- fill_lineage(tax_df)

          mag_mat <- as.matrix(wide[, levels_present, drop = FALSE])
          mag_mat[!is.finite(mag_mat)] <- 0
          colnames(mag_mat) <- levels_present

          # Per-level percent normalization so each dataset is a 0–100
          # composition (counts mode keeps raw per-level sums).
          if (val_col == "percent") {
            cs <- colSums(mag_mat, na.rm = TRUE)
            for (j in seq_along(cs)) {
              if (is.finite(cs[j]) && cs[j] > 0) mag_mat[, j] <- mag_mat[, j] / cs[j] * 100
            }
          }

          krona_opts <- list(root_label = paste0(krona_grp, " (by level)"))

          if (input$krona_use_abundance) {
            plot_krona(tax_df, mag_mat, dataset_group = levels_present, opts = krona_opts)
          } else {
            # Equal-weight wedges, but keep the per-level matrix so the dataset
            # dropdown still lists every level (presence per level).
            pres <- mag_mat
            pres[pres > 0] <- 1
            plot_krona(tax_df, pres, dataset_group = levels_present, opts = krona_opts)
          }
        }

      } else {
        # Fallback: no lineage columns (non-Kraken2 format)
        # Use simple flat approach
        rank_map <- c(
          "D" = "Domain", "K" = "Kingdom", "P" = "Phylum",
          "C" = "Class",  "O" = "Order",   "F" = "Family",
          "G" = "Genus",  "S" = "Species"
        )
        df_clean <- df[df$rank %in% names(rank_map) &
                        df$name != "unclassified", , drop = FALSE]
        req(nrow(df_clean) > 0)
        agg <- aggregate(
          as.formula(paste(val_col, "~ name + rank")),
          data = df_clean, FUN = sum
        )
        colnames(agg)[3] <- "magnitude"
        agg$rank_label <- rank_map[agg$rank]
        agg <- agg[!is.na(agg$rank_label), , drop = FALSE]

        tax_cols <- c("Domain", "Kingdom", "Phylum", "Class",
                      "Order", "Family", "Genus", "Species")
        tax_df <- data.frame(
          matrix(NA_character_, nrow = nrow(agg), ncol = length(tax_cols)),
          stringsAsFactors = FALSE
        )
        colnames(tax_df) <- tax_cols
        for (i in seq_len(nrow(agg))) {
          col <- agg$rank_label[i]
          if (col %in% tax_cols) tax_df[i, col] <- agg$name[i]
        }
        non_empty <- sapply(tax_df, function(x) !all(is.na(x)))
        tax_df <- tax_df[, non_empty, drop = FALSE]
        shiny::validate(need(nrow(tax_df) > 0, "No data."))

        tax_df[] <- lapply(tax_df, as.character)
        tax_df[is.na(tax_df)] <- "_unclassified_"

        sample_label <- if (input$krona_sample == "__ALL__") "All samples" else input$krona_sample
        krona_opts <- list(root_label = sample_label)

        if (input$krona_use_abundance) {
          mag <- as.numeric(agg$magnitude)
          mag[!is.finite(mag)] <- 0
          plot_krona(tax_df, mag, dataset_group = sample_label, opts = krona_opts)
        } else {
          plot_krona(tax_df, dataset_group = sample_label, opts = krona_opts)
        }
      }
    })
  }

  # ═══════════════════════════════════════════════════════════════════════
  # TAB 4: ALPHA DIVERSITY
  # ═══════════════════════════════════════════════════════════════════════

  output$diversity_color_selector <- renderUI({
    if (!is.null(rv$metadata)) {
      group_cols <- setdiff(colnames(rv$metadata), "Sample")
      selectInput("diversity_color", "Color By:",
                  choices = c("None", group_cols),
                  selected = restored_input("diversity_color", "None",
                                            valid = c("None", group_cols)))
    }
  })

  output$diversity_plot <- renderPlotly({
    df <- filtered_data()
    shiny::validate(shiny::need(nrow(df) > 0, empty_filter_message()))

    # Only species-level (or selected rank) and exclude unclassified
    clean <- df[df$name != "unclassified", ]
    shiny::validate(shiny::need(nrow(clean) > 0, empty_filter_message()))

    # Compute per-sample diversity
    samples <- unique(clean$Sample)
    div_list <- lapply(samples, function(s) {
      sample_data <- clean[clean$Sample == s, ]
      counts <- sample_data$reads_clade
      if (is.null(counts) || all(is.na(counts)) || sum(counts, na.rm = TRUE) == 0) {
        counts <- sample_data$percent
      }
      counts <- counts[!is.na(counts) & counts > 0]
      if (length(counts) == 0) {
        return(data.frame(Sample = s, Shannon = 0, Simpson = 0, InvSimpson = 0))
      }
      div <- compute_diversity(counts)
      div$Sample <- s
      div
    })

    div_df <- bind_rows(div_list)

    # Merge with metadata for coloring
    if (!is.null(rv$metadata)) {
      div_df <- merge(div_df, rv$metadata, by = "Sample", all.x = TRUE)
    }

    index <- input$diversity_index
    y_col <- switch(index,
      "Shannon" = "Shannon",
      "Simpson" = "Simpson",
      "Inverse Simpson" = "InvSimpson"
    )

    div_df <- div_df[order(div_df[[y_col]]), ]
    div_df$x_order <- seq_len(nrow(div_df))

    color_var <- if (!is.null(input$diversity_color) && input$diversity_color != "None" &&
                     input$diversity_color %in% colnames(div_df)) {
      input$diversity_color
    } else {
      NULL
    }

    div_df$y_value <- div_df[[y_col]]

    if (!is.null(color_var)) {
      div_df$color_val <- div_df[[color_var]]
      p <- ggplot(div_df, aes(x = reorder(Sample, y_value),
                               y = y_value,
                               fill = color_val,
                               text = paste0("Sample: ", Sample,
                                             "\n", index, ": ", sprintf("%.3f", y_value),
                                             "\n", color_var, ": ", color_val)))
      pal_colors <- get_palette_colors(length(unique(div_df[[color_var]])), input$color_palette)
      p <- p + scale_fill_manual(values = pal_colors)
    } else {
      p <- ggplot(div_df, aes(x = reorder(Sample, y_value),
                               y = y_value,
                               text = paste0("Sample: ", Sample,
                                             "\n", index, ": ", sprintf("%.3f", y_value))))
      p <- p + geom_col(fill = "#2980b9", alpha = 0.85)
    }

    p <- p +
      geom_col(width = 0.75, alpha = 0.85) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(face = "bold", size = 16),
        panel.grid.major.x = element_blank()
      ) +
      labs(
        title = paste(index, "Diversity Index"),
        x = NULL,
        y = index
      )

    ggplotly(p, tooltip = "text") %>%
      layout(margin = list(b = 120))
  })

  # ═══════════════════════════════════════════════════════════════════════
  # TAB: RAREFACTION CURVES
  # ═══════════════════════════════════════════════════════════════════════

  output$rarefaction_color_selector <- renderUI({
    if (!is.null(rv$metadata)) {
      group_cols <- setdiff(colnames(rv$metadata), "Sample")
      selectInput("rarefy_color", "Color By:",
                  choices = c("None", group_cols),
                  selected = restored_input("rarefy_color", "None",
                                            valid = c("None", group_cols)))
    }
  })

  #' Integer sample-by-taxa counts matrix at the selected rank.
  #' Rarefaction is only meaningful on raw read counts, so this returns NULL
  #' when the active dataset carries relative abundances only (e.g. MetaPhlAn,
  #' where reads_clade is 0). Empty (all-zero) samples are dropped because
  #' rarecurve() errors on them.
  rarefaction_matrix <- reactive({
    df <- filtered_data()
    shiny::validate(shiny::need(nrow(df) > 0, empty_filter_message()))
    if (!"reads_clade" %in% colnames(df)) return(NULL)

    clean <- df[df$name != "unclassified" & df$name != "Other", ]
    if (nrow(clean) == 0) return(NULL)
    if (sum(clean$reads_clade, na.rm = TRUE) <= 0) return(NULL)

    wide <- clean %>%
      select(Sample, name, reads_clade) %>%
      group_by(Sample, name) %>%
      summarise(reads_clade = sum(reads_clade, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = name, values_from = reads_clade, values_fill = 0)

    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$Sample
    mat <- round(mat)
    storage.mode(mat) <- "integer"
    mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
    mat
  })

  #' Heavy core: run rarecurve() once and return the tidy curve data plus the
  #' plotting metadata. Shared by the interactive render and the PDF/PNG exports.
  rarefaction_core <- reactive({
    # NB: use shiny::validate explicitly — if jsonlite is attached later its
    # validate() masks shiny's and treats the need() result as a JSON string
    # (stopifnot(is.character(txt)) → "is.character(txt) is not TRUE").
    shiny::validate(need(HAS_VEGAN,
                  "The 'vegan' package is required for rarefaction curves."))
    mat <- rarefaction_matrix()
    shiny::validate(need(!is.null(mat) && nrow(mat) > 0,
                  paste("Rarefaction requires read counts. The current dataset",
                        "has relative abundances only (e.g. MetaPhlAn profiles),",
                        "so rarefaction curves are not available.")))

    # Shallowest sample = the conventional rarefaction cutoff: rarefying every
    # sample to this depth keeps all samples comparable without discarding any.
    row_tot <- rowSums(mat, na.rm = TRUE)
    min_depth <- min(row_tot)
    normalize <- isTRUE(input$rarefy_normalize)

    # rarecurve() evaluates richness at seq(1, total_reads, by = step). Kraken/
    # Bracken read counts run to tens of millions, so a fixed small step would
    # generate hundreds of thousands of points per sample and hang the app.
    # Derive the step from the depth we actually plot to so every curve gets
    # ~`points` evaluations: the shallowest sample when capping, else the
    # largest library.
    points <- if (!is.null(input$rarefy_points) && input$rarefy_points >= 2) input$rarefy_points else 150
    plot_max <- if (normalize) min_depth else max(row_tot)
    step <- max(1, ceiling(plot_max / points))

    # rarecurve() draws to the active graphics device as a side effect; send
    # that to a null device so nothing leaks onto the real one (cf. dendro_gg).
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    rare <- vegan::rarecurve(mat, step = step, label = FALSE)

    samp_names <- rownames(mat)
    curve_df <- bind_rows(lapply(seq_along(rare), function(i) {
      data.frame(
        Sample   = samp_names[i],
        Depth    = as.numeric(attr(rare[[i]], "Subsample")),
        Richness = as.numeric(rare[[i]]),
        stringsAsFactors = FALSE
      )
    }))

    # Cap every curve at the shallowest sample's depth so all samples are
    # compared on a common, fully-supported x-range (the dashed line marks it).
    if (normalize) {
      curve_df <- curve_df[curve_df$Depth <= min_depth, , drop = FALSE]
    }

    rank_label <- names(rank_choices)[match(input$rank_select, rank_choices)]
    if (length(rank_label) == 0 || is.na(rank_label)) rank_label <- "taxa"
    palette_name <- if (!is.null(input$rarefy_palette)) input$rarefy_palette else "Set3"

    color_var <- if (!is.null(input$rarefy_color) && input$rarefy_color != "None" &&
                     !is.null(rv$metadata) && input$rarefy_color %in% colnames(rv$metadata)) {
      input$rarefy_color
    } else {
      NULL
    }

    if (!is.null(color_var)) {
      curve_df <- merge(curve_df, rv$metadata[, c("Sample", color_var)],
                        by = "Sample", all.x = TRUE)
      curve_df$grp <- curve_df[[color_var]]
      pal_colors <- get_palette_colors(length(unique(na.omit(curve_df$grp))), palette_name)
    } else {
      pal_colors <- get_palette_colors(length(unique(curve_df$Sample)), palette_name)
    }

    list(curve_df = curve_df, color_var = color_var, pal_colors = pal_colors,
         min_depth = min_depth, normalize = normalize, rank_label = rank_label)
  })

  #' Base ggplot of just the curves. `with_ref` adds the dashed min-depth line
  #' and its label as ggplot layers — used ONLY for the static PDF/PNG exports.
  #' The interactive plot deliberately keeps the reference line OUT of the
  #' ggplotly() conversion: a `y = Inf` annotation and an aesthetic-less
  #' geom_vline break older plotly builds with an opaque jsonlite error
  #' ("is.character(txt) is not TRUE"). The interactive render draws the line
  #' with plotly layout shapes/annotations instead (see output below).
  build_rarefaction_gg <- function(core, with_ref) {
    cd <- core$curve_df
    if (!is.null(core$color_var)) {
      p <- ggplot(cd, aes(x = Depth, y = Richness, group = Sample, color = grp,
            text = paste0("Sample: ", Sample, "\nDepth: ", round(Depth),
                          "\nRichness: ", round(Richness, 1),
                          "\n", core$color_var, ": ", grp))) +
        scale_color_manual(values = core$pal_colors) + labs(color = core$color_var)
    } else {
      p <- ggplot(cd, aes(x = Depth, y = Richness, group = Sample, color = Sample,
            text = paste0("Sample: ", Sample, "\nDepth: ", round(Depth),
                          "\nRichness: ", round(Richness, 1)))) +
        scale_color_manual(values = core$pal_colors)
    }
    p <- p + geom_line(linewidth = 0.7, alpha = 0.85)
    if (with_ref) {
      ymax <- max(cd$Richness, na.rm = TRUE)
      p <- p +
        geom_vline(xintercept = core$min_depth, linetype = "dashed",
                   color = "#7d6608", linewidth = 0.5) +
        annotate("text", x = core$min_depth, y = ymax,
                 label = paste0("min depth = ", format(round(core$min_depth), big.mark = ",")),
                 hjust = -0.05, vjust = 1, size = 3.2, color = "#7d6608")
    }
    p +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold", size = 16)) +
      labs(title = if (core$normalize) "Rarefaction Curves (capped at shallowest sample)" else "Rarefaction Curves",
           x = "Sequencing depth (reads sampled)",
           y = paste0("Observed ", tolower(core$rank_label), " richness"))
  }

  # Static (PDF/PNG): full ggplot including the reference line.
  rarefaction_gg <- reactive(build_rarefaction_gg(rarefaction_core(), with_ref = TRUE))

  output$rarefaction_plot <- renderPlotly({
    core <- rarefaction_core()
    gp <- ggplotly(build_rarefaction_gg(core, with_ref = FALSE), tooltip = "text")
    # Reference line + label added via plotly so they never enter the ggplotly
    # JSON path (see note on build_rarefaction_gg). yref="paper" spans the panel
    # height without any Inf coordinate.
    gp <- gp %>%
      layout(
        margin = list(b = 80),
        annotations = list(list(x = core$min_depth, yref = "paper", y = 1,
                                text = paste0("min depth = ",
                                              format(round(core$min_depth), big.mark = ",")),
                                showarrow = FALSE, xanchor = "left", yanchor = "top",
                                font = list(color = "#7d6608", size = 11)))
      )
    # layout(shapes=) is silently dropped by ggplotly's build, so assign the
    # dashed reference line directly on the widget, where it persists.
    gp$x$layout$shapes <- list(list(type = "line",
                                    x0 = core$min_depth, x1 = core$min_depth,
                                    yref = "paper", y0 = 0, y1 = 1,
                                    line = list(color = "#7d6608", dash = "dash", width = 1)))
    gp
  })

  output$download_rarefaction_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_rarefaction_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) save_gg(file, rarefaction_gg(), "pdf")
  )
  output$download_rarefaction_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_rarefaction_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_gg(file, rarefaction_gg(), "png")
  )
  output$download_rarefaction_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_rarefaction_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) save_gg(file, rarefaction_gg(), "svg")
  )

  # ═══════════════════════════════════════════════════════════════════════
  # TAB 5: PCA / BETA DIVERSITY
  # ═══════════════════════════════════════════════════════════════════════

  output$pca_color_selector <- renderUI({
    if (!is.null(rv$metadata)) {
      group_cols <- setdiff(colnames(rv$metadata), "Sample")
      selectInput("pca_color", "Color By:",
                  choices = c("None", group_cols),
                  selected = restored_input("pca_color", "None",
                                            valid = c("None", group_cols)))
    }
  })

  # Build sample-by-taxa matrix
  taxa_matrix <- reactive({
    df <- filtered_data()
    shiny::validate(shiny::need(nrow(df) > 0, empty_filter_message()))

    clean <- df[df$name != "unclassified" & df$name != "Other", ]
    shiny::validate(shiny::need(nrow(clean) > 0, empty_filter_message()))

    wide <- clean %>%
      select(Sample, name, percent) %>%
      group_by(Sample, name) %>%
      summarise(percent = sum(percent, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = name, values_from = percent, values_fill = 0)

    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$Sample
    mat
  })

  # PCA plot
  output$pca_plot <- renderPlotly({
    mat <- taxa_matrix()
    req(nrow(mat) >= 3, ncol(mat) >= 2)

    # Remove zero-variance columns
    col_var <- apply(mat, 2, var)
    mat_pca <- mat[, col_var > 0, drop = FALSE]
    req(ncol(mat_pca) >= 2)

    pca_res <- prcomp(mat_pca, center = TRUE, scale. = TRUE)
    var_exp <- round(summary(pca_res)$importance[2, ] * 100, 1)

    pca_df <- data.frame(
      PC1 = pca_res$x[, 1],
      PC2 = pca_res$x[, 2],
      Sample = rownames(mat)
    )

    # Merge metadata
    if (!is.null(rv$metadata)) {
      pca_df <- merge(pca_df, rv$metadata, by = "Sample", all.x = TRUE)
    }

    color_var <- if (!is.null(input$pca_color) && input$pca_color != "None" &&
                     input$pca_color %in% colnames(pca_df)) {
      input$pca_color
    } else {
      NULL
    }

    if (!is.null(color_var)) {
      pca_df$color_val <- pca_df[[color_var]]
      p <- ggplot(pca_df, aes(x = PC1, y = PC2,
                               color = color_val,
                               text = paste0("Sample: ", Sample,
                                             "\n", color_var, ": ", color_val,
                                             "\nPC1: ", round(PC1, 2),
                                             "\nPC2: ", round(PC2, 2)))) +
        scale_color_manual(values = get_palette_colors(
          length(unique(pca_df[[color_var]])), input$color_palette))
    } else {
      p <- ggplot(pca_df, aes(x = PC1, y = PC2,
                               text = paste0("Sample: ", Sample,
                                             "\nPC1: ", round(PC1, 2),
                                             "\nPC2: ", round(PC2, 2))))
    }

    p <- p +
      geom_point(size = 4, alpha = 0.8) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(fill = NA, color = "#d5dbdb")
      ) +
      labs(
        title = "PCA Ordination",
        x = paste0("PC1 (", var_exp[1], "%)"),
        y = paste0("PC2 (", var_exp[2], "%)")
      )

    ggplotly(p, tooltip = "text")
  })

  # Dendrogram
  output$dendro_plot <- renderPlotly({
    mat <- taxa_matrix()
    req(nrow(mat) >= 3)

    dist_method <- input$pca_dist_method
    if (dist_method == "bray" || dist_method == "jaccard") {
      if (HAS_VEGAN) {
        d <- vegdist(mat, method = dist_method)
      } else {
        d <- dist(mat, method = "euclidean")
      }
    } else {
      d <- dist(mat, method = "euclidean")
    }

    hc <- hclust(d, method = "ward.D2")

    # Prevent ggplot from trying to open a graphics device file
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)

    dend_data <- dendro_data(as.dendrogram(hc))
    seg <- dend_data$segments
    lab <- dend_data$labels
    lab$label <- as.character(lab$label)

    # Merge metadata for coloring leaf nodes
    if (!is.null(rv$metadata)) {
      lab <- merge(lab, rv$metadata, by.x = "label", by.y = "Sample", all.x = TRUE)
    }

    color_var <- if (!is.null(input$pca_color) && input$pca_color != "None" &&
                     input$pca_color %in% colnames(lab)) {
      input$pca_color
    } else {
      NULL
    }

    p <- ggplot() +
      geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend),
                   color = "#5d6d7e", linewidth = 0.5)

    if (!is.null(color_var)) {
      lab$color_val <- lab[[color_var]]
      p <- p + geom_point(data = lab, aes(x = x, y = y, color = color_val,
                                           text = label), size = 3) +
        scale_color_manual(values = get_palette_colors(
          length(unique(lab[[color_var]])), input$color_palette))
    } else {
      p <- p + geom_point(data = lab, aes(x = x, y = y, text = label),
                          size = 3, color = "#2980b9")
    }

    p <- p +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 14)
      ) +
      labs(title = paste("Hierarchical Clustering (", dist_method, ")"),
           x = "", y = "Distance")

    ggplotly(p, tooltip = "text")
  })

  # ═══════════════════════════════════════════════════════════════════════
  # TAB 6: LIFEMAPR TREE
  # ═══════════════════════════════════════════════════════════════════════

  output$lifemap_warning <- renderUI({
    df <- filtered_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    taxids <- unique(df$taxid)
    taxids <- taxids[!is.na(taxids) & taxids > 0]
    if (length(taxids) == 0) {
      div(class = "status-text status-warning",
          icon("exclamation-triangle"),
          " No valid NCBI TaxIDs found. LifemapR requires NCBI TaxIDs ",
          "(not GTDB). Please check your Kraken2 database.")
    }
  })

  observeEvent(input$run_lifemap, {
    req(HAS_LIFEMAPR)
    df <- filtered_data()
    if (nrow(df) == 0) {
      showNotification(empty_filter_message(), type = "warning", duration = 8)
      return()
    }

    clean <- df[df$name != "unclassified" & !is.na(df$taxid) & df$taxid > 0, ]
    if (nrow(clean) == 0) {
      showNotification(
        paste0("No taxa with a valid TaxID at rank '",
               names(rank_choices)[match(input$rank_select, rank_choices)] %||% input$rank_select,
               "'. LifemapR needs NCBI TaxIDs."),
        type = "warning", duration = 8
      )
      return()
    }

    taxid_df <- aggregate(percent ~ taxid + name, data = clean, FUN = mean)
    colnames(taxid_df) <- c("taxid", "name", "abundance")
    taxid_df$taxid <- as.integer(taxid_df$taxid)
    taxid_df <- taxid_df[!is.na(taxid_df$name) & nchar(taxid_df$name) > 0, , drop = FALSE]
    req(nrow(taxid_df) > 0)

    withProgress(message = "Building Lifemap tree...", value = 0.3, {
      tryCatch({
        build_output <- capture.output({
          lm_obj <- build_Lifemap(df = taxid_df, verbose = TRUE)
        }, type = "message")
        incProgress(0.5, detail = "Rendering...")

        not_found_lines <- grep("could not be found", build_output, value = TRUE)
        n_not_found <- length(not_found_lines)
        n_total <- nrow(taxid_df)
        n_found <- n_total - n_not_found

        if (n_found < 2) {
          showNotification(
            paste0("Only ", n_found, " of ", n_total,
                   " TaxIDs were resolved. Not enough for a tree."),
            type = "error", duration = 15
          )
          rv$lifemap_obj <- NULL
          return()
        }

        rv$lifemap_obj <- lm_obj
        rv$lifemap_df <- taxid_df

        warn_msg <- if (n_not_found > 0) {
          paste0(" (", n_not_found, " TaxIDs not found in NCBI)")
        } else ""

        showNotification(
          paste0("LifemapR tree built with ", n_found,
                 " of ", n_total, " taxa", warn_msg),
          type = "message", duration = 8
        )
      }, error = function(e) {
        showNotification(
          paste0("LifemapR error: ", e$message),
          type = "error", duration = 10
        )
        rv$lifemap_obj <- NULL
      })
    })
  })


  # Status text + embedded map for lifemap
  output$lifemap_status <- renderUI({
    req(rv$lifemap_obj)

    lm <- rv$lifemap_obj
    lm_df <- lm$df

    # Use display_map() which returns a standard leaflet htmlwidget
    # (unlike lifemap() which returns a full Shiny app that can't be serialized)
    vis <- tryCatch({
      m <- display_map(df = lm_df)

      # Filter to only requested taxa (not ancestors)
      req_df <- lm_df[lm_df$type == "requested", , drop = FALSE]
      if (nrow(req_df) == 0) req_df <- lm_df

      # Scale abundance for marker radius
      if ("abundance" %in% colnames(req_df) && any(!is.na(req_df$abundance))) {
        ab <- req_df$abundance
        ab[is.na(ab)] <- 0
        req_df$marker_radius <- scales::rescale(ab, to = c(5, 25))
      } else {
        req_df$marker_radius <- 10
      }

      # Build popup text
      req_df$popup_text <- paste0(
        "<b>", req_df$sci_name, "</b>",
        if ("name" %in% colnames(req_df)) paste0("<br>Name: ", req_df$name) else "",
        if ("abundance" %in% colnames(req_df)) paste0("<br>Abundance: ", round(req_df$abundance, 3), "%") else "",
        "<br>TaxID: ", req_df$taxid
      )

      # Add circle markers via leaflet
      leaflet::addCircleMarkers(
        m,
        lng = req_df$lon,
        lat = req_df$lat,
        radius = req_df$marker_radius,
        fillColor = "steelblue",
        fillOpacity = 0.7,
        stroke = TRUE,
        color = "#2c3e50",
        weight = 1,
        popup = req_df$popup_text,
        label = req_df$sci_name
      )
    }, error = function(e) {
      message("LifemapR display_map error: ", e$message)
      NULL
    })

    if (is.null(vis)) {
      return(div(class = "status-text status-warning",
                 icon("exclamation-triangle"),
                 " Could not render Lifemap widget."))
    }

    # Save the leaflet widget to temp HTML
    tmp_dir <- file.path(tempdir(), "lifemap_widget")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    tmp_file <- file.path(tmp_dir, "lifemap.html")

    save_ok <- tryCatch({
      htmlwidgets::saveWidget(vis, tmp_file, selfcontained = FALSE,
                               libdir = file.path(tmp_dir, "lib"))
      addResourcePath("lifemap_tmp", tmp_dir)
      TRUE
    }, error = function(e) {
      message("saveWidget error: ", e$message)
      FALSE
    })

    if (!save_ok) {
      return(div(class = "status-text status-warning",
                 icon("exclamation-triangle"),
                 " Could not save Lifemap widget. Check container permissions."))
    }

    # Show which samples were included
    sample_names <- unique(filtered_data()$Sample)
    sample_info <- if (length(sample_names) > 1) {
      paste0("Averaged across ", length(sample_names), " samples: ",
             paste(sample_names, collapse = ", "))
    } else {
      paste0("Sample: ", sample_names)
    }

    tagList(
      div(class = "status-text",
          icon("check-circle"),
          " Tree generated on Lifemap."),
      tags$small(style = "color: #666;", sample_info),
      br(), br(),
      fluidRow(
        column(4, downloadButton("download_lifemap_html", "HTML",
                                  icon = icon("download"))),
        column(4, downloadButton("download_lifemap_csv", "CSV",
                                  icon = icon("table"))),
        column(4, downloadButton("download_lifemap_script", "R Script",
                                  icon = icon("file-code")))
      ),
      br(),
      tags$iframe(
        src = "lifemap_tmp/lifemap.html",
        width = "100%", height = "700px",
        frameborder = "0",
        style = "border: 1px solid #ddd; border-radius: 4px;"
      )
    )
  })

  # Download handler: export HTML as zip (HTML + lib folder)
  output$download_lifemap_html <- downloadHandler(
    filename = function() {
      paste0("lifemap_tree_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      tmp_dir <- file.path(tempdir(), "lifemap_widget")
      req(file.exists(file.path(tmp_dir, "lifemap.html")))
      owd <- setwd(tmp_dir)
      on.exit(setwd(owd))
      files_to_zip <- list.files(".", recursive = TRUE)
      zip(file, files = files_to_zip)
    }
  )

  # Download handler: export CSV data
  output$download_lifemap_csv <- downloadHandler(
    filename = function() {
      paste0("lifemap_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(rv$lifemap_obj)
      write.csv(rv$lifemap_obj$df, file, row.names = FALSE)
    }
  )

  # Download handler: export R script for full interactive lifemap
  output$download_lifemap_script <- downloadHandler(
    filename = function() {
      paste0("lifemap_script_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R")
    },
    content = function(file) {
      req(rv$lifemap_obj)
      writeLines(c(
        "# LifemapR visualisation script",
        "# Generated by exploreMetaTax",
        paste0("# Date: ", Sys.time()),
        "",
        "# First download the CSV using the CSV button and place it",
        "# in the same directory as this script, renamed to lifemap_data.csv",
        "",
        "library(LifemapR)",
        "",
        'taxid_df <- read.csv("lifemap_data.csv")',
        'lm_obj <- build_Lifemap(df = taxid_df)',
        '',
        '# Interactive visualisation (opens in browser)',
        'lifemap(lm_obj) +',
        '  lm_markers(',
        '    radius = "abundance",',
        '    var_fillColor = "abundance",',
        '    fillColor = "PiYG",',
        '    popup = "name"',
        '  )'
      ), file)
    }
  )

  # ═══════════════════════════════════════════════════════════════════════
  # TAB 7: DATA TABLE
  # ═══════════════════════════════════════════════════════════════════════

  output$full_data_table <- renderDT({
    df <- merged_data()
    req(nrow(df) > 0)

    datatable(
      df,
      class = "cell-border stripe hover compact",
      filter = "top",
      extensions = "Buttons",
      options = list(
        pageLength = 20,
        lengthMenu = list(c(10, 20, 50, 100, 250, 500), c('10', '20', '50', '100', '250', '500')),
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "Blfrtip",
        buttons = list(
          list(extend = "copyHtml5",  exportOptions = dt_export_options),
          list(extend = "csvHtml5",   exportOptions = dt_export_options),
          list(extend = "excelHtml5", exportOptions = dt_export_options)
        )
      )
    )
  })

  # ═══════════════════════════════════════════════════════════════════════
  # DOWNLOADS
  # ═══════════════════════════════════════════════════════════════════════

  # Active plot for download (from composition tab)
  current_gg_plot <- reactive({
    comp <- composition_reactive()
    plot_df <- comp$plot_df
    all_colors <- comp$all_colors
    group_col <- comp$group_col
    req(nrow(plot_df) > 0)

    colors_to_use <- all_colors[names(all_colors) %in% levels(plot_df$name_grouped)]

    p <- ggplot(plot_df, aes(x = Sample, y = percent, fill = name_grouped)) +
      geom_bar(stat = "identity", width = 0.8) +
      scale_fill_manual(values = colors_to_use, drop = FALSE) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Taxonomic Composition (Rank:",
                         names(rank_choices[rank_choices == input$rank_select]), ")"),
           x = NULL, y = "% Reads", fill = "Taxon")

    if (input$show_labels) {
      p <- p + geom_text(data = subset(plot_df, percent > 1),
                          aes(label = sprintf("%.1f", percent)),
                          position = position_stack(vjust = 0.5),
                          size = 2.8, color = "black")
    }

    if (!is.null(group_col)) {
      p <- p + facet_grid(cols = vars(!!sym(group_col)), scales = "free_x", space = "free_x")
    }

    p
  })

  output$download_pdf <- downloadHandler(
    filename = function() {
      paste0("exploreMetaTax_composition_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = current_gg_plot(),
             width = input$plot_width, height = input$plot_height, device = "pdf")
    }
  )

  output$download_png <- downloadHandler(
    filename = function() {
      paste0("exploreMetaTax_composition_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      ggsave(file, plot = current_gg_plot(),
             width = input$plot_width, height = input$plot_height,
             dpi = 300, device = "png")
    }
  )

  output$download_svg <- downloadHandler(
    filename = function() {
      paste0("exploreMetaTax_composition_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg")
    },
    content = function(file) {
      save_gg(file, current_gg_plot(), "svg",
              width = input$plot_width, height = input$plot_height)
    }
  )

  output$download_csv <- downloadHandler(
    filename = function() {
      paste0("exploreMetaTax_data_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(merged_data(), file, row.names = FALSE)
    }
  )

  output$download_table_csv <- downloadHandler(
    filename = function() {
      paste0("exploreMetaTax_full_table_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(merged_data(), file, row.names = FALSE)
    }
  )

  # ── Helper: save a ggplot object as PDF, PNG, or SVG ──
  # SVG uses the svglite backend when available (vector text stays crisp);
  # ggsave falls back to grDevices::svg otherwise. dpi only matters for PNG.
  save_gg <- function(file, plot_expr, device, width = 10, height = 7) {
    p <- tryCatch(plot_expr, error = function(e) NULL)
    req(p)
    dev <- if (device == "svg" && requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite
    } else {
      device
    }
    ggsave(file, plot = p, width = width, height = height,
           dpi = if (device == "png") 300 else 72, device = dev)
  }

  # ── Heatmap downloads (native plotly → use ggplot heatmap recreation) ──
  # Mirrors the on-screen view, including the optional metadata facet, by
  # reusing heatmap_prep() and the shared build_heatmap_gg() builder.
  heatmap_gg <- reactive({
    build_heatmap_gg(heatmap_prep())
  })
  output$download_heatmap_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) save_gg(file, heatmap_gg(), "pdf")
  )
  output$download_heatmap_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_gg(file, heatmap_gg(), "png")
  )
  output$download_heatmap_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) save_gg(file, heatmap_gg(), "svg")
  )

  # ── Distribution downloads (mirrors the actual displayed plot) ──
  distrib_gg <- reactive({
    df <- filtered_data()
    req(nrow(df) > 0)
    # filtered_data() already carries rv$metadata; re-merging would suffix the
    # shared columns (.x/.y) and break group_var matching (see violin_plot).
    clean <- df[df$name != "unclassified", , drop = FALSE]
    req(nrow(clean) > 0)

    view_mode <- if (!is.null(input$violin_view_mode)) input$violin_view_mode else "single"
    plot_type <- input$violin_plot_type
    group_var <- if (!is.null(input$violin_group) && input$violin_group != "None" &&
                     input$violin_group %in% colnames(df)) input$violin_group else NULL

    if (view_mode == "faceted") {
      taxa_means <- aggregate(percent ~ name, data = clean, mean)
      taxa_means <- taxa_means[order(taxa_means$percent, decreasing = TRUE), ]
      top_n_val <- if (!is.null(input$top_n)) input$top_n else 10
      top_taxa <- head(taxa_means$name, top_n_val)
      plot_data <- df[df$name %in% top_taxa, , drop = FALSE]
      req(nrow(plot_data) > 0)
      plot_data$name <- factor(plot_data$name, levels = top_taxa)

      if (!is.null(group_var)) {
        p <- ggplot(plot_data, aes(x = .data[[group_var]], y = percent, fill = .data[[group_var]]))
      } else {
        plot_data$All <- "All"
        p <- ggplot(plot_data, aes(x = All, y = percent))
      }

      if (plot_type == "violin") {
        p <- p + geom_violin(alpha = 0.6) + geom_jitter(width = 0.15, size = 1.5, alpha = 0.6)
      } else if (plot_type == "box") {
        p <- p + geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          geom_jitter(width = 0.15, size = 1.5, alpha = 0.6)
      } else {
        p <- p + geom_jitter(width = 0.2, size = 2, alpha = 0.7)
      }

      p + facet_wrap(~name, scales = "free_y", ncol = 3) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = if (is.null(group_var)) "none" else "right") +
        labs(title = paste("Distribution of Top", top_n_val, "Taxa"),
             x = if (!is.null(group_var)) group_var else "", y = "Abundance (%)")

    } else {
      taxon <- input$violin_taxon
      req(taxon)
      taxon_data <- df[df$name == taxon, , drop = FALSE]
      req(nrow(taxon_data) > 0)

      if (!is.null(group_var)) {
        p <- ggplot(taxon_data, aes(x = .data[[group_var]], y = percent, fill = .data[[group_var]]))
      } else {
        taxon_data$All <- "All Samples"
        p <- ggplot(taxon_data, aes(x = All, y = percent))
      }

      if (plot_type == "violin") {
        p <- p + geom_violin(alpha = 0.6) +
          geom_boxplot(width = 0.15, alpha = 0.4) +
          geom_jitter(width = 0.1, size = 2, alpha = 0.6)
      } else if (plot_type == "box") {
        p <- p + geom_boxplot(alpha = 0.6, outlier.shape = NA) +
          geom_jitter(width = 0.15, size = 2, alpha = 0.6)
      } else {
        p <- p + geom_jitter(width = 0.15, size = 2, alpha = 0.7)
      }

      p + theme_minimal(base_size = 14) +
        theme(legend.position = if (is.null(group_var)) "none" else "right") +
        labs(title = paste("Distribution:", taxon),
             y = "Abundance (%)",
             x = if (!is.null(group_var)) group_var else "")
    }
  })
  output$download_distrib_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) save_gg(file, distrib_gg(), "pdf")
  )
  output$download_distrib_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_gg(file, distrib_gg(), "png")
  )
  output$download_distrib_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_distribution_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) save_gg(file, distrib_gg(), "svg")
  )

  # ── Alpha Diversity downloads (ggplotly based) ──
  alpha_gg <- reactive({
    df <- filtered_data()
    req(nrow(df) > 0)
    clean <- df[df$name != "unclassified", ]
    req(nrow(clean) > 0)
    samples <- unique(clean$Sample)
    div_list <- lapply(samples, function(s) {
      sample_data <- clean[clean$Sample == s, ]
      counts <- sample_data$reads_clade
      if (is.null(counts) || all(is.na(counts)) || sum(counts, na.rm = TRUE) == 0) {
        counts <- sample_data$percent
      }
      counts <- counts[!is.na(counts) & counts > 0]
      if (length(counts) == 0) return(data.frame(Sample = s, Shannon = 0, Simpson = 0, InvSimpson = 0))
      div <- compute_diversity(counts)
      div$Sample <- s
      div
    })
    div_df <- bind_rows(div_list)
    index <- input$diversity_index
    y_col <- switch(index, "Shannon" = "Shannon", "Simpson" = "Simpson", "Inverse Simpson" = "InvSimpson")
    div_df <- div_df[order(div_df[[y_col]]), ]
    div_df$y_value <- div_df[[y_col]]
    ggplot(div_df, aes(x = reorder(Sample, y_value), y = y_value)) +
      geom_col(fill = "#2980b9", alpha = 0.85) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Alpha Diversity:", index), x = NULL, y = index)
  })
  output$download_alpha_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_alpha_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) save_gg(file, alpha_gg(), "pdf")
  )
  output$download_alpha_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_alpha_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_gg(file, alpha_gg(), "png")
  )
  output$download_alpha_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_alpha_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) save_gg(file, alpha_gg(), "svg")
  )

  # ── PCA downloads (ggplotly based) ──
  pca_gg <- reactive({
    mat <- taxa_matrix()
    req(nrow(mat) >= 3)
    pca_res <- prcomp(mat, scale. = TRUE, center = TRUE)
    pc_df <- as.data.frame(pca_res$x[, 1:min(2, ncol(pca_res$x))])
    pc_df$Sample <- rownames(pc_df)
    var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
    ggplot(pc_df, aes(x = PC1, y = PC2, label = Sample)) +
      geom_point(size = 3, color = "#2980b9") +
      geom_text(vjust = -0.5, size = 3) +
      theme_minimal(base_size = 14) +
      labs(title = "PCA Ordination",
           x = paste0("PC1 (", var_exp[1], "%)"),
           y = paste0("PC2 (", var_exp[2], "%)"))
  })
  output$download_pca_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_pca_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) save_gg(file, pca_gg(), "pdf")
  )
  output$download_pca_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_pca_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_gg(file, pca_gg(), "png")
  )
  output$download_pca_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_pca_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) save_gg(file, pca_gg(), "svg")
  )

  # ── Dendro downloads (ggplotly based) ──
  dendro_gg <- reactive({
    mat <- taxa_matrix()
    req(nrow(mat) >= 3)
    dist_method <- input$pca_dist_method
    if (dist_method %in% c("bray", "jaccard") && HAS_VEGAN) {
      d <- vegdist(mat, method = dist_method)
    } else {
      d <- dist(mat, method = "euclidean")
    }
    hc <- hclust(d, method = "ward.D2")
    pdf(NULL); on.exit(dev.off(), add = TRUE)
    dend_data <- dendro_data(as.dendrogram(hc))
    seg <- dend_data$segments
    lab <- dend_data$labels
    lab$label <- as.character(lab$label)
    ggplot() +
      geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend),
                   color = "#5d6d7e", linewidth = 0.5) +
      geom_point(data = lab, aes(x = x, y = y), size = 3, color = "#2980b9") +
      geom_text(data = lab, aes(x = x, y = y, label = label), vjust = 2, size = 3) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            panel.grid = element_blank()) +
      labs(title = paste("Hierarchical Clustering (", dist_method, ")"),
           x = "", y = "Distance")
  })
  output$download_dendro_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_dendro_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) save_gg(file, dendro_gg(), "pdf")
  )
  output$download_dendro_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_dendro_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_gg(file, dendro_gg(), "png")
  )
  output$download_dendro_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_dendro_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) save_gg(file, dendro_gg(), "svg")
  )

  # ═══════════════════════════════════════════════════════════════════════
  # HUMAnN MODE — Info / per-table sub-tabs / Method
  # ═══════════════════════════════════════════════════════════════════════

  # ── Plot-only sample filter + rename controls ────────────────────────
  # These feed humann_slot_gg / humann_strat_gg via humann_effective_samples()
  # and humann_display(). The Full-table DTs deliberately do NOT consume
  # these — they always show the raw joined data.
  humann_all_samples <- reactive({
    req(rv$app_mode == "humann", rv$humann)
    unique(unlist(lapply(rv$humann, function(df) setdiff(colnames(df), "Feature")),
                  use.names = FALSE))
  })

  output$humann_sample_filter_ui <- renderUI({
    samples <- humann_all_samples()
    if (length(samples) == 0) return(NULL)
    # Preserve any current visible/hidden split across re-renders. Priority:
    # .rds restore state (isolated so its post-flush clear doesn't
    # re-invalidate this UI) → last-reported widget value → all-visible
    # default. New samples default to visible, appended at the end.
    pending <- isolate(rv$pending_ui_state)
    visible_saved <- pending$humann_visible_samples
    hidden_saved  <- pending$humann_hidden_samples
    if (is.null(visible_saved) && is.null(hidden_saved)) {
      visible_saved <- isolate(input$humann_visible_samples)
      hidden_saved  <- isolate(input$humann_hidden_samples)
    }
    known <- c(visible_saved, hidden_saved)
    if (length(intersect(known, samples)) > 0) {
      visible <- intersect(visible_saved, samples)
      hidden  <- intersect(hidden_saved,  samples)
      visible <- c(visible, setdiff(samples, c(visible, hidden)))
    } else {
      visible <- samples
      hidden  <- character(0)
    }
    if (HAS_SORTABLE) {
      sortable::bucket_list(
        header = "Drag to reorder / hide samples — order applies to plot bars:",
        group_name = "humann_sample_buckets",
        orientation = "horizontal",
        sortable::add_rank_list(
          text = "Visible (order applies to plots)",
          labels = visible, input_id = "humann_visible_samples",
          options = .multidrag_opts
        ),
        sortable::add_rank_list(
          text = "Hidden",
          labels = hidden, input_id = "humann_hidden_samples",
          options = .multidrag_opts
        )
      )
    } else {
      selectizeInput("humann_visible_samples",
                     label = "Include samples (leave empty for all):",
                     choices = samples, selected = visible,
                     multiple = TRUE, width = "100%",
                     options = list(plugins = list("remove_button")))
    }
  })

  # Group-by dropdown — splits the stacked bar into faceted panels, one
  # per level of the chosen metadata column. "None" = single-panel view.
  output$humann_group_by_ui <- renderUI({
    req(rv$app_mode == "humann")
    meta <- rv$metadata
    if (is.null(meta) || !("Sample" %in% colnames(meta))) return(NULL)
    cols <- setdiff(colnames(meta), "Sample")
    if (length(cols) == 0) return(NULL)
    selectInput("humann_group_by",
                label = "Group by metadata (splits into panels):",
                choices = c("None", cols),
                selected = restored_input("humann_group_by", "None",
                                          valid = c("None", cols)))
  })

  output$humann_meta_filter_ui <- renderUI({
    req(rv$app_mode == "humann")
    meta <- rv$metadata
    if (is.null(meta) || !("Sample" %in% colnames(meta))) {
      return(helpText(tags$em("No metadata columns found in the dataset.")))
    }
    cols <- setdiff(colnames(meta), "Sample")
    if (length(cols) == 0) {
      return(helpText(tags$em("No metadata columns found in the dataset.")))
    }
    tagList(
      tags$b("Filter by metadata (leave empty for all):"),
      lapply(cols, function(col) {
        vals <- unique(as.character(meta[[col]]))
        vals <- vals[!is.na(vals) & nzchar(vals)]
        if (length(vals) == 0) return(NULL)
        selectizeInput(paste0("humann_meta_", col),
                       label = col, choices = vals,
                       selected = restored_input(paste0("humann_meta_", col),
                                                 character(0), valid = vals),
                       multiple = TRUE,
                       options = list(placeholder = "(all)"))
      })
    )
  })

  humann_rename_df <- reactive({
    samples <- humann_all_samples()
    renames <- rv$sample_renames
    disp <- if (is.null(renames)) samples else {
      raw <- unname(renames[samples])
      ifelse(!is.na(raw) & nzchar(raw), raw, samples)
    }
    data.frame(Original = samples, Display = disp,
               stringsAsFactors = FALSE)
  })

  output$humann_rename_dt <- renderDT({
    datatable(
      humann_rename_df(),
      editable = list(target = "cell", disable = list(columns = 0)),
      rownames = FALSE,
      filter   = "none",
      options  = list(pageLength = 25, dom = "tp", scrollX = TRUE,
                      scrollY = "240px", scroller = TRUE)
    )
  }, server = FALSE)

  observeEvent(input$humann_rename_dt_cell_edit, {
    info <- input$humann_rename_dt_cell_edit
    df   <- isolate(humann_rename_df())
    if (info$row < 1 || info$row > nrow(df)) return()
    orig <- df$Original[info$row]
    new  <- trimws(as.character(info$value))
    if (!nzchar(new)) new <- orig
    cur <- rv$sample_renames
    if (is.null(cur)) cur <- setNames(df$Original, df$Original)
    cur[orig] <- new
    rv$sample_renames <- cur
  })

  # Effective sample list for HUMAnN plots. Honors the drag-bucket "Visible"
  # order first, then intersects with metadata filters (preserving order).
  humann_effective_samples <- reactive({
    all_samples <- humann_all_samples()
    keep <- input$humann_visible_samples
    if (is.null(keep) || length(keep) == 0) keep <- all_samples
    keep <- keep[keep %in% all_samples]        # preserve bucket order
    meta <- rv$metadata
    if (!is.null(meta) && "Sample" %in% colnames(meta)) {
      cols <- setdiff(colnames(meta), "Sample")
      for (col in cols) {
        sel <- input[[paste0("humann_meta_", col)]]
        if (!is.null(sel) && length(sel) > 0) {
          matched <- meta$Sample[as.character(meta[[col]]) %in% sel]
          keep <- keep[keep %in% matched]      # preserve order
        }
      }
    }
    keep
  })

  # Map original sample names → user-chosen display labels (falls back to
  # the original when no rename is set).
  humann_display <- function(orig) {
    m <- rv$sample_renames
    if (is.null(m)) return(orig)
    out <- ifelse(orig %in% names(m) & nzchar(m[orig]), m[orig], orig)
    unname(out)
  }

  # Look up the group-by metadata value for each original sample name.
  # Returns NULL if no grouping is active; otherwise a character vector
  # aligned to `orig` (missing → "(unknown)").
  humann_group_for <- function(orig) {
    col <- input$humann_group_by %||% "None"
    if (identical(col, "None")) return(NULL)
    meta <- rv$metadata
    if (is.null(meta) || !("Sample" %in% colnames(meta)) ||
        !(col %in% colnames(meta))) return(NULL)
    g <- setNames(as.character(meta[[col]]), meta$Sample)[orig]
    ifelse(is.na(g) | !nzchar(g), "(unknown)", g)
  }

  # Render one functional sub-tab: top-N stacked barplot + DT (full table,
  # including stratified rows when present). Plot drops specials + stratified
  # rows; the table keeps everything.
  render_humann_subtab <- function(slot_name, label, color_pal = "Set3",
                                   page_length = 10) {
    function() {
      df <- rv$humann[[slot_name]]
      if (is.null(df) || nrow(df) == 0) {
        return(div(class = "status-text status-warning",
                   icon("exclamation-triangle"),
                   sprintf(" No %s data found in this dataset.", label)))
      }
      tagList(
        h4(sprintf("Top features — %s (community-level, CPM)", label)),
        helpText(sprintf(
          "Stacked bar shows the 20 most abundant %s across samples. ",
          tolower(label)),
          "Stratified per-organism rows and HUMAnN special features ",
          "(UNMAPPED / READS_UNMAPPED / UNINTEGRATED / UNGROUPED) are ",
          "excluded from this view but kept in the table below."),
        div(style = "width: 100%;",
            plotlyOutput(paste0("humann_", slot_name, "_bar"),
                         width = "100%", height = "560px")),
        br(),
        div(
          downloadButton(paste0("download_humann_", slot_name, "_pdf"), "PDF"),
          downloadButton(paste0("download_humann_", slot_name, "_png"), "PNG"),
          downloadButton(paste0("download_humann_", slot_name, "_svg"), "SVG")
        ),
        br(),
        h4("Full table"),
        # R-side CSV/TSV downloads bypass DataTables Buttons — the JS export
        # was corrupting cells that contain a `"` (e.g. KEGG K00984's name),
        # and its "modifier.page = all" was being ignored on this DT build.
        div(
          downloadButton(paste0("download_humann_", slot_name, "_full_csv"),
                         "Download full table (CSV)"),
          downloadButton(paste0("download_humann_", slot_name, "_full_tsv"),
                         "Download full table (TSV)")
        ),
        br(),
        DTOutput(paste0("humann_", slot_name, "_dt"))
      )
    }
  }

  output$humann_pathways_ui      <- renderUI(render_humann_subtab(
    "pathways", "Pathways")())
  output$humann_reactions_ui     <- renderUI(render_humann_subtab(
    "reactions", "Reactions")())
  output$humann_kegg_kos_ui      <- renderUI(render_humann_subtab(
    "kegg_kos", "KEGG KOs")())
  output$humann_gene_families_ui <- renderUI(render_humann_subtab(
    "gene_families", "Gene families")())

  # Shared per-slot ggplot builder. Feeds both the on-screen plotly render
  # and the PDF/PNG/SVG download handlers so they stay in sync.
  humann_slot_gg <- function(slot_name) {
    df <- rv$humann[[slot_name]]
    req(df, nrow(df) > 0)
    keep_samples <- humann_effective_samples()
    shiny::validate(shiny::need(length(keep_samples) > 0,
      "No samples selected — expand 'Filters & renaming' above and pick at least one sample."))
    all_cols <- setdiff(colnames(df), "Feature")
    # Preserve the drag-bucket order (keep_samples already ordered).
    kept_cols <- keep_samples[keep_samples %in% all_cols]
    shiny::validate(shiny::need(length(kept_cols) > 0,
      sprintf("None of the selected samples appear in the %s table.", slot_name)))
    df <- df[, c("Feature", kept_cols), drop = FALSE]
    d <- humann_top_n(humann_drop_specials(humann_community_only(df)), n = 20)
    shiny::validate(shiny::need(!is.null(d) && nrow(d) > 0,
                                sprintf("No %s features to plot.", slot_name)))
    sample_cols <- setdiff(colnames(d), "Feature")
    # Reapply the intended order after top_n (which preserves it, but be
    # explicit) so the x-axis honours the Visible bucket.
    sample_cols <- kept_cols[kept_cols %in% sample_cols]
    d <- d[, c("Feature", sample_cols), drop = FALSE]
    display_samples <- humann_display(sample_cols)
    group_lut       <- humann_group_for(sample_cols)
    long <- data.frame(
      Feature = rep(d$Feature, times = length(sample_cols)),
      Sample  = rep(display_samples, each = nrow(d)),
      Value   = as.numeric(unlist(d[, sample_cols, drop = FALSE])),
      stringsAsFactors = FALSE
    )
    if (!is.null(group_lut)) {
      long$Group <- rep(group_lut, each = nrow(d))
    }
    # Keep Top-N order on the y-axis (highest mean at top) + drag order on x.
    long$Feature <- factor(long$Feature, levels = rev(d$Feature))
    long$Sample  <- factor(long$Sample, levels = display_samples)
    p <- ggplot(long, aes(x = Sample, y = Value, fill = Feature)) +
      geom_col(position = "stack") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            legend.text = element_text(size = 8)) +
      labs(x = NULL, y = "CPM (community-level)", fill = NULL)
    if (!is.null(group_lut)) {
      p <- p + facet_grid(cols = vars(Group), scales = "free_x",
                          space = "free_x")
    }
    p
  }

  # Generic per-slot renderers (plotly + DT). One function generates both
  # outputs so we don't repeat the filtering logic four times.
  attach_humann_outputs <- function(slot_name) {
    output[[paste0("humann_", slot_name, "_bar")]] <- renderPlotly({
      p <- humann_slot_gg(slot_name)
      ggplotly(p, tooltip = c("x", "y", "fill")) %>%
        layout(
          autosize = TRUE,
          legend   = list(font = list(size = 9), tracegroupgap = 2),
          margin   = list(l = 60, r = 10, t = 30, b = 110)
        ) %>%
        plotly::config(responsive = TRUE, displaylogo = FALSE)
    })

    local({
      s <- slot_name
      output[[paste0("download_humann_", s, "_pdf")]] <- downloadHandler(
        filename = function() paste0("exploreMetaTax_humann_", s, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
        content = function(file) save_gg(file, humann_slot_gg(s), "pdf")
      )
      output[[paste0("download_humann_", s, "_png")]] <- downloadHandler(
        filename = function() paste0("exploreMetaTax_humann_", s, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
        content = function(file) save_gg(file, humann_slot_gg(s), "png")
      )
      output[[paste0("download_humann_", s, "_svg")]] <- downloadHandler(
        filename = function() paste0("exploreMetaTax_humann_", s, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
        content = function(file) save_gg(file, humann_slot_gg(s), "svg")
      )
    })

    # R-side CSV / TSV downloads for the full slot table (bypasses the DT
    # Buttons pipeline entirely — some deployed builds mangle CSV cells on
    # rows whose plain data contains a `"`).
    local({
      s <- slot_name
      output[[paste0("download_humann_", s, "_full_csv")]] <- downloadHandler(
        filename = function() paste0("exploreMetaTax_humann_", s,
                                     "_", format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     ".csv"),
        content  = function(file) {
          d <- humann_slot_enrich(rv$humann[[s]], s)
          if (is.null(d)) d <- data.frame()
          write.csv(d, file, row.names = FALSE)
        }
      )
      output[[paste0("download_humann_", s, "_full_tsv")]] <- downloadHandler(
        filename = function() paste0("exploreMetaTax_humann_", s,
                                     "_", format(Sys.time(), "%Y%m%d_%H%M%S"),
                                     ".tsv"),
        content  = function(file) {
          d <- humann_slot_enrich(rv$humann[[s]], s)
          if (is.null(d)) d <- data.frame()
          write.table(d, file, sep = "\t", row.names = FALSE, quote = FALSE)
        }
      )
    })

    output[[paste0("humann_", slot_name, "_dt")]] <- renderDT({
      df <- humann_slot_enrich(rv$humann[[slot_name]], slot_name)
      req(df, nrow(df) > 0)
      # createdCell (below) injects the link into the DOM only — the
      # DataTables data store keeps the raw id, so anything reading it
      # (filter, sort, on-screen text) stays clean.
      spec <- humann_link_specs[[slot_name]]
      non_num  <- c("Feature", "Description",
                    vapply(humann_link_specs, `[[`, character(1), "col"))
      num_cols <- setdiff(colnames(df), non_num)
      opts <- list(
        pageLength = 10,
        lengthMenu = c(10, 25, 50, 100),
        scrollX    = TRUE,
        dom        = "Blfrtip",
        buttons    = list(
          list(extend = "copyHtml5",  exportOptions = dt_export_options),
          list(extend = "csvHtml5",   exportOptions = dt_export_options),
          list(extend = "excelHtml5", exportOptions = dt_export_options)
        )
      )
      if (!is.null(spec)) {
        opts$columnDefs <- list(list(
          targets = which(colnames(df) == spec$col) - 1L,
          createdCell = DT::JS(
            "function(td, cellData){",
            paste0("  var blank = '", spec$blank, "';"),
            "  if (!cellData || cellData === blank) return;",
            "  var a = document.createElement('a');",
            paste0("  a.href = '", spec$url, "' + encodeURIComponent(cellData);"),
            "  a.target = '_blank';",
            "  a.rel = 'noopener';",
            "  a.textContent = cellData;",
            "  td.textContent = '';",
            "  td.appendChild(a);",
            "}"
          )
        ))
      }
      dt <- datatable(
        df,
        rownames   = FALSE,
        filter     = "top",
        extensions = "Buttons",
        options    = opts
      )
      if (length(num_cols) > 0) dt <- formatRound(dt, columns = num_cols, digits = 2)
      dt
    })
  }
  for (s in c("pathways", "reactions", "kegg_kos", "gene_families")) {
    attach_humann_outputs(s)
  }

  # ── Lazy + async loader for the large gene-families table ────────────────
  # gene_families (~1.9M rows x N samples) is NOT read at initial dataset load
  # (load_humann_tables is called there without it). It is read the first time
  # the user opens a tab that needs it (Gene families or Stratified), IN A
  # BACKGROUND R PROCESS via a Shiny ExtendedTask, so the main session stays
  # responsive — the user can keep switching tabs while it reads. The
  # Gene-families/Stratified outputs stay empty (their req() halts) until the
  # data arrives, then fill in. A persistent bottom-right notification signals
  # progress instead of the blocking full-screen overlay.
  gf_invoked <- reactiveVal(FALSE)   # guards against a double-invoke on fast clicks

  # The ExtendedTask body runs in the main process but returns immediately with
  # a promise; the actual read+join happens in a callr child process. Helper
  # functions are passed in as args because that child is a fresh R session.
  gf_task <- ExtendedTask$new(
    function(paths, sample_names, read_fn, join_fn, quote, colClasses, nThread) {
      callr_promise(
        func = gf_read_join_bg,
        args = list(paths = paths, sample_names = sample_names,
                    read_fn = read_fn, join_fn = join_fn,
                    quote = quote, colClasses = colClasses, nThread = nThread)
      )
    }
  )

  observeEvent(input$humann_subtab, {
    if (!input$humann_subtab %in% c("Gene families", "Stratified by organisms"))
      return()
    if (isTRUE(gf_invoked())) return()                     # already loading/loaded
    if (!is.null(rv$humann[["gene_families"]])) return()   # already present
    pd <- rv$humann_parsed_ds
    if (is.null(pd)) return()                              # upload mode: already in
    col <- pd$humann_cols[["gene_families"]]
    if (is.null(col) || is.na(col)) return()              # no gene-family column
    ds           <- pd$raw
    sample_names <- as.character(ds[[pd$name_col]])
    rels         <- as.character(ds[[col]])
    paths        <- vapply(rels, resolve_gstore_file,
                           FUN.VALUE = NA_character_, USE.NAMES = FALSE)
    gf_invoked(TRUE)
    showNotification(
      tagList(tags$b("Loading gene families…"), tags$br(),
              "Reading in the background — you can keep using the other tabs."),
      id = "gf_loading", duration = NULL, closeButton = FALSE, type = "message")
    gf_task$invoke(paths, sample_names, read_humann_tsv, join_humann_wide,
                   "", c("character", "numeric"), HUMANN_GF_READ_THREADS)
  }, ignoreInit = TRUE)

  observeEvent(gf_task$status(), {
    st <- gf_task$status()
    if (identical(st, "success")) {
      removeNotification("gf_loading")
      gf_df <- gf_task$result()
      if (!is.null(gf_df)) {
        h <- isolate(rv$humann)
        h[["gene_families"]] <- gf_df
        rv$humann <- h
        showNotification(
          sprintf("Gene families loaded: %d features x %d sample(s)",
                  nrow(gf_df), ncol(gf_df) - 1L),
          type = "message", duration = 4)
      } else {
        showNotification("Gene families: no data could be loaded.",
                         type = "warning", duration = 6)
      }
    } else if (identical(st, "error")) {
      removeNotification("gf_loading")
      gf_invoked(FALSE)   # allow a retry on the next Gene-families/Stratified open
      msg <- tryCatch(gf_task$result(), error = function(e) conditionMessage(e))
      showNotification(paste0("Gene families failed to load: ", msg),
                       type = "error", duration = 8)
    }
  })

  # ─── Stratified per-species sub-tab ──────────────────────────────────
  # Users search for one or more species (partial substring, case-insensitive
  # on the s__Genus_species tag after "|") and get a stacked barplot +
  # DT of the matching stratified rows. Mirrors the "Organisms of interest"
  # search flow in the taxonomy mode.

  # Human-readable species label from a HUMAnN stratified suffix.
  humann_species_label <- function(x) {
    x <- sub("^[a-z]__", "", x, perl = TRUE)
    gsub("_", " ", x)
  }

  # Normalise a query / species string for substring matching.
  humann_species_norm <- function(x) {
    x <- humann_species_label(x)
    x <- tolower(x)
    x <- gsub("[^a-z0-9 ]", " ", x, perl = TRUE)
    trimws(gsub("\\s+", " ", x, perl = TRUE))
  }

  # Table-slot choices restricted to those with stratified rows.
  humann_strat_slots <- reactive({
    req(rv$app_mode == "humann", rv$humann)
    labels <- c(pathways = "Pathways", reactions = "Reactions",
                kegg_kos = "KEGG KOs", gene_families = "Gene families")
    keep <- vapply(names(rv$humann), function(s) {
      df <- rv$humann[[s]]
      !is.null(df) && nrow(df) > 0 && any(grepl("|", df$Feature, fixed = TRUE))
    }, logical(1))
    slots <- names(rv$humann)[keep]
    setNames(slots, labels[slots])
  })

  output$humann_stratified_ui <- renderUI({
    req(rv$app_mode == "humann")
    slots <- humann_strat_slots()
    if (length(slots) == 0) {
      return(div(class = "status-text status-warning",
                 icon("info-circle"),
                 " This HUMAnN run has no stratified per-organism rows",
                 " since it was run only in translated search mode with proteins.",
                 " A full run (tier-1 nucleotide + tier-2 translated) would be",
                 " needed to get insights beyond the community level."))
    }
    tagList(
      h4("Stratified per organism"),
      helpText(
        "Search for one or more organisms; matches are substring, ",
        "case-insensitive on the taxon tag after the ", tags$code("|"),
        " separator (typically ", tags$code("s__Genus_species"),
        " but any rank HUMAnN reported — e.g. ", tags$code("g__Bacteroides"),
        " or ", tags$code("unclassified"), " — is honoured). ",
        "The barplot stacks matched organisms per sample; the table lists ",
        "every matched stratified row."
      ),
      fluidRow(
        column(4,
          selectInput("humann_strat_slot", "Table:",
                      choices = slots,
                      selected = restored_input("humann_strat_slot", slots[[1]],
                                                valid = slots))
        ),
        column(8,
          textAreaInput("humann_strat_query",
                        "Organisms (one per line or comma-separated; leave empty to show all):",
                        value = "", width = "100%", height = "90px",
                        placeholder = "e.g.\nBacteroides fragilis\nEscherichia coli\nBacteroides"),
          actionButton("humann_strat_search_btn", "Search", icon = icon("search"),
                       class = "btn-primary"),
          actionButton("humann_strat_clear_btn", "Clear")
        )
      ),
      uiOutput("humann_strat_summary"),
      div(style = "width: 100%;",
          plotlyOutput("humann_strat_bar", width = "100%", height = "560px")),
      br(),
      div(
        downloadButton("download_humann_strat_pdf", "PDF"),
        downloadButton("download_humann_strat_png", "PNG"),
        downloadButton("download_humann_strat_svg", "SVG")
      ),
      br(),
      h4("Matched stratified rows"),
      DTOutput("humann_strat_dt"),
      br(),
      downloadButton("humann_strat_download_csv",
                     "Download matched rows (CSV)")
    )
  })

  observeEvent(input$humann_strat_clear_btn, {
    updateTextAreaInput(session, "humann_strat_query", value = "")
    rv$humann_strat_query_committed <- character(0)
  })

  # Commit queries only on Search click (mirrors organism_search_btn).
  rv$humann_strat_query_committed <- character(0)
  observeEvent(input$humann_strat_search_btn, {
    txt <- input$humann_strat_query %||% ""
    qs <- unlist(strsplit(txt, "[,\n]", perl = TRUE), use.names = FALSE)
    qs <- trimws(qs)
    rv$humann_strat_query_committed <- qs[nzchar(qs)]
  })

  # Stratified rows for the selected table, decomposed into feature + species.
  humann_strat_rows <- reactive({
    req(rv$app_mode == "humann", rv$humann, input$humann_strat_slot)
    df <- rv$humann[[input$humann_strat_slot]]
    req(!is.null(df), nrow(df) > 0)
    strat <- df[grepl("|", df$Feature, fixed = TRUE), , drop = FALSE]
    if (nrow(strat) == 0) return(strat)
    # Split "feature|taxon" with two vectorised sub() calls instead of
    # strsplit() + two vapply loops over ~1M rows (HUMAnN uses exactly one "|").
    strat$BaseFeature  <- sub("\\|.*$", "", strat$Feature)
    strat$Species      <- sub("^[^|]*\\|", "", strat$Feature)
    strat$SpeciesLabel <- humann_species_label(strat$Species)
    strat
  })

  # Rows matching the committed query set (or all stratified rows if empty).
  humann_strat_matched <- reactive({
    strat <- humann_strat_rows()
    req(nrow(strat) > 0)
    qs <- rv$humann_strat_query_committed
    if (length(qs) == 0) return(strat)
    species_norm <- humann_species_norm(strat$Species)
    qn <- humann_species_norm(qs)
    qn <- qn[nzchar(qn)]
    if (length(qn) == 0) return(strat)
    hit <- rep(FALSE, length(species_norm))
    for (q in qn) hit <- hit | grepl(q, species_norm, fixed = TRUE)
    strat[hit, , drop = FALSE]
  })

  # Capped view for the plot/table. With no query (or a very broad one) the
  # matched set can be ~1M rows x N samples, which would blow up the long-form
  # reshape in humann_strat_gg and the DT. Keep the top HUMANN_STRAT_ROW_CAP
  # rows by summed CPM and flag the cap so the summary can warn. A specific
  # organism query normally narrows the set well below the cap.
  humann_strat_view <- reactive({
    matched <- humann_strat_matched()
    total   <- nrow(matched)
    if (total <= HUMANN_STRAT_ROW_CAP)
      return(list(df = matched, capped = FALSE, total = total))
    sample_cols <- setdiff(colnames(matched),
                           c("Feature", "BaseFeature", "Species", "SpeciesLabel"))
    tot  <- rowSums(matched[, sample_cols, drop = FALSE], na.rm = TRUE)
    keep <- order(tot, decreasing = TRUE)[seq_len(HUMANN_STRAT_ROW_CAP)]
    list(df = matched[keep, , drop = FALSE], capped = TRUE, total = total)
  })

  output$humann_strat_summary <- renderUI({
    req(rv$app_mode == "humann")
    strat <- humann_strat_rows()
    if (nrow(strat) == 0)
      return(tags$p(tags$em("No stratified rows in this table.")))
    view    <- humann_strat_view()
    matched <- view$df
    qs      <- rv$humann_strat_query_committed
    n_sp    <- length(unique(matched$Species))
    cap_note <- if (isTRUE(view$capped))
      tags$p(class = "status-text status-warning",
             icon("exclamation-triangle"),
             sprintf(paste0(" Showing the top %s stratified rows by summed CPM ",
                            "(of %s matched). Enter an organism query and click ",
                            "Search to view a specific subset in full; the CSV ",
                            "download always contains every matched row."),
                     format(HUMANN_STRAT_ROW_CAP, big.mark = ","),
                     format(view$total, big.mark = ","))) else NULL
    body <- if (length(qs) == 0) {
      tags$p(sprintf(
        "Showing %d stratified row%s across %d organism%s (no query — enter one and click Search to filter).",
        nrow(matched), if (nrow(matched) == 1) "" else "s",
        n_sp, if (n_sp == 1) "" else "s"))
    } else {
      tags$p(sprintf(
        "%d row%s shown across %d organism%s for %d quer%s: %s",
        nrow(matched), if (nrow(matched) == 1) "" else "s",
        n_sp, if (n_sp == 1) "" else "s",
        length(qs), if (length(qs) == 1) "y" else "ies",
        paste(shQuote(qs, type = "sh"), collapse = ", ")))
    }
    tagList(cap_note, body)
  })

  # Shared ggplot builder for the stratified per-organism barplot. Used by
  # the on-screen renderPlotly and the PDF/PNG/SVG download handlers.
  humann_strat_gg <- reactive({
    matched <- humann_strat_view()$df
    shiny::validate(shiny::need(nrow(matched) > 0,
                                "No stratified rows matched — try a shorter substring, or clear the query to show all."))
    keep_samples <- humann_effective_samples()
    shiny::validate(shiny::need(length(keep_samples) > 0,
      "No samples selected — expand 'Filters & renaming' above and pick at least one sample."))

    all_cols <- setdiff(colnames(matched),
                        c("Feature", "BaseFeature", "Species", "SpeciesLabel"))
    # Preserve drag-bucket order.
    sample_cols <- keep_samples[keep_samples %in% all_cols]
    shiny::validate(shiny::need(length(sample_cols) > 0,
      "None of the selected samples appear in the stratified table."))
    group_lut <- humann_group_for(sample_cols)
    long <- do.call(rbind, lapply(seq_along(sample_cols), function(i) {
      s <- sample_cols[i]
      data.frame(Sample = humann_display(s),
                 Organism = matched$SpeciesLabel,
                 Value = as.numeric(matched[[s]]),
                 Group = if (!is.null(group_lut)) group_lut[i] else NA_character_,
                 stringsAsFactors = FALSE)
    }))
    long$Value[is.na(long$Value)] <- 0
    # Sample → group lookup (unique display sample → group value); reattached
    # after aggregate() drops non-grouping-by columns.
    sample_group <- if (!is.null(group_lut)) {
      setNames(long$Group[!duplicated(long$Sample)],
               long$Sample[!duplicated(long$Sample)])
    } else NULL
    display_samples <- humann_display(sample_cols)

    agg <- aggregate(Value ~ Sample + Organism, data = long, FUN = sum)

    sp_means <- aggregate(Value ~ Organism, data = agg, FUN = mean)
    sp_means <- sp_means[order(sp_means$Value, decreasing = TRUE), ]
    top_sp <- head(sp_means$Organism, 20)
    agg$Organism <- ifelse(agg$Organism %in% top_sp, agg$Organism, "Other")
    agg <- aggregate(Value ~ Sample + Organism, data = agg, FUN = sum)
    ordering <- aggregate(Value ~ Organism, data = agg, FUN = sum)
    ordering <- ordering[order(ordering$Value), , drop = FALSE]
    lvls <- as.character(ordering$Organism)
    if ("Other" %in% lvls) lvls <- c("Other", setdiff(lvls, "Other"))
    agg$Organism <- factor(agg$Organism, levels = lvls)
    # Honour the drag-bucket order on the x-axis.
    agg$Sample <- factor(agg$Sample, levels = display_samples)
    if (!is.null(sample_group)) {
      agg$Group <- unname(sample_group[as.character(agg$Sample)])
    }

    p <- ggplot(agg, aes(x = Sample, y = Value, fill = Organism,
                         text = paste0("Sample: ", Sample,
                                       "\nOrganism: ", Organism,
                                       "\nCPM (summed): ", sprintf("%.2f", Value)))) +
      geom_col(position = "stack") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            legend.text = element_text(size = 9)) +
      labs(x = NULL, y = "Summed CPM (matched stratified rows)", fill = "Organism")
    if (!is.null(sample_group)) {
      p <- p + facet_grid(cols = vars(Group), scales = "free_x",
                          space = "free_x")
    }
    p
  })

  output$humann_strat_bar <- renderPlotly({
    ggplotly(humann_strat_gg(), tooltip = "text") %>%
      layout(
        autosize = TRUE,
        legend   = list(font = list(size = 9), tracegroupgap = 2),
        margin   = list(l = 60, r = 10, t = 30, b = 110)
      ) %>%
      plotly::config(responsive = TRUE, displaylogo = FALSE)
  })

  output$download_humann_strat_pdf <- downloadHandler(
    filename = function() paste0("exploreMetaTax_humann_stratified_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) save_gg(file, humann_strat_gg(), "pdf")
  )
  output$download_humann_strat_png <- downloadHandler(
    filename = function() paste0("exploreMetaTax_humann_stratified_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) save_gg(file, humann_strat_gg(), "png")
  )
  output$download_humann_strat_svg <- downloadHandler(
    filename = function() paste0("exploreMetaTax_humann_stratified_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content = function(file) save_gg(file, humann_strat_gg(), "svg")
  )

  output$humann_strat_dt <- renderDT({
    matched <- humann_strat_view()$df
    shiny::validate(shiny::need(nrow(matched) > 0,
                                "No stratified rows matched. Adjust the query and click Search."))
    sample_cols <- setdiff(colnames(matched),
                           c("Feature", "BaseFeature", "Species", "SpeciesLabel"))
    show <- matched[, c("BaseFeature", "SpeciesLabel", sample_cols), drop = FALSE]
    colnames(show)[1:2] <- c("Feature", "Organism")
    dt <- datatable(
      show, rownames = FALSE, filter = "top",
      extensions = "Buttons",
      options = list(pageLength = 15, lengthMenu = c(10, 15, 25, 50, 100),
                     scrollX = TRUE, dom = "Blfrtip",
                     buttons = list(
                       list(extend = "copyHtml5",  exportOptions = dt_export_options),
                       list(extend = "csvHtml5",   exportOptions = dt_export_options),
                       list(extend = "excelHtml5", exportOptions = dt_export_options)
                     ))
    )
    if (length(sample_cols) > 0) dt <- formatRound(dt, columns = sample_cols, digits = 2)
    dt
  })

  output$humann_strat_download_csv <- downloadHandler(
    filename = function()
      paste0("humann_stratified_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      matched <- humann_strat_matched()
      write.csv(matched, file, row.names = FALSE)
    }
  )

  humann_summary_df <- reactive({
    req(rv$app_mode == "humann", rv$humann)
    samples <- names(rv$humann_modes %||% character(0))
    if (length(samples) == 0) {
      any_df <- rv$humann[[which(!vapply(rv$humann, is.null, logical(1)))[1]]]
      samples <- setdiff(colnames(any_df), "Feature")
    }
    modes_disp <- if (!is.null(rv$humann_modes)) {
      ifelse(is.na(rv$humann_modes) | !nzchar(rv$humann_modes),
             "unknown", rv$humann_modes)
    } else {
      rep("unknown", length(samples))
    }
    summary_df <- data.frame(
      Sample = samples,
      Mode   = modes_disp[seq_along(samples)],
      stringsAsFactors = FALSE
    )
    for (slot in names(rv$humann)) {
      df <- rv$humann[[slot]]
      if (is.null(df)) next
      community_n  <- sum(!grepl("|", df$Feature, fixed = TRUE))
      stratified_n <- nrow(df) - community_n
      label <- switch(slot,
                      pathways = "Pathways",
                      reactions = "Reactions",
                      kegg_kos = "KEGG KOs",
                      gene_families = "Gene families",
                      slot)
      summary_df[[label]] <- sprintf("%d (%d strat.)",
                                     community_n, stratified_n)
    }
    summary_df
  })

  output$humann_summary_dt <- renderDT({
    datatable(humann_summary_df(), rownames = FALSE,
              options = list(dom = "t", pageLength = 50, scrollX = TRUE))
  })

  output$humann_info_ui <- renderUI({
    if (!isTRUE(rv$app_mode == "humann")) {
      return(div("No HUMAnN dataset loaded."))
    }
    tagList(
      h4("HUMAnN run summary"),
      helpText(sprintf("Source folder: %s",
                       rv$humann_source_dir %||% "(unknown)")),
      tags$ul(
        tags$li(tags$b("full"), " — HUMAnN ran with both tier-1 (CHOCOPhlAn ",
                "nucleotide search) and tier-2 (UniRef translated search). ",
                "Output includes per-SGB stratified rows."),
        tags$li(tags$b("translated"), " — HUMAnN ran with translated search ",
                "only (tier-2). Output is community-level; no stratification.")
      ),
      DTOutput("humann_summary_dt")
    )
  })

  output$humann_method_ui <- renderUI({
    tagList(
      h4("Method notes"),
      tags$ul(
        tags$li("Gene family / reaction abundances are renormalized to CPM ",
                "(humann_renorm_table -u cpm) after stripping the ",
                code("READS_UNMAPPED"), " row (HUMAnN 4 alpha renamed ",
                code("UNMAPPED"), " but humann_renorm_table's --special ",
                "filter does not yet recognise it)."),
        tags$li("Pathway abundances are the raw ", code("_4_pathabundance.tsv"),
                " from HUMAnN (already in copies-per-million-derived units)."),
        tags$li("KEGG KO is regrouped from the CPM gene families via ",
                code("humann_regroup_table -g uniref90_ko"),
                ". HUMAnN 4 alpha's utility_mapping only covers ",
                "EC-annotated KOs, so coverage is partial."),
        tags$li("Stratified rows (", code("feature|species"), ") are kept in ",
                "the full table but excluded from top-N / heatmap views.")
      ),
      h4("References"),
      tags$ul(
        tags$li(tags$a(href = "https://github.com/biobakery/humann",
                       "HUMAnN GitHub repository")),
        tags$li(tags$a(href = "https://forum.biobakery.org/t/metaphlan-4-humann-4-compatibility/8523",
                       "HUMAnN 4 / MetaPhlAn 4 compatibility (bioBakery forum)"))
      )
    )
  })

  # ═══════════════════════════════════════════════════════════════════════
  # LOAD & EXPORT SETTINGS
  # ─────────────────────────────────────────────────────────────────────
  # Serialize the whole app state (rv$ data + input$ widget snapshot) to a
  # single RDS bundle the user downloads locally, and restore it from the
  # same file. Uploaded file paths, background processes and computed plot
  # objects are excluded — they can't survive a round-trip.
  # ═══════════════════════════════════════════════════════════════════════

  # Inputs we never restore: file uploads (stale datapath), action buttons
  # (would re-fire observers), DT / plotly UI-state echoes, dynamic
  # per-sample rename_* boxes (rebuilt from rv$sample_renames), and the
  # settings tab's own widgets. Regex is matched against each input name.
  ## NB: build with paste0 and keep every "|" INSIDE a fragment. An earlier
  ## paste(sep="|") inserted a stray "|" between the 2nd and 3rd fragments
  ## ("_row_last_clicked" + "|" + "|_state" -> "||"), and that empty alternative
  ## made the trailing (...)$ group match the empty string at the end of EVERY
  ## id — so grepl() skipped every input and the restore map was always empty
  ## (only the sortable buckets, which bypass this regex, restored).
  .settings_skip_pat <- paste0(
    "^(load_settings_file|data_files|metadata_file|rename_)",
    "|(_btn|_cell_edit|_cell_clicked|_rows_selected|_row_last_clicked",
    "|_state|_search|_hover|_click|_relayout|_brush|_selected)$")

  # sortable::bucket_list widgets have no update_*() function. Instead of
  # sendInputMessage (which is lost when the widget doesn't yet exist), we
  # stash their values in rv$pending_ui_state and let the sample_checkboxes
  # / humann_sample_filter_ui renderUIs pick them up on next render.
  .settings_sortable_ids <- c(
    "samples_visible", "samples_hidden",
    "humann_visible_samples", "humann_hidden_samples",
    "sample_filter"
  )

  # rv fields we persist. Skips lifemap_bg (OS process handle), lifemap_port
  # (bound to a running sub-process), lifemap_obj (bulky; rebuilt via the
  # "Run LifemapR" button), and the pending_* helpers themselves.
  .settings_rv_fields <- c(
    "raw_data", "metadata", "sample_renames", "loaded", "load_time",
    "app_mode", "humann", "humann_modes", "humann_source_dir",
    "discard_committed"
  )

  output$download_settings <- downloadHandler(
    filename = function() {
      sprintf("exploreMetaTax_settings_%s.rds",
              format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      rv_snap <- setNames(
        lapply(.settings_rv_fields, function(nm) rv[[nm]]),
        .settings_rv_fields
      )
      input_snap <- reactiveValuesToList(input, all.names = TRUE)
      bundle <- list(
        meta = list(
          app       = "exploreMetaTax",
          version   = 2L,
          saved_at  = Sys.time(),
          app_mode  = rv$app_mode,
          n_samples = if (!is.null(rv$raw_data))
                        length(unique(rv$raw_data$sample)) else 0L
        ),
        rv     = rv_snap,
        inputs = input_snap
      )
      saveRDS(bundle, file)
    }
  )

  output$settings_save_status <- renderText({
    if (isTRUE(rv$loaded) || isTRUE(rv$app_mode == "humann")) {
      sprintf("Ready to save — %d rv fields + %d input widgets.",
              length(.settings_rv_fields),
              length(reactiveValuesToList(input, all.names = TRUE)))
    } else {
      "No data loaded yet. You can still export current widget values."
    }
  })

  # Retry-until-applied loop for restored input values. sendInputMessage() is the
  # ONLY reliable way to overwrite an input that already exists on the client:
  # re-rendering a renderUI widget with the same inputId PRESERVES its current
  # value and ignores the new `selected=`, so a `selected = restored_input(...)`
  # only lands on a widget's FIRST render. Widgets that were already on screen at
  # load time (e.g. Composition's Group By on the default tab) therefore need
  # this message instead.
  #
  # Widgets on not-yet-opened tabs aren't in the DOM yet, so a message sent now
  # is lost. Rather than give up after a fixed window (and drop the queue), we
  # run a short burst, then go DORMANT while KEEPING the queue, and re-arm the
  # burst every time the user opens a tab (see the input$tabs observer below) —
  # that's when a tab's widgets mount. Each id is dropped once input$<id> matches
  # its target, so the queue drains as widgets are visited; it is cleared wholesale
  # only when a fresh dataset is loaded.
  .restore_ticks <- reactiveVal(0L)
  observe({
    pending <- rv$pending_inputs
    ticks <- .restore_ticks()
    if (length(pending) == 0) return()
    if (ticks >= 30L) return()   # burst done — stay dormant, KEEP the queue
    invalidateLater(200, session)
    still <- list()
    for (id in names(pending)) {
      v <- pending[[id]]
      if (is.null(v)) next
      cur <- isolate(input[[id]])
      if (!is.null(cur) && identical(cur, v)) next   # already applied
      tryCatch(session$sendInputMessage(id, list(value = v)),
               error = function(e) NULL)
      still[[id]] <- v
    }
    if (length(still) != length(pending)) rv$pending_inputs <- still
    .restore_ticks(ticks + 1L)
  })

  # Opening a tab mounts its dynamic widgets; re-arm the retry burst so any
  # still-unapplied restored values get pushed to them now.
  observeEvent(input$tabs, {
    if (length(rv$pending_inputs) > 0) .restore_ticks(0L)
  }, ignoreInit = TRUE)

  observeEvent(input$load_settings_file, {
    f <- input$load_settings_file
    req(f)
    bundle <- tryCatch(readRDS(f$datapath),
                       error = function(e) list(.error = conditionMessage(e)))
    if (!is.null(bundle$.error) ||
        !is.list(bundle) ||
        !identical(bundle$meta$app, "exploreMetaTax")) {
      msg <- paste0("✗ Not a valid exploreMetaTax settings file",
                    if (!is.null(bundle$.error))
                      paste0(" (", bundle$.error, ")") else "",
                    ".")
      output$load_settings_status <- renderText(msg)
      showNotification(msg, type = "error", duration = 8)
      return()
    }

    inputs_bundle <- if (is.list(bundle$inputs)) bundle$inputs else list()

    # Flat map of every restorable input value (skip file uploads, action
    # buttons, DT/plotly echoes, and the sortable ids handled separately). Feeds
    # BOTH restore paths: the sendInputMessage retry loop (rv$pending_inputs, for
    # static widgets) and the seed-on-render store (rv$restore_inputs, for
    # renderUI widgets). Built up-front and stored NOW so the store is populated
    # before the rv fields below are restored — restoring rv$metadata etc.
    # re-renders the dynamic selectors, which then read the store to seed their
    # saved value.
    restore <- list()
    for (id in names(inputs_bundle)) {
      if (grepl(.settings_skip_pat, id)) next
      if (id %in% .settings_sortable_ids)      next
      v <- inputs_bundle[[id]]
      if (is.null(v)) next
      restore[[id]] <- v
    }
    rv$restore_inputs <- restore

    ## Diagnostic (container stdout log): confirm the dropdown values actually
    ## made it into the saved file. If a selector is missing here, the problem is
    ## on the SAVE side, not the restore side.
    .dyn_ids <- c("group_by", "organism_group_by", "organism_distrib_group",
                  "heatmap_group_by", "violin_group", "diversity_color",
                  "rarefy_color", "pca_color", "krona_sample", "krona_group_by",
                  "humann_group_by")
    .present <- intersect(.dyn_ids, names(restore))
    message("[restore] ", length(restore), " input(s) queued; dynamic selectors saved: ",
            if (length(.present)) paste0(vapply(.present,
              function(k) paste0(k, "=", paste(restore[[k]], collapse = "|")),
              character(1)), collapse = ", ") else "(none)")

    # ── (1) Seed sortable widget state — renderUIs consult this on next
    # render (see output$sample_checkboxes and output$humann_sample_filter_ui)
    # so the initial DOM order matches the saved state on first paint.
    # We ALSO ship the order to the client-side reorder handler (below),
    # which corrects the DOM after the widget mounts — belt-and-suspenders
    # because add_rank_list's initial-labels-order doesn't always survive
    # sortable's client-side re-init after uiOutput swaps the fragment in.
    sortable_state <- list()
    for (id in .settings_sortable_ids) {
      v <- inputs_bundle[[id]]
      if (!is.null(v)) sortable_state[[id]] <- v
    }
    rv$pending_ui_state <- sortable_state

    send_bucket_pair <- function(vis_id, hid_id) {
      vis <- inputs_bundle[[vis_id]]
      hid <- inputs_bundle[[hid_id]]
      if (is.null(vis) && is.null(hid)) return(invisible())
      session$sendCustomMessage("restore_sortable_buckets", list(
        visible_id    = vis_id,
        hidden_id     = hid_id,
        visible_order = as.list(if (is.null(vis)) character(0) else vis),
        hidden_order  = as.list(if (is.null(hid)) character(0) else hid)
      ))
    }
    send_bucket_pair("samples_visible",        "samples_hidden")
    send_bucket_pair("humann_visible_samples", "humann_hidden_samples")

    # ── (2) Restore rv fields. Setting raw_data / metadata / sample_renames
    # invalidates every downstream renderUI, so dynamic widgets (group_by,
    # humann_meta_*, colour selectors, …) will be re-created — the retry
    # observer above then pushes their saved values as they materialise.
    rv_bundle <- if (is.list(bundle$rv)) bundle$rv else list()
    for (nm in intersect(.settings_rv_fields, names(rv_bundle))) {
      rv[[nm]] <- rv_bundle[[nm]]
    }

    # ── (3) Queue the same map for the sendInputMessage retry loop, which drives
    # the STATIC (always-present) widgets. Dynamic renderUI widgets are handled
    # by the rv$restore_inputs seed set above.
    rv$pending_inputs <- restore
    .restore_ticks(0L)

    # After the renderUIs have consumed the sortable snapshot on the next
    # flush, drop it — otherwise a later re-render (e.g. new data load) would
    # revert the user's post-restore drags to the restored order.
    session$onFlushed(function() { rv$pending_ui_state <- list() },
                      once = TRUE)

    saved_at <- if (inherits(bundle$meta$saved_at, "POSIXct"))
                  format(bundle$meta$saved_at, "%Y-%m-%d %H:%M:%S") else "unknown"
    msg <- sprintf(
      paste0("✓ Restored from '%s' (saved %s). ",
             "%d rv fields, %d sortable buckets, %d widgets queued."),
      f$name, saved_at,
      length(intersect(.settings_rv_fields, names(rv_bundle))),
      length(sortable_state),
      length(restore))
    output$load_settings_status <- renderText(msg)
    showNotification(msg, type = "message", duration = 8)
  })
}

# ═══════════════════════════════════════════════════════════════════════════
# RUN
# ═══════════════════════════════════════════════════════════════════════════
shinyApp(ui = ui, server = server)

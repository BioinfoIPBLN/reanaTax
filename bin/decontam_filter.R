#!/usr/bin/env Rscript
# decontam_filter.R -- identify reagent contaminants from negative controls.
#
# Every filter this pipeline has until now judges a taxon by its own evidence:
# how abundant it is, how much of its reference its reads cover, whether its
# counts scale across samples. None of them can answer the question a blank
# answers directly - is this taxon in my samples, or in my kit?
#
# decontam (Davis et al., Microbiome 2018) settles it statistically. Its
# prevalence method compares how often a taxon appears in true samples against
# how often it appears in negative controls; a taxon commoner in the blanks is
# scored as a contaminant. Its frequency method needs no blanks but does need a
# post-PCR DNA concentration per sample, and tests whether a taxon's relative
# abundance falls as total DNA rises - the signature of a constant amount of
# contaminating template diluted by varying amounts of real sample.
#
# This is optional, and stays optional, because most public datasets have
# neither blanks nor concentrations. A reanalysis pipeline that required them
# could not run on the archives it exists to mine. What it does instead is fail
# loudly when asked for a method whose inputs are absent, rather than quietly
# reporting that nothing was contaminated.
#
# Like the minimizer and host-k-mer filters, nothing is removed here: the
# taxids go to the abundance filter, so one step owns every removal.
#
# Arguments are parsed by hand, matching bin/differential_abundance.R: the
# biocontainers carry their Bioconductor stack and little else, and a missing
# CRAN package would fail the task for no good reason.

.defaults <- list(
  counts = NA, metadata = NA, neg_column = NA, neg_value = "true",
  conc_column = "", method = "prevalence", threshold = 0.1,
  batch_column = "", batch_combine = "minimum", prefix = "decontam"
)
.numeric <- c("threshold")

parseArgs <- function(argv, defaults) {
  opt <- defaults
  i <- 1
  while (i <= length(argv)) {
    key <- argv[i]
    if (!startsWith(key, "--")) stop("Unexpected argument '", key, "'.")
    name <- gsub("-", "_", sub("^--", "", key))
    if (!name %in% names(defaults)) stop("Unknown option '", key, "'.")
    if (i + 1 > length(argv)) stop("Option '", key, "' needs a value.")
    opt[[name]] <- if (name %in% .numeric) as.numeric(argv[i + 1]) else argv[i + 1]
    i <- i + 2
  }
  opt
}

opt <- parseArgs(commandArgs(trailingOnly = TRUE), .defaults)
for (required in c("counts", "metadata")) {
  if (is.na(opt[[required]])) stop("--", gsub("_", "-", required), " is required.")
}

needsNeg  <- opt$method %in% c("prevalence", "combined", "either", "minimum")
needsConc <- opt$method %in% c("frequency", "combined", "either", "minimum")

## ---- counts -------------------------------------------------------------
# The combined Bracken table, same contract as differential_abundance.R:
# `_num` columns are integer counts, `_frac` are relative abundances. decontam's
# frequency method models abundance against DNA concentration, so it needs the
# counts, not a profile that has already been closed to 1.
tab <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
numCols <- grep("_num$", colnames(tab), value = TRUE)
if (length(numCols) == 0) {
  stop("No `<sample>_num` columns in '", opt$counts, "'. This needs the combined ",
       "Bracken table (bracken_combined_<level>.txt).")
}
if (!"name" %in% colnames(tab)) stop("No `name` column in '", opt$counts, "'.")

counts <- as.matrix(tab[, numCols, drop = FALSE])
mode(counts) <- "numeric"
rownames(counts) <- make.unique(as.character(tab$name))
colnames(counts) <- sub("_num$", "", numCols)

## ---- metadata -----------------------------------------------------------
meta <- read.delim(opt$metadata, check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- as.character(meta[[1]])
shared <- intersect(colnames(counts), rownames(meta))
if (length(shared) == 0) {
  stop("No sample in '", opt$metadata, "' matches a column of '", opt$counts,
       "'. The first metadata column must hold the sample ids the pipeline used.")
}
missingFromMeta <- setdiff(colnames(counts), shared)
if (length(missingFromMeta)) {
  message("[decontam_filter] not in the metadata, dropped from the test: ",
          paste(missingFromMeta, collapse = ", "))
}
counts <- counts[, shared, drop = FALSE]
meta <- meta[shared, , drop = FALSE]

## ---- which samples are blanks -------------------------------------------
isNeg <- NULL
if (needsNeg) {
  if (is.na(opt$neg_column)) {
    stop("--neg-column is required for method '", opt$method, "'. It names the ",
         "metadata column marking negative controls.")
  }
  if (!opt$neg_column %in% colnames(meta)) {
    stop("Column '", opt$neg_column, "' is not in '", opt$metadata, "'. Columns present: ",
         paste(colnames(meta), collapse = ", "))
  }
  isNeg <- tolower(trimws(as.character(meta[[opt$neg_column]]))) == tolower(opt$neg_value)
  isNeg[is.na(isNeg)] <- FALSE
  if (sum(isNeg) == 0) {
    stop("No sample has '", opt$neg_column, "' == '", opt$neg_value, "'. The prevalence ",
         "method cannot run without blanks - either name the right column and value, or ",
         "use --decontam_method frequency with a DNA-concentration column, or turn ",
         "--decontam off. It is optional precisely because most public datasets have neither.")
  }
  if (sum(!isNeg) < 2) {
    stop("Only ", sum(!isNeg), " true sample(s) against ", sum(isNeg), " control(s). ",
         "The prevalence test compares prevalence between the two groups and says nothing useful here.")
  }
  message("[decontam_filter] ", sum(isNeg), " negative control(s), ", sum(!isNeg), " true sample(s).")
}

conc <- NULL
if (needsConc) {
  if (!nzchar(opt$conc_column)) {
    stop("--conc-column is required for method '", opt$method,
         "'. It names the metadata column holding post-PCR DNA concentration.")
  }
  if (!opt$conc_column %in% colnames(meta)) {
    stop("Column '", opt$conc_column, "' is not in '", opt$metadata, "'.")
  }
  conc <- suppressWarnings(as.numeric(meta[[opt$conc_column]]))
  if (anyNA(conc) || any(conc <= 0, na.rm = TRUE)) {
    stop("Column '", opt$conc_column, "' must be positive numbers - decontam's frequency ",
         "model regresses abundance on concentration, so a zero or a missing value has no meaning.")
  }
}

## ---- batches ------------------------------------------------------------
# A contaminant is a property of a KIT, not of a study. Pool two sequencing
# runs, two extraction days or two centres and the same taxon can be heavy in
# one batch and absent from the other; scored across the pool it looks like a
# taxon that varies between samples, which is exactly what a real organism
# looks like. decontam's own answer is to identify contaminants independently
# within each batch and combine the verdicts, so a taxon has to be judged
# against the blanks that were processed alongside it.
#
# The combination rule matters as much as the split. "minimum" (decontam's
# default) takes the smallest p across batches - a taxon condemned in any one
# batch is condemned overall, which is the right stance when a contaminant may
# only have been present in one run. "product" and "fisher" require agreement
# across batches and are correspondingly conservative.
batch <- NULL
if (nzchar(opt$batch_column)) {
  if (!opt$batch_column %in% colnames(meta)) {
    stop("Column '", opt$batch_column, "' is not in '", opt$metadata, "'. Columns present: ",
         paste(colnames(meta), collapse = ", "))
  }
  if (!opt$batch_combine %in% c("minimum", "product", "fisher")) {
    stop("--batch-combine must be one of minimum, product, fisher; got '", opt$batch_combine, "'.")
  }
  batch <- factor(trimws(as.character(meta[[opt$batch_column]])))
  if (anyNA(batch)) {
    stop("Column '", opt$batch_column, "' has missing values. Every sample must be assigned ",
         "to a batch, or the samples with no batch are silently judged against blanks that ",
         "never touched them.")
  }
  if (nlevels(batch) < 2) {
    stop("Column '", opt$batch_column, "' has one level ('", levels(batch)[1], "'). Batch-wise ",
         "identification of a single batch is the same run as without --decontam-batch-column, ",
         "so name the column that actually separates the runs, or drop the option.")
  }

  # Fail here rather than let decontam return NA for a whole batch. A batch
  # without its own control contributes nothing to the prevalence test, and a
  # taxon absent from every OTHER batch would then never be scored at all.
  sizes <- table(batch)
  if (needsNeg) {
    perBatch <- table(batch, ifelse(isNeg, "control", "sample"))
    starved <- rownames(perBatch)[perBatch[, "control"] == 0]
    if (length(starved)) {
      stop("Batch(es) with no negative control: ", paste(starved, collapse = ", "),
           ". The prevalence method compares each taxon against the blanks of its OWN batch, ",
           "so a batch without one cannot be scored. Either merge it into a batch that shares ",
           "its controls, or drop --decontam-batch-column and accept the pooled test.")
    }
    thin <- rownames(perBatch)[perBatch[, "sample"] < 2]
    if (length(thin)) {
      stop("Batch(es) with fewer than two true samples: ", paste(thin, collapse = ", "),
           ". Prevalence within such a batch is 0 or 1 and carries no information.")
    }
  }
  message("[decontam_filter] ", nlevels(batch), " batch(es) from '", opt$batch_column, "': ",
          paste(paste0(names(sizes), " (n=", as.integer(sizes), ")"), collapse = ", "),
          "; combined by '", opt$batch_combine, "'.")
}

## ---- decontam -----------------------------------------------------------
suppressPackageStartupMessages(library(decontam))

# decontam wants samples as ROWS; the Bracken table is taxa by sample.
seqtab <- t(counts)

res <- isContaminant(
  seqtab,
  conc = conc,
  neg = isNeg,
  method = opt$method,
  batch = batch,
  batch.combine = opt$batch_combine,
  threshold = opt$threshold,
  normalize = TRUE
)

res$name <- rownames(res)
res$taxid <- if ("taxonomy_id" %in% colnames(tab)) {
  tab$taxonomy_id[match(res$name, make.unique(as.character(tab$name)))]
} else {
  NA
}
res$reads <- rowSums(counts)[res$name]
res$contaminant[is.na(res$contaminant)] <- FALSE

# Which batches condemned it, not just whether the combination did. With
# `batch.combine = "minimum"` a taxon flagged in one batch of six is flagged
# overall, and that is a materially different claim from one flagged in all
# six: the first is a contaminant of one run, the second is a contaminant of
# the kit. decontam returns only the combined verdict, so the per-batch scores
# are recovered by scoring each batch on its own. Cheap - the table is taxa by
# a handful of samples - and it is the whole reason for batching.
if (!is.null(batch)) {
  flagged <- matrix(FALSE, nrow = nrow(res), ncol = nlevels(batch),
                    dimnames = list(rownames(res), levels(batch)))
  for (level in levels(batch)) {
    inBatch <- batch == level
    one <- tryCatch(
      isContaminant(
        seqtab[inBatch, , drop = FALSE],
        conc = if (is.null(conc)) NULL else conc[inBatch],
        neg = if (is.null(isNeg)) NULL else isNeg[inBatch],
        method = opt$method,
        threshold = opt$threshold,
        normalize = TRUE
      ),
      error = function(e) {
        message("[decontam_filter] batch '", level, "' could not be scored on its own (",
                conditionMessage(e), "); its column is reported as FALSE throughout.")
        NULL
      }
    )
    if (!is.null(one)) {
      call <- one[rownames(res), "contaminant"]
      call[is.na(call)] <- FALSE
      flagged[, level] <- call
    }
  }
  res$n_batches_flagged <- rowSums(flagged)
  res$batches_flagged <- apply(flagged, 1, function(row) {
    hits <- colnames(flagged)[row]
    if (length(hits)) paste(hits, collapse = ";") else ""
  })
}

ordered <- res[order(res$p, -res$reads), ]
keep <- intersect(
  c("taxid", "name", "reads", "freq", "prev", "p.freq", "p.prev", "p",
    "n_batches_flagged", "batches_flagged", "contaminant"),
  colnames(ordered)
)
write.table(ordered[, keep, drop = FALSE],
            paste0(opt$prefix, ".decontam_evidence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Only taxa carrying a taxid can be handed to the abundance filter, which keys
# on the taxid column. A contaminant without one is reported and said out loud
# rather than dropped silently from the drop list.
contaminants <- ordered[ordered$contaminant %in% TRUE, , drop = FALSE]
withTaxid <- contaminants[!is.na(contaminants$taxid), , drop = FALSE]
writeLines(as.character(withTaxid$taxid), paste0(opt$prefix, ".decontam_drop.txt"))
if (nrow(contaminants) > nrow(withTaxid)) {
  message("[decontam_filter] ", nrow(contaminants) - nrow(withTaxid),
          " contaminant(s) had no taxonomy_id and could not be put on the drop list.")
}

mqc <- paste0(opt$prefix, "_decontam_mqc.tsv")
writeLines(c(
  "# id: 'reanatax_decontam'",
  "# section_name: 'Contaminant identification'",
  "# description: 'Taxa scored by decontam against the negative controls. A taxon commoner",
  "#     in the blanks than in the true samples is called a contaminant. This is the only",
  "#     filter here that can distinguish a reagent contaminant from a genuinely rare",
  "#     organism, because it is the only one given an external measurement of the kit.'",
  "# plot_type: 'bargraph'",
  "# pconfig:",
  "#     id: 'reanatax_decontam_plot'",
  "#     title: 'reanaTax: decontam'",
  "#     ylab: 'Taxa'",
  "Sample\tContaminant\tNot contaminant",
  paste0("all taxa\t", nrow(contaminants), "\t", nrow(ordered) - nrow(contaminants))
), mqc)

message("[decontam_filter] method '", opt$method, "', threshold ", opt$threshold, ": ",
        nrow(contaminants), "/", nrow(ordered), " taxa called contaminants; ",
        nrow(withTaxid), " taxid(s) listed for removal.")

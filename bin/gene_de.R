#!/usr/bin/env Rscript
#
# Differential expression of GENE counts between two groups of samples.
#
# This is deliberately not the same machinery as bin/differential_abundance.R.
# ALDEx2 and ANCOM-BC2 model compositional data: a microbial profile carries no
# absolute scale, only proportions, and their whole purpose is to recover
# something interpretable from that. Gene counts are not compositional in that
# sense - library size is a nuisance that DESeq2's median-of-ratios and edgeR's
# TMM estimate directly - so the right models here are negative binomial.
#
# The contrast comes from the same --da_* metadata as the microbial side, so one
# metadata file and one design describe both halves of the same library.

# Arguments are parsed by hand rather than with optparse, following
# bin/differential_abundance.R: the Bioconductor containers carry their own
# stack and little else, and a missing CRAN package would fail the task for no
# good reason.
.defaults <- list(
  counts = NA, metadata = NA, method = "deseq2", grouping = NA,
  sample_group = NA, ref_group = NA, covariates = "", samples_to_drop = "",
  min_count = 10, min_samples = 2, p_threshold = 0.05, lfc_threshold = 1.0,
  output = "gene_de.results.tsv", mqc = ""
)
.numeric <- c("min_count", "min_samples", "p_threshold", "lfc_threshold")

parseArgs <- function(argv, defaults) {
  opt <- defaults
  i <- 1
  while (i <= length(argv)) {
    key <- argv[i]
    if (!startsWith(key, "--")) stop("Unexpected argument '", key, "'.")
    name <- gsub("-", "_", sub("^--", "", key))
    if (!name %in% names(defaults)) stop("Unknown option '", key, "'.")
    if (i + 1 > length(argv)) stop("Option '", key, "' needs a value.")
    value <- argv[i + 1]
    opt[[name]] <- if (name %in% .numeric) as.numeric(value) else value
    i <- i + 2
  }
  opt
}

opt <- parseArgs(commandArgs(trailingOnly = TRUE), .defaults)
for (required in c("counts", "metadata", "grouping", "sample_group", "ref_group")) {
  if (is.na(opt[[required]])) stop("--", gsub("_", "-", required), " is required.")
}

fail <- function(...) stop(paste0("[gene_de] ", ...), call. = FALSE)

counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
meta <- read.delim(opt$metadata, check.names = FALSE)
rownames(meta) <- as.character(meta[[1]])

if (nzchar(opt$samples_to_drop)) {
  drop <- trimws(strsplit(opt$samples_to_drop, ",")[[1]])
  counts <- counts[, !(colnames(counts) %in% drop), drop = FALSE]
}

shared <- intersect(colnames(counts), rownames(meta))
if (length(shared) < 2) {
  fail("only ", length(shared), " sample(s) are in both the count matrix and ",
       opt$metadata, ". The metadata's first column must hold the same sample ids ",
       "the pipeline uses (", paste(head(colnames(counts), 3), collapse = ", "), ", ...).")
}
counts <- counts[, shared, drop = FALSE]
meta <- meta[shared, , drop = FALSE]

if (!opt$grouping %in% colnames(meta)) {
  fail("--da_grouping '", opt$grouping, "' is not a column of ", opt$metadata,
       " (have: ", paste(colnames(meta), collapse = ", "), ")")
}

group <- as.character(meta[[opt$grouping]])
keep <- group %in% c(opt$sample_group, opt$ref_group)
counts <- counts[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
group <- factor(group[keep], levels = c(opt$ref_group, opt$sample_group))

n_ref <- sum(group == opt$ref_group)
n_test <- sum(group == opt$sample_group)
if (n_ref < 2 || n_test < 2) {
  fail("differential expression needs at least two samples per group; got ",
       n_test, " '", opt$sample_group, "' and ", n_ref, " '", opt$ref_group,
       "'. With one replicate the dispersion cannot be estimated and any p-value ",
       "would be an artefact of the model's fallback, not evidence.")
}

covariates <- character(0)
if (nzchar(opt$covariates)) {
  covariates <- trimws(strsplit(opt$covariates, ",")[[1]])
  missing <- setdiff(covariates, colnames(meta))
  if (length(missing)) fail("--da_covariates not in metadata: ", paste(missing, collapse = ", "))
}

# Genes with almost no reads carry no information and cost multiple-testing
# power, so they are removed before the test rather than filtered from the
# results afterwards - filtering after the fact would leave the adjusted
# p-values computed against a gene set that was never really tested.
expressed <- rowSums(counts >= opt$min_count) >= opt$min_samples
message("[gene_de] ", sum(expressed), " of ", nrow(counts), " genes pass the expression filter")
counts <- counts[expressed, , drop = FALSE]
if (nrow(counts) < 2) fail("fewer than two genes survived the expression filter")

design_data <- data.frame(group = group, row.names = colnames(counts))
for (cv in covariates) design_data[[cv]] <- factor(as.character(meta[[cv]]))
design_formula <- as.formula(paste("~", paste(c(covariates, "group"), collapse = " + ")))

if (opt$method == "deseq2") {
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), design_data, design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", opt$sample_group, opt$ref_group))
  out <- data.frame(
    gene_id = rownames(res),
    base_mean = res$baseMean,
    log2fc = res$log2FoldChange,
    lfc_se = res$lfcSE,
    stat = res$stat,
    p_value = res$pvalue,
    q_value = res$padj,
    stringsAsFactors = FALSE
  )
} else if (opt$method == "edger") {
  suppressPackageStartupMessages(library(edgeR))
  dge <- DGEList(counts = as.matrix(counts), group = group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(design_formula, data = design_data)
  dge <- estimateDisp(dge, design)
  # Quasi-likelihood F, not the exact test: it is the better-calibrated choice
  # at low replication, and it is the only one of the two that admits the
  # covariates in the design above.
  fit <- glmQLFit(dge, design)
  coef_name <- paste0("group", opt$sample_group)
  if (!coef_name %in% colnames(design)) fail("contrast coefficient '", coef_name, "' not in the design")
  qlf <- glmQLFTest(fit, coef = coef_name)
  tab <- topTags(qlf, n = Inf, sort.by = "none")$table
  out <- data.frame(
    gene_id = rownames(tab),
    base_mean = 2^tab$logCPM,
    log2fc = tab$logFC,
    lfc_se = NA_real_,
    stat = tab$F,
    p_value = tab$PValue,
    q_value = tab$FDR,
    stringsAsFactors = FALSE
  )
} else {
  fail("unknown --method '", opt$method, "'")
}

out$significant <- !is.na(out$q_value) & out$q_value < opt$p_threshold &
  !is.na(out$log2fc) & abs(out$log2fc) >= opt$lfc_threshold
out <- out[order(out$q_value, na.last = TRUE), ]
write.table(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

n_sig <- sum(out$significant)
message("[gene_de] ", opt$method, ": ", n_sig, " genes at q<", opt$p_threshold,
        " and |log2FC|>=", opt$lfc_threshold)

if (nzchar(opt$mqc)) {
  up <- sum(out$significant & out$log2fc > 0)
  down <- sum(out$significant & out$log2fc < 0)
  writeLines(c(
    paste0("# id: 'host_de_", opt$method, "'"),
    paste0("# section_name: 'Host gene DE (", opt$method, ")'"),
    "# format: 'tsv'",
    "# plot_type: 'bargraph'",
    paste0("# description: '", opt$sample_group, " vs ", opt$ref_group,
           ": genes at q<", opt$p_threshold, " with |log2FC|>=", opt$lfc_threshold,
           ", from ", nrow(out), " tested.'"),
    "Sample\tUp\tDown",
    paste0(opt$method, "\t", up, "\t", down)
  ), opt$mqc)
}

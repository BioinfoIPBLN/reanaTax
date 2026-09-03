#!/usr/bin/env Rscript
# differential_abundance.R -- test which taxa differ between two groups.
#
# Modelled on FGCZ's EzAppDiffShot (ezRun/R/app-diffShot.R and its per-method
# child Rmds), reduced to what this pipeline already has. DiffShot stages
# per-sample Bracken reports, builds a BIOM with kraken-biom and embeds the
# metadata through the biom Python API; here the combined Bracken table is
# already a taxa-by-sample count matrix, so the BIOM round-trip buys nothing
# and is skipped. The filters, the grouping/sampleGroup/refGroup contrast
# convention and the covariate handling are DiffShot's.
#
# Both supported methods need INTEGER counts: ALDEx2 draws Monte-Carlo
# instances from a Dirichlet posterior, ANCOM-BC2 estimates a per-sample
# sampling fraction from library size. So the `_num` columns are used, never
# `_frac`, and MetaPhlAn relative-abundance profiles are rejected rather than
# silently rescaled.

# Arguments are parsed by hand rather than with optparse: the ALDEx2 and
# ANCOMBC biocontainers carry their Bioconductor stack and little else, and a
# missing CRAN package would fail the task for no good reason.
.defaults <- list(
  counts = NA, metadata = NA, method = "ancombc", grouping = NA,
  sample_group = NA, ref_group = NA, covariates = "",
  prevalence_min = 0.10, relabund_min = 0.10, lib_cut = 1000,
  samples_to_drop = "", p_threshold = 0.05, lfc_threshold = 1.0,
  prefix = "differential_abundance"
)
.numeric <- c("prevalence_min", "relabund_min", "lib_cut", "p_threshold", "lfc_threshold")

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

splitList <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

## ---- counts -------------------------------------------------------------
tab <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
numCols <- grep("_num$", colnames(tab), value = TRUE)
if (length(numCols) == 0) {
  stop("No `<sample>_num` columns in '", opt$counts, "'. This needs the combined ",
       "Bracken table (bracken_combined_<level>.txt). MetaPhlAn profiles carry ",
       "relative abundances, which neither ALDEx2 nor ANCOM-BC2 can model.")
}
if (!"name" %in% colnames(tab)) stop("No `name` column in '", opt$counts, "'.")

counts <- as.matrix(tab[, numCols, drop = FALSE])
mode(counts) <- "numeric"
rownames(counts) <- make.unique(as.character(tab$name))
colnames(counts) <- sub("_num$", "", numCols)
counts <- round(counts)
storage.mode(counts) <- "integer"

## ---- metadata -----------------------------------------------------------
meta <- read.delim(opt$metadata, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(meta) < 2) stop("--metadata needs a sample-id column plus at least the grouping column.")
rownames(meta) <- as.character(meta[[1]])

drop <- splitList(opt$samples_to_drop)
keep <- setdiff(intersect(colnames(counts), rownames(meta)), drop)
missing <- setdiff(colnames(counts), c(rownames(meta), drop))
if (length(missing)) {
  message("[differential_abundance] no metadata row for: ", paste(missing, collapse = ", "),
          " -- excluded.")
}
if (length(keep) < 3) stop("Only ", length(keep), " sample(s) have both counts and metadata; need at least 3.")
counts <- counts[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

if (!opt$grouping %in% colnames(meta)) {
  stop("--grouping '", opt$grouping, "' is not a column of --metadata. Available: ",
       paste(colnames(meta)[-1], collapse = ", "))
}
covariates <- splitList(opt$covariates)
badCov <- setdiff(covariates, colnames(meta))
if (length(badCov)) stop("--covariates not in --metadata: ", paste(badCov, collapse = ", "))

group <- as.character(meta[[opt$grouping]])
if (is.null(opt$sample_group) || is.null(opt$ref_group)) {
  stop("--sample-group and --ref-group are required.")
}
if (identical(opt$sample_group, opt$ref_group)) stop("--sample-group and --ref-group must differ.")
inContrast <- group %in% c(opt$sample_group, opt$ref_group)
if (sum(group == opt$sample_group) < 2 || sum(group == opt$ref_group) < 2) {
  stop("Each side of the contrast needs at least 2 samples; got ",
       sum(group == opt$sample_group), " vs ", sum(group == opt$ref_group), ".")
}
counts <- counts[, inContrast, drop = FALSE]
meta <- meta[inContrast, , drop = FALSE]
# refGroup first so every model reads sampleGroup-over-refGroup.
meta[[opt$grouping]] <- factor(as.character(meta[[opt$grouping]]),
                               levels = c(opt$ref_group, opt$sample_group))

## ---- filters (DiffShot's) ----------------------------------------------
libSize <- colSums(counts)
prevalence <- rowSums(counts > 1) / ncol(counts)
pooled <- rowSums(counts) / max(sum(counts), 1) * 100
keepTaxa <- prevalence >= opt$prevalence_min & pooled >= opt$relabund_min
message(sprintf("[differential_abundance] %d/%d taxa pass prevalence >= %.2f and pooled abundance >= %.2f%%.",
                sum(keepTaxa), length(keepTaxa), opt$prevalence_min, opt$relabund_min))
counts <- counts[keepTaxa, , drop = FALSE]
if (nrow(counts) < 2) stop("Fewer than 2 taxa survive the filters; loosen --da_prevalence_min / --da_relabund_min.")

method <- tolower(gsub("[^A-Za-z0-9]", "", opt$method))

## ---- test ---------------------------------------------------------------
if (method %in% c("aldex2", "aldex")) {
  suppressPackageStartupMessages(library(ALDEx2))
  covFormula <- length(covariates) > 0
  if (covFormula) {
    mm <- model.matrix(as.formula(paste("~", paste(c(opt$grouping, covariates), collapse = " + "))), data = meta)
    clr <- aldex.clr(counts, mm, mc.samples = 128, denom = "all", verbose = FALSE)
    glmRes <- aldex.glm(clr, mm)
    coefName <- grep(paste0("^", opt$grouping, opt$sample_group), colnames(glmRes), value = TRUE)
    # ALDEx2 names the glm columns "<term>:Est", "<term>:pval", "<term>:pval.padj"
    # (older releases used "<term>Estimate"/"pval"/"pval.holm") - accept both.
    est <- grep("(Estimate|:Est)$", coefName, value = TRUE)
    pv  <- grep("(^|:)pval$", coefName, value = TRUE)
    qv  <- grep("(pval\\.(padj|holm)|pval\\.fdr)$", coefName, value = TRUE)
    if (length(est) == 0 || length(pv) == 0) {
      stop("ALDEx2 returned no coefficient for '", opt$grouping, "' (columns: ",
           paste(colnames(glmRes), collapse = ", "), ").")
    }
    res <- data.frame(taxon = rownames(glmRes),
                      lfc = glmRes[[est[1]]],
                      pvalue = glmRes[[pv[1]]],
                      qvalue = if (length(qv)) glmRes[[qv[1]]] else p.adjust(glmRes[[pv[1]]], "BH"),
                      stringsAsFactors = FALSE)
  } else {
    clr <- aldex.clr(counts, as.character(meta[[opt$grouping]]), mc.samples = 128, denom = "all", verbose = FALSE)
    tt <- aldex.ttest(clr)
    ef <- aldex.effect(clr)
    res <- data.frame(taxon = rownames(tt),
                      lfc = ef$diff.btw,
                      pvalue = tt$we.ep,
                      qvalue = tt$we.eBH,
                      stringsAsFactors = FALSE)
  }
} else if (method %in% c("ancombc", "ancombc2")) {
  suppressPackageStartupMessages(library(ANCOMBC))
  # The ancombc biocontainer ships ANCOMBC without TreeSummarizedExperiment
  # (a Suggests-only dependency, so the Galaxy build omits it). ancombc2()
  # accepts a plain count matrix with tax_level = "none" + a metadata data.frame
  # whose first column is the sample id and whose rownames are the sample ids -
  # exactly what `meta` already is, so no TSE object is needed.
  fix <- paste(c(opt$grouping, covariates), collapse = " + ")
  out <- ancombc2(data = counts, tax_level = "none", fix_formula = fix,
                  p_adj_method = "BH", prv_cut = 0, lib_cut = opt$lib_cut,
                  group = opt$grouping, struc_zero = TRUE, alpha = opt$p_threshold,
                  n_cl = 1, verbose = FALSE, meta_data = meta)
  primary <- out$res
  coef <- grep(paste0("^lfc_", opt$grouping), colnames(primary), value = TRUE)[1]
  if (is.na(coef)) stop("ANCOM-BC2 returned no coefficient for '", opt$grouping, "'.")
  suffix <- sub("^lfc_", "", coef)
  res <- data.frame(taxon = primary$taxon,
                    lfc = primary[[coef]],
                    pvalue = primary[[paste0("p_", suffix)]],
                    qvalue = primary[[paste0("q_", suffix)]],
                    stringsAsFactors = FALSE)
} else {
  stop("Unknown --method '", opt$method, "'; expected aldex2 or ancombc.")
}

res <- res[order(res$qvalue, res$pvalue), ]
res$significant <- !is.na(res$qvalue) & res$qvalue < opt$p_threshold &
                   !is.na(res$lfc) & abs(res$lfc) >= opt$lfc_threshold
comparison <- paste0(opt$sample_group, "_over_", opt$ref_group)
res$comparison <- comparison

write.table(res, sprintf("%s.%s.results.tsv", opt$prefix, method),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

## ---- volcano ------------------------------------------------------------
png(sprintf("%s.%s.volcano.png", opt$prefix, method), width = 1400, height = 1100, res = 150)
ok <- !is.na(res$lfc) & !is.na(res$pvalue) & res$pvalue > 0
plot(res$lfc[ok], -log10(res$pvalue[ok]),
     pch = 19, cex = 0.7,
     col = ifelse(res$significant[ok], "#c0392b", "#95a5a6"),
     xlab = sprintf("log2 fold change (%s over %s)", opt$sample_group, opt$ref_group),
     ylab = "-log10(p)",
     main = sprintf("%s - %s", toupper(method), comparison))
abline(h = -log10(opt$p_threshold), lty = 2, col = "grey40")
abline(v = c(-opt$lfc_threshold, opt$lfc_threshold), lty = 2, col = "grey40")
top <- head(res[ok & res$significant, , drop = FALSE], 15)
if (nrow(top)) text(top$lfc, -log10(top$pvalue), labels = top$taxon, pos = 3, cex = 0.55)
invisible(dev.off())

## ---- MultiQC section ----------------------------------------------------
mqc <- sprintf("%s.%s_mqc.tsv", opt$prefix, method)
writeLines(c(
  sprintf("# id: 'reanatax_da_%s'", method),
  sprintf("# section_name: 'Differential abundance (%s)'", toupper(method)),
  sprintf("# description: 'Taxa differing between %s and %s, %s on the combined Bracken counts. Ranked by adjusted p-value; the full table and a volcano plot are in the differential_abundance/ output directory.'",
          opt$sample_group, opt$ref_group, toupper(method)),
  "# plot_type: 'table'",
  sprintf("# pconfig:\n#     id: 'reanatax_da_%s_table'\n#     title: 'reanaTax: differential abundance'", method),
  paste(c("Taxon", "log2FC", "p", "adj.p"), collapse = "\t")
), mqc)
topN <- head(res, 25)
write.table(data.frame(topN$taxon, signif(topN$lfc, 4), signif(topN$pvalue, 4), signif(topN$qvalue, 4)),
            mqc, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE, na = "")

message(sprintf("[differential_abundance] %s: %d taxa tested, %d significant (adj.p < %.3g, |log2FC| >= %.3g).",
                method, nrow(res), sum(res$significant, na.rm = TRUE), opt$p_threshold, opt$lfc_threshold))

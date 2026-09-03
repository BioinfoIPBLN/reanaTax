#!/usr/bin/env Rscript
# diversity.R -- alpha diversity, and what explains the variation between samples.
#
# The pipeline already answers "which taxa differ between two groups" with
# ALDEx2 and ANCOM-BC2. It has never answered the question that usually comes
# first when reanalysing someone else's dataset: how much of the variation in
# this cohort does each recorded variable actually explain, and which variables
# are redundant with each other? A dataset where 30% of the community variation
# tracks sequencing batch and 2% tracks the disease of interest is telling you
# something no per-taxon test will.
#
# The design follows Lloréns-Rico et al. (Nat Commun 2021), whose analysis code
# is public (raeslab/covid19_respiratory_microbiome):
#
#   * Counts are CLR-transformed and compared by Euclidean distance. That pair
#     is the Aitchison distance, and it is used rather than Bray-Curtis because
#     sequencing yields compositions: only ratios between taxa carry
#     information, and CLR is the transform that makes a distance respect that.
#   * Each variable is tested on its own with a distance-based RDA
#     (vegan::capscale), giving an adjusted R^2 - the share of variation it
#     explains - and a permutation p-value, BH-corrected across variables.
#   * A stepwise model (vegan::ordiR2step) then adds variables in order of the
#     variation they explain that the ones already in the model do not. The gap
#     between a variable's univariate R^2 and its contribution here IS the
#     confounding: two variables can each explain 20% and jointly explain 21%.
#   * PERMANOVA (vegan::adonis2) on the same distance, for the categorical
#     variables, since that is what most readers will expect to see quoted.
#
# Metadata is optional. Without it there is still alpha diversity and an
# unconstrained ordination; the variance partitioning simply has nothing to
# partition by and is skipped rather than faked.

.defaults <- list(
  counts = NA, metadata = "", prefix = "diversity",
  permutations = 999, min_samples = 3, max_levels = 0
)
.numeric <- c("permutations", "min_samples", "max_levels")

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
if (is.na(opt$counts)) stop("--counts is required.")

suppressPackageStartupMessages(library(vegan))

## ---- counts -------------------------------------------------------------
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
counts[is.na(counts)] <- 0

# Samples as rows, which is what vegan expects throughout.
mat <- t(counts)
mat <- mat[, colSums(mat) > 0, drop = FALSE]
if (ncol(mat) < 2) stop("Fewer than two taxa with any reads; nothing to measure.")

## ---- alpha --------------------------------------------------------------
# On raw counts, deliberately: Shannon and Simpson are defined on the
# composition and a CLR-transformed matrix has negative entries.
observed <- rowSums(mat > 0)
shannon <- diversity(mat, index = "shannon")
simpson <- diversity(mat, index = "simpson")
invsimpson <- diversity(mat, index = "invsimpson")
# Pielou's evenness: observed diversity over the maximum possible at that
# richness. Separates "many taxa" from "no taxon dominating".
evenness <- ifelse(observed > 1, shannon / log(observed), NA_real_)

alpha <- data.frame(
  sample = rownames(mat), reads = rowSums(mat), observed = observed,
  shannon = shannon, simpson = simpson, invsimpson = invsimpson,
  evenness = evenness, row.names = NULL, stringsAsFactors = FALSE
)
write.table(alpha, paste0(opt$prefix, ".alpha_diversity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

## ---- Aitchison distance -------------------------------------------------
# Zeros have to go before the log. Each zero is replaced by half the smallest
# non-zero count in its own sample - the usual multiplicative-simple
# replacement. It is a simplification: a principled treatment (zCompositions'
# Bayesian-multiplicative replacement) models the zeros rather than imputing
# them, and matters more as sparsity rises. Stated here rather than buried,
# because CLR on sparse data is sensitive to this choice.
replaceZeros <- function(row) {
  nonzero <- row[row > 0]
  if (!length(nonzero)) return(row + 1)
  row[row == 0] <- min(nonzero) / 2
  row
}
positive <- t(apply(mat, 1, replaceZeros))
logged <- log(positive)
clr <- logged - rowMeans(logged)
aitchison <- dist(clr, method = "euclidean")

write.table(as.matrix(aitchison), paste0(opt$prefix, ".aitchison_distance.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

## ---- unconstrained ordination ------------------------------------------
ordination <- NULL
if (nrow(mat) >= 3) {
  pcoa <- cmdscale(aitchison, k = min(2, nrow(mat) - 1), eig = TRUE)
  variance <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])
  ordination <- data.frame(
    sample = rownames(mat),
    PCo1 = pcoa$points[, 1],
    PCo2 = if (ncol(pcoa$points) > 1) pcoa$points[, 2] else NA_real_,
    row.names = NULL, stringsAsFactors = FALSE
  )
  attr(ordination, "variance") <- variance[1:2]
  write.table(ordination, paste0(opt$prefix, ".ordination.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

## ---- variance partitioning ---------------------------------------------
betaRows <- NULL
if (nzchar(opt$metadata) && file.exists(opt$metadata)) {
  meta <- read.delim(opt$metadata, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(meta) <- as.character(meta[[1]])
  shared <- intersect(rownames(mat), rownames(meta))

  if (length(shared) < opt$min_samples) {
    message("[diversity] only ", length(shared), " sample(s) shared with the metadata; ",
            "variance partitioning needs at least ", opt$min_samples, " and is skipped.")
  } else {
    metaShared <- meta[shared, , drop = FALSE]
    distShared <- dist(clr[shared, , drop = FALSE], method = "euclidean")
    n <- length(shared)
    # A variable with as many distinct values as samples saturates the model:
    # it would explain everything and mean nothing. --max-levels 0 sets the
    # cap from the sample count instead of fixing it.
    maxLevels <- if (opt$max_levels > 0) opt$max_levels else max(2, n - 1)

    usable <- character(0)
    for (column in setdiff(colnames(metaShared), colnames(metaShared)[1])) {
      values <- metaShared[[column]]
      if (all(is.na(values))) next
      distinct <- length(unique(values[!is.na(values)]))
      if (distinct < 2) next
      if (!is.numeric(values) && distinct > maxLevels) {
        message("[diversity] '", column, "' has ", distinct, " levels over ", n,
                " samples; skipped as saturated.")
        next
      }
      if (anyNA(values)) {
        message("[diversity] '", column, "' has missing values; skipped.")
        next
      }
      usable <- c(usable, column)
    }

    if (!length(usable)) {
      message("[diversity] no metadata column varies usefully across these samples.")
    } else {
      rows <- list()
      for (column in usable) {
        frame <- data.frame(value = metaShared[[column]], row.names = shared)
        if (!is.numeric(frame$value)) frame$value <- factor(frame$value)
        fit <- try(capscale(distShared ~ value, data = frame), silent = TRUE)
        if (inherits(fit, "try-error")) {
          message("[diversity] capscale failed for '", column, "'; skipped.")
          next
        }
        test <- try(anova.cca(fit, permutations = opt$permutations), silent = TRUE)
        perm <- try(adonis2(distShared ~ value, data = frame,
                            permutations = opt$permutations), silent = TRUE)
        rows[[column]] <- data.frame(
          variable = column,
          levels = if (is.factor(frame$value)) nlevels(frame$value) else NA_integer_,
          r2_adj = RsquareAdj(fit)$adj.r.squared,
          F = if (inherits(test, "try-error")) NA_real_ else test$F[1],
          p = if (inherits(test, "try-error")) NA_real_ else test$`Pr(>F)`[1],
          permanova_r2 = if (inherits(perm, "try-error")) NA_real_ else perm$R2[1],
          permanova_p = if (inherits(perm, "try-error")) NA_real_ else perm$`Pr(>F)`[1],
          stringsAsFactors = FALSE
        )
      }
      betaRows <- do.call(rbind, rows)

      if (!is.null(betaRows)) {
        betaRows$p_adj <- p.adjust(betaRows$p, method = "BH")
        betaRows$cumulative_r2_adj <- NA_real_

        # Stepwise, for the non-redundant share. Only attempted when there are
        # enough samples to fit the full model; ordiR2step on a saturated one
        # returns the full model and says nothing.
        if (length(usable) > 1 && n > length(usable) + 1) {
          frame <- metaShared[, usable, drop = FALSE]
          for (column in usable) {
            if (!is.numeric(frame[[column]])) frame[[column]] <- factor(frame[[column]])
          }
          null <- try(capscale(distShared ~ 1, data = frame), silent = TRUE)
          full <- try(capscale(distShared ~ ., data = frame), silent = TRUE)
          if (!inherits(null, "try-error") && !inherits(full, "try-error")) {
            step <- try(ordiR2step(null, scope = formula(full), direction = "forward",
                                  permutations = opt$permutations, trace = FALSE),
                        silent = TRUE)
            if (!inherits(step, "try-error") && !is.null(step$anova)) {
              chosen <- gsub("^\\+ ", "", rownames(step$anova))
              cumulative <- step$anova$R2.adj
              hit <- match(betaRows$variable, chosen)
              betaRows$cumulative_r2_adj <- cumulative[hit]
            }
          }
        }
        betaRows <- betaRows[order(-betaRows$r2_adj), ]
        write.table(betaRows, paste0(opt$prefix, ".beta_variance.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
      }
    }
  }
} else if (nzchar(opt$metadata)) {
  message("[diversity] metadata '", opt$metadata, "' not found; alpha diversity only.")
}

## ---- MultiQC ------------------------------------------------------------
alphaOut <- alpha
alphaOut[] <- lapply(alphaOut, function(column) {
  if (is.numeric(column)) round(column, 4) else column
})
writeLines(c(
  "# id: 'reanatax_alpha_diversity'",
  "# section_name: 'Alpha diversity'",
  "# description: 'Within-sample diversity of the Bracken profile. Observed is the number of",
  "#     taxa seen; Shannon and Simpson weight that by how evenly reads are spread over them;",
  "#     evenness (Pielou) is Shannon divided by its maximum at that richness, which separates",
  "#     >many taxa< from >no taxon dominating<. All are sensitive to sequencing depth, so read",
  "#     them next to the reads column.'",
  "# plot_type: 'table'",
  "# pconfig:",
  "#     id: 'reanatax_alpha_diversity_table'",
  "#     title: 'reanaTax: alpha diversity'",
  paste(colnames(alphaOut), collapse = "\t"),
  apply(alphaOut, 1, function(row) paste(row, collapse = "\t"))
), paste0(opt$prefix, "_alpha_diversity_mqc.tsv"))

if (!is.null(betaRows)) {
  writeLines(c(
    "# id: 'reanatax_beta_variance'",
    "# section_name: 'Variance explained'",
    "# description: 'Share of between-sample variation (Aitchison distance) that each metadata",
    "#     variable explains on its own, by distance-based RDA. Read it before any per-taxon",
    "#     test: a variable that explains more than the one you care about is either a",
    "#     confounder or the actual story. Values are adjusted R-squared, so they are",
    "#     comparable between variables with different numbers of levels and can go negative",
    "#     when a variable explains less than chance.'",
    "# plot_type: 'bargraph'",
    "# pconfig:",
    "#     id: 'reanatax_beta_variance_plot'",
    "#     title: 'reanaTax: variance explained (adjusted R2)'",
    "#     ylab: 'Adjusted R2'",
    "Sample\tAdjusted R2",
    apply(betaRows, 1, function(row) {
      paste(row[["variable"]], max(0, round(as.numeric(row[["r2_adj"]]), 5)), sep = "\t")
    })
  ), paste0(opt$prefix, "_beta_variance_mqc.tsv"))
}

## ---- plots --------------------------------------------------------------
# Wrapped: a plotting failure must not lose the tables, which are the result.
try({
  png(paste0(opt$prefix, ".alpha_diversity.png"), width = 1200, height = 800, res = 130)
  op <- par(mar = c(9, 4, 3, 1))
  barplot(alpha$shannon, names.arg = alpha$sample, las = 2, ylab = "Shannon",
          main = "Alpha diversity (Shannon)", col = "#4C72B0", border = NA)
  par(op)
  dev.off()
}, silent = TRUE)

if (!is.null(ordination)) {
  try({
    variance <- attr(ordination, "variance")
    png(paste0(opt$prefix, ".ordination.png"), width = 1000, height = 900, res = 130)
    plot(ordination$PCo1, ordination$PCo2, pch = 19, col = "#4C72B0",
         xlab = sprintf("PCo1 (%.1f%%)", 100 * variance[1]),
         ylab = sprintf("PCo2 (%.1f%%)", 100 * variance[2]),
         main = "Aitchison PCoA")
    text(ordination$PCo1, ordination$PCo2, labels = ordination$sample, pos = 3, cex = 0.7)
    dev.off()
  }, silent = TRUE)
}

message("[diversity] ", nrow(alpha), " samples, ", ncol(mat), " taxa; ",
        if (is.null(betaRows)) "no variance partitioning"
        else paste0(nrow(betaRows), " variable(s) tested"), ".")

# ==============================================================================
# MODULE 2 · SCRIPT 4 of 7 — TRANSFORM, PCA, AND SAMPLE-DISTANCE EXPLORATION
# ==============================================================================
# ONE NEW IDEA on top of Script 3: Script 3 decided which genes to keep. This
# script looks at the data BEFORE fitting any model, to catch problems a
# p-value would hide - mislabeled samples, an unexpected outlier, or a batch
# effect that dominates the biological signal.
#
# TEACHING NOTE — why raw counts are never plotted directly:
#   Raw counts have variance that scales with the mean (a highly expressed
#   gene naturally has more absolute variability than a lowly expressed one).
#   Plotting raw counts on a PCA or distance heatmap lets the few
#   highest-count genes dominate the plot for a purely technical reason. The
#   Variance Stabilizing Transformation (VST) makes variance roughly constant
#   across the expression range so that PCA/clustering reflects genuine
#   sample-to-sample differences, not just which genes happen to be highly
#   expressed. blind=FALSE means the transformation is allowed to use the
#   known design when estimating dispersion (appropriate once you already
#   trust the design from Script 2) - use blind=TRUE only for a fully
#   unsupervised first look with no assumed design.
#
# Requires: DESeq2, ggplot2, pheatmap, RColorBrewer
# Input:    counts_filtered.tsv (from Script 3), metadata.csv
# Output:   pca_plot.pdf, sample_distance_heatmap.pdf; prints any samples
#           flagged as PCA outliers.
# ==============================================================================

suppressMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
counts_path   <- if (length(args) >= 1) args[1] else "counts_filtered.tsv"
metadata_path <- if (length(args) >= 2) args[2] else "metadata.csv"

cat("==================================================================\n")
cat(" MODULE 2 / SCRIPT 4 — Exploratory analysis (VST + PCA + distances)\n")
cat("==================================================================\n")

counts <- read.table(counts_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta   <- read.csv(metadata_path, stringsAsFactors = FALSE)
meta   <- meta[match(colnames(counts), meta$sample), ]
rownames(meta) <- meta$sample

has_batch <- "batch" %in% colnames(meta)
design_formula <- if (has_batch) ~ batch + condition else ~ condition

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)

# blind=FALSE: use the known design when estimating dispersions for the
# transform (see teaching note above). Appropriate here because Script 2
# already validated this design is not confounded.
#
# TEACHING NOTE - a real gotcha this script guards against: vst()'s default
# fast path fits its dispersion trend on a random SUBSET of genes (nsub=1000
# by default) and errors out if fewer than 1000 genes are available -
# something that happens routinely after aggressive filtering (Script 3),
# on a small custom gene panel, or on a bacterial/viral genome with only a
# few hundred genes total. Capping nsub at the actual gene count keeps this
# script working across genome sizes instead of failing only on small ones.
n_genes_available <- nrow(dds)
if (n_genes_available < 1000) {
  cat(sprintf(
    "[INFO] only %d genes available (< vst()'s default nsub=1000) - capping nsub to %d.\n",
    n_genes_available, n_genes_available))
  vsd <- vst(dds, blind = FALSE, nsub = n_genes_available)
} else {
  vsd <- vst(dds, blind = FALSE)
}

# --- PCA -----------------------------------------------------------------
pca_data <- plotPCA(vsd, intgroup = if (has_batch) c("condition", "batch") else "condition",
                     returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# NOTE: the shape aesthetic is built as a separate mapping and only added
# when a batch column exists. Putting a conditional expression like
# `if (has_batch) batch else NULL` directly inside aes() technically works,
# but ggplot2 deparses that expression as the legend title verbatim -
# producing a legend labeled "if (has_batch) batch else NULL" instead of
# "batch". Building the aes() list programmatically avoids that.
point_aes <- aes(x = PC1, y = PC2, color = condition)
if (has_batch) point_aes <- modifyList(point_aes, aes(shape = batch))

p <- ggplot(pca_data, point_aes) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#01696F", "#DA7101", "#7A39BB", "#A13544")) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  labs(shape = "batch") +
  ggtitle("PCA of variance-stabilized counts") +
  theme_bw()

ggsave("pca_plot.pdf", plot = p, width = 6, height = 5)
cat("Wrote pca_plot.pdf\n")

# --- Outlier flag: any sample > 3 SD from its group's PC1/PC2 centroid ----
outliers <- character(0)
for (cond in unique(pca_data$condition)) {
  grp <- pca_data[pca_data$condition == cond, ]
  if (nrow(grp) < 2) next
  centroid <- colMeans(grp[, c("PC1", "PC2")])
  dists <- sqrt((grp$PC1 - centroid[1])^2 + (grp$PC2 - centroid[2])^2)
  thresh <- mean(dists) + 3 * sd(dists)
  flagged <- grp$name[dists > thresh]
  outliers <- c(outliers, flagged)
}
if (length(outliers) > 0) {
  cat("[WARN] possible PCA outlier(s), inspect before excluding anything:",
      paste(outliers, collapse = ", "), "\n")
} else {
  cat("[PASS] no sample sits more than 3 SD from its group's PCA centroid\n")
}

# --- Sample-to-sample distance heatmap ------------------------------------
sample_dists <- dist(t(assay(vsd)))
dist_matrix <- as.matrix(sample_dists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf("sample_distance_heatmap.pdf", width = 6, height = 5)
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colors,
         main = "Sample-to-sample distances (VST)")
dev.off()
cat("Wrote sample_distance_heatmap.pdf\n")

cat("\n==================================================================\n")
cat(" EXPLORATION COMPLETE - inspect both PDFs for unexpected clustering\n")
cat(" before moving to Script 5 (design formula and contrasts).\n")
cat("==================================================================\n")

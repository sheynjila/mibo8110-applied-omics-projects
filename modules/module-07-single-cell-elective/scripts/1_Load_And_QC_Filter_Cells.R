###############################################################################
# MODULE 7 · SCRIPT 1 of 5 — LOAD AND QC-FILTER CELLS
# ==============================================================================
# ONE NEW IDEA on top of master_scrnaseq_pipeline.sh: that script's STARsolo
# "filtered" matrix already applied STARsolo's OWN cell-calling (real cell
# vs. empty droplet), but it has no concept of a LOW-QUALITY real cell -
# a droplet that did contain a cell, but one that was dying, ruptured, or
# otherwise not representative of real biology. This script is the first
# actual QC gate in this module's own pipeline, and it is deliberately
# separate from Script 2 (doublet detection) - a low-quality cell and a
# doublet are different failure modes with different fixes, and conflating
# them into one filter would make it impossible to tell which one is
# discarding a given cell and why.
#
# TEACHING NOTE — why mitochondrial percentage is the primary QC metric,
# not just a nice-to-have:
#   When a cell's membrane ruptures during dissociation, cytoplasmic mRNA
#   leaks out into the droplet buffer while mitochondria (larger, membrane-
#   bound organelles) are retained longer. A ruptured cell's captured RNA is
#   therefore disproportionately mitochondrial relative to a healthy cell's
#   - a high percent.mt is evidence of exactly this failure mode, not an
#   arbitrary threshold. nFeature_RNA (genes detected) and nCount_RNA (UMIs)
#   floors catch the complementary failure mode: a droplet with too little
#   captured material to say anything meaningful about that cell's biology,
#   whether or not it was ever a real cell.
#
# Input:  STARsolo's filtered matrix.mtx/features.tsv/barcodes.tsv directory
#         (master_scrnaseq_pipeline.sh's own output).
# Output: prints a before/after QC filtering report; writes qc_filtered.rds
#         (a Seurat object) and qc_summary.tsv.
###############################################################################

# Install once if needed (Bioconductor, not CRAN):
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("Seurat"))

library(Seurat)
library(ggplot2)   # needed for ggsave() on VlnPlot's ggplot object below

# ==============================================================================
# CONFIGURATION
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 1_Load_And_QC_Filter_Cells.R <path to Solo.out/Gene/filtered/> [output_dir]")
}
MATRIX_DIR <- args[1]
OUTDIR <- if (length(args) >= 2) args[2] else "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# QC thresholds. These are commonly-cited starting points for human tissue,
# NOT universal constants — always inspect this script's own QC plots for
# your actual data before trusting the defaults blindly (same caution
# amplicon_dada2_analysis.R gives for its own truncLen defaults).
MIN_FEATURES <- 200    # floor: droplets with fewer genes detected are likely empty/ambient, not low-quality cells
MIN_COUNTS <- 500       # floor: droplets with fewer UMIs are likely empty/ambient
MAX_PCT_MT <- 15        # ceiling: cells above this mitochondrial fraction are treated as ruptured/dying

cat("==================================================================\n")
cat(" MODULE 7 / SCRIPT 1 — Load and QC-filter cells\n")
cat("==================================================================\n")
cat(sprintf("Matrix directory: %s\n\n", MATRIX_DIR))

if (!dir.exists(MATRIX_DIR)) {
    stop(sprintf("Matrix directory not found: %s — did master_scrnaseq_pipeline.sh finish successfully?", MATRIX_DIR))
}

# ==============================================================================
# PHASE A: LOAD THE STARSOLO MATRIX
# ==============================================================================
counts <- Read10X(data.dir = MATRIX_DIR)
seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 0)
n_cells_raw <- ncol(seurat_obj)
cat(sprintf("Loaded %d cell barcodes x %d genes from STARsolo's filtered matrix.\n", n_cells_raw, nrow(seurat_obj)))

# ==============================================================================
# PHASE B: COMPUTE PER-CELL QC METRICS
# ==============================================================================
# "^MT-" matches human mitochondrial gene symbols (MT-ND1, MT-CO1, etc.) as
# annotated by Ensembl/GENCODE gene naming — confirm this pattern matches
# your reference's actual gene symbol convention (e.g. "^mt-" for mouse)
# before trusting percent.mt on non-human data.
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

qc_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(OUTDIR, "qc_metrics_violin.pdf"), plot = qc_plot, width = 10, height = 5)
cat(sprintf("Wrote %s — inspect before trusting the thresholds below on new data.\n",
            file.path(OUTDIR, "qc_metrics_violin.pdf")))

# ==============================================================================
# PHASE C: APPLY QC FILTERS
# ==============================================================================
pass_mask <- seurat_obj$nFeature_RNA >= MIN_FEATURES &
             seurat_obj$nCount_RNA >= MIN_COUNTS &
             seurat_obj$percent.mt <= MAX_PCT_MT

n_fail_features <- sum(seurat_obj$nFeature_RNA < MIN_FEATURES)
n_fail_counts <- sum(seurat_obj$nCount_RNA < MIN_COUNTS)
n_fail_mt <- sum(seurat_obj$percent.mt > MAX_PCT_MT)

cat(sprintf("\nQC filter results (thresholds: nFeature>=%d, nCount>=%d, percent.mt<=%d%%):\n",
            MIN_FEATURES, MIN_COUNTS, MAX_PCT_MT))
cat(sprintf("  failed nFeature_RNA floor : %d cells\n", n_fail_features))
cat(sprintf("  failed nCount_RNA floor   : %d cells\n", n_fail_counts))
cat(sprintf("  failed percent.mt ceiling : %d cells\n", n_fail_mt))

seurat_obj_filtered <- subset(seurat_obj, cells = colnames(seurat_obj)[pass_mask])
n_cells_pass <- ncol(seurat_obj_filtered)
pct_retained <- 100 * n_cells_pass / n_cells_raw

cat(sprintf("\nRetained %d/%d cells (%.1f%%) after QC filtering.\n", n_cells_pass, n_cells_raw, pct_retained))
if (pct_retained < 50) {
    cat("[WARN] fewer than half the input cells passed QC — check qc_metrics_violin.pdf before proceeding;\n")
    cat("       this can mean the sample itself was poor quality, or that the thresholds above don't fit it.\n")
} else {
    cat("[PASS] a majority of input cells passed QC.\n")
}

# ==============================================================================
# PHASE D: WRITE OUTPUT
# ==============================================================================
saveRDS(seurat_obj_filtered, file.path(OUTDIR, "qc_filtered.rds"))

summary_df <- data.frame(
    metric = c("cells_before_qc", "cells_after_qc", "pct_retained",
               "failed_min_features", "failed_min_counts", "failed_max_pct_mt"),
    value = c(n_cells_raw, n_cells_pass, round(pct_retained, 1),
              n_fail_features, n_fail_counts, n_fail_mt)
)
write.table(summary_df, file.path(OUTDIR, "qc_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\nWrote %s and %s\n", file.path(OUTDIR, "qc_filtered.rds"), file.path(OUTDIR, "qc_summary.tsv")))
cat("==================================================================\n")

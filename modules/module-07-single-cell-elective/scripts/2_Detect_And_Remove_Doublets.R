###############################################################################
# MODULE 7 · SCRIPT 2 of 5 — DETECT AND REMOVE DOUBLETS
# ==============================================================================
# ONE NEW IDEA on top of Script 1: Script 1 removed low-quality DROPLETS
# (too little material, too much mitochondrial signal). This script removes
# a structurally different problem that Script 1's per-cell QC metrics
# cannot catch on their own: a droplet that captured TWO cells at once. A
# doublet often looks like a perfectly healthy, high-count, low-mito cell —
# high nCount/nFeature is actually the expected signature of a doublet, not
# a red flag by Script 1's own filters — but its transcriptome is a mixture
# of two cells' biology, and left in place it can look like a novel,
# intermediate, or hybrid cell type in Script 3's clustering that does not
# correspond to any real cell.
#
# TEACHING NOTE — why this uses scDblFinder instead of a hand-rolled
# high-count heuristic:
#   A naive "flag the top N% by nCount" rule conflates doublets with real,
#   biologically large/transcriptionally active cells (which absolutely
#   exist and are not doublets). scDblFinder instead simulates realistic
#   doublets by summing pairs of this dataset's OWN real cells' UMI counts,
#   then trains a classifier to distinguish real cells from those synthetic
#   doublets by their k-nearest-neighbor expression profile, not by count
#   depth alone — the field-standard approach (Bioconductor), not a
#   heuristic reinvented here.
#
# TEACHING NOTE — this dataset's own ground truth, and why this script does
# not use it:
#   This module's worked-example cohort (BioProject PRJNA593571 — see
#   smoke_test_prjna593571.sh's header) is a deliberate species-mixing
#   ("barnyard") experiment: human PBMC + mouse colon + three cell lines in
#   one library. A droplet containing both human and mouse transcripts is an
#   unambiguous, ground-truth doublet — the classic experimental doublet-rate
#   validation method. This script does NOT use that ground truth, because
#   master_scrnaseq_pipeline.sh (Tier A, unchanged) aligns only against the
#   human reference by design; a droplet's mouse content is simply absent
#   from this pipeline's count matrix, not visible to check against.
#   Extending Tier A to a combined human+mouse reference specifically to
#   cross-validate scDblFinder's calls against this dataset's real barnyard
#   ground truth is a genuine, not-yet-built extension — see this module's
#   README "still to develop".
#
# Input:  qc_filtered.rds (Script 1).
# Output: prints a doublet-rate report; writes doublets_removed.rds (a
#         Seurat object) and doublet_summary.tsv.
###############################################################################

# Install once if needed (Bioconductor, not CRAN):
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("scDblFinder", "SingleCellExperiment"))

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 2_Detect_And_Remove_Doublets.R <qc_filtered.rds> [output_dir]")
}
QC_RDS_PATH <- args[1]
OUTDIR <- if (length(args) >= 2) args[2] else "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("==================================================================\n")
cat(" MODULE 7 / SCRIPT 2 — Doublet detection (scDblFinder)\n")
cat("==================================================================\n")

if (!file.exists(QC_RDS_PATH)) {
    stop(sprintf("QC-filtered object not found: %s — run Script 1 first.", QC_RDS_PATH))
}
seurat_obj <- readRDS(QC_RDS_PATH)
n_cells_in <- ncol(seurat_obj)
cat(sprintf("Loaded %d QC-passed cells from Script 1.\n\n", n_cells_in))

# ==============================================================================
# PHASE A: RUN scDblFinder
# ==============================================================================
set.seed(1)
sce <- as.SingleCellExperiment(seurat_obj)
sce <- scDblFinder(sce)

seurat_obj$scDblFinder.class <- sce$scDblFinder.class
seurat_obj$scDblFinder.score <- sce$scDblFinder.score

doublet_table <- table(seurat_obj$scDblFinder.class)
n_doublets <- if ("doublet" %in% names(doublet_table)) doublet_table[["doublet"]] else 0
pct_doublets <- 100 * n_doublets / n_cells_in

cat(sprintf("scDblFinder classified %d/%d cells (%.1f%%) as doublets.\n", n_doublets, n_cells_in, pct_doublets))
if (pct_doublets > 20) {
    cat("[WARN] doublet rate exceeds 20% — this is high even for a densely-loaded 10x run;\n")
    cat("       confirm the actual cell-loading concentration for this sample before treating\n")
    cat("       this number as expected rather than a sign something upstream went wrong\n")
    cat("       (e.g. ambient-RNA-heavy droplets misclassified as doublets).\n")
} else {
    cat("[PASS] doublet rate is within a typical range for standard 10x loading densities.\n")
}

# ==============================================================================
# PHASE B: REMOVE DOUBLETS
# ==============================================================================
seurat_obj_singlets <- subset(seurat_obj, subset = scDblFinder.class == "singlet")
n_cells_out <- ncol(seurat_obj_singlets)
cat(sprintf("\nRetained %d singlet cells for downstream clustering.\n", n_cells_out))

# ==============================================================================
# PHASE C: WRITE OUTPUT
# ==============================================================================
saveRDS(seurat_obj_singlets, file.path(OUTDIR, "doublets_removed.rds"))

summary_df <- data.frame(
    metric = c("cells_in", "cells_flagged_doublet", "pct_doublet", "cells_out"),
    value = c(n_cells_in, n_doublets, round(pct_doublets, 1), n_cells_out)
)
write.table(summary_df, file.path(OUTDIR, "doublet_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\nWrote %s and %s\n", file.path(OUTDIR, "doublets_removed.rds"), file.path(OUTDIR, "doublet_summary.tsv")))
cat("==================================================================\n")

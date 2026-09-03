###############################################################################
# MODULE 7 · SCRIPT 4 of 5 — IDENTIFY CLUSTER MARKERS
# ==============================================================================
# ONE NEW IDEA on top of Script 3: Script 3 produced clusters, which are
# purely a statement about transcriptional SIMILARITY — it makes no claim
# about what any cluster actually IS biologically. This script computes,
# for every cluster, which genes are differentially expressed in it versus
# every other cluster — the evidence a human would use to assign a cell-type
# label, not a label itself. This module deliberately does NOT assign
# cell-type labels automatically (see Script 5's Section 5) — a marker gene
# list is evidence for a judgment call, not a substitute for one.
#
# TEACHING NOTE — why min.pct and logfc.threshold matter as much as which
# test is used:
#   FindAllMarkers' pre-filters (only test genes expressed in at least
#   min.pct of cells in either group, only report genes above
#   logfc.threshold) exist because testing every one of ~20,000 genes for
#   every cluster is both slow and statistically wasteful — a gene expressed
#   in 2% of cells in both groups cannot produce a meaningful marker call
#   regardless of its p-value. These are not just performance knobs: setting
#   them too loose lets noisy, low-expression genes flood the marker table;
#   too strict and real but modest markers get filtered out before testing
#   even happens. The values below (Seurat's own documented defaults) are a
#   reasonable starting point, not a claim they suit every tissue.
#
# TEACHING NOTE — p_val_adj, not p_val, is the column to actually read:
#   FindAllMarkers tests every gene against every cluster, so the raw p_val
#   column is uncorrected for the resulting multiple-comparisons problem.
#   p_val_adj (Bonferroni-corrected across all genes tested) is what should
#   actually be used to call a gene significant — this script sorts and
#   reports by p_val_adj for exactly that reason, not p_val.
#
# Input:  clustered.rds (Script 3).
# Output: prints a per-cluster top-marker report; writes all_markers.tsv
#         (every significant marker, every cluster) and
#         top_markers_per_cluster.tsv (top 10 by p_val_adj, per cluster).
###############################################################################

library(Seurat)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 4_Identify_Cluster_Markers.R <clustered.rds> [output_dir]")
}
INPUT_RDS_PATH <- args[1]
OUTDIR <- if (length(args) >= 2) args[2] else "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MIN_PCT <- 0.25          # only test genes expressed in >=25% of cells in either group
LOGFC_THRESHOLD <- 0.25  # only report genes with at least this log-fold-change
TOP_N_PER_CLUSTER <- 10
PADJ_THRESHOLD <- 0.05

cat("==================================================================\n")
cat(" MODULE 7 / SCRIPT 4 — Cluster marker identification\n")
cat("==================================================================\n")

if (!file.exists(INPUT_RDS_PATH)) {
    stop(sprintf("Clustered object not found: %s — run Script 3 first.", INPUT_RDS_PATH))
}
seurat_obj <- readRDS(INPUT_RDS_PATH)
n_clusters <- length(unique(Idents(seurat_obj)))
cat(sprintf("Loaded %d cells across %d clusters from Script 3.\n\n", ncol(seurat_obj), n_clusters))

# ==============================================================================
# PHASE A: FIND MARKERS FOR EVERY CLUSTER
# ==============================================================================
all_markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,   # report genes UP in each cluster vs. the rest — the standard "what defines this cluster" question
    min.pct = MIN_PCT,
    logfc.threshold = LOGFC_THRESHOLD,
    verbose = FALSE)

all_markers_sig <- all_markers[all_markers$p_val_adj < PADJ_THRESHOLD, ]
cat(sprintf("%d/%d candidate marker calls pass p_val_adj < %.2f across all clusters.\n",
            nrow(all_markers_sig), nrow(all_markers), PADJ_THRESHOLD))

n_clusters_no_markers <- sum(!levels(Idents(seurat_obj)) %in% unique(all_markers_sig$cluster))
if (n_clusters_no_markers > 0) {
    cat(sprintf("[WARN] %d cluster(s) have NO significant marker gene at these thresholds —\n", n_clusters_no_markers))
    cat("       this usually means that cluster is not well-separated from its neighbors at the\n")
    cat("       current resolution (Script 3), not that the cluster has no real biology; consider\n")
    cat("       re-running Script 3 at a lower resolution if this happens for several clusters.\n")
} else {
    cat("[PASS] every cluster has at least one significant marker gene.\n")
}

# ==============================================================================
# PHASE B: TOP N MARKERS PER CLUSTER
# ==============================================================================
if (nrow(all_markers_sig) == 0) {
    cat("\nNo cluster has any significant marker gene at these thresholds — a legitimate\n")
    cat("(if unusual) outcome, not a script error. top_markers_per_cluster.tsv will be empty.\n")
    top_markers <- all_markers_sig[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
} else {
    top_markers <- do.call(rbind, lapply(split(all_markers_sig, all_markers_sig$cluster), function(df) {
        df[order(df$p_val_adj, -df$avg_log2FC), ][seq_len(min(TOP_N_PER_CLUSTER, nrow(df))), , drop = FALSE]
    }))
    cat(sprintf("\nTop %d marker(s) per cluster:\n", TOP_N_PER_CLUSTER))
    print(top_markers[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])
}

# ==============================================================================
# PHASE C: WRITE OUTPUT
# ==============================================================================
write.table(all_markers_sig, file.path(OUTDIR, "all_markers.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(top_markers, file.path(OUTDIR, "top_markers_per_cluster.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\nWrote %s and %s\n",
            file.path(OUTDIR, "all_markers.tsv"),
            file.path(OUTDIR, "top_markers_per_cluster.tsv")))
cat("Cell-type identity is NOT assigned here — cross-reference top_markers_per_cluster.tsv\n")
cat("against known marker genes for the expected tissue by hand (Script 5's Section 5).\n")
cat("==================================================================\n")

###############################################################################
# MODULE 7 · SCRIPT 3 of 5 — CLUSTER CELLS
# ==============================================================================
# ONE NEW IDEA on top of Scripts 1-2: those scripts decided WHICH cells to
# trust. This script is the first one that looks at what those cells' gene
# expression actually says — grouping cells by transcriptional similarity
# into clusters, the input every downstream biological question (Script 4's
# marker genes, and the interpretation this module's report leaves to the
# analyst) depends on.
#
# TEACHING NOTE — why PCA runs on SCALED data, and why clustering runs on
# PCA space, not the raw gene matrix directly:
#   ScaleData z-scores each gene across cells so that highly-expressed genes
#   (which have naturally larger raw variance) do not dominate PCA purely by
#   expression magnitude rather than by how informative their variation
#   actually is across cells. FindVariableFeatures narrows PCA's input
#   further to genes that vary meaningfully across THIS dataset specifically
#   — most genes are uninformative "housekeeping" noise for distinguishing
#   cell types. Clustering (FindNeighbors + FindClusters, a graph-based
#   Louvain/Leiden method) then runs on a handful of PCA dimensions, not the
#   full ~20,000-gene expression matrix directly — PCA denoises and
#   compresses the signal that matters into far fewer dimensions than genes,
#   which is what makes a k-nearest-neighbor graph over cells tractable and
#   less dominated by single-gene noise.
#
# TEACHING NOTE — resolution is a knob, not a fact about the data:
#   FindClusters' `resolution` parameter directly trades off cluster COUNT
#   against cluster GRANULARITY — a higher resolution finds more, finer
#   clusters; a lower one finds fewer, coarser ones. There is no single
#   "correct" resolution independent of the biological question being asked
#   (e.g. "immune cell types" vs. "T cell subtypes" are different, both
#   valid, answers at different resolutions on the SAME data) — the value
#   below is a reasonable default for a broad first pass, not a claim that
#   it is the right number of clusters for this tissue.
#
# Input:  doublets_removed.rds (Script 2).
# Output: prints a clustering summary; writes clustered.rds (a Seurat
#         object with cluster identities), umap_by_cluster.pdf, and
#         cluster_sizes.tsv.
###############################################################################

library(Seurat)
library(ggplot2)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 3_Cluster_Cells.R <doublets_removed.rds> [output_dir]")
}
INPUT_RDS_PATH <- args[1]
OUTDIR <- if (length(args) >= 2) args[2] else "."
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

N_PCS <- 20          # number of PCA dimensions carried into neighbor-graph construction and UMAP
CLUSTER_RESOLUTION <- 0.8   # standard default first-pass resolution (see TEACHING NOTE above)

cat("==================================================================\n")
cat(" MODULE 7 / SCRIPT 3 — Clustering\n")
cat("==================================================================\n")

if (!file.exists(INPUT_RDS_PATH)) {
    stop(sprintf("Input object not found: %s — run Script 2 first.", INPUT_RDS_PATH))
}
seurat_obj <- readRDS(INPUT_RDS_PATH)
cat(sprintf("Loaded %d cells x %d genes from Script 2.\n\n", ncol(seurat_obj), nrow(seurat_obj)))

# ==============================================================================
# PHASE A: NORMALIZE, FIND VARIABLE FEATURES, SCALE
# ==============================================================================
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# ==============================================================================
# PHASE B: PCA
# ==============================================================================
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = N_PCS, verbose = FALSE)
cat(sprintf("Ran PCA on %d variable genes, retaining %d principal components.\n", length(VariableFeatures(seurat_obj)), N_PCS))

# ==============================================================================
# PHASE C: CLUSTER (GRAPH-BASED, LOUVAIN)
# ==============================================================================
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:N_PCS, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = CLUSTER_RESOLUTION, verbose = FALSE)

n_clusters <- length(unique(Idents(seurat_obj)))
cat(sprintf("Found %d clusters at resolution=%.1f.\n", n_clusters, CLUSTER_RESOLUTION))

cluster_sizes <- table(Idents(seurat_obj))
cat("Cluster sizes:\n")
print(cluster_sizes)

n_singleton_clusters <- sum(cluster_sizes < 10)
if (n_singleton_clusters > 0) {
    cat(sprintf("[WARN] %d cluster(s) have fewer than 10 cells — check whether these are real rare\n", n_singleton_clusters))
    cat("       populations or an artifact of over-clustering at this resolution before trusting\n")
    cat("       Script 4's marker genes for them; a marker list from 6 cells is far less reliable\n")
    cat("       than one from 600.\n")
} else {
    cat("[PASS] no cluster has fewer than 10 cells.\n")
}

# ==============================================================================
# PHASE D: UMAP (VISUALIZATION ONLY — NOT USED FOR CLUSTERING ITSELF)
# ==============================================================================
# UMAP is computed here purely to visualize the clusters found above in 2D;
# clustering itself (Phase C) already happened in PCA space, not UMAP space
# — UMAP coordinates are a projection for human eyes, not the basis for
# any of this script's own decisions.
seurat_obj <- RunUMAP(seurat_obj, dims = 1:N_PCS, verbose = FALSE)
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggtitle(sprintf("%d clusters, resolution=%.1f", n_clusters, CLUSTER_RESOLUTION))
ggsave(file.path(OUTDIR, "umap_by_cluster.pdf"), plot = umap_plot, width = 8, height = 6)

# ==============================================================================
# PHASE E: WRITE OUTPUT
# ==============================================================================
saveRDS(seurat_obj, file.path(OUTDIR, "clustered.rds"))

cluster_df <- as.data.frame(cluster_sizes)
colnames(cluster_df) <- c("cluster", "n_cells")
write.table(cluster_df, file.path(OUTDIR, "cluster_sizes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\nWrote %s, %s, and %s\n",
            file.path(OUTDIR, "clustered.rds"),
            file.path(OUTDIR, "umap_by_cluster.pdf"),
            file.path(OUTDIR, "cluster_sizes.tsv")))
cat("==================================================================\n")

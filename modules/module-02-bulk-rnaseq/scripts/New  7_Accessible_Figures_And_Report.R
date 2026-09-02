# ==============================================================================
# MODULE 2 · SCRIPT 7 of 7 — ACCESSIBLE FIGURES + THE PORTFOLIO REPORT
# ==============================================================================
# ONE NEW IDEA on top of Script 6: Script 6 produced the numbers. This script
# turns those numbers into the curriculum's required portfolio artifact - "a
# reproducible report containing metadata validation, exploratory plots,
# model specification, differential-expression results, effect-size
# visualization, interpretation, and limitations" - and makes the figures
# themselves accessible, not just correct.
#
# TEACHING NOTE — accessible figure design (mirrors Module 1 Script 7's
# report-generation habit, one level up: figures, not just text):
#   - Never use red/green alone to mean up/down - ~8% of men have red-green
#     color vision deficiency. This script uses a colorblind-safe blue/orange
#     diverging scheme (matches the Nexus data-viz palette) instead.
#   - Label the few genes that matter directly on the plot rather than
#     forcing the reader to cross-reference a separate legend or table.
#   - State the finding in the title ("14 genes pass FDR < 0.05"), not a
#     generic label like "Volcano Plot" - a chart should carry its own
#     conclusion.
#
# Requires: DESeq2 (for results class), ggplot2, ggrepel (optional, falls
#           back to geom_text if unavailable)
# Input:    deseq2_results_shrunk.csv (from Script 6)
# Output:   volcano_plot.pdf, ma_plot.pdf, top_genes.csv, DE_REPORT.md
# ==============================================================================

suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
results_path <- if (length(args) >= 1) args[1] else "deseq2_results_shrunk.csv"
alpha        <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
top_n        <- if (length(args) >= 3) as.integer(args[3]) else 10

cat("==================================================================\n")
cat(" MODULE 2 / SCRIPT 7 — Accessible figures + portfolio report\n")
cat("==================================================================\n")

res <- read.csv(results_path, row.names = 1)
res$gene <- rownames(res)
res$sig  <- ifelse(!is.na(res$padj) & res$padj < alpha, "FDR < threshold", "not significant")

n_sig <- sum(res$sig == "FDR < threshold")
cat(sprintf("Loaded %d genes; %d pass FDR < %.2f\n", nrow(res), n_sig, alpha))

# --- 1. Volcano plot - colorblind-safe blue/orange, top genes labeled -------
top_hits <- res[order(res$padj), ][seq_len(min(top_n, nrow(res))), ]

volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("FDR < threshold" = "#DA7101", "not significant" = "#5591C7")) +
  geom_text(data = top_hits, aes(label = gene), size = 2.5, vjust = -0.6,
            color = "black", check_overlap = TRUE) +
  labs(
    title = sprintf("%d of %d genes pass FDR < %.2f", n_sig, nrow(res), alpha),
    subtitle = "Color encodes FDR-adjusted significance, not raw p-value",
    x = "log2 fold change (shrunk)", y = expression(-log[10](p~value)), color = NULL
  ) +
  theme_bw()
ggsave("volcano_plot.pdf", plot = volcano, width = 7, height = 5.5)
cat("Wrote volcano_plot.pdf\n")

# --- 2. MA plot - same color scheme for consistency -------------------------
ma_plot <- ggplot(res, aes(x = baseMean, y = log2FoldChange, color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_x_log10() +
  scale_color_manual(values = c("FDR < threshold" = "#DA7101", "not significant" = "#5591C7")) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
  labs(title = "MA plot: effect size vs. mean expression",
       x = "Mean of normalized counts (log10)", y = "log2 fold change (shrunk)", color = NULL) +
  theme_bw()
ggsave("ma_plot.pdf", plot = ma_plot, width = 7, height = 5.5)
cat("Wrote ma_plot.pdf\n")

# --- 3. Top-gene table -------------------------------------------------------
write.csv(top_hits[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")],
          "top_genes.csv", row.names = FALSE)
cat("Wrote top_genes.csv (", nrow(top_hits), " genes)\n", sep = "")

# --- 4. Auto-generate the portfolio report -----------------------------------
report <- c(
  "# Differential Expression Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  "## 1. Metadata validation",
  "- See Script 1 console output for the full PASS/WARN/FAIL log (count matrix",
  "  integrity, sample-name matching, library-size checks).",
  "",
  "## 2. Design and confounding check",
  "- See Script 2 console output. Paste the replicate-count table and the",
  "  condition-vs-batch cross-tabulation here, and note the PASS/FAIL verdict.",
  "",
  "## 3. Gene filtering rule applied",
  "- Documented rule from Script 3: keep genes with >= 10 reads in at least",
  "  as many samples as the smallest condition group (not a naive",
  "  rowSums-only threshold). Paste the before/after gene counts here.",
  "",
  "## 4. Exploratory analysis",
  "- See pca_plot.pdf and sample_distance_heatmap.pdf (Script 4). Note any",
  "  flagged outliers and whether clustering matches the expected condition",
  "  grouping.",
  "",
  "## 5. Model specification",
  "- See Script 5 console output for the exact design formula, reference",
  "  level, and contrast vector used below.",
  "",
  "## 6. Differential expression results",
  sprintf("- %d of %d testable genes pass FDR (padj) < %.2f.", n_sig, nrow(res), alpha),
  "- Ranking and fold-change values use the apeglm-SHRUNK log2FoldChange",
  "  (deseq2_results_shrunk.csv), not the raw DESeq2 output, to avoid",
  "  overstating effect sizes for low-count genes.",
  "- See top_genes.csv for the full top-hit table.",
  "",
  "## 7. Effect-size visualization",
  "- volcano_plot.pdf and ma_plot.pdf (this script). Both use a",
  "  colorblind-safe blue/orange scheme and label top hits directly.",
  "",
  "## 8. Interpretation (fill in before submitting)",
  "- What do the top genes do biologically? Do they form a coherent",
  "  pathway/process, or look like scattered noise?",
  "- Does the direction of change (up/down) make biological sense given the",
  "  experimental condition?",
  "",
  "## 9. Limitations (fill in before submitting)",
  "- State the actual replicate count per group and whether it meets the",
  "  curriculum's >= 3 recommendation (Script 2).",
  "- Note any flagged PCA outliers (Script 4) and whether they were",
  "  excluded, and why.",
  "- Note any significant-but-trivial-effect genes flagged by Script 6",
  "  (|log2FC| < 0.5) that were excluded from the biological narrative.",
  ""
)
writeLines(report, "DE_REPORT.md")
cat("Wrote DE_REPORT.md (fill in Sections 8-9 by hand before submitting as the portfolio artifact)\n")

cat("\n==================================================================\n")
cat(" PIPELINE COMPLETE (Scripts 1-7). DE_REPORT.md is the curriculum's\n")
cat(" required portfolio artifact for Module 2 once Sections 8-9 are filled in.\n")
cat("==================================================================\n")

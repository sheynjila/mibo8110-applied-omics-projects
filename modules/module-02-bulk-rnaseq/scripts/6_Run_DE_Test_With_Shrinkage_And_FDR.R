# ==============================================================================
# MODULE 2 · SCRIPT 6 of 7 — RUN THE DE TEST, APPLY SHRINKAGE, INTERPRET FDR
# ==============================================================================
# ONE NEW IDEA on top of Script 5: Script 5 finalized the design and named the
# contrast. This script actually fits the model, and - the curriculum's
# "Interpret effect sizes, uncertainty, p-values, and false-discovery-rate
# adjustment" learning outcome - makes each of those four things explicit
# and separately checked, rather than reading only the p-value column.
#
# TEACHING NOTE — the four numbers in a DESeq2 results table, and why each
# matters on its own:
#   1. log2FoldChange (EFFECT SIZE): the magnitude of change. A gene can be
#      "significant" with a trivial 1.05-fold change if sample size is large
#      - statistical significance is not the same as biological importance.
#   2. lfcSE (UNCERTAINTY): the standard error on that fold change. A large
#      log2FoldChange with a huge lfcSE is not a reliable estimate, however
#      small its p-value.
#   3. pvalue (RAW SIGNIFICANCE): the probability of seeing this result by
#      chance for ONE gene. Never use this directly for gene-level calls -
#      testing thousands of genes at p<0.05 guarantees thousands of false
#      positives by chance alone.
#   4. padj (FDR-ADJUSTED, Benjamini-Hochberg by default in DESeq2): the
#      number that actually controls the false discovery rate across all
#      genes tested. THIS is the column to threshold on - conventionally
#      padj < 0.05 - never raw pvalue.
#
# lfcShrink() (apeglm method) additionally shrinks log2FoldChange estimates
# for low-count / high-uncertainty genes toward zero, which stabilizes
# ranking and visualization without changing which genes are called
# significant (that still comes from the Wald test's padj).
#
# Requires: DESeq2, apeglm
# Input:    dds_prepared.rds (from Script 5)
# Output:   deseq2_results.csv (raw), deseq2_results_shrunk.csv (shrunk LFC);
#           prints a plain-language summary of how many genes passed FDR.
# ==============================================================================

suppressMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
prepared_path <- if (length(args) >= 1) args[1] else "dds_prepared.rds"
alpha         <- if (length(args) >= 2) as.numeric(args[2]) else 0.05

cat("==================================================================\n")
cat(" MODULE 2 / SCRIPT 6 — Differential expression test + FDR interpretation\n")
cat("==================================================================\n")

prepared <- readRDS(prepared_path)
dds      <- prepared$dds
contrast <- prepared$contrast
cat("Loaded design:", deparse(design(dds)), "\n")
cat("Contrast:", paste(contrast, collapse = " / "), "\n")
cat("FDR threshold (alpha):", alpha, "\n\n")

# --- 1. Fit the model ----------------------------------------------------
dds <- DESeq(dds)

# --- 2. Extract results using the EXPLICIT contrast, and the SAME alpha ----
# Passing alpha here matters: DESeq2's internal default filtering threshold
# for the "optimal" p-value cutoff should match the alpha you will actually
# report against, or the printed summary and your own padj<alpha count can
# disagree.
res <- results(dds, contrast = contrast, alpha = alpha)

cat("DESeq2 built-in summary:\n")
summary(res)

# --- 3. Shrink log2FoldChange for stable ranking/visualization -------------
# coef name must match the fitted model's coefficient - for a two-level
# factor this is "<column>_<comparison>_vs_<reference>".
coef_name <- paste0(contrast[1], "_", contrast[2], "_vs_", contrast[3])
res_shrunk <- tryCatch({
  suppressMessages(library(apeglm))
  lfcShrink(dds, coef = coef_name, type = "apeglm")
}, error = function(e) {
  cat("[WARN] apeglm shrinkage unavailable (", conditionMessage(e), "). Falling back to normal shrinkage.\n", sep = "")
  lfcShrink(dds, contrast = contrast, type = "normal")
})

# --- 4. Explicit FDR-based gene call, using padj (never raw pvalue) ---------
n_total <- sum(!is.na(res$padj))
n_sig   <- sum(res$padj < alpha, na.rm = TRUE)
n_sig_raw_p_only <- sum(res$pvalue < alpha, na.rm = TRUE)

cat(sprintf(
  "\n%d of %d testable genes (%.1f%%) pass padj < %.2f - THIS is the number to report.\n",
  n_sig, n_total, 100 * n_sig / n_total, alpha))
if (n_sig > 0) {
  cat(sprintf(
    "%d genes have raw pvalue < %.2f (uncorrected - do NOT report this number; it overstates significance by %.1fx).\n",
    n_sig_raw_p_only, alpha, n_sig_raw_p_only / n_sig))
} else {
  cat(sprintf(
    "%d genes have raw pvalue < %.2f (uncorrected - do NOT report this number; it is not FDR-controlled).\n",
    n_sig_raw_p_only, alpha))
}

# --- 5. Flag "significant but trivial effect" genes for the interpretation step
sig_idx <- which(res$padj < alpha & !is.na(res$padj))
trivial_effect <- sum(abs(res$log2FoldChange[sig_idx]) < 0.5)
if (trivial_effect > 0) {
  cat(sprintf(
    "\n[WARN] %d of the %d significant genes have |log2FoldChange| < 0.5 (less than\n",
    trivial_effect, length(sig_idx)))
  cat("a 1.4-fold change) - statistically significant but a small effect. Flag\n")
  cat("these for the interpretation step (Script 7) rather than presenting them\n")
  cat("as top hits without qualification.\n")
}

write.csv(as.data.frame(res), "deseq2_results.csv")
write.csv(as.data.frame(res_shrunk), "deseq2_results_shrunk.csv")
cat("\nWrote deseq2_results.csv (raw LFC + padj) and deseq2_results_shrunk.csv (shrunk LFC).\n")

cat("\n==================================================================\n")
cat(" DE TEST COMPLETE - use deseq2_results_shrunk.csv for ranking/plots and\n")
cat(" the padj column (not pvalue) for any significance claim in Script 7.\n")
cat("==================================================================\n")

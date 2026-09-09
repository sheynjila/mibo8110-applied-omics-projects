# ==============================================================================
# MODULE 2 · SCRIPT 5 of 7 — SPECIFY THE DESIGN FORMULA AND CONTRASTS EXPLICITLY
# ==============================================================================
# ONE NEW IDEA on top of Script 4: Script 4 looked at the data. This script
# turns that understanding into an explicit, written-down statistical design
# - the single most consequential decision in the whole analysis, and the
# curriculum's "Specify designs and contrasts in DESeq2" learning outcome.
#
# TEACHING NOTE — two silent failure modes this script exists to prevent:
#   1. DEFAULT REFERENCE LEVEL: by default, R factors a character column
#      alphabetically, so "treatment" vs "control" becomes the reference
#      level purely because "c" < "t". If your labels were instead "KO" and
#      "WT", DESeq2 would treat "KO" as the reference and report every log2
#      fold change with the SIGN FLIPPED relative to what you probably
#      intended (WT vs KO). relevel() makes the reference level an explicit,
#      reviewable choice instead of an alphabetical accident.
#   2. IMPLICIT VS NAMED CONTRASTS: results(dds) with no arguments returns
#      the LAST coefficient in the design formula by default. In a design
#      with a batch term (~ batch + condition), that is usually what you
#      want - but in a design with 3+ condition levels, or a designer who
#      later adds a covariate, "the last coefficient" silently becomes the
#      wrong comparison. Always pass an explicit contrast=c("condition",
#      "treatment", "control") so the comparison being tested is visible in
#      the code, not implied by column order.
#
# Requires: DESeq2
# Input:    counts_filtered.tsv, metadata.csv
# Output:   prints the design formula, reference level, and the exact contrast
#           vector that Script 6 will use; saves dds_prepared.rds for Script 6.
# ==============================================================================

suppressMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
counts_path     <- if (length(args) >= 1) args[1] else "counts_filtered.tsv"
metadata_path   <- if (length(args) >= 2) args[2] else "metadata.csv"
reference_level <- if (length(args) >= 3) args[3] else NULL  # e.g. "control"

cat("==================================================================\n")
cat(" MODULE 2 / SCRIPT 5 — Design formula and contrast specification\n")
cat("==================================================================\n")

counts <- read.table(counts_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta   <- read.csv(metadata_path, stringsAsFactors = FALSE)
meta   <- meta[match(colnames(counts), meta$sample), ]
rownames(meta) <- meta$sample

# --- 1. Choose the design formula based on what Script 2 found ---------------
has_batch <- "batch" %in% colnames(meta)
design_formula <- if (has_batch) ~ batch + condition else ~ condition
cat("Design formula:", deparse(design_formula), "\n")
if (has_batch) {
  cat("  'batch' is included as a covariate (Script 2 confirmed it is not\n")
  cat("  perfectly confounded with condition) to absorb technical variance\n")
  cat("  before estimating the condition effect. 'condition' is listed LAST\n")
  cat("  so it is the term results() reports on by default - but this script\n")
  cat("  still names the contrast explicitly (see below) rather than relying\n")
  cat("  on that default.\n")
}

# --- 2. Set an EXPLICIT reference level, never rely on alphabetical order ----
meta$condition <- factor(meta$condition)
levels_found <- levels(meta$condition)
if (is.null(reference_level)) {
  reference_level <- levels_found[1]
  cat(sprintf(
    "\n[WARN] no reference level was specified - defaulting to '%s' because it\n",
    reference_level))
  cat("       sorts first alphabetically. Pass the intended control/baseline\n")
  cat("       label explicitly as the 3rd script argument in real use.\n")
} else {
  if (!(reference_level %in% levels_found)) {
    stop(sprintf("requested reference level '%s' is not among the condition levels found: %s",
                 reference_level, paste(levels_found, collapse = ", ")), call. = FALSE)
  }
  cat(sprintf("\nReference level explicitly set to '%s'.\n", reference_level))
}
meta$condition <- relevel(meta$condition, ref = reference_level)
cat("Condition factor levels (reference first):", paste(levels(meta$condition), collapse = ", "), "\n")

# --- 3. Build the DESeqDataSet with the finalized design ----------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)

# --- 4. Name the contrast EXPLICITLY - never rely on default results() -------
non_reference_levels <- setdiff(levels(meta$condition), reference_level)
if (length(non_reference_levels) > 1) {
  cat(sprintf(
    "\n%d non-reference condition levels found (%s). This script demonstrates\n",
    length(non_reference_levels), paste(non_reference_levels, collapse = ", ")))
  cat("the FIRST one below; repeat with a different contrast vector for each\n")
  cat("additional group-vs-reference comparison you need.\n")
}
comparison_level <- non_reference_levels[1]
contrast_vec <- c("condition", comparison_level, reference_level)
cat(sprintf(
  "\nExplicit contrast for Script 6: contrast = c(%s)\n",
  paste(sprintf('"%s"', contrast_vec), collapse = ", ")))
cat(sprintf("  -> tests: %s vs %s (positive log2FoldChange = higher in %s)\n",
            comparison_level, reference_level, comparison_level))

saveRDS(list(dds = dds, contrast = contrast_vec), "dds_prepared.rds")
cat("\nWrote dds_prepared.rds (contains the un-fit DESeqDataSet + contrast vector)\n")

cat("\n==================================================================\n")
cat(" DESIGN SPECIFICATION COMPLETE - the printed design formula and contrast\n")
cat(" vector both belong verbatim in the portfolio report (Script 7).\n")
cat("==================================================================\n")

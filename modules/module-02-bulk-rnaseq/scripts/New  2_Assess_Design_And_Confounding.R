# ==============================================================================
# MODULE 2 · SCRIPT 2 of 7 — ASSESS EXPERIMENTAL DESIGN AND CONFOUNDING
# ==============================================================================
# ONE NEW IDEA on top of Script 1: Script 1 confirmed the data are readable
# and aligned. This script asks a harder question BEFORE any statistics run:
# is the experimental design even capable of answering the biological
# question, or is a technical variable confounded with the condition of
# interest such that no statistical method can separate them?
#
# TEACHING NOTE — replication vs. confounding (curriculum learning outcome:
# "Explain biological replication, confounding, covariates, and batch
# effects"):
#   - Biological replicates = independent animals/cultures/patients. Technical
#     replicates (re-sequencing the same library) are NOT a substitute for
#     biological replicates and must never be pooled as if they were.
#   - A design is CONFOUNDED when a covariate (batch, sequencing lane,
#     extraction date, sex, etc.) perfectly tracks the condition of interest
#     - e.g., every control was sequenced in batch 1 and every treatment in
#     batch 2. No amount of modeling can then distinguish "treatment effect"
#     from "batch effect" - the design itself has to change (re-randomize,
#     add samples), not the statistics.
#
# Input:  metadata.csv - sample,condition[,batch,sex,...]
# Output: a cross-tabulation of condition against every other metadata column,
#         and an explicit PASS/WARN/FAIL confounding verdict per column.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
metadata_path <- if (length(args) >= 1) args[1] else "metadata.csv"

cat("==================================================================\n")
cat(" MODULE 2 / SCRIPT 2 — Design and confounding assessment\n")
cat("==================================================================\n")

meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
if (!"condition" %in% colnames(meta)) stop("metadata.csv must contain a 'condition' column", call. = FALSE)

# --- 1. Replicate count per condition -----------------------------------------
rep_counts <- table(meta$condition)
cat("\nReplicates per condition:\n")
print(rep_counts)

min_n <- min(rep_counts)
if (min_n < 2) {
  cat("\n[FAIL-STYLE WARNING] At least one condition has only", min_n,
      "replicate(s). DESeq2 will RUN with n=1 per group, but dispersion\n")
  cat("cannot be estimated from real biological variance, so p-values are not\n")
  cat("trustworthy. Treat any results from this design as exploratory only -\n")
  cat("do not report them as statistically significant findings.\n")
} else if (min_n < 3) {
  cat("\n[WARN] Minimum replication is", min_n, "per group. This RUNS, but the\n")
  cat("curriculum recommends >= 3 biological replicates per condition for\n")
  cat("reliable dispersion estimates. Treat effect sizes as provisional.\n")
} else {
  cat("\n[PASS] Minimum replication is", min_n, "per group (>= 3) - adequate for\n")
  cat("stable dispersion estimation.\n")
}

# --- 2. Confounding check: cross-tabulate condition against every other column
other_cols <- setdiff(colnames(meta), c("sample", "condition"))

if (length(other_cols) == 0) {
  cat("\nNo additional metadata columns (batch, sex, etc.) were provided.\n")
  cat("Design formula will be the simple case: ~ condition (see Script 5).\n")
} else {
  cat("\nChecking each additional column against 'condition' for confounding:\n")
  for (col in other_cols) {
    cat("\n---", col, "vs condition ---\n")
    tab <- table(meta$condition, meta[[col]])
    print(tab)

    # A column is PERFECTLY confounded with condition if every condition group
    # maps to exactly one level of that column (i.e., each row of the table
    # has exactly one non-zero cell).
    nonzero_per_row <- apply(tab, 1, function(r) sum(r > 0))
    if (all(nonzero_per_row == 1)) {
      cat(sprintf(
        "[FAIL] '%s' is PERFECTLY CONFOUNDED with condition - every condition group falls in exactly one '%s' level.\n",
        col, col))
      cat("        This cannot be fixed by adding '", col, "' to the design formula.\n", sep = "")
      cat("        The treatment effect and the '", col, "' effect are mathematically\n", sep = "")
      cat("        indistinguishable. The only fix is a better-randomized design\n")
      cat("        (re-sequence with conditions balanced across '", col, "').\n", sep = "")
    } else {
      cat(sprintf(
        "[PASS] '%s' is NOT perfectly confounded with condition - it can be added\n",
        col))
      cat("        to the design formula as a covariate (see Script 5) to remove its\n")
      cat("        effect and reduce residual variance.\n")
    }
  }
}

cat("\n==================================================================\n")
cat(" DESIGN ASSESSMENT COMPLETE - resolve any [FAIL] above before Script 3.\n")
cat(" PASS/WARN findings are informational; they do not block the pipeline.\n")
cat("==================================================================\n")

# ==============================================================================
# MODULE 2 · SCRIPT 3 of 7 — FILTER LOW-INFORMATION GENES WITH A DOCUMENTED RULE
# ==============================================================================
# ONE NEW IDEA on top of Script 2: Script 2 confirmed the design can support a
# statistical test. This script decides WHICH GENES are even worth testing -
# and, critically, forces that decision to be explicit and documented instead
# of an unexamined magic number.
#
# TEACHING NOTE — why "rowSums(counts) >= 10" (the rule used in the original
# rnaseq_deseq2_analysis.R) is not good enough on its own:
#   A gene with 40 reads in ONE sample and 0 in every other sample passes
#   "total >= 10" but is not "detected" in any meaningful biological sense -
#   it is noise or a single-sample artifact. A better rule ties the count
#   threshold to GROUP SIZE: require a minimum count in AT LEAST as many
#   samples as the smallest group, so a gene must be plausibly detectable
#   within at least one whole condition group, not just scattered across
#   unrelated samples.
#
# This script does NOT replace DESeq2's own independent filtering (which
# happens later, on p-values, inside results()). This is a separate,
# earlier, count-based filter — the curriculum's "filter low-information
# features using a documented rule" learning outcome.
#
# Input:  counts.txt, metadata.csv (same as Scripts 1-2)
# Output: prints before/after gene counts and the exact rule used; writes
#         counts_filtered.tsv for use by Script 4 onward.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
counts_path   <- if (length(args) >= 1) args[1] else "counts.txt"
metadata_path <- if (length(args) >= 2) args[2] else "metadata.csv"
min_count     <- if (length(args) >= 3) as.numeric(args[3]) else 10

cat("==================================================================\n")
cat(" MODULE 2 / SCRIPT 3 — Low-count gene filtering\n")
cat("==================================================================\n")

raw <- read.table(counts_path, header = TRUE, row.names = 1,
                   sep = "\t", check.names = FALSE, comment.char = "#")
known_annotation_cols <- c("Chr", "Start", "End", "Strand", "Length")
counts <- raw[, !(colnames(raw) %in% known_annotation_cols), drop = FALSE]

meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
meta <- meta[match(colnames(counts), meta$sample), ]

# --- The documented rule ------------------------------------------------------
# Keep a gene if it has >= min_count reads in at least as many samples as the
# smallest condition group. Rationale: a real biological signal should be
# detectable across (at least) one whole group's worth of replicates, not
# just scattered single samples.
group_sizes  <- table(meta$condition)
min_group_n  <- min(group_sizes)

cat(sprintf(
  "\nRULE: keep genes with >= %d reads in at least %d sample(s)\n",
  min_count, min_group_n))
cat(sprintf(
  "(%d = size of the smallest condition group: %s)\n\n",
  min_group_n, paste(names(group_sizes)[group_sizes == min_group_n], collapse = ", ")))

samples_passing <- rowSums(counts >= min_count)
keep_documented <- samples_passing >= min_group_n

# --- For comparison: the naive rule from the original script ------------------
keep_naive <- rowSums(counts) >= min_count

cat("Comparison of the two rules:\n")
cat(sprintf("  Naive rule       (rowSums(counts) >= %d):        %d / %d genes kept\n",
            min_count, sum(keep_naive), nrow(counts)))
cat(sprintf("  Documented rule  (>= %d reads in >= %d samples): %d / %d genes kept\n",
            min_count, min_group_n, sum(keep_documented), nrow(counts)))

only_naive_keeps <- sum(keep_naive & !keep_documented)
cat(sprintf(
  "\n%d gene(s) pass the naive rule ONLY because reads are concentrated in a\n",
  only_naive_keeps))
cat("single sample rather than spread across a group - the documented rule\n")
cat("correctly excludes these from downstream testing.\n")

counts_filtered <- counts[keep_documented, , drop = FALSE]
cat(sprintf(
  "\nFinal filtered matrix: %d genes x %d samples (%.1f%% of genes retained)\n",
  nrow(counts_filtered), ncol(counts_filtered),
  100 * nrow(counts_filtered) / nrow(counts)))

out_path <- "counts_filtered.tsv"
write.table(counts_filtered, out_path, sep = "\t", quote = FALSE, col.names = NA)
cat("\nWrote", out_path, "for use by Script 4 (exploratory analysis) onward.\n")

cat("\n==================================================================\n")
cat(" FILTERING COMPLETE - the rule above belongs in the portfolio report\n")
cat(" (Script 7) verbatim, so a reader can see exactly what was excluded and why.\n")
cat("==================================================================\n")

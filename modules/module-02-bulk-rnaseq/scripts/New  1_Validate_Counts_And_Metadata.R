# ==============================================================================
# MODULE 2 · SCRIPT 1 of 7 — VALIDATE THE COUNT MATRIX AND SAMPLE METADATA
# ==============================================================================
# ONE NEW IDEA on top of Module 1: Module 1 validated *reads* before alignment.
# This script validates the *count matrix* before any statistical test — the
# same "check inputs before trusting outputs" habit, one level up the pipeline.
#
# TEACHING NOTE — why this has to be its own step, not folded into DESeq2:
#   DESeqDataSetFromMatrix() will happily accept a count matrix and metadata
#   that are silently misaligned, non-integer, or missing samples, and will
#   not error until much later (or never — it may just give you wrong
#   p-values). By the time a wrong result surfaces, hours of compute and a
#   plausible-looking volcano plot have already been produced. Validating
#   BEFORE the test is the only way to catch these problems while they are
#   still cheap to fix.
#
# Input:  counts.txt   — featureCounts output from
#                         master_rnaseq_pipeline_consolidated.sh (Module 2 wet
#                         -lab script), OR any gene x sample integer matrix.
#         metadata.csv — one row per sample: sample,condition[,batch,...]
# Output: prints a PASS/FAIL validation report to the console; stops with a
#         clear error on the first hard failure (never on a soft warning).
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
counts_path   <- if (length(args) >= 1) args[1] else "counts.txt"
metadata_path <- if (length(args) >= 2) args[2] else "metadata.csv"

cat("==================================================================\n")
cat(" MODULE 2 / SCRIPT 1 — Count matrix + metadata validation\n")
cat("==================================================================\n")
cat("Counts file:   ", counts_path, "\n")
cat("Metadata file: ", metadata_path, "\n\n")

fail <- function(msg) {
  cat("[FAIL]", msg, "\n")
  stop(msg, call. = FALSE)
}
warn_flag <- function(msg) cat("[WARN]", msg, "\n")
pass <- function(msg) cat("[PASS]", msg, "\n")

# --- 1. Files must exist -----------------------------------------------------
if (!file.exists(counts_path))   fail(paste("counts file not found:", counts_path))
if (!file.exists(metadata_path)) fail(paste("metadata file not found:", metadata_path))
pass("both input files exist")

# --- 2. Read the count matrix -------------------------------------------------
# featureCounts output has 6 leading annotation columns (Chr, Start, End,
# Strand, Length, ...) before the per-sample count columns. Detect and strip
# them automatically so this script works on raw featureCounts output OR an
# already-cleaned matrix without the caller having to know which one it has.
raw <- read.table(counts_path, header = TRUE, row.names = 1,
                   sep = "\t", check.names = FALSE, comment.char = "#")

known_annotation_cols <- c("Chr", "Start", "End", "Strand", "Length")
annotation_present <- intersect(known_annotation_cols, colnames(raw))
if (length(annotation_present) > 0) {
  cat("Detected raw featureCounts annotation columns:",
      paste(annotation_present, collapse = ", "), "- stripping them.\n")
  counts <- raw[, !(colnames(raw) %in% known_annotation_cols), drop = FALSE]
} else {
  counts <- raw
}

# --- 3. Counts must be non-negative integers ----------------------------------
# TEACHING NOTE: DESeq2's negative-binomial model assumes raw integer counts.
# TPM, FPKM, or already-normalized values will run through DESeq2 without
# error and produce meaningless statistics. This check exists because that
# mistake produces no crash — only a wrong answer.
all_numeric <- all(sapply(counts, is.numeric))
if (!all_numeric) fail("count matrix contains non-numeric columns after stripping annotation - check the input file format")

non_integer <- sum(abs(as.matrix(counts) - round(as.matrix(counts))) > 1e-6)
if (non_integer > 0) {
  fail(sprintf(
    "%d values are non-integer. This looks like normalized data (TPM/FPKM/CPM), not raw counts. DESeq2 requires raw integer counts.",
    non_integer))
}
if (any(as.matrix(counts) < 0)) fail("count matrix contains negative values")
pass(sprintf("count matrix is non-negative integers (%d genes x %d samples)",
             nrow(counts), ncol(counts)))

# --- 4. Read and validate metadata --------------------------------------------
meta <- read.csv(metadata_path, stringsAsFactors = FALSE)
if (!all(c("sample", "condition") %in% colnames(meta))) {
  fail("metadata.csv must contain at least 'sample' and 'condition' columns")
}
if (any(duplicated(meta$sample))) {
  fail(paste("duplicate sample name(s) in metadata.csv:",
             paste(meta$sample[duplicated(meta$sample)], collapse = ", ")))
}
pass(sprintf("metadata.csv has %d unique sample rows", nrow(meta)))

# --- 5. Metadata <-> counts must match BY NAME, not by position ---------------
# TEACHING NOTE: this is the exact bug documented in rnaseq_deseq2_analysis.R's
# corrected-version header — the original script hardcoded condition labels in
# column order instead of matching by sample name. This script makes that
# check a named, reusable, first-class step instead of something buried
# inside the analysis script.
missing_in_meta   <- setdiff(colnames(counts), meta$sample)
missing_in_counts <- setdiff(meta$sample, colnames(counts))
if (length(missing_in_meta) > 0) {
  fail(paste("count columns with no metadata row:", paste(missing_in_meta, collapse = ", ")))
}
if (length(missing_in_counts) > 0) {
  fail(paste("metadata rows with no matching count column:", paste(missing_in_counts, collapse = ", ")))
}
pass("every count column has exactly one matching metadata row, matched by sample name")

# --- 6. Library size and gene-detection sanity checks -------------------------
lib_sizes <- colSums(counts)
detected  <- colSums(counts > 0)
report <- data.frame(
  sample          = colnames(counts),
  library_size    = lib_sizes,
  genes_detected  = detected,
  pct_of_genes    = round(100 * detected / nrow(counts), 1)
)
cat("\nPer-sample summary:\n")
print(report, row.names = FALSE)

med_lib <- median(lib_sizes)
low_lib <- report$sample[lib_sizes < 0.3 * med_lib]
if (length(low_lib) > 0) {
  warn_flag(paste0(
    "sample(s) with library size < 30% of the cohort median (possible failed ",
    "library, not necessarily fatal - inspect before excluding): ",
    paste(low_lib, collapse = ", ")))
} else {
  pass("no sample's library size is below 30% of the cohort median")
}

zero_count_genes <- sum(rowSums(counts) == 0)
cat(sprintf("\n%d of %d genes (%.1f%%) have zero counts across ALL samples ",
            zero_count_genes, nrow(counts), 100 * zero_count_genes / nrow(counts)))
cat("- expected, and handled by the filtering step (Script 3), not this one.\n")

cat("\n==================================================================\n")
cat(" VALIDATION COMPLETE - safe to proceed to Script 2 (design/confounding)\n")
cat("==================================================================\n")

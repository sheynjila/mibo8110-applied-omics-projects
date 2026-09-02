# ==============================================================================
# MODULE 4 · SCRIPT 5 of 7 — COHORT VCF QC SUMMARY
# ==============================================================================
# ONE NEW IDEA on top of Scripts 1-4: those scripts validated and improved
# each SAMPLE's own calls. This script looks at the merged COHORT VCF as a
# whole, the way Module 2 Script 4 looked at the whole count matrix with
# PCA before trusting any DE test — a first read of overall call-set
# trustworthiness before variant_analysis.R's SNP-distance clustering
# (Tier A) is asked to mean anything biologically.
#
# TEACHING NOTE — why Ti/Tv is the headline number here:
#   Transitions (A<->G, C<->T) are chemically favored over transversions
#   (the other four substitution pairs) by real mutational processes, so
#   genuine biological variant sets consistently show MORE transitions
#   than transversions. Sequencing/alignment noise has no such preference:
#   of the 6 possible substitution pairs, 2 are transitions and 4 are
#   transversions, so a call set dominated by artifacts drifts toward
#   Ti/Tv = 2/4 = 0.5 (the count expected from picking substitution pairs
#   at random), not above it. A Ti/Tv ratio sitting AT or BELOW ~0.5 is a
#   cohort-level red flag that the call set skews toward noise, no matter
#   how clean any single sample's QUAL/DP looked in Script 3. There is no
#   single universal "correct" Ti/Tv value across all organisms and
#   library preps — this script reports the number and the direction of
#   the flag, not a pass/fail verdict, the same caveat-not-verdict pattern
#   Module 3 used for N50.
#
# TEACHING NOTE — why heterozygosity ties back to Script 4:
#   This reference genome is haploid (see Script 4). If this cohort VCF was
#   built from *_raw_haploid.vcf (ploidy-corrected) inputs, per-sample
#   heterozygosity should be at or near 0%. Non-trivial heterozygosity here
#   is the cohort-level echo of the diploid-default bug Script 4 exists to
#   catch — this script flags it again at the merged level as a final
#   check that the correction actually propagated forward.
#
# Input:  merged_cohort.vcf (from master_snp_pipeline.sh Phase 5, or your own
#         re-merge of the *_raw_haploid.vcf / *_filtered_v2.vcf files per
#         this module's README "How to run").
# Output: prints a PASS/WARN QC report; writes cohort_vcf_qc_summary.tsv.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
vcf_path <- if (length(args) >= 1) args[1] else "merged_cohort.vcf"

cat("==================================================================\n")
cat(" MODULE 4 / SCRIPT 5 — Cohort VCF QC summary\n")
cat("==================================================================\n")
cat("Cohort VCF: ", vcf_path, "\n\n")

pass <- function(msg) cat("[PASS]", msg, "\n")
warn_flag <- function(msg) cat("[WARN]", msg, "\n")
fail <- function(msg) { cat("[FAIL]", msg, "\n"); stop(msg, call. = FALSE) }

if (!file.exists(vcf_path)) fail(paste("cohort VCF not found:", vcf_path))

suppressPackageStartupMessages(library(vcfR))
vcf <- read.vcfR(vcf_path, verbose = FALSE)

fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
gt  <- extract.gt(vcf, element = "GT")
samples <- colnames(gt)
n_sites <- nrow(fix)
pass(sprintf("loaded %d sites x %d samples from %s", n_sites, length(samples), vcf_path))

# --- 1. Classify each ALT allele as SNP or indel (length-of-REF/ALT test,
#        splitting multi-allelic ALT lists) -----------------------------------
is_transition <- function(ref, alt) {
  pair <- toupper(paste0(ref, alt))
  pair %in% c("AG", "GA", "CT", "TC")
}

n_snp <- 0; n_indel <- 0; n_ti <- 0; n_tv <- 0
for (i in seq_len(n_sites)) {
  ref <- fix$REF[i]
  alts <- strsplit(fix$ALT[i], ",", fixed = TRUE)[[1]]
  for (alt in alts) {
    if (nchar(ref) == 1 && nchar(alt) == 1 && alt != "*") {
      n_snp <- n_snp + 1
      if (is_transition(ref, alt)) n_ti <- n_ti + 1 else n_tv <- n_tv + 1
    } else {
      n_indel <- n_indel + 1
    }
  }
}
titv <- if (n_tv > 0) round(n_ti / n_tv, 3) else NA

cat(sprintf("\nSNPs: %d   Indels: %d   Transitions: %d   Transversions: %d\n",
            n_snp, n_indel, n_ti, n_tv))
if (is.na(titv)) {
  warn_flag("no transversions observed — Ti/Tv is undefined; too few SNPs to interpret")
} else if (titv <= 0.6) {
  warn_flag(sprintf(
    "Ti/Tv = %.3f, at or near the 0.5 value expected from random substitutions — this cohort's SNP calls skew toward noise, not biology",
    titv))
} else {
  pass(sprintf("Ti/Tv = %.3f, above the 0.5 random-noise floor", titv))
}

# --- 2. Per-sample missingness and heterozygosity -----------------------------
cat("\nPer-sample summary:\n")
summary_rows <- data.frame(sample = character(), pct_missing = numeric(),
                            pct_heterozygous = numeric(), stringsAsFactors = FALSE)

split_gt <- function(g) strsplit(g, "[/|]")

any_het <- FALSE
for (s in samples) {
  calls <- gt[, s]
  n_missing <- sum(is.na(calls) | calls %in% c(".", "./.", ".|."))
  called <- calls[!(is.na(calls) | calls %in% c(".", "./.", ".|."))]
  n_het <- 0
  if (length(called) > 0) {
    alleles <- split_gt(called)
    n_het <- sum(vapply(alleles, function(a) length(unique(a)) > 1, logical(1)))
  }
  pct_missing <- round(100 * n_missing / n_sites, 1)
  pct_het <- if (length(called) > 0) round(100 * n_het / length(called), 1) else NA

  summary_rows <- rbind(summary_rows, data.frame(
    sample = s, pct_missing = pct_missing, pct_heterozygous = pct_het))
  if (!is.na(pct_het) && pct_het > 0) any_het <- TRUE
}
print(summary_rows, row.names = FALSE)

if (any_het) {
  warn_flag(paste(
    "one or more samples show non-zero heterozygosity against a haploid reference.",
    "If this cohort VCF was built from the original diploid-default *_raw.vcf files",
    "instead of Script 4's *_raw_haploid.vcf, re-merge from the haploid-corrected",
    "calls and re-run this script before trusting these numbers."))
} else {
  pass("no heterozygous calls found — consistent with a haploid reference")
}

out_path <- "cohort_vcf_qc_summary.tsv"
write.table(
  data.frame(metric = c("n_sites", "n_snp", "n_indel", "n_transitions", "n_transversions", "titv"),
             value  = c(n_sites, n_snp, n_indel, n_ti, n_tv, ifelse(is.na(titv), "NA", titv))),
  out_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summary_rows, "cohort_vcf_qc_per_sample.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nWrote %s and cohort_vcf_qc_per_sample.tsv\n", out_path))

cat("\n==================================================================\n")
cat(" COHORT QC COMPLETE — safe to proceed to Script 6 (limitations)\n")
cat("==================================================================\n")

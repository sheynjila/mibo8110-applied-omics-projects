# ==============================================================================
# MODULE 4 · SCRIPT 7 of 7 — ASSEMBLE THE VARIANT REPORT
# ==============================================================================
# ONE NEW IDEA on top of Scripts 1-6: each prior script computed ONE fact
# about this call set in isolation (input validity, mapping/reference bias,
# a defensible filter, a ploidy correction, cohort-level QC, computable
# limitations). This script assembles those facts into the curriculum's
# required portfolio artifact — "a filtered VCF + a QC summary + an explicit
# statement of the call set's limitations" — the same auto-fill-the-objective-
# sections/leave-interpretation-blank pattern Module 2's DE_REPORT.md and
# Module 3's ASSEMBLY_REPORT.md both use. No script in this module, including
# this one, auto-generates a verdict on whether a variant call is real; that
# is the analyst's judgment call, informed by — not replaced by — the six
# prior scripts' output.
#
# Input:  WORKDIR — reads mapping_bias_summary.tsv (Script 2),
#         cohort_vcf_qc_summary.tsv / cohort_vcf_qc_per_sample.tsv (Script
#         5), LIMITATIONS.md (Script 6), and outbreak_clusters.tsv
#         (variant_analysis.R, Tier A) if present. Missing inputs are noted
#         explicitly rather than silently skipped.
# Output: VARIANT_REPORT.md in WORKDIR.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[1] else "."

cat("==================================================================\n")
cat(" MODULE 4 / SCRIPT 7 — Assembling VARIANT_REPORT.md\n")
cat("==================================================================\n")
cat("Working directory: ", workdir, "\n\n")

read_text_or_note <- function(path, missing_note) {
  if (file.exists(path)) paste(readLines(path, warn = FALSE), collapse = "\n")
  else missing_note
}

read_tsv_as_md_table <- function(path, missing_note) {
  if (!file.exists(path)) return(missing_note)
  df <- tryCatch(read.delim(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(missing_note)
  header <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  sep    <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows   <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

mapping_md  <- read_tsv_as_md_table(file.path(workdir, "mapping_bias_summary.tsv"),
                                     "_Not available — run Script 2 (`2_Assess_Mapping_And_Reference_Bias.sh`)._")
qc_summary  <- read_tsv_as_md_table(file.path(workdir, "cohort_vcf_qc_summary.tsv"),
                                     "_Not available — run Script 5 (`5_Cohort_VCF_QC_Summary.R`) on the merged cohort VCF._")
qc_persamp  <- read_tsv_as_md_table(file.path(workdir, "cohort_vcf_qc_per_sample.tsv"),
                                     "_Not available — run Script 5._")
limitations <- read_text_or_note(file.path(workdir, "LIMITATIONS.md"),
                                  "_Not available — run Script 6 (`6_Coverage_Gaps_And_Limitations.sh`)._")

hap_files <- Sys.glob(file.path(workdir, "variants", "*_raw_haploid.vcf"))
ploidy_status <- if (length(hap_files) > 0) {
  sprintf("Script 4 was run — ploidy-corrected calls available for %d sample(s).", length(hap_files))
} else {
  "**Script 4 was not run.** Calls referenced below may still reflect bcftools call's diploid default against this haploid bacterial reference (see Script 4)."
}

clusters_path <- file.path(workdir, "variants", "outbreak_clusters.tsv")
clusters_md <- if (file.exists(clusters_path)) {
  read_tsv_as_md_table(clusters_path, "_present but empty — no isolate pairs within the SNP threshold._")
} else {
  "_Not available — run `variant_analysis.R` (Tier A) on the cohort VCF. No `outbreak_clusters.tsv` means either it hasn't been run, or no isolate pairs fell within `SNP_THRESHOLD`._"
}

report <- c(
  "# Variant Calling Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  "Auto-assembled by `7_Assemble_Variant_Report.R` from Scripts 1-6's own",
  "output plus Tier A's `variant_analysis.R`. Sections marked **[FILL IN]**",
  "are intentionally left blank — no script in this module should generate",
  "a verdict on whether a variant is biologically real on its own.",
  "",
  "---",
  "",
  "## 1. Input validation (Script 1)",
  "- See Script 1 console output for the full per-sample PASS/FAIL log:",
  "  BAM/VCF structural integrity, and BAM `@RG SM:` / VCF sample-column",
  "  identity matched **by name**, not by filename assumption.",
  "",
  "## 2. Mapping and reference-bias assessment (Script 2)",
  mapping_md,
  "",
  "## 3. Variant filter design (Script 3)",
  "- Naive fixed filter (Tier A, `master_snp_pipeline.sh`): QUAL >= 20, DP >= 10, MQ >= 30.",
  "- Redesigned filter (`*_filtered_v2.vcf`): per-sample thresholds at or",
  "  above the naive values, driven by this cohort's own 10th-percentile",
  "  QUAL/DP/MQ, plus an INFO/DP4-derived strand-bias exclusion the naive",
  "  filter never applied. See Script 3 console output for the exact",
  "  per-sample thresholds used and the naive-vs-v2 PASS-count comparison.",
  "",
  "## 4. Ploidy assumption check (Script 4)",
  paste0("- ", ploidy_status),
  "- Reference GCF_000006945.2 (*S.* Typhimurium LT2) is haploid; any",
  "  heterozygous genotype call in this cohort should be treated as an",
  "  artifact of the wrong diploid default until proven otherwise by other",
  "  evidence (e.g. a genuinely mixed culture), not as ordinary variation.",
  "",
  "## 5. Cohort VCF QC summary (Script 5)",
  "**Cohort-level metrics:**",
  "",
  qc_summary,
  "",
  "**Per-sample missingness / heterozygosity:**",
  "",
  qc_persamp,
  "",
  "## 6. Outbreak / relatedness signal (Tier A: `variant_analysis.R`)",
  "- SNP/Hamming-distance clusters within `SNP_THRESHOLD` SNPs of each other:",
  "",
  clusters_md,
  "",
  "- See `phylogenetic_tree.pdf` and `variant_quality_qc.pdf` for the",
  "  corresponding visualizations.",
  "- A shared SNP cluster is a hypothesis, not proof of an epidemiological",
  "  link — interpretation must also weigh reference choice, filtering",
  "  criteria, missing data, recombination, and sampling context (see",
  "  `variant_analysis.R`'s own closing note).",
  "",
  "## 7. Limitations (Script 6, computed)",
  "",
  limitations,
  "",
  "## 8. Interpretation **[FILL IN]**",
  "- Do the flagged samples in Section 2 (mapping/coverage) change how much",
  "  you trust their variant calls specifically, or the cohort as a whole?",
  "- Does Section 5's Ti/Tv ratio support treating this SNP set as",
  "  biological signal rather than noise?",
  "- Do the Section 6 clusters correspond to any known epidemiological",
  "  context you have for these isolates (shared source, timeframe,",
  "  geography), or are they SNP-distance alone?",
  "",
  "(blank — to be completed by the analyst)",
  ""
)

out_path <- file.path(workdir, "VARIANT_REPORT.md")
writeLines(report, out_path)
cat(sprintf("Wrote %s (%d lines)\n", out_path, length(report)))
cat("Section 8 is intentionally blank — fill it in before treating this as a\n")
cat("finished portfolio artifact.\n")

cat("\n==================================================================\n")
cat(" PIPELINE COMPLETE (Scripts 1-7). VARIANT_REPORT.md is this module's\n")
cat(" required portfolio artifact once Section 8 is filled in.\n")
cat("==================================================================\n")

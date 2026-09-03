###############################################################################
# CAPSTONE · SCRIPT 1 of 3 — ACQUIRE TCGA-BRCA BULK RNA-SEQ + CLINICAL DATA
# ==============================================================================
# NEW SCRIPT, not part of the delivered material. PRJNA482620_TCGA_BRCA_downstream.R
# (this project's Tier A) assumes three files already exist in its working
# directory — TCGA_BRCA_raw_counts.tsv, TCGA_BRCA_tpm.tsv, TCGA_BRCA_clinical.csv
# — without saying how to obtain them. TCGA-BRCA is NOT retrievable as raw SRA
# FASTQ the way every other module in this repository downloads its data
# (dbGaP-controlled access for the raw reads); the standard, legitimate route
# is the Genomic Data Commons (GDC) portal's already-aligned, already-quantified
# gene-level output, via the TCGAbiolinks Bioconductor package. This script is
# that acquisition step — "Phase 2: Data and compute plan" of this capstone's
# 8 phases, made concrete and reproducible rather than a manual download.
#
# TEACHING NOTE — why STAR-Counts, and why "unstranded" vs. "tpm_unstrand":
#   GDC's current harmonized RNA-seq workflow (STAR alignment + STAR's own
#   gene-level quantification, replacing the older HTSeq-based pipeline) ships
#   several quantification columns in the same file: raw counts (three
#   strandedness variants) and derived TPM/FPKM values from the same counts.
#   The source paper's own methods (see this project's README "Data provenance")
#   state raw counts were used for DESeq2 (which requires actual integer counts
#   — it models count noise directly, not a pre-normalized value) and TPM,
#   log2-transformed, for everything else (ssGSEA, WGCNA, LASSO) — this script
#   pulls both from the SAME GDCprepare() call so they are guaranteed to be the
#   same samples, same gene set, not two separate downloads that could drift.
#
# Install once if needed (Bioconductor, not CRAN):
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))
###############################################################################

library(TCGAbiolinks)
library(SummarizedExperiment)

# ==============================================================================
# PHASE A: QUERY AND DOWNLOAD GENE EXPRESSION DATA
# ==============================================================================
cat("Querying GDC for TCGA-BRCA STAR-Counts gene expression data...\n")
query_expr <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts")

# GDCdownload caches to GDCdata/ in the current working directory - re-running
# this script will NOT re-download files already present, unlike this
# repository's shell-pipeline SRA downloads which check for output files, not
# GDC's own cache directly; confirm GDCdata/ has real content before assuming
# a fast re-run means nothing downloaded.
GDCdownload(query_expr, method = "api", files.per.chunk = 50)
brca_se <- GDCprepare(query_expr)

cat(sprintf("Prepared expression data: %d genes x %d samples.\n",
            nrow(brca_se), ncol(brca_se)))

# ==============================================================================
# PHASE B: EXTRACT RAW COUNTS AND TPM (SAME SAMPLES, SAME GENE SET)
# ==============================================================================
raw_counts <- assay(brca_se, "unstranded")
tpm <- assay(brca_se, "tpm_unstrand")

# The source paper's training cohort is 1,066 tumor + 99 normal (n=1,165) -
# TCGA-BRCA's actual current sample count on GDC may differ (additional
# samples deposited, or withdrawn/redacted ones removed, since that paper's
# own download date) - this is expected drift in a living database, not a
# script error; do not force the count to match 1,165 by dropping or
# duplicating samples.
n_tumor <- sum(substr(colnames(brca_se), 14, 15) == "01")
n_normal <- sum(substr(colnames(brca_se), 14, 15) == "11")
cat(sprintf("Sample composition (by TCGA barcode sample-type code): %d tumor (01), %d normal (11), %d other.\n",
            n_tumor, n_normal, ncol(brca_se) - n_tumor - n_normal))

# ==============================================================================
# PHASE C: QUERY CLINICAL DATA AND DERIVE OVERALL-SURVIVAL FIELDS
# ==============================================================================
cat("\nQuerying GDC for TCGA-BRCA clinical data...\n")
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# TCGAbiolinks' clinical table has vital_status / days_to_death /
# days_to_last_follow_up, not the plain time/status columns
# PRJNA482620_TCGA_BRCA_downstream.R's Surv(time=, event=) call expects -
# this is the standard TCGA overall-survival derivation (event=1 if
# deceased, using days_to_death; event=0 if alive, using
# days_to_last_follow_up as the censoring time), not a shortcut invented
# for this script.
clinical$status <- ifelse(clinical$vital_status == "Dead", 1, 0)
clinical$time <- ifelse(clinical$vital_status == "Dead",
                         clinical$days_to_death,
                         clinical$days_to_last_follow_up)

n_missing_time <- sum(is.na(clinical$time))
if (n_missing_time > 0) {
    cat(sprintf("[WARN] %d/%d patients have no usable survival time (missing both days_to_death\n",
                n_missing_time, nrow(clinical)))
    cat("       and days_to_last_follow_up) - these will be dropped by any Surv()-based analysis\n")
    cat("       downstream, not silently imputed. This is a real, expected gap in TCGA's clinical\n")
    cat("       annotation, not a bug in this script.\n")
}

# ==============================================================================
# PHASE D: WRITE OUTPUT (matching PRJNA482620_TCGA_BRCA_downstream.R's expected filenames)
# ==============================================================================
write.table(raw_counts, "TCGA_BRCA_raw_counts.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(tpm, "TCGA_BRCA_tpm.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.csv(clinical, "TCGA_BRCA_clinical.csv", row.names = FALSE)

cat(sprintf("\nWrote TCGA_BRCA_raw_counts.tsv, TCGA_BRCA_tpm.tsv (%d genes x %d samples each), and TCGA_BRCA_clinical.csv (%d patients).\n",
            nrow(raw_counts), ncol(raw_counts), nrow(clinical)))
cat("NOTE: PRJNA482620_TCGA_BRCA_downstream.R reads TCGA_BRCA_clinical.csv with row.names=1 -\n")
cat("      confirm its first column (patient barcode) is suitable as a row name, or adjust that\n")
cat("      read.csv() call, before running Script 3.\n")

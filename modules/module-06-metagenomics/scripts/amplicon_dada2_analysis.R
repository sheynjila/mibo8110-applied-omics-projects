###############################################################################
# NEW PIPELINE — amplicon_dada2_analysis.R
# Downstream half of the 16S/ITS amplicon workflow (see
# master_amplicon_pipeline.sh for the upstream download/QC/primer-removal
# stage this script picks up from). Fills the amplicon gap identified in
# domain_coverage_audit_2026-08-31.md ("Microbiome (plant & animal)").
#
# Target input: trimmed_reads/<SRR>_R1_trimmed.fastq / _R2_trimmed.fastq
# produced by master_amplicon_pipeline.sh.
#
# Install once if needed (Bioconductor, not CRAN):
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("dada2")
###############################################################################

library(dada2)
library(ggplot2)

# ==============================================================================
# CONFIGURATION — must match the AMPLICON_TYPE used in master_amplicon_pipeline.sh
# ==============================================================================
AMPLICON_TYPE <- "16S"   # "16S" or "ITS" - drives truncation length and reference DB below

# Path to the primer-trimmed reads directory produced by the bash pipeline.
READS_DIR <- "trimmed_reads"

# Reference database for taxonomy assignment (assignTaxonomy expects a
# training-set FASTA in the DADA2-specific format - NOT the raw SILVA/UNITE
# release files). Download the DADA2-formatted training sets from the
# official DADA2 tutorial page (https://benjjneb.github.io/dada2/training.html):
#   16S -> SILVA v138.1 (e.g. "silva_nr99_v138.1_train_set.fa.gz")
#   ITS -> UNITE general FASTA release (e.g. "sh_general_release_dynamic.fasta")
REF_DB_PATH <- switch(AMPLICON_TYPE,
                       "16S" = "/path/to/silva_nr99_v138.1_train_set.fa.gz",
                       "ITS" = "/path/to/unite_general_release_dynamic.fasta",
                       stop("Unknown AMPLICON_TYPE - use '16S' or 'ITS'."))

# ==============================================================================
# PHASE A: LOCATE SAMPLES
# ==============================================================================
fnFs <- sort(list.files(READS_DIR, pattern = "_R1_trimmed.fastq$", full.names = TRUE))
fnRs <- sort(list.files(READS_DIR, pattern = "_R2_trimmed.fastq$", full.names = TRUE))
sample_names <- sub("_R1_trimmed.fastq$", "", basename(fnFs))

if (length(fnFs) == 0 || length(fnFs) != length(fnRs)) {
    stop("No matched forward/reverse trimmed FASTQ pairs found in ", READS_DIR,
         " - confirm master_amplicon_pipeline.sh completed successfully.")
}
cat(sprintf("Found %d samples.\n", length(fnFs)))

# Inspect quality profiles before choosing truncation length - always look at
# these plots for your own data rather than trusting the defaults blindly.
pdf("quality_profiles.pdf", width = 8, height = 6)
print(plotQualityProfile(fnFs[1:min(2, length(fnFs))]))
print(plotQualityProfile(fnRs[1:min(2, length(fnRs))]))
dev.off()

# ==============================================================================
# PHASE B: FILTER & TRIM
# ==============================================================================
filtFs <- file.path("filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

if (AMPLICON_TYPE == "16S") {
    # 16S amplifies a fixed-length hypervariable region (e.g. V4), so a
    # fixed truncation length is appropriate and standard. 240/160 are
    # typical starting points for 2x250bp MiSeq V4 data after primer
    # removal - adjust based on quality_profiles.pdf for your actual run.
    truncLen_F <- 240
    truncLen_R <- 160
} else {
    # ITS region length varies genuinely across fungal taxa. Applying a
    # fixed truncLen here would truncate real amplicons of different
    # (correct) lengths and bias the resulting ASVs. DADA2's official
    # ITS workflow recommendation is truncLen = 0 (no truncation) -
    # quality filtering alone (via truncQ / maxEE below) is used instead.
    truncLen_F <- 0
    truncLen_R <- 0
}

filt_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                           truncLen = c(truncLen_F, truncLen_R),
                           maxN = 0, maxEE = c(2, 2), truncQ = 2,
                           rm.phix = TRUE, compress = TRUE, multithread = TRUE)
cat("Filter/trim summary:\n")
print(filt_out)

# ==============================================================================
# PHASE C: DENOISE & MERGE (ASV inference)
# ==============================================================================
cat("Learning error rates...\n")
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

cat("Running DADA2 core denoising algorithm...\n")
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

cat("Merging paired reads...\n")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
cat(sprintf("Sequence table: %d samples x %d ASVs (pre-chimera-removal)\n",
            nrow(seqtab), ncol(seqtab)))

# Remove chimeric sequences (PCR artifacts, not real biological ASVs)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                     multithread = TRUE, verbose = TRUE)
cat(sprintf("Retained %.1f%% of reads after chimera removal.\n",
            100 * sum(seqtab_nochim) / sum(seqtab)))

# ==============================================================================
# PHASE D: TAXONOMY ASSIGNMENT
# ==============================================================================
cat("Assigning taxonomy against reference database...\n")
if (!file.exists(REF_DB_PATH)) {
    stop("REF_DB_PATH not found at '", REF_DB_PATH, "'. Download the ",
         "DADA2-formatted training set for ", AMPLICON_TYPE,
         " before running this script (see comment above REF_DB_PATH).")
}
taxa <- assignTaxonomy(seqtab_nochim, REF_DB_PATH, multithread = TRUE)

# HOST/ORGANELLE CONTAMINATION FILTER (16S only):
# 16S primers targeting conserved bacterial rRNA regions can co-amplify
# plant chloroplast and animal/plant mitochondrial rRNA-like sequences,
# since both organelles retain bacterial-like rRNA genes from their
# endosymbiotic origin. For plant/animal host-associated 16S samples,
# these should be removed before computing relative abundances, or they
# will inflate false "bacterial" signal. This step is a no-op for ITS
# (fungal-specific primers do not have this failure mode).
if (AMPLICON_TYPE == "16S") {
    is_organelle <- !is.na(taxa[, "Order"]) & taxa[, "Order"] %in% c("Chloroplast") |
                    (!is.na(taxa[, "Family"]) & taxa[, "Family"] %in% c("Mitochondria"))
    n_organelle <- sum(is_organelle)
    if (n_organelle > 0) {
        cat(sprintf(">> Removing %d ASV(s) classified as Chloroplast/Mitochondria (host organelle contamination).\n",
                    n_organelle))
        seqtab_nochim <- seqtab_nochim[, !is_organelle, drop = FALSE]
        taxa <- taxa[!is_organelle, , drop = FALSE]
    } else {
        cat(">> No Chloroplast/Mitochondria ASVs detected.\n")
    }
}

# ==============================================================================
# PHASE E: EXPORT ASV TABLE + TAXONOMY (R-ready outputs)
# ==============================================================================
asv_ids <- paste0("ASV", seq_len(ncol(seqtab_nochim)))
asv_table <- t(seqtab_nochim)
rownames(asv_table) <- asv_ids
write.table(cbind(ASV_ID = asv_ids, asv_table), "ASV_table.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

taxa_out <- cbind(ASV_ID = asv_ids, taxa)
write.table(taxa_out, "taxonomy.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Long-format table matching the ALL_SAMPLES_MICROBIOME.tsv convention used
# by the shotgun pipeline's master_microbiome_pipeline.R, so the same
# downstream visualization pattern applies to amplicon data.
sample_totals <- rowSums(seqtab_nochim)
long_rows <- list()
for (s in seq_len(nrow(seqtab_nochim))) {
    sample_id <- rownames(seqtab_nochim)[s]
    for (a in seq_len(ncol(seqtab_nochim))) {
        reads <- seqtab_nochim[s, a]
        if (reads > 0) {
            genus <- taxa[a, "Genus"]
            label <- if (is.na(genus)) asv_ids[a] else genus
            long_rows[[length(long_rows) + 1]] <- data.frame(
                Sample = sample_id, Taxon = label,
                Estimated_Reads = reads,
                Fraction_Total = reads / sample_totals[s]
            )
        }
    }
}
all_samples_long <- do.call(rbind, long_rows)
write.table(all_samples_long, "ALL_SAMPLES_AMPLICON.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# PHASE F: COMPOSITION PLOT
# ==============================================================================
plot_df <- all_samples_long[all_samples_long$Fraction_Total > 0.01, ]
p <- ggplot(plot_df, aes(x = Sample, y = Fraction_Total, fill = Taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = sprintf("%s Amplicon Community Composition (Genus level, >1%% abundance)", AMPLICON_TYPE),
         x = "Sample", y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("amplicon_composition.pdf", plot = p, width = 10, height = 6)

cat("Analysis complete. Outputs: ASV_table.tsv, taxonomy.tsv, ALL_SAMPLES_AMPLICON.tsv, amplicon_composition.pdf, quality_profiles.pdf\n")

###############################################################################
# METHODOLOGY NOTES
# ---------------------------------------------------------------------------
# 1. ASVs vs. OTUs: DADA2 infers Amplicon Sequence Variants (exact sequence
#    resolution, denoised for sequencing error) rather than clustering into
#    97%-similarity OTUs. ASVs are the current best practice for amplicon
#    studies - they are directly comparable across studies/runs without
#    re-clustering, unlike OTUs.
# 2. Why truncLen differs by AMPLICON_TYPE: this is the single most common
#    correctness mistake in ITS pipelines built by copying 16S defaults -
#    truncating ITS reads to a fixed length silently discards true-length
#    ASVs and biases community composition toward taxa with shorter ITS
#    regions. Always confirm truncLen=0 for ITS unless you have a specific,
#    validated reason to do otherwise.
# 3. Chimera removal (removeBimeraDenovo) is essential for amplicon data:
#    PCR of a mixed community routinely produces chimeric amplicons (a
#    single molecule stitched from two different template sequences) at
#    non-trivial rates. Skipping this step inflates apparent diversity with
#    artifacts that were never present in the original sample.
# 4. This script assigns genus-level taxonomy for the composition plot;
#    assignTaxonomy() returns the full rank hierarchy (Kingdom...Genus) in
#    taxonomy.tsv if you need finer control (e.g. species-level calls via
#    addSpecies(), not included here since it requires an additional
#    reference file).
###############################################################################

###############################################################################
# CORRECTED VERSION of variant_analysis.R
# See domain_coverage_audit_2026-08-31.md ("Pathogen tracing / molecular
# epidemiology" section) for the full review. Fixes applied:
#
#   1. Removed the stray literal "1." that prefixed the line-1 comment
#      block in the original file (cosmetic artifact, harmless to execution
#      but confusing to read).
#   2. Replaced the Euclidean dist() on the genind object with a proper
#      SNP/Hamming distance (ape::dist.gene, pairwise deletion of missing
#      data) - the standard metric for bacterial outbreak phylogenetics,
#      where "genetic distance" should mean "number of differing SNP
#      sites", not a Euclidean norm over allele-count vectors.
#   3. Added explicit SNP-distance-threshold cluster calling: pairs/groups
#      of isolates within SNP_THRESHOLD SNPs of each other are now flagged
#      as a probable epidemiological cluster and written to
#      outbreak_clusters.tsv, instead of leaving cluster interpretation
#      entirely to eyeballing the tree image.
#
# Target Input: merged_cohort.vcf (generated via bcftools merge - see the
# BCFtools merge instructions retained below).
###############################################################################

# THIS IS TO BE COPIED LINE BY LINE AND EXECUTED BEFORE USING R
#  module load BCFtools/1.21-GCC-13.3.0
#  cd /scratch/$(whoami)/master_snp_pipeline/variants
#  Index individual files
#  for vcf in *_filtered.vcf; do bcftools index $vcf; done
#  Merge them into a single multi-sample file
#  bcftools merge *_filtered.vcf -O v -o merged_cohort.vcf

# Transfer merged_cohort.vcf to your local computer to use in RStudio.


# 2. The R Script: variant_analysis.R
# Save this in RStudio. You may need to run
# install.packages(c("vcfR", "ape", "adegenet")) in your R console first if you do not already have these genomics packages installed.


# ==============================================================================
# EPIDEMIOLOGICAL VARIANT ANALYSIS & PHYLOGENETICS
# ==============================================================================
# Target Input: merged_cohort.vcf (Generated via bcftools merge)
# ==============================================================================

# 1. Load Required Genomics Libraries
library(vcfR)       # For parsing VCF files
library(adegenet)    # For genetic data manipulation
library(ape)         # For phylogenetic tree construction

# ------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------
# SNP-distance threshold (in number of differing SNP sites) below which a
# pair of isolates is flagged as a probable epidemiological cluster.
# This is pathogen- and context-dependent (published outbreak thresholds
# range roughly 0-20 SNPs depending on organism mutation rate and outbreak
# timescale) - 10 is a commonly cited general-purpose starting point for
# bacterial foodborne outbreak work. Adjust based on your organism and the
# literature for your specific study system before drawing conclusions.
SNP_THRESHOLD <- 10

# ==============================================================================
# PHASE A: DATA IMPORT & QUALITY CONTROL
# ==============================================================================
cat("Loading multi-sample VCF file...\n")
vcf <- read.vcfR("merged_cohort.vcf")

# Extract and visualize overall quality scores to confirm Bash filtering worked
qual_scores <- getQUAL(vcf)

pdf("variant_quality_qc.pdf", width=8, height=6)
hist(qual_scores,
     breaks=50,
     col="steelblue",
     main="Distribution of Variant Quality Scores",
     xlab="Phred-scaled Quality (QUAL)",
     ylab="Frequency")
abline(v=20, col="red", lwd=2, lty=2) # Our bash pipeline threshold
legend("topright", legend=c("Filtering Threshold (20)"), col="red", lty=2, lwd=2)
dev.off()

# ==============================================================================
# PHASE B: GENETIC CONVERSION
# ==============================================================================
cat("Converting VCF to genind object for distance calculations...\n")
# Convert the vcfR object into a 'genind' object (used by adegenet for population genetics)
genind_obj <- vcfR2genind(vcf)

# ==============================================================================
# PHASE C: SNP DISTANCE CALCULATION (CORRECTED)
# ==============================================================================
cat("Calculating pairwise SNP (Hamming) distances...\n")

# FIX #2: dist() on a genind object computes Euclidean distance over
# allele-count vectors, which does not correspond to "number of SNPs that
# differ between two isolates" - the quantity outbreak investigators
# actually mean by "genetic distance". ape::dist.gene() instead counts, for
# each pair of isolates, the number of loci at which their genotypes
# differ (pairwise deletion of missing calls), which is the standard
# SNP/Hamming distance used in bacterial genomic epidemiology.
genotype_matrix <- tab(genind_obj)
genetic_dist <- dist.gene(genotype_matrix, method = "pairwise", pairwise.deletion = TRUE)

cat("SNP distance matrix (first few pairs):\n")
print(round(as.matrix(genetic_dist)[1:min(5, nrow(as.matrix(genetic_dist))),
                                     1:min(5, nrow(as.matrix(genetic_dist)))], 1))

# ==============================================================================
# PHASE D: OUTBREAK CLUSTER CALLING (NEW)
# ==============================================================================
# FIX #3: Explicitly report which isolate pairs fall within SNP_THRESHOLD
# SNPs of each other, rather than requiring visual inspection of the tree.
cat(sprintf("Identifying isolate pairs within %d SNPs (probable transmission links)...\n", SNP_THRESHOLD))

dist_matrix <- as.matrix(genetic_dist)
sample_names <- rownames(dist_matrix)
cluster_pairs <- data.frame(Sample_A = character(), Sample_B = character(),
                             SNP_Distance = numeric(), stringsAsFactors = FALSE)

if (!is.null(sample_names) && length(sample_names) > 1) {
    for (i in 1:(length(sample_names) - 1)) {
        for (j in (i + 1):length(sample_names)) {
            d <- dist_matrix[i, j]
            if (!is.na(d) && d <= SNP_THRESHOLD) {
                cluster_pairs <- rbind(cluster_pairs, data.frame(
                    Sample_A = sample_names[i],
                    Sample_B = sample_names[j],
                    SNP_Distance = d
                ))
            }
        }
    }
}

if (nrow(cluster_pairs) > 0) {
    cluster_pairs <- cluster_pairs[order(cluster_pairs$SNP_Distance), ]
    write.table(cluster_pairs, "outbreak_clusters.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("Found %d isolate pair(s) within the %d-SNP threshold - written to outbreak_clusters.tsv\n",
                nrow(cluster_pairs), SNP_THRESHOLD))
} else {
    cat(sprintf("No isolate pairs found within the %d-SNP threshold. No outbreak_clusters.tsv written.\n", SNP_THRESHOLD))
}

# ==============================================================================
# PHASE E: PHYLOGENETIC TREE CONSTRUCTION
# ==============================================================================
cat("Building Neighbor-Joining tree from SNP distances...\n")

# Construct a Neighbor-Joining (NJ) phylogenetic tree
# NJ is a standard method for bacterial outbreak clustering
nj_tree <- nj(genetic_dist)

# Export the tree visualization
pdf("phylogenetic_tree.pdf", width=8, height=6)
plot(nj_tree,
     main="Neighbor-Joining Tree of Bacterial Isolates (SNP distance)",
     cex=1.2,           # Label size
     type="unrooted",   # Standard for outbreak clusters without a known ancestor
     edge.width=2)
add.scale.bar(length=10, col="red") # Adds a scale bar for genetic distance (in SNPs)
dev.off()

cat("Analysis complete! Check your working directory for the PDFs and outbreak_clusters.tsv.\n")



# What This Script Accomplishes

# Quality Verification: It extracts the QUAL scores from your VCF and generates a histogram (variant_quality_qc.pdf). You will see a red dashed line at 20, proving visually that your bash script successfully dropped all low-quality garbage mutations.
# Data Conversion: The vcfR2genind() function acts as a bridge, translating raw genomic alleles into a statistical format that R can use for matrix mathematics.
# SNP Distance & Cluster Calling: It calculates a true SNP/Hamming distance (how many mutations separate Sample A from Sample B) rather than a generic Euclidean distance, flags isolate pairs within SNP_THRESHOLD SNPs as probable transmission links (outbreak_clusters.tsv), and plots a Neighbor-Joining tree (phylogenetic_tree.pdf). In a real outbreak investigation, isolates that cluster tightly together (both on the tree and in outbreak_clusters.tsv) share a recent common ancestor (e.g., the contaminated food source under investigation).

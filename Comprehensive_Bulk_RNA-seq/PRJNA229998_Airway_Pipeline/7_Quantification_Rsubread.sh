#!/bin/bash
#SBATCH --job-name=qc07_quant
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00

###############################################################################
# 7_Quantification_Rsubread.sh
#
# NARRATIVE:
# This script transitions the data from genomic coordinates to a discrete gene 
# count matrix. Because the standalone C-version of Subread was unavailable, 
# this script dynamically constructs and executes a transient R script. It uses 
# `Rsubread::featureCounts` to summarize BAM alignments against the GTF annotation, 
# outputting the numeric matrix required for DESeq2 statistical modeling.
###############################################################################
set -e
set -o pipefail

module load R

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep"
REFDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_ref"
mkdir -p "${WORKDIR}/counts"
cd "${WORKDIR}"

SRR="SRR1039508"
BAM="alignments/${SRR}_Aligned.sortedByCoord.out.bam"

cat << EOF > counts/run_quant.R
suppressMessages(library(Rsubread))
fc <- featureCounts(files = "${BAM}",
                    annot.ext = "${REFDIR}/GRCh38.gtf",
                    isGTFAnnotationFile = TRUE,
                    isPairedEnd = TRUE,
                    nthreads = 4)
write.table(fc\$counts, file="counts/${SRR}_counts.txt", sep="\t", quote=FALSE, col.names=NA)
write.table(fc\$stat, file="counts/${SRR}_counts.txt.summary", sep="\t", quote=FALSE, row.names=FALSE)
EOF

echo "Quantifying ${SRR} via Rscript..."
Rscript counts/run_quant.R
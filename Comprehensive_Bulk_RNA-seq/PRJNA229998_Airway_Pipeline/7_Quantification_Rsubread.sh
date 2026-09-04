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
# Quantifies gene expression using Rsubread::featureCounts. Includes an 
# automated bootstrap check to set up a local R library and install Rsubread 
# via BiocManager if missing from the system module.
###############################################################################
set -e
set -o pipefail

module purge
module load R

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep"
REFDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_ref"
mkdir -p "${WORKDIR}"/{counts,R_libs}
cd "${WORKDIR}"

export R_LIBS_USER="${WORKDIR}/R_libs"
SRR="SRR1039508"
BAM="alignments/${SRR}_Aligned.sortedByCoord.out.bam"

cat << EOF > counts/run_quant.R
user_lib <- "${WORKDIR}/R_libs"
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = user_lib)
}
if (!requireNamespace("Rsubread", quietly = TRUE)) {
    BiocManager::install("Rsubread", update = FALSE, ask = FALSE, lib = user_lib)
}

suppressMessages(library(Rsubread, lib.loc = user_lib))

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
#!/bin/bash
#SBATCH --job-name=qc06_align
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

###############################################################################
# 6_STAR_Alignment.sh
#
# NARRATIVE:
# This script maps the surviving, paired-verified reads to the human genome. 
# It first downloads the GRCh38 FASTA and GTF into a *shared* reference 
# directory to avoid redundant re-indexing across scripts. It uses STAR to 
# perform splice-aware alignment (crucial for eukaryotic RNA-seq) and directly 
# outputs a coordinate-sorted BAM file ready for feature counting.
###############################################################################
set -e
set -o pipefail

module load STAR/2.7.10b-GCC-11.3.0

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep"
REFDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_ref"
mkdir -p "${WORKDIR}/alignments" "${REFDIR}"
cd "${WORKDIR}"

SRR="SRR1039508"

FASTA_URL="http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

if [ ! -f "${REFDIR}/GRCh38.fa" ]; then
    wget -qO- "${FASTA_URL}" | gunzip > "${REFDIR}/GRCh38.fa"
fi

if [ ! -f "${REFDIR}/GRCh38.gtf" ]; then
    wget -qO- "${GTF_URL}" | gunzip > "${REFDIR}/GRCh38.gtf"
fi

if [ ! -d "${REFDIR}/star_index" ] || [ -z "$(ls -A ${REFDIR}/star_index)" ]; then
    mkdir -p "${REFDIR}/star_index"
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir "${REFDIR}/star_index" \
         --genomeFastaFiles "${REFDIR}/GRCh38.fa" --sjdbGTFfile "${REFDIR}/GRCh38.gtf" --sjdbOverhang 74
fi

echo "Aligning ${SRR}..."
STAR --runThreadN 8 --genomeDir "${REFDIR}/star_index" \
     --readFilesIn "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix "alignments/${SRR}_"
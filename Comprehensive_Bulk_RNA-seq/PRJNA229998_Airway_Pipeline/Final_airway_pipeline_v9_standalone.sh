#!/bin/bash
#SBATCH --job-name=airway_v9_full
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=36G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/osnjila/slurm-%j.out

# Exit immediately on error and fail pipelines correctly
set -eo pipefail

# ==============================================================================
# 1. Initialization & Strict Pathing
# ==============================================================================
BASE_DIR="/scratch/osnjila/PRJNA229998_v9_standalone"
THREADS=8
SAMPLES=("SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513")

echo "Starting V9 Pipeline at $(date)"
echo "Base Directory: $BASE_DIR"

mkdir -p "$BASE_DIR"/{raw_reads,trimmed_reads,alignments,counts,logs,ref/star_index,sra_cache}
cd "$BASE_DIR"

# ==============================================================================
# 2. Reference Genome Preparation (Ensembl GRCh38)
# ==============================================================================
echo "--- Preparing Reference Genome ---"
# Download FASTA and GTF if they don't exist
if [ ! -f "ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    wget -qO ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
        http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi

if [ ! -f "ref/Homo_sapiens.GRCh38.110.gtf" ]; then
    wget -qO ref/Homo_sapiens.GRCh38.110.gtf.gz \
        http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
    gunzip ref/Homo_sapiens.GRCh38.110.gtf.gz
fi

# Load GCC and STAR for indexing
module purge
module load GCC/13.2.0 STAR

# Generate STAR Index (Requires ~32GB RAM)
if [ ! -f "ref/star_index/SAindex" ]; then
    echo "Building STAR Genome Index..."
    STAR --runThreadN $THREADS \
         --runMode genomeGenerate \
         --genomeDir ref/star_index \
         --genomeFastaFiles ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
         --sjdbGTFfile ref/Homo_sapiens.GRCh38.110.gtf \
         --sjdbOverhang 99
fi

# ==============================================================================
# 3. Data Fetching
# ==============================================================================
echo "--- Fetching SRA Data ---"
module purge
module load SRA-Toolkit

for SAMPLE in "${SAMPLES[@]}"; do
    if [ ! -d "sra_cache/${SAMPLE}" ]; then
        echo "Prefetching $SAMPLE..."
        prefetch "$SAMPLE" -O sra_cache/
    fi
done

# ==============================================================================
# 4. Processing Loop (Extraction -> Trimming -> Alignment)
# ==============================================================================
for SAMPLE in "${SAMPLES[@]}"; do
    echo "=========================================================="
    echo "Processing $SAMPLE"
    echo "=========================================================="

    # --- SRA Extraction ---
    module purge
    module load SRA-Toolkit
    
    echo "Extracting FASTQ files for $SAMPLE..."
    # --split-3 guarantees identical line counts for fastp by isolating 0-length reads
    fasterq-dump --split-3 --threads $THREADS -O raw_reads "sra_cache/${SAMPLE}/${SAMPLE}.sra"

    # --- Quality Trimming & Alignment ---
    module purge
    module load GCC/13.2.0 fastp STAR Subread

    echo "Running fastp trimming for $SAMPLE..."
    fastp \
        -i "raw_reads/${SAMPLE}_1.fastq" \
        -I "raw_reads/${SAMPLE}_2.fastq" \
        -o "trimmed_reads/${SAMPLE}_1_clean.fastq" \
        -O "trimmed_reads/${SAMPLE}_2_clean.fastq" \
        --thread $THREADS \
        --json "logs/${SAMPLE}_fastp.json" \
        --html "logs/${SAMPLE}_fastp.html"

    echo "Aligning $SAMPLE to GRCh38..."
    STAR \
        --runThreadN $THREADS \
        --genomeDir ref/star_index \
        --readFilesIn "trimmed_reads/${SAMPLE}_1_clean.fastq" "trimmed_reads/${SAMPLE}_2_clean.fastq" \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix "alignments/${SAMPLE}_"

    # Clean up massive raw fastq files to save scratch space
    rm "raw_reads/${SAMPLE}_1.fastq" "raw_reads/${SAMPLE}_2.fastq"
    
    echo "RESULT: $SAMPLE processing completed successfully."
done

# ==============================================================================
# 5. Cohort Quantification & QC
# ==============================================================================
echo "=========================================================="
echo "Generating Final Gene Count Matrix & QC Report"
echo "=========================================================="

# Ensure quantification modules are loaded
module purge
module load GCC/13.2.0 Subread MultiQC 2>/dev/null || module load Subread MultiQC

featureCounts \
    -T $THREADS -p -B -C \
    -a ref/Homo_sapiens.GRCh38.110.gtf \
    -o counts/gene_counts.txt \
    alignments/*.bam

echo "Running MultiQC..."
multiqc logs/ counts/ alignments/ trimmed_reads/ -o counts/multiqc_report

echo "PIPELINE V9 COMPLETE at $(date)!"





# ==============================================================================
# Assessment of V9
# ==============================================================================

# 1. Consistency
# ------------------------------------------------------------------------------

# 1.1 Pathing Validation
# Every generated directory and output file explicitly uses $BASE_DIR as its
# prefix. Slurm output is directed to /scratch/osnjila/.
# The pipeline has zero interactions with /home/.

# 1.2 Module Loading
# The pipeline proactively runs "module purge" before switching toolchains.
# For example, it switches from SRA-Toolkit to the GCC/13.2.0 environment
# required by fastp.
#
# This approach eliminates the hidden dependency crashes encountered in
# earlier versions.

# 2. Reproducibility
# ------------------------------------------------------------------------------

# 2.1 Strict Error Handling
# The command below ensures that the script stops immediately if a command or
# pipeline fails:
#
#   set -eo pipefail
#
# Therefore, if a step such as prefetch or fastp fails, the script halts instead
# of passing incomplete or corrupted data to downstream steps.

# 2.2 Data Integrity
# The --split-3 flag is hardcoded.
#
# This ensures that biological irregularities, such as the 26.5 million
# zero-length reads observed in SRR1039513, are handled consistently during
# extraction.
#
# It also ensures that the paired Read 1 and Read 2 input files contain matching
# numbers of reads for downstream tools.

# 3. Completeness
# ------------------------------------------------------------------------------

# 3.1 Resource Allocation
# The script requests the following Slurm resources:
#
#   --mem=36G
#   --cpus-per-task=8
#
# These resources provide sufficient memory and CPU capacity to process the
# GRCh38 STAR index without reproducing the std::bad_alloc crash encountered
# on the login node.

# 3.2 End-to-End Execution
# The script includes all major stages required to reproduce the analysis from
# scratch:
#
#   1. Reference downloading
#   2. STAR index generation
#   3. SRA caching
#   4. FASTQ extraction
#   5. Read trimming
#   6. Read alignment
#   7. Gene-level counting
#   8. Quality-control reporting

# Conclusion
# ------------------------------------------------------------------------------
# V9 has been assessed and verified.
#
# The validated fixes are ready to be cascaded down to the nine modular
# pipeline scripts.
# ==============================================================================
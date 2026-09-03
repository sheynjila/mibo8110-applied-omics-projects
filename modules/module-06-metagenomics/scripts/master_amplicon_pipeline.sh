#!/bin/bash
#SBATCH --job-name=master_amplicon
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=master_amplicon_%j.out
#SBATCH --error=master_amplicon_%j.err

# ==============================================================================
# NEW PIPELINE — 16S / ITS AMPLICON SEQUENCING (primer trimming stage)
# See domain_coverage_audit_2026-08-31.md, "Microbiome (plant & animal)" gap:
# the shotgun pipeline (master_microbiome_pipeline.sh, Kraken2+Bracken) had
# no amplicon (16S rRNA / ITS) workflow. This script is the upstream half of
# that workflow: download -> QC -> primer removal. The downstream half
# (denoising, ASV inference, taxonomy) is amplicon_dada2_analysis.R.
#
# WHY THIS IS A SEPARATE PIPELINE FROM master_microbiome_pipeline.sh:
#   Shotgun metagenomics classifies whole-genome fragments directly with
#   Kraken2/Bracken. Amplicon sequencing instead sequences one targeted PCR
#   product per sample (16S rRNA gene for bacteria/archaea, or the ITS
#   region for fungi) and needs a completely different analytical path:
#   primer removal -> denoising into Amplicon Sequence Variants (ASVs) ->
#   reference-based taxonomy assignment. Kraken2/Bracken are not used at all
#   here.
#
# ONE SCRIPT, TWO MODES (AMPLICON_TYPE):
#   16S — bacterial/archaeal community profiling. Fixed-length amplicon
#         region (e.g. the V4 hypervariable region), so downstream DADA2
#         can safely truncate reads to a fixed length.
#   ITS — fungal community profiling. The ITS region has HIGHLY VARIABLE
#         length across taxa - fixed-length truncation in DADA2 would
#         discard real biological length variation, so ITS reads must NOT
#         be truncated downstream (this script only removes primers here;
#         see amplicon_dada2_analysis.R's TRUNC_LEN=0 handling for ITS).
#
# IMPORTANT CAVEAT for plant/animal host-associated samples: 16S primers
# (which target conserved bacterial rRNA regions) can still co-amplify
# host mitochondrial and, for plant samples, chloroplast rRNA-like
# sequences. This script does not filter those out (it happens after
# taxonomy assignment, in amplicon_dada2_analysis.R) - see that script's
# Phase D for the chloroplast/mitochondria removal step.
# ==============================================================================

set -e
set -o pipefail

# ==========================================
# CONFIGURATION
# ==========================================
AMPLICON_TYPE="${AMPLICON_TYPE:-16S}"   # 16S | ITS

case "${AMPLICON_TYPE}" in
    16S)
        # Earth Microbiome Project V4 primers (515F/806R) - the most common
        # bacterial/archaeal 16S primer pair in current use. Swap these if
        # your study used a different hypervariable region (e.g. V3-V4).
        FWD_PRIMER="${FWD_PRIMER:-GTGYCAGCMGCCGCGGTAA}"
        REV_PRIMER="${REV_PRIMER:-GGACTACNVGGGTWTCTAAT}"
        ;;
    ITS)
        # ITS1F / ITS2 fungal primers (targets the ITS1 region).
        FWD_PRIMER="${FWD_PRIMER:-CTTGGTCATTTAGAGGAAGTAA}"
        REV_PRIMER="${REV_PRIMER:-GCTGCGTTCTTCATCGATGC}"
        ;;
    *)
        echo "CRITICAL ERROR: Unknown AMPLICON_TYPE '${AMPLICON_TYPE}'. Use 16S or ITS."
        exit 1
        ;;
esac

echo "=========================================================="
echo " INITIATING AMPLICON PIPELINE (${AMPLICON_TYPE}) "
echo " Forward primer: ${FWD_PRIMER}"
echo " Reverse primer: ${REV_PRIMER}"
echo "=========================================================="

# ==========================================
# PHASE 1: DIRECTORY SETUP
# ==========================================
WORKDIR="/scratch/$(whoami)/master_amplicon_pipeline_${AMPLICON_TYPE}"
mkdir -p ${WORKDIR}/{raw_reads,trimmed_reads,qc,logs}
cd ${WORKDIR}

if [ ! -f "srr_list.txt" ]; then
    echo "WARNING: srr_list.txt not found. Generating default test list..."
    echo -e "SRR11092056\nSRR11092057" > srr_list.txt
fi
sed -i 's/\r$//' srr_list.txt

# ==========================================
# PHASE 2: DOWNLOAD, QC & PRIMER REMOVAL LOOP
# ==========================================
echo -e "\n[PHASE 2] EXECUTING PER-SAMPLE PRIMER TRIMMING..."

FAILED_SAMPLES=()

for SRR in $(cat srr_list.txt); do
    echo "----------------------------------------------------------"
    echo " PROCESSING SAMPLE: ${SRR}"
    echo "----------------------------------------------------------"

    set +e
    (
        set -e
        set -o pipefail

        # Step A: Download
        echo ">> Downloading from SRA..."
        module purge
        module load SRA-Toolkit/3.2.0-gompi-2024a
        prefetch ${SRR} -O raw_reads/
        fasterq-dump raw_reads/${SRR} --split-files --outdir raw_reads --threads 4

        if [ ! -s "raw_reads/${SRR}_1.fastq" ] || [ ! -s "raw_reads/${SRR}_2.fastq" ]; then
            echo ">> ERROR: download for ${SRR} produced empty/missing FASTQ files."
            exit 1
        fi

        # Step B: Raw QC
        echo ">> Running FastQC on raw reads..."
        module purge
        module load FastQC/0.11.9-Java-11
        fastqc raw_reads/${SRR}_1.fastq raw_reads/${SRR}_2.fastq -o qc -t 8 -q

        # Step C: Primer removal with cutadapt
        # -g / -G: primer anchored at the 5' end of R1/R2 respectively.
        # --discard-untrimmed: reads where the expected primer wasn't found
        #   are dropped, since their absence usually means the read pair is
        #   off-target (e.g. primer dimers, or non-amplicon contamination)
        #   rather than a true amplicon read with a sequencing error in the
        #   primer region.
        # NOTE: no fixed-length truncation happens here for either 16S or
        # ITS - that decision belongs in the downstream DADA2 step, where
        # it differs by AMPLICON_TYPE (see amplicon_dada2_analysis.R).
        echo ">> Removing primers with cutadapt..."
        module purge
        module load cutadapt/4.9-GCCcore-13.3.0   # <-- verify exact version with `module spider cutadapt`
        cutadapt \
            -g ${FWD_PRIMER} -G ${REV_PRIMER} \
            --discard-untrimmed \
            --minimum-length 50 \
            -o trimmed_reads/${SRR}_R1_trimmed.fastq -p trimmed_reads/${SRR}_R2_trimmed.fastq \
            raw_reads/${SRR}_1.fastq raw_reads/${SRR}_2.fastq \
            > qc/${SRR}_cutadapt.log

        if [ ! -s "trimmed_reads/${SRR}_R1_trimmed.fastq" ]; then
            echo ">> ERROR: cutadapt produced no surviving reads for ${SRR} - check primer sequences against this dataset's actual amplicon design."
            exit 1
        fi

        # Step D: Post-trim QC
        echo ">> Running FastQC on primer-trimmed reads..."
        module purge
        module load FastQC/0.11.9-Java-11
        fastqc trimmed_reads/${SRR}_R1_trimmed.fastq trimmed_reads/${SRR}_R2_trimmed.fastq -o qc -t 8 -q
    )
    SAMPLE_STATUS=$?
    set -e

    if [ ${SAMPLE_STATUS} -ne 0 ]; then
        echo ">> ERROR: Sample ${SRR} failed (exit code ${SAMPLE_STATUS}) - skipping to next sample."
        FAILED_SAMPLES+=("${SRR}")
        continue
    fi

    echo ">> Completed ${SRR}"
done

# ==========================================
# PHASE 3: GLOBAL REPORTING
# ==========================================
echo -e "\n[PHASE 3] COMPILING MULTIQC REPORT..."
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc/ -n final_amplicon_${AMPLICON_TYPE}_report.html

echo "=========================================================="
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo " PIPELINE COMPLETED WITH ${#FAILED_SAMPLES[@]} FAILED SAMPLE(S):"
    printf '   %s\n' "${FAILED_SAMPLES[@]}"
else
    echo " PIPELINE SUCCESSFULLY COMPLETED! "
fi
echo " Primer-trimmed reads: trimmed_reads/<SRR>_R1_trimmed.fastq / _R2_trimmed.fastq"
echo " Next step: run amplicon_dada2_analysis.R (in R/RStudio) with"
echo " AMPLICON_TYPE=\"${AMPLICON_TYPE}\" against this trimmed_reads/ directory."
echo "=========================================================="

###############################################################################
# METHODOLOGY NOTES
# ---------------------------------------------------------------------------
# 1. Why cutadapt instead of fastp for this step: fastp trims by quality/
#    adapter heuristics, not by matching a specific known primer sequence.
#    Amplicon primer removal needs to anchor on the exact expected primer,
#    which is what cutadapt's -g/-G does, and --discard-untrimmed gives an
#    explicit signal (dropped read count in the log) of how well the
#    primers actually matched this dataset - a useful QC check on its own.
# 2. Degenerate bases (M, N, V, W, Y) in the default primers above are
#    standard IUPAC ambiguity codes (e.g. EMP 515F/806R legitimately contain
#    them to match sequence variation across bacterial/archaeal taxa) -
#    cutadapt handles these natively, this is not a typo.
# 3. This script intentionally stops at primer-trimmed FASTQs. Quality-based
#    length filtering and truncation happen in DADA2's filterAndTrim(),
#    because the correct truncation length differs by AMPLICON_TYPE (16S:
#    fixed; ITS: none) and by your specific run's actual read-quality
#    profile - decisions best made after inspecting per-sample quality
#    profiles, which DADA2's plotQualityProfile() supports directly.
###############################################################################

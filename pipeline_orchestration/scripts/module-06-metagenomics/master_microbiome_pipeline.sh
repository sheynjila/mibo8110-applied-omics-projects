#!/bin/bash
#SBATCH --job-name=master_microbiome
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=master_microbiome_%j.out
#SBATCH --error=master_microbiome_%j.err

# ==============================================================================
# CORRECTED VERSION of master_microbiome_pipeline.sh
# See pipeline_appraisal_2026-08-31.md for the full review. Fixes applied:
#
#   1. Bracken's -r (read length) is no longer hardcoded to 150. It's now
#      measured directly from each sample's own trimmed FASTQ (fastp can
#      change read lengths via quality/adapter trimming, so a fixed 150
#      assumption drifts from reality per-sample and biases abundance
#      re-estimation). Falls back to 150 only if length can't be measured.
#   2. The per-sample loop is fault-tolerant: a failed download, Kraken2 run,
#      or Bracken run is logged and skipped instead of aborting the whole
#      job, and the aggregation step only includes samples that succeeded.
#   3. [Added per domain_coverage_audit_2026-08-31.md] Optional host-read
#      depletion step. Plant- and animal-associated microbiome samples
#      (gut, tissue, rhizosphere, etc.) typically contain substantial host
#      DNA. Left unfiltered, host reads waste Kraken2 compute and can
#      distort relative-abundance estimates if any host sequence happens to
#      match a microbial reference. Set HOST_REFERENCE to a host genome
#      FASTA to align trimmed reads with Bowtie2 and classify only the
#      reads that do NOT map to the host. Leave HOST_REFERENCE empty (the
#      default) to reproduce the original pipeline's behavior exactly -
#      this feature is fully backward-compatible and off by default.
# ==============================================================================

# ORCHESTRATION COPY -- see pipeline_orchestration/README.md. Differs from
# the original module script only by RUN_TAG-suffixed WORKDIR and
# PIPELINE_SCRIPT_DIR (so it reuses that module's load_modules.sh).

set -e
set -o pipefail

RUN_TAG="${RUN_TAG:-$(date -u +%Y%m%dT%H%M%SZ)}"
export PIPELINE_SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../../../modules/module-06-metagenomics/scripts" && pwd)}"
source "${PIPELINE_SCRIPT_DIR}/load_modules.sh"

# Path to a host reference genome FASTA (e.g. the plant or animal host of
# your microbiome samples). Leave empty to skip host depletion entirely.
HOST_REFERENCE="${HOST_REFERENCE:-}"

echo "=========================================================="
echo " INITIATING MICROBIOME PIPELINE (Kraken2 + Bracken) "
echo "=========================================================="

# ==========================================
# PHASE 1: DIRECTORY SETUP & DB CONFIG
# ==========================================
WORKDIR="/scratch/$(whoami)/master_microbiome_pipeline_${RUN_TAG}"
mkdir -p ${WORKDIR}/{raw_reads,trimmed_reads,host_depleted_reads,kraken_output,bracken_output,qc,logs}
cd ${WORKDIR}

if [ ! -f "srr_list.txt" ]; then
    echo "WARNING: srr_list.txt not found. Generating default test list..."
    echo -e "SRR11092056\nSRR11092057" > srr_list.txt
fi
sed -i 's/\r$//' srr_list.txt

# IMPORTANT: Kraken2 requires a massive pre-built database (>50GB).
# Update this path to wherever your HPC admins store the shared Kraken2 database.
KRAKEN_DB="/usr/local/apps/kraken2_db/standard"

if [ ! -d "${KRAKEN_DB}" ]; then
    echo "CRITICAL ERROR: Kraken2 Database not found at ${KRAKEN_DB}."
    echo "Please check your HPC documentation for the correct path."
    exit 1
fi

# FIX #3: If host depletion is requested, build the Bowtie2 index ONCE up
# front (not per-sample - indexing a host genome is expensive and identical
# across all samples in a run).
HOST_INDEX=""
if [ -n "${HOST_REFERENCE}" ]; then
    if [ ! -s "${HOST_REFERENCE}" ]; then
        echo "CRITICAL ERROR: HOST_REFERENCE set to '${HOST_REFERENCE}' but file not found."
        exit 1
    fi
    echo "Host depletion enabled. Building Bowtie2 index for host reference..."
    module purge
    load_bowtie2 || exit 1
    mkdir -p host_index
    HOST_INDEX="host_index/host_ref"
    if [ ! -f "${HOST_INDEX}.1.bt2" ]; then
        bowtie2-build --threads 8 ${HOST_REFERENCE} ${HOST_INDEX}
    fi
else
    echo "Host depletion disabled (HOST_REFERENCE not set) - classifying all trimmed reads directly, as in the original pipeline."
fi

# ==========================================
# PHASE 2: TAXONOMIC CLASSIFICATION LOOP
# ==========================================
echo -e "\n[PHASE 2] EXECUTING METAGENOMIC PROFILING..."

# FIX #2: track failures instead of a hard abort on the first bad sample.
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
        load_sra_toolkit || exit 1
        prefetch ${SRR} -O raw_reads/
        fasterq-dump raw_reads/${SRR} --split-files --outdir raw_reads --threads 4

        # Step B: Quality Control & Trimming
        echo ">> Trimming with fastp..."
        module purge
        load_fastp || exit 1
        fastp -i raw_reads/${SRR}_1.fastq -I raw_reads/${SRR}_2.fastq \
              -o trimmed_reads/${SRR}_1_clean.fastq -O trimmed_reads/${SRR}_2_clean.fastq \
              --thread 8 --html qc/${SRR}_fastp.html 2> qc/${SRR}_fastp.log

        # FIX #1: measure this sample's actual post-trim read length instead
        # of assuming 150bp. Average length over the first 40,000 reads
        # (line 2 of every 4-line FASTQ record) is enough to characterize a
        # trimmed run without scanning the entire file.
        READ_LEN=$(awk 'NR%4==2 {sum+=length($0); n++} n>=40000{exit} END{if (n>0) print int(sum/n); else print 150}' \
                   trimmed_reads/${SRR}_1_clean.fastq)
        if [ -z "${READ_LEN}" ] || [ "${READ_LEN}" -le 0 ]; then
            echo ">> WARNING: could not measure read length for ${SRR}, defaulting to 150."
            READ_LEN=150
        fi
        echo ">> Using measured average read length ${READ_LEN}bp for Bracken (-r)."

        # Step B.5 (NEW): Optional host-read depletion.
        # Aligns trimmed reads to the host genome with Bowtie2 and retains
        # only read pairs where NEITHER mate maps to the host (--un-conc),
        # before those reads ever reach Kraken2. Skipped entirely (reads
        # passed through unchanged) when HOST_REFERENCE is not set.
        CLASSIFY_R1="trimmed_reads/${SRR}_1_clean.fastq"
        CLASSIFY_R2="trimmed_reads/${SRR}_2_clean.fastq"
        if [ -n "${HOST_INDEX}" ]; then
            echo ">> Depleting host reads with Bowtie2..."
            module purge
            load_bowtie2 || exit 1
            bowtie2 -x ${HOST_INDEX} \
                    -1 trimmed_reads/${SRR}_1_clean.fastq -2 trimmed_reads/${SRR}_2_clean.fastq \
                    --threads 16 \
                    --un-conc host_depleted_reads/${SRR}_nonhost_%.fastq \
                    -S /dev/null \
                    2> qc/${SRR}_host_depletion.log
            # Bowtie2's --un-conc names mates %1/%2; normalize to _1/_2 for
            # consistency with the rest of the pipeline's naming.
            mv host_depleted_reads/${SRR}_nonhost_1.fastq host_depleted_reads/${SRR}_nonhost_R1.fastq
            mv host_depleted_reads/${SRR}_nonhost_2.fastq host_depleted_reads/${SRR}_nonhost_R2.fastq
            CLASSIFY_R1="host_depleted_reads/${SRR}_nonhost_R1.fastq"
            CLASSIFY_R2="host_depleted_reads/${SRR}_nonhost_R2.fastq"
            HOST_PCT=$(awk '/overall alignment rate/{print $1}' qc/${SRR}_host_depletion.log)
            echo ">> Host alignment rate for ${SRR}: ${HOST_PCT:-unknown} (these reads were removed before classification)."
        fi

        # Step C: Taxonomic Classification (Kraken2)
        echo ">> Running Kraken2 Taxonomic Classification..."
        module purge
        load_kraken2 || exit 1
        kraken2 --db ${KRAKEN_DB} \
                --paired ${CLASSIFY_R1} ${CLASSIFY_R2} \
                --threads 16 \
                --report kraken_output/${SRR}_report.txt \
                --output kraken_output/${SRR}_kraken.out

        # Step D: Abundance Re-estimation (Bracken)
        echo ">> Running Bracken for Species Abundance Estimation..."
        module purge
        load_bracken || exit 1
        # -r now uses this sample's measured read length; -l S = Species level.
        bracken -d ${KRAKEN_DB} \
                -i kraken_output/${SRR}_report.txt \
                -o bracken_output/${SRR}_bracken.tsv \
                -r ${READ_LEN} -l S
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
# PHASE 3: AGGREGATE RESULTS FOR R
# ==========================================
echo -e "\n[PHASE 3] AGGREGATING R-READY FINAL REPORT..."
cd bracken_output

BRACKEN_TSVS=$(ls *_bracken.tsv 2>/dev/null || true)

if [ -z "${BRACKEN_TSVS}" ]; then
    echo "WARNING: No Bracken results were produced - skipping aggregation."
else
    echo -e "Sample\tTaxon\tEstimated_Reads\tFraction_Total" > ALL_SAMPLES_MICROBIOME.tsv
    for TSV in ${BRACKEN_TSVS}; do
        BASENAME=$(basename ${TSV} _bracken.tsv)
        awk -v sample="${BASENAME}" 'FS="\t" {if (NR!=1) print sample "\t" $1 "\t" $6 "\t" $7}' ${TSV} >> ALL_SAMPLES_MICROBIOME.tsv
    done
fi
cd ${WORKDIR}

# ==========================================
# PHASE 4: GLOBAL REPORTING
# ==========================================
echo -e "\n[PHASE 4] COMPILING MULTIQC REPORT..."
module purge
load_python  || exit 1   # must be loaded before MultiQC -- see load_modules.sh
load_multiqc || exit 1
multiqc qc/ kraken_output/ -n final_microbiome_report.html

echo "=========================================================="
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo " PIPELINE COMPLETED WITH ${#FAILED_SAMPLES[@]} FAILED SAMPLE(S):"
    printf '   %s\n' "${FAILED_SAMPLES[@]}"
else
    echo " PIPELINE SUCCESSFULLY COMPLETED! "
fi
echo " R-ready Microbiome Data: bracken_output/ALL_SAMPLES_MICROBIOME.tsv"
echo "=========================================================="

###############################################################################
# EDUCATIONAL COMMENTARY & METHODOLOGY (carried over)
#
# 1. Metagenomics vs. Isolate Sequencing
#    Unlike pipelines that map reads from a pure, isolated bacterium to a
#    single reference genome, this pipeline analyzes mixed communities
#    (e.g., soil, gut, or water samples).
#
# 2. Kraken2 (Taxonomic Assignment)
#    Kraken2 breaks reads into k-mers and compares them against a database
#    of known genomes to assign a taxonomic label to each read. Accuracy
#    depends heavily on database completeness and compatibility.
#
# 3. Bracken (Abundance Estimation)
#    Bracken mathematically re-estimates Kraken2's classifications to give a
#    more accurate relative abundance per sample. Its -r parameter should
#    match the actual (post-trim) read length used for classification.
#
# 4. R Integration & Limitations
#    ALL_SAMPLES_MICROBIOME.tsv is formatted for R (e.g. stacked bar charts
#    of species composition). Relative read abundance can still be biased by
#    extraction efficiency and genome size, and is not always exactly equal
#    to true biological abundance.
#
# 5. Host-read depletion (optional, HOST_REFERENCE)
#    For plant- or animal-associated samples, host DNA can be a large share
#    of total reads. Removing it before classification (a) prevents wasted
#    Kraken2 compute on non-microbial reads and (b) avoids any host
#    sequence spuriously matching a database entry. The reported
#    "overall alignment rate" in qc/<SRR>_host_depletion.log is the
#    fraction of reads identified as host and removed - a useful QC metric
#    to report alongside your abundance results (e.g. "82% of reads were
#    host-derived and excluded prior to taxonomic classification").
###############################################################################

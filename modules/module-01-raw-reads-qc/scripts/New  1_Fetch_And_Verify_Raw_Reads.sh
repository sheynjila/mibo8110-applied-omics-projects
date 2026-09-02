#!/bin/bash
#SBATCH --job-name=qc01_fetch
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00

###############################################################################
# 1_Fetch_And_Verify_Raw_Reads.sh
# Module 1 (Raw Sequencing Reads & QC) — Script 1 of 7
#
# ONE NEW IDEA over a blank slate: before you can ask "are these reads any
# good?" you first have to prove you actually HAVE the reads you think you
# have. This script retrieves ONE hardcoded paired-end run, validates the
# SRA cache before extracting it, and records checksums of the resulting
# FASTQ files — the first entry in the audit trail the module's portfolio
# artifact (a QC report) will eventually require.
#
# Study used throughout this module (same one used in Module 0 and Module 2,
# for continuity): PRJNA1518998, "Effect of adipocyte depletion on the
# retina," Mus musculus, Vanderbilt University Medical Center
# (GEO accession GSE345194). Real paired-end run used here: SRR40359383.
###############################################################################

set -e
set -o pipefail

module load SRA-Toolkit/3.0.3-gompi-2022a

cd /scratch/$(whoami)
mkdir -p module1_qc/{raw_reads,checksums,logs}
cd module1_qc

SRR="SRR40359383"
LOG="logs/${SRR}_fetch.log"

echo "==========================================================" | tee "${LOG}"
echo " Fetching and verifying: ${SRR}" | tee -a "${LOG}"
echo "==========================================================" | tee -a "${LOG}"

# ------------------------------------------------------------------------------
# STEP 1: prefetch into the SRA cache, THEN validate the cache, BEFORE paying
# the cost of fasterq-dump extracting it to FASTQ.
#
# TEACHING NOTE: prefetch can complete "successfully" while still leaving a
# truncated or corrupted .sra file behind (interrupted connection, disk
# quota hit mid-download, etc.). vdb-validate catches that corruption at the
# CACHE stage, which is cheap to re-run. If you skip straight to
# fasterq-dump on a bad cache, you either get a cryptic extraction error, or
# — worse — a truncated FASTQ that LOOKS complete and silently poisons every
# downstream QC metric in this module.
# ------------------------------------------------------------------------------
prefetch "${SRR}" --output-directory sra_cache 2>&1 | tee -a "${LOG}"

if ! vdb-validate "sra_cache/${SRR}/${SRR}.sra" 2>&1 | tee -a "${LOG}" | grep -q "is consistent"; then
    echo "ERROR: vdb-validate did not confirm ${SRR}.sra is consistent. Stopping before extraction." | tee -a "${LOG}"
    exit 1
fi
echo "OK: SRA cache for ${SRR} validated." | tee -a "${LOG}"

# ------------------------------------------------------------------------------
# STEP 2: extract to paired FASTQ, then confirm BOTH mate files exist and are
# non-empty. A partial extraction (disk full, killed job) can leave one mate
# file present and the other missing or zero-byte — checking file EXISTENCE
# alone is not enough, hence the `-s` (non-empty) test on both mates.
# ------------------------------------------------------------------------------
fasterq-dump "sra_cache/${SRR}/${SRR}.sra" --split-files --outdir raw_reads --threads 4 2>&1 | tee -a "${LOG}"

if [ ! -s "raw_reads/${SRR}_1.fastq" ] || [ ! -s "raw_reads/${SRR}_2.fastq" ]; then
    echo "ERROR: Extraction for ${SRR} did not produce two non-empty mate files. Stopping." | tee -a "${LOG}"
    exit 1
fi
echo "OK: raw_reads/${SRR}_1.fastq and _2.fastq both present and non-empty." | tee -a "${LOG}"

# ------------------------------------------------------------------------------
# STEP 3: checksum the FASTQ files. This is the "data source and checksums"
# line item the curriculum's portfolio artifact explicitly requires — write it
# now, at the moment of retrieval, not retroactively from memory later.
# ------------------------------------------------------------------------------
md5sum "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" > "checksums/${SRR}_raw.md5"
echo "OK: checksums written to checksums/${SRR}_raw.md5" | tee -a "${LOG}"
cat "checksums/${SRR}_raw.md5" | tee -a "${LOG}"

echo "Done: ${SRR} fetched and verified." | tee -a "${LOG}"


# ==============================================================================
# FASTQ record anatomy (read this before script 2)
# ==============================================================================
# Every read in a FASTQ file is exactly 4 lines:
#   1. @<read ID> <optional metadata>       — starts with "@"
#   2. <the sequence itself>                — A/C/G/T/N
#   3. "+" (optionally repeating the read ID again)
#   4. <quality string>                     — SAME LENGTH as line 2
#
# Line 4 is not raw probability — it is a single ASCII character per base,
# offset by 33 (Phred+33 encoding). "!" (ASCII 33) = Phred quality 0 (worst),
# "I" (ASCII 73) = Phred quality 40 (excellent). Phred quality Q relates to
# the PROBABILITY that base call is wrong: P(error) = 10^(-Q/10). So Q30
# means a 1-in-1000 chance that base is wrong; Q20 means 1-in-100. This is
# exactly what FastQC's "per base sequence quality" plot in script 2 is
# showing you, one position at a time, averaged across every read.
#
# WHY VERIFY BEFORE INSPECTING (nuances & pitfalls):
#   - `prefetch` retries transient network failures internally, but a
#     completed download is not the same as a VALID one — always validate.
#   - Checking file existence (`-f`) is not enough after `fasterq-dump`;
#     a killed job mid-extraction can leave a 0-byte file that still "exists".
#     Use `-s` (exists AND non-empty) as the minimum bar.
#   - Checksums recorded at fetch time are cheap insurance: if a later script
#     in this module reports a metric you don't trust, comparing the current
#     file's md5 against checksums/${SRR}_raw.md5 tells you in one command
#     whether the underlying data changed since retrieval.
# ==============================================================================

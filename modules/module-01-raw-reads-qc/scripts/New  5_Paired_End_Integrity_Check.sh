#!/bin/bash
#SBATCH --job-name=qc05_pairing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00

###############################################################################
# 5_Paired_End_Integrity_Check.sh
# Module 1 (Raw Sequencing Reads & QC) — Script 5 of 7
#
# ONE NEW IDEA on top of script 4: matching READ COUNTS between mate files is
# necessary but NOT sufficient. Two mate files can have identical read counts
# while being desynchronized — read 5000 in file 1 no longer corresponds to
# read 5000 in file 2 — if any tool in the chain drops or reorders reads
# asymmetrically. This script checks BOTH conditions, on BOTH the raw and the
# trimmed read sets, directly answering the curriculum's learning outcome:
# "process paired-end reads without breaking read pairing."
###############################################################################

set -e
set -o pipefail

cd /scratch/$(whoami)/module1_qc

SRR="SRR40359383"

check_pair() {
    local LABEL="$1"
    local MATE1="$2"
    local MATE2="$3"

    if [ ! -s "${MATE1}" ] || [ ! -s "${MATE2}" ]; then
        echo "  ERROR [${LABEL}]: one or both mate files missing/empty (${MATE1}, ${MATE2})."
        return 1
    fi

    # ------------------------------------------------------------------------
    # CHECK 1: read counts match. FASTQ records are exactly 4 lines each, so
    # (line count / 4) gives the read count without needing a FASTQ parser.
    # ------------------------------------------------------------------------
    local COUNT1 COUNT2
    COUNT1=$(( $(wc -l < "${MATE1}") / 4 ))
    COUNT2=$(( $(wc -l < "${MATE2}") / 4 ))

    if [ "${COUNT1}" -ne "${COUNT2}" ]; then
        echo "  FAIL [${LABEL}]: read count mismatch — mate 1 has ${COUNT1} reads, mate 2 has ${COUNT2}."
        return 1
    fi
    echo "  OK [${LABEL}]: read counts match (${COUNT1} reads per mate)."

    # ------------------------------------------------------------------------
    # CHECK 2: read ID ORDER matches, not just the count. Extract every read
    # ID (line 1 of every 4-line record), strip the trailing /1 or /2 (or
    # " 1:N:..." style suffix) mate marker, and diff the two ID lists
    # positionally. A mismatch here means the mates are the same SIZE but no
    # longer point at the same underlying fragments — a subtler, more
    # dangerous failure than a count mismatch because nothing about the file
    # sizes looks wrong.
    # ------------------------------------------------------------------------
    local IDS1 IDS2
    IDS1=$(mktemp)
    IDS2=$(mktemp)
    awk 'NR%4==1' "${MATE1}" | sed -E 's#[ /][12]?(:.*)?$##' > "${IDS1}"
    awk 'NR%4==1' "${MATE2}" | sed -E 's#[ /][12]?(:.*)?$##' > "${IDS2}"

    if ! diff -q "${IDS1}" "${IDS2}" > /dev/null; then
        echo "  FAIL [${LABEL}]: read ID order does not match between mates (counts matched, but pairing is broken)."
        rm -f "${IDS1}" "${IDS2}"
        return 1
    fi
    echo "  OK [${LABEL}]: read IDs match in order between mates — pairing intact."
    rm -f "${IDS1}" "${IDS2}"
    return 0
}

FAILED=0

echo "Checking raw reads for ${SRR}..."
check_pair "raw" "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" || FAILED=1

echo "Checking trimmed reads for ${SRR}..."
if [ -s "trimmed_reads/${SRR}_1_clean.fastq" ]; then
    check_pair "trimmed" "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq" || FAILED=1
else
    echo "  SKIPPED: no trimmed output found yet — run script 4 first if you want this check."
fi

if [ "${FAILED}" -eq 1 ]; then
    echo ""
    echo "RESULT: one or more pairing checks FAILED for ${SRR}. Do not proceed to alignment on this sample."
    exit 1
fi

echo ""
echo "RESULT: ${SRR} passed all pairing checks (raw and trimmed)."


# ==============================================================================
# Nuances & Pitfalls
# ==============================================================================
# - fastp (used correctly, with -i/-I and -o/-O both mates passed together as
#   in script 4) preserves pairing by design — it drops or keeps BOTH mates of
#   a pair together. Pairing breaks are much more common when reads are
#   filtered with single-end tools applied independently to each mate file,
#   or when mate files are concatenated/subset with plain `grep`/`head`
#   without accounting for the 4-line record boundary.
# - The read-ID suffix stripping regex here (`[ /][12]?(:.*)?$`) handles both
#   common SRA-style ID conventions (`.../1` `.../2` suffix, and modern
#   Illumina-style ` 1:N:...` / ` 2:N:...` space-separated suffix). If your
#   own data uses a different ID convention, verify this pattern actually
#   strips ONLY the mate marker — an overly greedy pattern can make two
#   genuinely different reads look identical, silently hiding a real pairing
#   break.
# - This check is fast (pure text processing, no aligner needed) precisely so
#   it can run as a cheap gate BEFORE any expensive alignment step — that is
#   exactly how it is used inside script 6's integrated pipeline.
# ==============================================================================

#!/bin/bash
#SBATCH --job-name=qc05_pairing
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00

###############################################################################
# 5_Paired_End_Integrity_Check.sh
#
# NARRATIVE:
# This script acts as a strict defensive gatekeeper before alignment[cite: 9]. Poorly 
# configured trimming can silently drop single reads, breaking the synchronization 
# of paired-end FASTQ mates[cite: 9]. This script extracts and diffs the actual read IDs 
# sequentially (ignoring mate suffixes) to prove that file 1 and file 2 remain 
# perfectly 1:1 synchronized, preventing cryptic aligner crashes downstream[cite: 9].
###############################################################################
set -e
set -o pipefail

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep"
cd "${WORKDIR}"

SRR="SRR1039508"

check_pair() {
    local MATE1="$1" MATE2="$2"
    local C1=$(( $(wc -l < "${MATE1}") / 4 ))
    local C2=$(( $(wc -l < "${MATE2}") / 4 ))

    if [ "${C1}" -ne "${C2}" ]; then
        echo "FAIL: Read count mismatch (${C1} vs ${C2})."
        return 1
    fi

    local IDS1=$(mktemp); local IDS2=$(mktemp)
    awk 'NR%4==1' "${MATE1}" | sed -E 's#[ /][12]?(:.*)?$##' > "${IDS1}"
    awk 'NR%4==1' "${MATE2}" | sed -E 's#[ /][12]?(:.*)?$##' > "${IDS2}"

    if ! diff -q "${IDS1}" "${IDS2}" > /dev/null; then
        echo "FAIL: Read ID order mismatch (pairing broken)."
        rm -f "${IDS1}" "${IDS2}"; return 1
    fi
    echo "OK: Pairing intact."
    rm -f "${IDS1}" "${IDS2}"; return 0
}

echo "Checking trimmed reads for ${SRR}..."
check_pair "trimmed_reads/${SRR}_1_clean.fastq" "trimmed_reads/${SRR}_2_clean.fastq"
#!/bin/bash
#SBATCH --job-name=qc01_fetch
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00

###############################################################################
# 1_Fetch_And_Verify_Raw_Reads.sh
# 
# NARRATIVE:
# This script initiates the data acquisition phase for a single sample. 
# To enforce data integrity, it does not blindly extract reads. Instead, it 
# first pre-fetches the SRA file, strictly validates the local cache to detect 
# network corruption, and only then executes `fasterq-dump` to generate 
# paired-end FASTQ files. Finally, it computes MD5 checksums to establish a 
# cryptographic baseline for downstream auditing.
###############################################################################
set -e
set -o pipefail

module load SRA-Toolkit/3.0.3-gompi-2022a

WORKDIR="/scratch/$(whoami)/PRJNA229998_airway_pipeline_stepbystep"
mkdir -p "${WORKDIR}"/{raw_reads,checksums,logs,sra_cache}
cd "${WORKDIR}"

SRR="SRR1039508"
LOG="logs/${SRR}_fetch.log"

echo "Fetching and verifying: ${SRR}..." | tee "${LOG}"

prefetch "${SRR}" --output-directory sra_cache 2>&1 | tee -a "${LOG}"

if ! vdb-validate "sra_cache/${SRR}/${SRR}.sra" 2>&1 | tee -a "${LOG}" | grep -q "is consistent"; then
    echo "ERROR: vdb-validate failed. Stopping before extraction." | tee -a "${LOG}"
    exit 1
fi

fasterq-dump "sra_cache/${SRR}/${SRR}.sra" --split-files --outdir raw_reads --threads 4 2>&1 | tee -a "${LOG}"

md5sum "raw_reads/${SRR}_1.fastq" "raw_reads/${SRR}_2.fastq" > "checksums/${SRR}_raw.md5"
echo "OK: ${SRR} fetched and verified. Checksums saved." | tee -a "${LOG}"
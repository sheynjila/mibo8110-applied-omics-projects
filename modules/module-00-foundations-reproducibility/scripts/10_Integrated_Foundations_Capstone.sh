#!/bin/bash
#SBATCH --job-name=module0_capstone
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=24:00:00

###############################################################################
# 10_Integrated_Foundations_Capstone.sh
#
# THIS is the Module 0 portfolio artifact the curriculum asks for: "a
# repository containing a small pipeline, environment specification, README,
# input checks, log, outputs, and a reproducibility reflection."
#
# NEW IDEA on top of scripts 1-9: nothing here is a new technique — every
# habit below was already taught individually. The new idea is INTEGRATION:
# a single run that visibly uses the dynamic storage guard (script 5), logs
# and predicts memory (script 6), records its software environment (script
# 7), stamps itself with the exact Git commit that produced it (script 8),
# and validates its own input AND output file formats (script 9) — then
# writes a reflection document summarizing all of it. This is the shape every
# capstone pipeline in this project (AMR, RNA-seq, SNP calling, microbiome)
# already takes; this script just makes that shape explicit for the first
# time, at a small scale, using the same PRJNA1518998 mouse retina samples
# referenced throughout this module.
#
# Prerequisites (run these once, on the login node, before submitting this
# job with sbatch):
#   bash 7_Environment_Reproducibility.sh create    # or export, if env exists
#   bash 8_Git_Reproducibility_Workflow.sh init
#   bash 8_Git_Reproducibility_Workflow.sh snapshot
#   # a working srr_list.txt in the working directory (see module README)
###############################################################################

SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/modules/module-00-foundations-reproducibility/scripts}}"
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
source "${SCRIPT_DIR}/load_modules.sh"

module purge
load_sra_toolkit || exit 1
load_fastqc      || exit 1
load_fastp       || exit 1
load_multiqc     || exit 1
load_samtools    || exit 1

cd /scratch/$(whoami)
mkdir -p capstone/module0/{raw_reads,qc_before,trimmed_reads,qc_after,fastp_reports,logs}
cd capstone/module0

RUN_LOG="logs/run_$(date -u '+%Y%m%dT%H%M%SZ').log"
MEM_LOG="mem_usage_log.txt"
SAFE_BUFFER_GB=20
SRR_LIST="srr_list.txt"

# ----------------------------------------------------------------------------
# Step 0: input checks (habit from scripts 2-6, applied first as always)
# ----------------------------------------------------------------------------
if [ ! -s "$SRR_LIST" ]; then
    echo "ERROR: ${SRR_LIST} not found or empty in $(pwd)." | tee -a "$RUN_LOG"
    echo "This module's srr_list.txt should contain PRJNA1518998 accessions," | tee -a "$RUN_LOG"
    echo "e.g. SRR40359383 and SRR40359384 (one per line)." | tee -a "$RUN_LOG"
    exit 1
fi

{
    echo "=============================================================="
    echo "Module 0 capstone run starting $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "=============================================================="
} | tee -a "$RUN_LOG"

# ----------------------------------------------------------------------------
# Step 1: provenance stamp — WHICH CODE (script 8) ran this job.
# Run in "stamp" mode from the parent scripts/ directory so it reports the
# state of the Git repository the scripts themselves live in.
# ----------------------------------------------------------------------------
GIT_STAMP_SCRIPT="../../../modules/module-00-foundations-reproducibility/scripts/8_Git_Reproducibility_Workflow.sh"
{
    echo "--- Provenance (script 8) ---"
    if [ -f "$GIT_STAMP_SCRIPT" ]; then
        bash "$GIT_STAMP_SCRIPT" stamp
    else
        echo "git_commit=UNAVAILABLE (8_Git_Reproducibility_Workflow.sh not found at expected path)"
    fi
} | tee -a "$RUN_LOG"

# ----------------------------------------------------------------------------
# Step 2: environment fingerprint — WHICH TOOLS (script 7) ran this job.
# We don't require a Conda environment to be active on the compute node
# (module-loaded tools are also valid), but we always record what actually
# ran, the same "measure, don't assume" instinct as everywhere else.
# ----------------------------------------------------------------------------
{
    echo "--- Environment fingerprint (script 7) ---"
    echo "fasterq-dump: $(fasterq-dump --version 2>&1 | head -1)"
    echo "fastqc:       $(fastqc --version 2>&1)"
    echo "fastp:        $(fastp --version 2>&1)"
    echo "samtools:     $(samtools --version 2>&1 | head -1)"
    [ -f ../../../modules/module-00-foundations-reproducibility/scripts/environment.yml ] && \
        echo "environment.yml present: yes" || echo "environment.yml present: no (run script 7 export/create)"
} | tee -a "$RUN_LOG"

# ----------------------------------------------------------------------------
# Step 3: the actual mini-pipeline — same fault-tolerant loop, storage guard,
# and format checks taught in scripts 2, 5, and 9.
# ----------------------------------------------------------------------------
FAILED_LOG="failed_downloads.txt"
> "$FAILED_LOG"

while IFS= read -r SRR || [ -n "$SRR" ]; do
    SRR="${SRR%$'\r'}"
    [ -z "$SRR" ] && continue

    AVAILABLE_GB=$(df -BG . | awk 'NR==2 {print $4}' | tr -d 'G')
    if [[ "$AVAILABLE_GB" =~ ^[0-9]+$ ]] && [ "$AVAILABLE_GB" -le "$SAFE_BUFFER_GB" ]; then
        echo "CRITICAL: only ${AVAILABLE_GB}GB left. Halting the loop." | tee -a "$RUN_LOG"
        break
    fi

    echo "Processing ${SRR} (${AVAILABLE_GB}GB free)..." | tee -a "$RUN_LOG"
    fasterq-dump "${SRR}" --split-files --outdir raw_reads --threads 4

    if [ ! -s "raw_reads/${SRR}_1.fastq" ] || [ ! -s "raw_reads/${SRR}_2.fastq" ]; then
        echo "WARNING: download for ${SRR} failed or incomplete. Skipping." | tee -a "$RUN_LOG"
        echo "${SRR}" >> "$FAILED_LOG"
        continue
    fi

    fastqc raw_reads/${SRR}_1.fastq raw_reads/${SRR}_2.fastq -o qc_before -t 4
    fastp -i raw_reads/${SRR}_1.fastq -I raw_reads/${SRR}_2.fastq \
          -o trimmed_reads/${SRR}_1_clean.fastq -O trimmed_reads/${SRR}_2_clean.fastq \
          --thread 4 --html fastp_reports/${SRR}_fastp.html
    fastqc trimmed_reads/${SRR}_1_clean.fastq trimmed_reads/${SRR}_2_clean.fastq -o qc_after -t 4
done < "$SRR_LIST"

multiqc qc_before/ -n multiqc_raw_final.html
multiqc qc_after/ -n multiqc_trimmed_final.html

# ----------------------------------------------------------------------------
# Step 4: output format validation (script 9) — trust the pipeline's own
# outputs the same way input files were checked, not just the input.
# ----------------------------------------------------------------------------
FORMAT_CHECK_SCRIPT="../../../modules/module-00-foundations-reproducibility/scripts/9_File_Format_Validation.sh"
{
    echo "--- Output format validation (script 9) ---"
    if [ -f "$FORMAT_CHECK_SCRIPT" ]; then
        bash "$FORMAT_CHECK_SCRIPT" trimmed_reads
    else
        echo "SKIPPED: 9_File_Format_Validation.sh not found at expected path."
    fi
} | tee -a "$RUN_LOG"

# ----------------------------------------------------------------------------
# Step 5: memory logging (script 6) — last step, same as script 6.
# ----------------------------------------------------------------------------
{
    echo "--- Memory usage (script 6) ---"
    echo "Job ${SLURM_JOBID} finished $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    seff "$SLURM_JOBID" 2>&1
} | tee -a "$MEM_LOG" | tee -a "$RUN_LOG"

# ----------------------------------------------------------------------------
# Step 6: generate the reproducibility reflection template. This is the
# curriculum's exact "reproducibility reflection" deliverable — pre-filled
# with everything this run can answer automatically, leaving only the
# genuinely reflective questions for the student to answer by hand.
# ----------------------------------------------------------------------------
REFLECTION="REPRODUCIBILITY_REFLECTION.md"
{
    echo "# Reproducibility Reflection — Module 0 Capstone"
    echo ""
    echo "**Run date (UTC):** $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "**Slurm job ID:** ${SLURM_JOBID:-N/A}"
    echo "**Code version:** $(bash "$GIT_STAMP_SCRIPT" stamp 2>/dev/null | grep git_commit || echo 'git_commit=UNAVAILABLE')"
    echo "**Samples processed:** $(cat "$SRR_LIST" 2>/dev/null | tr '\n' ' ')"
    echo "**Failed downloads:** $([ -s "$FAILED_LOG" ] && cat "$FAILED_LOG" | tr '\n' ' ' || echo 'none')"
    echo ""
    echo "## Auto-recorded facts"
    echo "- Environment spec present: $([ -f ../../../modules/module-00-foundations-reproducibility/scripts/environment.yml ] && echo yes || echo no)"
    echo "- Run log: \`${RUN_LOG}\`"
    echo "- Memory log: \`${MEM_LOG}\`"
    echo "- Output format validation: see 'Output format validation' section of \`${RUN_LOG}\`"
    echo ""
    echo "## Questions to answer by hand"
    echo "1. If a labmate cloned this repository today with no other context,"
    echo "   what is the SINGLE command they would run first, and how would they"
    echo "   know it worked?"
    echo "2. Which of this run's steps depended on something NOT captured in this"
    echo "   repository (e.g. a specific HPC's default \`module load\` version,"
    echo "   a manually-downloaded reference file, a value typed by hand)?"
    echo "3. Responsible practice: is PRJNA1518998 (or your chosen dataset) public"
    echo "   and licensed for reuse in coursework? Does this run output, log, or"
    echo "   commit any personally identifiable or otherwise sensitive information?"
    echo "4. If this run's numbers were challenged by a reviewer, which single"
    echo "   artifact in this repository would you point to first, and why?"
} > "$REFLECTION"

echo "Wrote reflection template to ${REFLECTION}." | tee -a "$RUN_LOG"
echo "Module 0 capstone run complete." | tee -a "$RUN_LOG"

###############################################################################
# NUANCES & PITFALLS
#
#   - This script does not introduce any single new technique — its entire
#     lesson is that reproducibility is the SUM of small habits, not one
#     big feature. Skipping any one piece (no commit stamp, no environment
#     record, no format check, no memory log) leaves a gap that quietly
#     breaks the "repeatable, portable, auditable" guiding question this
#     whole module is built around.
#
#   - The reflection template's auto-recorded facts are only as trustworthy
#     as the scripts that produced them — if script 8's commit was made
#     with UNCOMMITTED changes still present (see its own "dirty" warning),
#     the "Code version" line above is misleading, not wrong; students
#     should check `git_dirty` before trusting it.
#   - Path assumptions: the relative paths to scripts 7/8/9 above assume
#     this script is run from `capstone/module0/` under `/scratch/$(whoami)`,
#     three directories below this module's `scripts/` folder, matching
#     this module's own repository layout. Adjust the paths if you copy
#     this script somewhere else.
###############################################################################

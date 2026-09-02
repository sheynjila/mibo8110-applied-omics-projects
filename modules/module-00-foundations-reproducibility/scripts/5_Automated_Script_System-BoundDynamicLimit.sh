#!/bin/bash
#SBATCH --job-name=qc_system_limit
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out

###############################################################################
# Automated Script (System-Bound Dynamic Limit)
#
# Instead of tracking how much storage has already been consumed, this workflow
# queries the GACRC file system directly to determine how much free space
# remains on the target filesystem[cite: 7]. The script continues downloading and
# processing sequencing runs until the available storage falls below a defined
# safety threshold (for example, 20 GB)[cite: 7].
#
# Once the free space drops below this buffer, the workflow exits gracefully,
# preventing disk exhaustion and ensuring that all completed outputs remain
# valid and intact[cite: 7].
###############################################################################

# Set up directory structure
cd /scratch/$(whoami)
mkdir -p automated_system_bound_dynamicLimit/{raw_reads,qc_before,trimmed_reads,qc_after}
cd automated_system_bound_dynamicLimit

# We keep a safety buffer so the job doesn't completely crash the drive[cite: 7]
SAFE_BUFFER_GB=20

for SRR in $(cat srr_list.txt); do
    
    # 1. Ask the HPC file system exactly how much space is left[cite: 7]
    # df -BG checks the drive, awk isolates the "Available" column, tr removes the 'G'[cite: 7]
    AVAILABLE_GB=$(df -BG . | awk 'NR==2 {print $4}' | tr -d 'G')
    
    # 2. Check if the system is effectively full[cite: 7]
    if [ "$AVAILABLE_GB" -le "$SAFE_BUFFER_GB" ]; then
        echo "CRITICAL: The GACRC system only has ${AVAILABLE_GB}GB left."
        echo "Halting the loop to prevent job failure and file corruption."
        break
    fi

    echo "System has ${AVAILABLE_GB}GB available. Processing ${SRR}..."
    
    # 1. Retrieve the reads (Isolated Environment + Prefetch Fix)
    module purge
    module load SRA-Toolkit/3.2.0-gompi-2024a
    prefetch ${SRR} -O raw_reads/
    fasterq-dump raw_reads/${SRR} --split-files --outdir raw_reads --threads 4
    
    # 2. Initial FastQC (Isolated Environment)
    module purge
    module load FastQC/0.12.1-Java-11
    fastqc raw_reads/${SRR}_1.fastq raw_reads/${SRR}_2.fastq -o qc_before -t 4
    
    # 3. Evidence-based trimming (Isolated Environment)
    module purge
    module load fastp/0.23.4-GCC-13.2.0
    fastp -i raw_reads/${SRR}_1.fastq -I raw_reads/${SRR}_2.fastq -o trimmed_reads/${SRR}_1_clean.fastq -O trimmed_reads/${SRR}_2_clean.fastq --thread 4 --html ${SRR}_fastp.html
    
    # 4. Post-trimming FastQC (Isolated Environment)
    module purge
    module load FastQC/0.12.1-Java-11
    fastqc trimmed_reads/${SRR}_1_clean.fastq trimmed_reads/${SRR}_2_clean.fastq -o qc_after -t 4
    
done

echo "Compiling final MultiQC reports for all processed samples..."
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc_before/ -n multiqc_raw_final.html
multiqc qc_after/ -n multiqc_trimmed_final.html
echo "Workflow finished!"

###############################################################################
# Why this approach is safe
#
# Relying on a hard system quota or allowing the filesystem to fill completely
# can cause downloads to fail in the middle of writing a FASTQ file[cite: 7]. Partial
# downloads produce incomplete files that cannot be processed correctly by
# downstream tools such as fastp[cite: 7].
#
# By monitoring available disk space with the `df` command and maintaining a
# 20 GB reserve, the script can stop before storage becomes critically low[cite: 7].
# This avoids incomplete downloads, prevents downstream processing failures,
# and allows completed analyses and MultiQC report generation to finish
# successfully[cite: 7].
###############################################################################

# CREATING SRR list
# # 1. Create the working directory
#   mkdir -p /scratch/$(whoami)/automated_system_bound_dynamicLimit

# 2. Enter the directory
#    cd /scratch/$(whoami)/automated_system_bound_dynamicLimit

# 3. Generate the SRR list
#    echo -e "SRR11092056\nSRR11092057" > srr_list.txt

# 4. Clean formatting
#    dos2unix srr_list.txt

# 5. Return home and submit
#     cd ~
#     batch 5_Automated_Script_System-BoundDynamicLimit.sh
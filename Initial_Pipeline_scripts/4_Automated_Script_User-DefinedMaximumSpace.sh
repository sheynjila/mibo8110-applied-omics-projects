#!/bin/bash
#SBATCH --job-name=qc_custom_limit
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out

# ==========================================
# USER SETTINGS
# ==========================================
# Define your maximum allowed space in Gigabytes here:
MAX_ALLOWED_GB=50
# ==========================================

# Set up directory structure
cd /scratch/$(whoami)
mkdir -p automated_user_definedmaxspace/{raw_reads,qc_before,trimmed_reads,qc_after}
cd automated_user_definedmaxspace

# ==========================================
# MISSING FILE SAFETY NET
# ==========================================
# If srr_list.txt doesn't exist, create a default one to prevent crashes!
if [ ! -f "srr_list.txt" ]; then
    echo "WARNING: srr_list.txt not found. Generating a test list..."
    echo -e "SRR11092056\nSRR11092057" > srr_list.txt
    dos2unix srr_list.txt
fi
# ==========================================

for SRR in $(cat srr_list.txt); do
    # Check current directory size in Gigabytes
    CURRENT_GB=$(du -s -BG . | awk '{print $1}' | tr -d 'G')
    
    # Stop if we hit or exceed the limit
    if [ "$CURRENT_GB" -ge "$MAX_ALLOWED_GB" ]; then
        echo "WARNING: Folder size is ${CURRENT_GB}GB. Your limit is ${MAX_ALLOWED_GB}GB."
        echo "Halting the loop to protect your quota."
        break
    fi

    echo "Space check passed (${CURRENT_GB}GB / ${MAX_ALLOWED_GB}GB). Starting ${SRR}..."
    
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

echo "Compiling MultiQC reports..."
module purge
module load MultiQC/1.28-foss-2024a
multiqc qc_before/ -n multiqc_raw_custom.html
multiqc qc_after/ -n multiqc_trimmed_custom.html
echo "Workflow finished!"

# ==============================================================================
# CONFIGURATION DESIGN & SAFETY NOTE
# ==============================================================================
# This version pulls the max space constraint into a highly visible variable 
# at the very top of the script. This is how you build a tool for other 
# people to use—they only have to change one number at the top of the file 
# without digging into the code logic.
#
# (Note: Because it checks the space before downloading a run, the final 
# downloaded sample might push the folder size slightly over your limit 
# by 1-2GB. If you need a strict hard stop, you would set your limit a bit 
# lower than your actual absolute quota).
# ==============================================================================

# Run Script in home Directory:
#  sbatch 4._Automated_Script_User-DefinedMaximumSpace.sh
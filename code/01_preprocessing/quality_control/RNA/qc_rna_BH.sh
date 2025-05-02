#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J fastqc_BH
#SBATCH -t 01:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# Loading FastQC module
module load bioinfo-tools FastQC/0.11.9

# === CONFIGURATION VARIABLES ===
export SRCDIR=/home/edman/genomeAnalysis/data/raw_data/illumina/RNA/seq_BH
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/01_preprocessing/rna/BH
mkdir -p "$RESULT_DIR"

# Copying compressed FASTQ files to local temporary storage
cp "$SRCDIR"/* "$SNIC_TMP/"

# Changing to local temporary directory
cd "$SNIC_TMP" || exit 1

# Running FastQC on all FASTQ files
fastqc * > "$SNIC_TMP/fastqc_BH.log" 2>&1

# Moving FastQC output and log to result directory
mv *.html *.zip fastqc_BH.log "$RESULT_DIR/"


#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J fastqc_Serum
#SBATCH -t 01:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# Loading FastQC module
module load bioinfo-tools FastQC/0.11.9

# === CONFIGURATION VARIABLES ===
export SRCDIR=/home/edman/genomeAnalysis/data/raw_data/illumina/RNA/seq_serum
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/01_preprocessing/rna/post_trim/serum
mkdir -p "$RESULT_DIR"

# Copying compressed FASTQ files to local temporary storage
cp "$SRCDIR"/* "$SNIC_TMP/"

# Changing to local temporary directory
cd "$SNIC_TMP" || exit 1

# Running FastQC on all FASTQ files
fastqc * > "$SNIC_TMP/fastqc_serum.log" 2>&1

# Moving FastQC output and log to result directory
mv *.html *.zip fastqc_serum.log "$RESULT_DIR/"


#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J fastqc_illumina
#SBATCH -t 00:40:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# Loading FastQC module
module load bioinfo-tools FastQC/0.11.9

# === CONFIGURATION VARIABLES ===
export SRCDIR=/home/edman/genomeAnalysis/data/trimmed_data/dna
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/01_preprocessing/dna/illumina/post_trim

# Creating destination directory (if it doesn't exist)
mkdir -p "$RESULT_DIR"

# Copying compressed FASTQ files to local temporary storage
cp "$SRCDIR"/*_paired.fq.gz cp "$SRCDIR"/*_unpaired.fq.gz "$SNIC_TMP/"

# Changing to local temporary directory
cd "$SNIC_TMP" || exit 1

# Unzipping FASTQ files
gunzip *.fq.gz

# Running FastQC on all uncompressed FASTQ files,
# redirecting stdout and stderr to the log file
fastqc *.fq > "$SNIC_TMP/fastqc.log" 2>&1

# Moving FastQC output files (HTML, ZIP reports, and log) to the result directory
mv *.html *.zip "$SNIC_TMP/fastqc.log" "$RESULT_DIR/"

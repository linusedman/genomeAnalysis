#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J trimmomatic_run
#SBATCH -t 02:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# Loading required modules
module load bioinfo-tools trimmomatic/0.39

# CONFIGURATION VARIABLES
export SRCDIR=/home/edman/genomeAnalysis/data/raw_data/illumina
export RESULT_DIR=/home/edman/genomeAnalysis/data/trimmed_data/dna
export LOG_DIR=/home/edman/genomeAnalysis/analyses/01_preprocessing/dna/illumina/trim

# Creating the destination directories
mkdir -p "$RESULT_DIR"
mkdir -p "$LOG_DIR"

# Copying the FASTQ files to local temporary storage
cp "$SRCDIR"/*.fq.gz "$SNIC_TMP/"

# Changing to local temporary directory
cd "$SNIC_TMP" || exit 1

# Unzipping FASTQ files
gunzip *.fq.gz

# Running Trimmomatic in paired-end mode without adapter clipping
trimmomatic PE -threads 4 \
    E745-1.L500_SZAXPI015146-56_1_clean.fq E745-1.L500_SZAXPI015146-56_2_clean.fq \
    E745-1.L500_SZAXPI015146-56_1_paired.fq E745-1.L500_SZAXPI015146-56_1_unpaired.fq \
    E745-1.L500_SZAXPI015146-56_2_paired.fq E745-1.L500_SZAXPI015146-56_2_unpaired.fq \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 > "$SNIC_TMP/trimmomatic.log" 2>&1

# Compressing the output FASTQ files
gzip *.fq

# Moving the compressed files and log to the result directory
mv *.fq.gz "$RESULT_DIR/"
mv "$SNIC_TMP/trimmomatic.log" "$LOG_DIR/"


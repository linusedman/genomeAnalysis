#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J trimmomatic_rna_batch
#SBATCH -t 02:00:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# Loading required modules
module load bioinfo-tools trimmomatic/0.39

# CONFIGURATION VARIABLES
SRCDIR=/home/edman/genomeAnalysis/data/raw_data/illumina/RNA/seq_BH
RESULT_DIR=/home/edman/genomeAnalysis/data/trimmed_data/rna/BH
LOG_DIR=/home/edman/genomeAnalysis/analyses/01_preprocessing/rna/trim
ADAPTERS=/sw/bioinfo/trimmomatic/0.39/rackham/adapters/TruSeq3-PE.fa

mkdir -p "$RESULT_DIR" "$LOG_DIR"

# Copying the FASTQ files to local temporary storage
cp "$SRCDIR"/trim_paired_ERR*.fastq.gz "$SNIC_TMP/"

# Changing to local temporary directory
cd "$SNIC_TMP" || exit 1

# Unzipping FASTQ files
gunzip *.fastq.gz

# SAMPLE LIST
SAMPLES=("ERR1797972" "ERR1797973" "ERR1797974")

# LOOP OVER SAMPLES
# Running Trimmomatic in paired-end mode with adapter clipping
for SAMPLE in "${SAMPLES[@]}"; do
    IN1="trim_paired_${SAMPLE}_pass_1.fastq"
    IN2="trim_paired_${SAMPLE}_pass_2.fastq"
    OUT1P="${SAMPLE}_1_paired.fq"
    OUT1U="${SAMPLE}_1_unpaired.fq"
    OUT2P="${SAMPLE}_2_paired.fq"
    OUT2U="${SAMPLE}_2_unpaired.fq"

    trimmomatic PE -threads 8 \
        "$IN1" "$IN2" \
        "$OUT1P" "$OUT1U" "$OUT2P" "$OUT2U" \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
        >> "$SNIC_TMP/trimmomatic.log" 2>&1
done

# Compressing the output FASTQ files
gzip *.fq

# Moving the compressed files and log to the result directory
mv *.fq.gz "$RESULT_DIR/"
mv trimmomatic.log "$LOG_DIR/"


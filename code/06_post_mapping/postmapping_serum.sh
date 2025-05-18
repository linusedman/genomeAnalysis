#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 06:00:00
#SBATCH -J bwa_serum
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bwa/0.7.18 samtools/1.20

# CONFIGURATION VARIABLES
ASSEMBLY_GZ=/home/edman/genomeAnalysis/data/assembly_data/pacBio/assembly.contigs.fasta.gz
READS_DIR=/home/edman/genomeAnalysis/data/raw_data/illumina/RNA/seq_serum
OUT_DIR=/home/edman/genomeAnalysis/analyses/06_post_mapping/serum
mkdir -p "$OUT_DIR"

# Copy and unzip the assembly
cp "$ASSEMBLY_GZ" "$SNIC_TMP/"
cd "$SNIC_TMP"
gunzip -f "$(basename "$ASSEMBLY_GZ")"
ASSEMBLY=${ASSEMBLY_GZ##*/}  # extract filename
ASSEMBLY="${ASSEMBLY%.gz}"   # remove .gz
bwa index "$ASSEMBLY"

# Sample List
SAMPLES=(ERR1797969 ERR1797970 ERR1797971)

for SAMPLE in "${SAMPLES[@]}"; do
    bwa mem -t 4 "$ASSEMBLY" "$READS_DIR/trim_paired_${SAMPLE}_pass_1.fastq.gz" "$READS_DIR/trim_paired_${SAMPLE}_pass_2.fastq.gz" |
    samtools view -bS - |
    samtools sort -o "$OUT_DIR/${SAMPLE}_sorted.bam"

    samtools index "$OUT_DIR/${SAMPLE}_sorted.bam"
done


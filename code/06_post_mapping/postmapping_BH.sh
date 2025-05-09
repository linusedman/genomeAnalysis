#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 06:00:00
#SBATCH -J bwa_brain_heart
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools bwa/0.7.18 samtools/1.20

# === CONFIGURATION ===
ASSEMBLY_GZ=/home/edman/genomeAnalysis/data/assembly_data/pacBio/run2/assembly.contigs.fasta.gz
READS_DIR=/home/edman/genomeAnalysis/data/raw_data/illumina/RNA/seq_BH
OUT_DIR=/home/edman/genomeAnalysis/analyses/06_post_mapping/BH
mkdir -p "$OUT_DIR"

# Copy and unzip the assembly
cp "$ASSEMBLY_GZ" "$SNIC_TMP/"
cd "$SNIC_TMP"
gunzip -f "$(basename "$ASSEMBLY_GZ")"
ASSEMBLY=${ASSEMBLY_GZ##*/}  # get filename
ASSEMBLY="${ASSEMBLY%.gz}"   # remove .gz
bwa index "$ASSEMBLY"

# === Sample List ===
SAMPLES=(ERR1797972 ERR1797973 ERR1797974)

for SAMPLE in "${SAMPLES[@]}"; do
    bwa mem -t 4 "$ASSEMBLY" "$READS_DIR/trim_paired_${SAMPLE}_pass_1.fastq.gz" "$READS_DIR/trim_paired_${SAMPLE}_pass_2.fastq.gz" |
    samtools view -bS - |
    samtools sort -o "$OUT_DIR/${SAMPLE}_sorted.bam"

    samtools index "$OUT_DIR/${SAMPLE}_sorted.bam"
done


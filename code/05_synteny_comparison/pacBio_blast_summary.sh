#!/bin/bash
#SBATCH -A uppmax2025-2-288
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J blast_summary
#SBATCH -t 00:15:00
#SBATCH --mail-user=linus.edman.8474@student.uu.se
#SBATCH --mail-type=ALL

# === CONFIGURATION VARIABLES ===
export BLAST_OUTPUT=/home/edman/genomeAnalysis/analyses/05_synteny_comparison/blast_output/comparison.tab
export RESULT_DIR=/home/edman/genomeAnalysis/analyses/05_synteny_comparison/blast_output
mkdir -p "$RESULT_DIR"

# Copy BLAST file to local temp directory
cp "$BLAST_OUTPUT" "$SNIC_TMP/"
cd "$SNIC_TMP" || exit 1


# Run AWK summary
awk '{
    contig=$1;
    ref=$2;
    identity=$3;
    length=$4;
    qstart=$7; qend=$8;
    strand = (qstart <= qend ? "+" : "-");

    hits[contig]++;
    total_len[contig] += length;
    total_id[contig] += identity;
    refs[contig][ref]++;
    strands[contig][strand]++;
}
END {
    printf "%-25s %-6s %-10s %-10s %-20s %-10s\n", "Contig", "Hits", "Total_bp", "Avg_ID(%)", "Refs", "Strands";
    for (c in hits) {
        ref_list=""; strand_list="";
        for (r in refs[c]) { ref_list = ref_list r "," }
        for (s in strands[c]) { strand_list = strand_list s "," }

        printf "%-25s %-6d %-10d %-10.2f %-20s %-10s\n", c, hits[c], total_len[c], total_id[c]/hits[c], ref_list, strand_list
    }
}' comparison.tab > contig_summary.tsv

# Move summary back to result directory
mv contig_summary.tsv "$RESULT_DIR/"


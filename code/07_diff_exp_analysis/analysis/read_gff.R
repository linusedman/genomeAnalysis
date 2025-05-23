# Read GFF lines
gff_lines <- readLines("pacbio.gff")

# Truncate after ##FASTA
fasta_line <- grep("##FASTA", gff_lines)
if (length(fasta_line) > 0) {
  gff_lines <- gff_lines[1:(fasta_line - 1)]
}

# Write only annotation part
writeLines(gff_lines, "annotation_only.gff")

# Load rtracklayer and import GFF
library(rtracklayer)
gff <- import.gff("annotation_only.gff")

# Convert to data frame
annotations <- as.data.frame(gff)

# Extract UniProt ID from inference field
annotations$uniprot_id <- ifelse(
  grepl("UniProtKB:", annotations$inference),
  sub(".*UniProtKB:([A-Z0-9]+).*", "\\1", annotations$inference),
  NA
)

# Create mapping data frame
gene_mapping <- data.frame(
  gene_id = annotations$locus_tag,
  gene_name = ifelse(is.na(annotations$gene), annotations$product, annotations$gene),
  product = annotations$product,
  uniprot_id = annotations$uniprot_id,
  stringsAsFactors = FALSE
)

# Preview and export
head(gene_mapping)
write.csv(gene_mapping, "gene_mapping_from_gff.csv", row.names = FALSE)

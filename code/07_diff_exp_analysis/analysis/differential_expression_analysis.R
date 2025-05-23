# Load libraries
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(pheatmap)

# Set working directory
setwd("/Users/linus/Library/CloudStorage/OneDrive-Uppsalauniversitet/Genome analysis/projekt/read_count_files")

# List and sort count files
files <- sort(list.files(pattern = "_counts.txt$"))
print(files)

# Read and merge count files
count_list <- lapply(files, function(f) {
  read.table(f, header = FALSE, row.names = 1, stringsAsFactors = FALSE)
})
counts <- do.call(cbind, count_list)
colnames(counts) <- sub("_counts.txt", "", files)
counts <- counts[!grepl("^__", rownames(counts)), ]

# Sample metadata
sample_conditions <- data.frame(
  row.names = colnames(counts),
  condition = factor(c("Serum", "Serum", "Serum", "BHI", "BHI", "BHI"), levels = c("BHI", "Serum"))
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = sample_conditions,
                              design = ~ condition)


# Summary statistics with consistent filtering threshold
gene_matrix <- counts(dds, normalized = FALSE)
total_genes <- nrow(gene_matrix)
expressed_genes_strict <- sum(rowSums(gene_matrix >=10) >= 3)

cat("Total number of genes (before filtering):", total_genes, "\n")
cat("Number of expressed genes (>=10 counts in ≥3 samples):", expressed_genes_strict, "\n")


# Filter low count genes (>=10)
keep <- rowSums(counts(dds) >=10) >= 3
dds <- dds[keep, ]

# Calculate total counts per gene across all samples
gene_counts <- rowSums(counts(dds))

# Create a histogram of total read counts per gene
hist_data <- data.frame(TotalCounts = gene_counts)

# Plot and save histogram (with 'final' in filename)
hist_plot <- ggplot(hist_data, aes(x = TotalCounts)) +
  geom_histogram(bins = 100, color = "black", fill = "steelblue") +
  scale_x_log10() +  # Log10 scale for better readability
  labs(title = "Distribution of Total Counts per Gene",
       x = "Total Read Counts per Gene (log10 scale)",
       y = "Number of Genes") +
  theme_minimal()

ggsave("histogram_gene_counts_final.png", plot = hist_plot, width = 6, height = 5)
print(hist_plot)

# Run DESeq
dds <- DESeq(dds)

# Shrink log2 fold changes
res <- lfcShrink(dds,
                 coef = "condition_serum_vs_BH",
                 type = "apeglm")

# Results table
res_df <- as.data.frame(res)
res_df <- cbind(`Gene ID` = rownames(res_df), res_df)
res_df <- res_df[!is.na(res_df$log2FoldChange), ]
res_df$fold_change <- 2^res_df$log2FoldChange

# Merge with gene mapping
gene_mapping <- read.csv("gene_mapping_from_gff.csv", stringsAsFactors = FALSE)
res_df <- merge(res_df, gene_mapping, by.x = "Gene ID", by.y = "gene_id", all.x = TRUE)

# Save results
write.csv(res_df, "deseq2_all_results_final.csv", row.names = FALSE)
write.xlsx(res_df, "deseq2_all_results_final.xlsx", rowNames = FALSE)

# Filter significant results
sig_res_df <- res_df[res_df$padj < 0.001, ]
write.csv(sig_res_df, "deseq2_significant_results_final.csv", row.names = FALSE)
write.xlsx(sig_res_df, "deseq2_significant_results_final.xlsx", rowNames = FALSE)

# MA plot
plotMA(res, main = "DESeq2 MA Plot", ylim = c(-2, 2))

# PCA plot with labels
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggpca <- ggplot(pcaData, aes(PC1, PC2, color = condition, label = rownames(pcaData))) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA Plot")
ggsave("pca_plot_final.png", plot = ggpca, width = 6, height = 5)
print(ggpca)

# Mark significance
res_df$significant <- ifelse(
  !is.na(res_df$padj) & res_df$padj < 0.001 & abs(res_df$log2FoldChange) > 1,
  "q<0.001 & \n |log2FC|>1", "Not significant"
)

# Top genes
top_up <- res_df[order(-res_df$log2FoldChange), ][1:10, ]
top_down <- res_df[order(res_df$log2FoldChange), ][1:10, ]
label_genes <- rbind(top_up, top_down)

# Save top 20
write.csv(label_genes, "top20_labeled_genes_final.csv", row.names = FALSE)
write.xlsx(label_genes, "top20_labeled_genes_final.xlsx", rowNames = FALSE)

# Top 5 for labeling
top_up_5 <- top_up[1:5, ]
top_down_5 <- top_down[1:5, ]
label_genes_10 <- rbind(top_up_5, top_down_5)

# log10(padj)
res_df$log10_padj <- -log10(res_df$padj)
label_genes_10$log10_padj <- -log10(label_genes_10$padj)

# Label manual positions
label_left_manual <- data.frame(
  gene_name = label_genes_10$gene_name[label_genes_10$log2FoldChange < 0],
  x = -7,
  y = seq(300, 260, length.out = 5)
)

label_right <- subset(label_genes_10, log2FoldChange >= 0)
label_right_manual <- data.frame(
  gene_name = label_right$gene_name,
  x = c(9.5, 9.5, 9.5, 9.5, 9.5),
  y = c(290, 50, 150, 100, 200)
)
label_left_actual <- merge(label_left_manual, label_genes_10, by = "gene_name")
label_right_actual <- merge(label_right_manual, label_genes_10, by = "gene_name")

# Volcano plot
ggvolcano <- ggplot(res_df, aes(x = log2FoldChange, y = log10_padj, color = significant)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "gray") +
  geom_segment(data = label_left_actual,
               aes(x = log2FoldChange, y = log10_padj, xend = x, yend = y),
               inherit.aes = FALSE, color = "gray40", linetype = "dotted", linewidth = 0.3) +
  geom_text(data = label_left_actual,
            aes(x = x, y = y, label = gene_name),
            inherit.aes = FALSE, hjust = 1, size = 3) +
  geom_segment(data = label_right_actual,
               aes(x = log2FoldChange, y = log10_padj, xend = x, yend = y),
               inherit.aes = FALSE, color = "gray40", linetype = "dotted", linewidth = 0.3) +
  geom_text(data = label_right_actual,
            aes(x = x, y = y, label = gene_name),
            inherit.aes = FALSE, hjust = 0, size = 3) +
  scale_color_manual(
    name = "Significance level",
    values = c("q<0.001 & \n |log2FC|>1" = "orange", "Not significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot  Serum vs. BHI",
       x = "log2 Fold Change",
       y = "-log10 q-value") +
  xlim(-11, 11)

ggsave("volcano_plot_final.png", plot = ggvolcano, width = 8, height = 6)
print(ggvolcano)

# Heatmap
heat_data <- assay(vsd)[label_genes$`Gene ID`, ]
heat_data_z <- t(scale(t(heat_data)))
rownames(heat_data_z) <- label_genes$gene_name

png("heatmap_top20_log2FC_genes_final.png", width = 2400, height = 1800, res = 300)
pheatmap(heat_data_z,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = sample_conditions,
         fontsize_row = 8,
         main = "Top 10 Up/Down Regulated Genes (Z-score scaled)")

dev.off()


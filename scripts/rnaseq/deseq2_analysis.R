# DESeq2 differential expression on the synthetic dataset (Week 4 frozen output).
# Run from repo root:  Rscript scripts/rnaseq/deseq2_analysis.R
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

counts <- read.csv("data/rnaseq/counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("data/rnaseq/coldata.csv", stringsAsFactors = TRUE)
stopifnot(identical(colnames(counts), as.character(coldata$sample_id)))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
res$gene <- rownames(res)

dir.create("scripts/rnaseq/output", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "scripts/rnaseq/output/deseq2_results.csv", row.names = FALSE)

sig <- subset(res, padj < 0.05)
cat("DESeq2: ", nrow(sig), " genes with padj < 0.05;",
    " ", sum(sig$log2FoldChange > 0, na.rm = TRUE), " up /",
    " ", sum(sig$log2FoldChange < 0, na.rm = TRUE), " down\n", sep = "")

planted <- sprintf("GENE%03d", 1:40)
recovered <- sum(planted %in% sig$gene)
cat("planted DE genes recovered:", recovered, "/ 40\n")

# ---- plots ----
png("scripts/rnaseq/output/deseq2_ma.png", width = 900, height = 600, res = 130)
plotMA(results(dds), main = "DESeq2 MA-plot (synthetic RNA-seq)")
dev.off()

rld <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_dat <- as.data.frame(plotPCA(rld, intgroup = "condition", returnData = TRUE))
pct_var <- round(100 * attr(pca_dat, "percentVar"), 1)

p <- ggplot(pca_dat, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA (vst), DESeq2") +
  theme_minimal()
ggsave("scripts/rnaseq/output/deseq2_pca.png", p, width = 6, height = 4.5, dpi = 130)

vol <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.4) +
  scale_color_manual(values = c("grey60", "firebrick"), name = "padj < 0.05") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  ggtitle("Volcano — DESeq2") +
  theme_minimal()
ggsave("scripts/rnaseq/output/deseq2_volcano.png", vol, width = 6, height = 4.5, dpi = 130)

cat("top 10 by padj:\n")
print(head(res[order(res$padj), c("gene", "baseMean", "log2FoldChange", "padj")], 10),
      row.names = FALSE)

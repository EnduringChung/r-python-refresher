# edgeR differential expression on the synthetic dataset (Week 4 frozen output).
# Same data, same questions, second dialect — results should echo DESeq2's.
# Run from repo root:  Rscript scripts/rnaseq/edger_analysis.R
suppressPackageStartupMessages({
  library(edgeR)
})

counts <- read.csv("data/rnaseq/counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("data/rnaseq/coldata.csv", stringsAsFactors = TRUE)

y <- DGEList(counts = counts, group = coldata$condition)
y <- calcNormFactors(y, method = "TMM")        # the edgeR analogue of size factors
design <- model.matrix(~ condition, data = coldata)
y <- estimateDisp(y, design)                    # NB dispersion
fit <- glmQLFit(y, design)                      # quasi-likelihood
qlf <- glmQLFTest(fit, coef = "conditiontreated")

res <- as.data.frame(topTags(qlf, n = Inf)$table)
res$gene <- rownames(res)
res$padj <- res$FDR

write.csv(res, "scripts/rnaseq/output/edger_results.csv", row.names = FALSE)

sig <- subset(res, FDR < 0.05)
cat("edgeR:  ", nrow(sig), " genes with FDR < 0.05\n", sep = "")

planted <- sprintf("GENE%03d", 1:40)
cat("planted DE genes recovered:", sum(planted %in% sig$gene), "/ 40\n")

cat("top 10 by FDR:\n")
print(head(res[order(res$FDR), c("gene", "logCPM", "logFC", "FDR")], 10),
      row.names = FALSE)

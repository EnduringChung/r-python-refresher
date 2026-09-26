# Trajectory — the slingshot mirror of Week 5 Day 2.
# Target: local RStudio/Positron.
#   BiocManager::install("slingshot")   # pulls SingleCellExperiment
# STATUS: reference implementation authored for this course; not executed
# in CI. Runs on slingshot >= 2.

library(SingleCellExperiment)
library(slingshot)

set.seed(9)
n <- 300

# ------------------------------------------------- seeded branched truth
ps <- rbeta(n, 1.2, 1.2)
branch <- ifelse(ps < 0.5, "trunk", sample(c("A", "B"), n, replace = TRUE))
prog <- 2 * ps + rnorm(n * 60, 0, .4)
brA <- 3 * (branch == "A" & ps > 0.5) + rnorm(n * 40, 0, .4)
brB <- 3 * (branch == "B" & ps > 0.5) + rnorm(n * 40, 0, .4)
hk <- rnorm(n * 200, 0, .6)
E <- cbind(matrix(prog, n, 60), matrix(brA, n, 40),
           matrix(brB, n, 40), matrix(hk, n, 200))
counts <- t(rpois(n * 440, 10 * exp(E / 2)))   # genes x cells

# ------------------------------------ normalize + log + PCA (manual, W5D2)
cell_tot <- colSums(counts)
Xn <- log1p(counts / cell_tot * median(cell_tot))
pca <- prcomp(t(Xn), rank. = 15)$x             # cells x PCs

sce <- SingleCellExperiment(
  assays  = list(logcounts = Xn),
  reducedDims = list(PCA = pca),
  colData = data.frame(true_ps = ps, branch = branch)
)

# slingshot needs cluster labels: coarse k-means along the lineage
set.seed(4)
sce$clusters <- factor(kmeans(pca[, 1:6], centers = 5)$cluster)

# root = lowest progressive-gene score (first 60 genes) — the biological claim
prog_score <- rowMeans(Xn[1:60, ])
start_cluster <- names(which.min(tapply(prog_score, sce$clusters, mean)))

slc <- slingshot(sce, clusterLabels = "clusters", reducedDim = "PCA",
                 start.clus = start_cluster)

pst <- slc$slingPseudotime_1
ok <- !is.na(pst)
cat("Spearman(truth, slingshot pseudotime):",
    round(cor(ps[ok], pst[ok], method = "spearman"), 3), "\n")

plot(pst[ok], ps[ok], xlab = "slingshot pseudotime", ylab = "true ps",
     pch = 16, col = factor(branch)[ok])
legend("topleft", legend = c("A", "B", "trunk"),
       col = 1:3, pch = 16)

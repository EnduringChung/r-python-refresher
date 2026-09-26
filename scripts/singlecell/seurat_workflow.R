# Single-cell core workflow — the Seurat mirror of Week 5 Day 1.
# Target: local RStudio/Positron (install.packages("Seurat")).
# STATUS: reference implementation authored for this course; not executed
# in CI. Runs on Seurat >= 5.
#
# scanpy idiom -> Seurat equivalent, same order as the live lesson:
#   normalize_total -> LogNormalize, HVG -> VariableFeatures,
#   PCA -> RunPCA, neighbors -> FindNeighbors, leiden -> FindClusters,
#   umap -> RunUMAP, rank_genes_groups -> FindAllMarkers.

library(Seurat)
library(dplyr)

set.seed(5)

# ---------------------------------------------------------------- data
# swap for real data: Read10X() / ReadH5AD() / LoadH5Seurat()
n_per <- 60; n_markers <- 6; n_hk <- 194
types <- rep(c("T-like", "B-like", "Mono-like"), each = n_per)
depth <- exp(rnorm(3 * n_per, 0, 0.4))
effects <- matrix(0, 3 * n_per, n_markers)
for (t in 1:3) effects[(t - 1) * n_per + 1:n_per, (t - 1) * 2 + 1:2] <- 2
lam <- depth * 5 * exp(effects)
counts <- cbind(matrix(rpois(3 * n_per * n_markers, as.vector(lam)), 3 * n_per),
                matrix(rpois(3 * n_per * n_hk, rep(depth * 5, n_hk)), 3 * n_per))
rownames(counts) <- paste0("g", seq_len(ncol(counts)))
colnames(counts) <- paste0("cell", seq_len(nrow(counts)))

obj <- CreateSeuratObject(counts = counts, meta.data = data.frame(true_type = types))
obj <- UpdateSeuratObject(obj)

# ------------------------------------------------------------- pipeline
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")   # QC metrics
obj <- subset(obj, subset = nFeature_RNA > 20)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- FindVariableFeatures(obj, nfeatures = 50)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 20)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:10)

# ------------------------------------------------------- evaluate + plot
ari_available <- requireNamespace("aricode", quietly = TRUE)
if (ari_available) {
  print(aricode::ARI(obj$true_type, obj$seurat_clusters))
}
DimPlot(obj, group.by = c("seurat_clusters", "true_type"))
obj <- FindAllMarkers(obj, only.pos = TRUE)
head(obj, 5)

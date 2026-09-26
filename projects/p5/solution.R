# Capstone — R cross-check of the DE section (reference solution)
# Run from the repo root:  Rscript projects/p5/solution.R
#
# Purpose: the parity thread — R must recover the same DE gene set as
# projects/p5/solution.py (Welch t-test + Benjamini-Hochberg, same data).

library(tidyverse)

meta <- read_csv("data/capstone_meta.csv", show_col_types = FALSE)
rna  <- read_csv("data/capstone_rna_counts.csv", show_col_types = FALSE)

## QC (same fence as the Python solution) ----
X <- as.matrix(rna[, -1])
totals <- rowSums(X)
q <- quantile(totals, c(0.25, 0.75))
iqr <- q[2] - q[1]
fence <- totals >= q[1] - 1.5 * iqr & totals <= q[2] + 1.5 * iqr
cat("flagged:", sum(!fence), "samples\n")

X <- X[fence, ]
meta_f <- meta[fence, ]
treated <- meta_f$condition == "treated"

## normalize + log (identical formula) ----
Xn <- log1p(X / rowSums(X) * median(rowSums(X)))

## vectorized Welch t-test per gene (base R, no loop) ----
nx <- sum(treated); ny <- sum(!treated)
mx <- colMeans(Xn[treated, ]);  my <- colMeans(Xn[!treated, ])
vx <- colSums((Xn[treated, ] - rep(mx, each = nx))^2) / (nx - 1)
vy <- colSums((Xn[!treated, ] - rep(my, each = ny))^2) / (ny - 1)

se2 <- vx / nx + vy / ny
tstat <- (mx - my) / sqrt(se2)
df <- se2^2 / ((vx/nx)^2 / (nx - 1) + (vy/ny)^2 / (ny - 1))
pvals <- 2 * pt(-abs(tstat), df)

## Benjamini-Hochberg ----
o <- order(pvals)
adj <- pvals[o] * length(pvals) / seq_along(pvals)
fdr <- numeric(length(pvals))
fdr[o] <- pmin(rev(cummin(rev(adj))), 1)   # cummin runs LARGEST -> smallest
sig <- which(fdr < 0.05)
cat("DE genes at FDR < 0.05:", length(sig), "\n")
cat("up-regulated:", sum((mx - my)[sig] > 0), "\n")

## log2FC for the top genes (sanity vs the Python solution) ----
top10 <- order(pvals)[1:10]
cat("top-10 genes (index, treated/control mean ratio):\n")
print(round((mx / my)[top10], 2))

## write the gene set for cross-language Jaccard (used by report) ----
writeLines(as.character(sig - 1L), "projects/p5/de_genes_r.txt")
cat("wrote projects/p5/de_genes_r.txt (0-based indices, R/Python comparable)\n")

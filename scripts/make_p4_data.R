# Project 4 dataset: synthetic bulk expression matrix with a planted tumor signature.
# Seed-fixed; run from repo root:  Rscript scripts/make_p4_data.R
#
# Calibrated for pedagogy: logistic CV ~0.85, forest ~0.78 — realistic
# difficulty (imperfect, fold-to-fold variance), not a 1.0 toy.
set.seed(2025)

n_genes <- 200
n_per_group <- 30
n_sig <- 10
depth_sd <- 0.25
fc_mu <- 1.6; fc_sd <- 0.35; size <- 10

genes <- sprintf("SIG%03d", seq_len(n_genes))
sig_genes <- genes[seq_len(n_sig)]                     # SIG001..SIG010

samples <- c(paste0("N", 1:n_per_group), paste0("T", 1:n_per_group))
condition <- factor(rep(c("normal", "tumor"), each = n_per_group),
                    levels = c("normal", "tumor"))

base <- rlnorm(n_genes, meanlog = log(120), sdlog = 0.8)
names(base) <- genes
depths <- rlnorm(length(samples), meanlog = 0, sdlog = depth_sd)

counts <- matrix(NA_integer_, n_genes, length(samples),
                 dimnames = list(genes, samples))
for (g in genes) {
  fc <- if (g %in% sig_genes) rlnorm(1, meanlog = log(fc_mu), sdlog = fc_sd) else 1
  for (j in seq_along(samples)) {
    mu <- base[g] * (if (condition[j] == "tumor") fc else 1) * depths[j]
    counts[g, j] <- rnbinom(1, size = size, mu = mu)
  }
}

dir.create("data/p4", showWarnings = FALSE, recursive = TRUE)
write.csv(counts, "data/p4/expr.csv")
write.csv(data.frame(sample_id = samples, condition = condition),
          "data/p4/labels.csv", row.names = FALSE)

score <- colMeans(counts[sig_genes, ])
cat("mean sig score normal:", round(mean(score[condition == "normal"]), 1),
    "| tumor:", round(mean(score[condition == "tumor"]), 1), "\n")
cat("wrote data/p4/expr.csv (", n_genes, "genes x", length(samples), "samples )\n",
    sep = "")

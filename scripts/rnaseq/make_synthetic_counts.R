# Generates a synthetic bulk RNA-seq dataset with planted differential expression.
# Seed-fixed: everyone who runs this gets byte-identical data.
set.seed(2024)

n_genes <- 500
n_per_group <- 3

genes <- sprintf("GENE%03d", seq_len(n_genes))
de_genes <- sprintf("GENE%03d", 1:40)   # 40 genes will be truly DE

coldata <- data.frame(
  sample_id = c(paste0("ctrl", 1:n_per_group), paste0("treat", 1:n_per_group)),
  condition = factor(rep(c("control", "treated"), each = n_per_group),
                     levels = c("control", "treated"))
)

# base expression: log-normal, most genes modest, a few highly expressed
base_expr <- rlnorm(n_genes, meanlog = log(150), sdlog = 1)
names(base_expr) <- genes

counts <- matrix(NA_integer_, nrow = n_genes, ncol = nrow(coldata),
                 dimnames = list(genes, coldata$sample_id))

for (g in genes) {
  is_de <- g %in% de_genes
  log_fc <- if (is_de) rnorm(1, mean = 1.8, sd = 0.5) else 0   # up-regulated
  size <- 40                                                   # NB dispersion
  mu_ctrl <- base_expr[g]
  mu_treat <- base_expr[g] * 2^log_fc
  for (i in 1:n_per_group) {
    counts[g, i] <- rnbinom(1, size = size, mu = mu_ctrl)
    counts[g, n_per_group + i] <- rnbinom(1, size = size, mu = mu_treat)
  }
}

dir.create("data/rnaseq", showWarnings = FALSE, recursive = TRUE)
write.csv(counts, "data/rnaseq/counts.csv")
write.csv(coldata, "data/rnaseq/coldata.csv", row.names = FALSE)

cat("wrote", n_genes, "genes x", ncol(counts), "samples;",
    length(de_genes), "planted DE genes\n")

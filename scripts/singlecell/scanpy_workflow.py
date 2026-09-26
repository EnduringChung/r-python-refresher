# Single-cell core workflow — the production version of Week 5 Day 1.
# Target: Google Colab (%pip install scanpy leidenalg) or a local env.
# STATUS: reference implementation authored for this course; not executed
# in CI (no browser/scanpy). Runs top-to-bottom on scanpy >= 1.10.
#
# The live lesson (w5/day1-singlecell.qmd) implements every step below in
# raw numpy — read the two side by side to map idioms to mechanics.

import numpy as np
import scanpy as sc

sc.settings.verbosity = 1
rng = np.random.default_rng(5)

# ---------------------------------------------------------------- data
# In Colab, swap this block for real data:
#   adata = sc.read_h5ad("your_data.h5ad")
#   # or sc.read_10x_mtx("filtered_feature_bc_matrix/", var_names="gene_symbols")
n_per, n_markers, n_hk = 60, 6, 194
types = np.repeat(["T-like", "B-like", "Mono-like"], n_per)
depth = rng.lognormal(0.0, 0.4, 3 * n_per)
effects = np.zeros((3 * n_per, n_markers))
for t in range(3):
    effects[t * n_per:(t + 1) * n_per, t * 2:(t + 1) * 2] = 2.0
lam = depth[:, None] * 5.0 * np.exp(effects)
counts = np.hstack([
    rng.poisson(lam),
    rng.poisson(depth[:, None] * 5.0, size=(3 * n_per, n_hk)),
])

adata = sc.AnnData(counts.astype(np.float32))
adata.obs["true_type"] = types

# ------------------------------------------------------------- pipeline
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False,
                           inplace=True)
sc.pp.filter_cells(adata, min_genes=50)
sc.pp.filter_genes(adata, min_cells=3)

adata.layers["counts"] = adata.X.copy()          # keep raw for DE
sc.pp.normalize_total(adata, target_sum=1e4)     # CPM-to-10k (W5D1 step 2)
sc.pp.log1p(adata)                               # step 2b
sc.pp.highly_variable_genes(adata, n_top_genes=50)   # step 3
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata, n_comps=20)                     # step 4
sc.pp.neighbors(adata, n_neighbors=10)           # step 5
sc.tl.leiden(adata, resolution=0.5, key_added="leiden")

# ------------------------------------------------------- evaluate + plot
from sklearn.metrics import adjusted_rand_score

ari = adjusted_rand_score(adata.obs["true_type"], adata.obs["leiden"])
print("ARI vs planted truth:", round(ari, 3))

sc.tl.umap(adata)
sc.pl.umap(adata, color=["leiden", "true_type"], save="_day1.png")
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")   # naming step
sc.pl.rank_genes_groups(adata, n_genes=5, save="_day1.png")

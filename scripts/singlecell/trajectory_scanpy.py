# Trajectory & pseudotime — the production version of Week 5 Day 2.
# Target: Google Colab (%pip install scanpy leidenalg) or a local env.
# STATUS: reference implementation authored for this course; not executed
# in CI. Runs top-to-bottom on scanpy >= 1.10.

import numpy as np
import scanpy as sc
from scipy.stats import spearmanr

rng = np.random.default_rng(9)
n = 300

# ------------------------------------------------- seeded branched truth
ps = rng.beta(1.2, 1.2, n)
branch = np.where(ps < 0.5, "trunk", rng.choice(["A", "B"], n))
prog = 2.0 * ps[:, None] + rng.normal(0, 0.4, (n, 60))
onA = ((branch == "A") & (ps > 0.5))[:, None]
onB = ((branch == "B") & (ps > 0.5))[:, None]
brA = 3.0 * onA + rng.normal(0, 0.4, (n, 40))
brB = 3.0 * onB + rng.normal(0, 0.4, (n, 40))
hk = rng.normal(0, 0.6, (n, 200))
counts = rng.poisson(10 * np.exp(np.hstack([prog, brA, brB, hk]) / 2))

adata = sc.AnnData(counts.astype(np.float32))
adata.obs["true_ps"] = ps
adata.obs["branch"] = branch

# ------------------------------------------------------- scanpy pipeline
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=60)
sc.tl.pca(adata, n_comps=15)
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.diffmap(adata)

# root choice is a BIOLOGICAL claim — here the planted early marker score
prog_genes = [f"{i}" for i in range(60)]
adata.obs["prog_score"] = adata[:, prog_genes].X.mean(axis=1)
adata.uns["iroot"] = int(np.argmin(adata.obs["prog_score"]))

sc.tl.dpt(adata)
rho = spearmanr(adata.obs["true_ps"], adata.obs["dpt_pseudotime"]).statistic
print("Spearman(truth, scanpy DPT):", round(rho, 3))   # ~0.7 like the live cell

# PAGA: which clusters are connected beyond chance?
sc.tl.leiden(adata, resolution=0.4, key_added="clusters")
sc.tl.paga(adata, groups="clusters")
sc.pl.paga(adata, save="_day2.png")

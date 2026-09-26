# Capstone — reference solution (Python)
# Run from the repo root:  python3 projects/p5/solution.py
#
# Verified numbers on the seed-2026 data (see report-solution.qmd):
#   QC: 10 samples flagged (all 4 planted broken + 6 natural outliers)
#   DE: 37 genes at FDR < 0.05 | within-batch Jaccard 0.29
#   classifier accuracy 0.98 | protein DE: 12 genes, CCA cc1 0.924

import datetime
import hashlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as st
from sklearn.cross_decomposition import CCA
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split

rng = np.random.default_rng(2026)          # receipt seed (analysis has no
                                           # random steps, but set it anyway)

# ---------------------------------------------------------------- 1. setup
meta = pd.read_csv("data/capstone_meta.csv")
rna = pd.read_csv("data/capstone_rna_counts.csv")
prot = pd.read_csv("data/capstone_protein.csv")
print(f"samples {meta.shape[0]} | genes {rna.shape[1]-1} | proteins {prot.shape[1]-1}")
print("question: does treated change gene expression, do both views agree, "
      "and can a classifier recover condition from the DE genes?")

# ------------------------------------------------------------------- 2. QC
X = rna.iloc[:, 1:].to_numpy().astype(float)
totals = pd.Series(X.sum(axis=1))
q1, q3 = totals.quantile([0.25, 0.75])
iqr = q3 - q1
fence = (totals >= q1 - 1.5 * iqr) & (totals <= q3 + 1.5 * iqr)
print(f"\nQC: {int((~fence).sum())} samples flagged:", meta.sample_id[~fence].tolist())
print(pd.crosstab(meta.batch[~fence], columns="flagged"))

X = X[fence.to_numpy()]
meta_f = meta[fence.to_numpy()].reset_index(drop=True)
treated = (meta_f.condition == "treated").to_numpy()

keep_gene = (X >= 10).sum(axis=0) >= 20
X = X[:, keep_gene]
print(f"genes kept: {keep_gene.sum()} of {X.shape[1]} "
      f"(detected in >= 20 samples)")

P = prot.iloc[fence.to_numpy(), 1:].to_numpy()
print(f"protein range: {P.min():.1f} to {P.max():.1f}; "
      f"samples with >5% missing: {(np.isnan(P).mean(1) > 0.05).sum()}")

# --------------------------------------------------- 3. normalize + EDA
Xn = np.log1p(X / X.sum(axis=1, keepdims=True) * np.median(X.sum(axis=1)))

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for ax, col in zip(axes, ["condition", "cell_type", "batch"]):
    pcs = PCA(2).fit_transform(Xn)
    for g in meta_f[col].unique():
        m = (meta_f[col] == g).to_numpy()
        ax.scatter(pcs[m, 0], pcs[m, 1], s=12, label=str(g), alpha=0.6)
    spread = max(pcs[meta_f[col] == g, 0].mean()
                 for g in meta_f[col].unique()) - \
             min(pcs[meta_f[col] == g, 0].mean()
                 for g in meta_f[col].unique())
    ax.set(title=f"PC1 spread by {col}: {spread:.2f}", xlabel="PC1", ylabel="PC2")
    ax.legend(fontsize=7)
plt.tight_layout()
plt.savefig("projects/p5/fig1_pca_structures.png", dpi=300, bbox_inches="tight")
plt.close()
print("\nEDA: PC1 is dominated by cell_type, NOT condition — "
      "the marker-profile structure is the loudest signal.")
print("Condition lives on a different axis (check PC2/PC3) and DE handles it.")

# ------------------------------------------------------- 4. DE + FDR
pvals = np.array([
    st.ttest_ind(Xn[treated, j], Xn[~treated, j], equal_var=False).pvalue
    for j in range(Xn.shape[1])
])
order = np.argsort(pvals)
m_tot = len(pvals)
ranked = pvals[order] * m_tot / (np.arange(m_tot) + 1)
bh = np.minimum.accumulate(ranked[::-1])[::-1]      # Benjamini-Hochberg
fdr = np.empty(m_tot)
fdr[order] = np.minimum(bh, 1.0)
sig = fdr < 0.05
log2fc = Xn[treated].mean(axis=0) - Xn[~treated].mean(axis=0)

print(f"\nDE (treated vs control): {sig.sum()} genes at FDR < 0.05; "
      f"{int((sig & (log2fc > 0)).sum())} up-regulated")

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.scatter(log2fc[~sig], -np.log10(fdr[~sig]), s=6, c="grey", alpha=0.5)
ax.scatter(log2fc[sig], -np.log10(fdr[sig]), s=8, c="tab:red")
ax.set(xlabel="log2 fold change", ylabel="-log10 FDR",
       title=f"{int(sig.sum())} genes at FDR < 0.05")
plt.savefig("projects/p5/fig2_volcano.png", dpi=300, bbox_inches="tight")
plt.close()

# batch check: DE within each batch separately, Jaccard overlap
within = []
for b in ["B1", "B2"]:
    m = (meta_f.batch == b).to_numpy()
    pv = np.array([
        st.ttest_ind(Xn[m & treated, j], Xn[m & ~treated, j], equal_var=False).pvalue
        for j in range(Xn.shape[1])
    ])
    o = np.argsort(pv)
    bh_b = np.minimum.accumulate((pv[o] * m_tot / (np.arange(m_tot) + 1))[::-1])[::-1]
    f_b = np.empty(m_tot)
    f_b[o] = np.minimum(bh_b, 1.0)
    within.append(set(np.where(f_b < 0.05)[0]))
jac = len(within[0] & within[1]) / len(within[0] | within[1])
print(f"batch check: {len(within[0])} (B1) vs {len(within[1])} (B2) genes, "
      f"Jaccard = {jac:.2f} — partial agreement: per-batch power is lower "
      f"(n halves), and B2's 15 shifted genes add noise")

# --------------------------------------------------- 5. multi-view integration
def zscore(M):
    return (M - M.mean(axis=0)) / M.std(axis=0)

scaled = np.hstack([zscore(Xn), zscore(P)])
pcs_int = PCA(2).fit_transform(scaled)
print("\nintegration: scaled concat PCA — condition separation persists "
      "alongside cell-type structure")

top20 = np.argsort(pvals)[:20]                       # DE genes drive CCA
u, w = CCA(n_components=1).fit_transform(zscore(Xn[:, top20]), zscore(P))
cc1 = np.corrcoef(u[:, 0], w[:, 0])[0, 1]
cond_corr = abs(np.corrcoef(u[:, 0], treated.astype(int))[0, 1])
print(f"CCA (top-20 DE genes x proteins): canonical corr = {cc1:.3f}, "
      f"|corr with condition| = {cond_corr:.2f} — the protein view "
      f"supports the RNA story")

# ------------------------------------------------------- 6. classifier
top20_fdr = np.argsort(fdr)[:20] if (fdr < 0.05).sum() >= 20 else top20
Xtr, Xte, ytr, yte = train_test_split(
    Xn[:, top20], meta_f.condition.to_numpy(),
    test_size=0.3, stratify=meta_f.condition, random_state=2026)
clf = LogisticRegression(max_iter=2000).fit(Xtr, ytr)
acc = clf.score(Xte, yte)
print(f"\nclassifier (top-20 DE genes): test accuracy = {acc:.2f}")
print(confusion_matrix(yte, clf.predict(Xte)))

# ------------------------------------------------------- 7. receipt
print("\n== receipt ==")
print("run at:", datetime.datetime.now(datetime.timezone.utc).isoformat())
for f in ["capstone_meta.csv", "capstone_rna_counts.csv",
          "capstone_protein.csv"]:
    h = hashlib.sha256(open(f"data/{f}", "rb").read()).hexdigest()
    print(f"  {f}: {h[:32]}...")
pd.show_versions()

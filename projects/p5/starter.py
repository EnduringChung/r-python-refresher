# Capstone — analysis starter (Python)
# Run from the repo root:  python3 projects/p5/starter.py
# Data: data/capstone_{meta,rna_counts,protein}.csv (seed 2026, regenerable)

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

rng = np.random.default_rng(2026)

# 1. Setup ----------------------------------------------------------------
meta = pd.read_csv("data/capstone_meta.csv")
rna = pd.read_csv("data/capstone_rna_counts.csv")
prot = pd.read_csv("data/capstone_protein.csv")

# print dimensions + your analysis question in one sentence

# 2. QC -------------------------------------------------------------------
X = rna.iloc[:, 1:].to_numpy().astype(float)
totals = pd.Series(X.sum(axis=1))

# a) IQR fence on totals: flag samples, report WHICH BATCH they're from
# b) keep genes detected (>= 10 counts) in >= 20 samples
# c) protein view: value range + samples with > 5% missing

# 3. normalize + EDA ------------------------------------------------------
# a) CPM-to-median + log1p  (log AFTER depth-scaling!)
# b) histogram of log-expression — Day 3 pre-flight rules
# c) PCA: one panel per coloring (condition / cell_type / batch) — say
#    honestly which structure dominates PC1

# 4. DE + FDR -------------------------------------------------------------
# a) Welch t-test per gene (equal_var=False), treated vs control
# b) Benjamini-Hochberg over all genes (no statsmodels shortcut needed:
#    sort, scale by m/i, cummin from the largest — W2D5 pattern)
# c) how many at FDR < 0.05; volcano plot (log2FC vs -log10 FDR)
# d) batch check: repeat DE within each batch; Jaccard overlap of the two
#    gene lists; interpret

# 5. multi-view integration ----------------------------------------------
# a) z-score BOTH views, concat, PCA
# b) CCA between the top-20 DE genes (RNA) and the protein matrix:
#    canonical correlation + correlation with condition

# 6. classifier -----------------------------------------------------------
# top-20 DE genes as features; stratified 70/30 split, random_state=2026;
# LogisticRegression(max_iter=2000); test accuracy + confusion matrix +
# one sentence on what it does and does not prove

# 7. receipt --------------------------------------------------------------
# date (UTC), SHA-256 of the three data files, pd.show_versions()

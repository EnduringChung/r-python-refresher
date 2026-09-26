# Project 2 — Peptide assay EDA + statistics (Python starter)
# Run from the repo root:  python3 projects/p2/starter.py
# Data: data/peptides_assay_clean.csv (your Project 1 output, committed)

import datetime
import hashlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as st
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

rng = np.random.default_rng(42)   # the rng object IS your seed receipt

dat = pd.read_csv("data/peptides_assay_clean.csv")

# 1. Structure: shape, dtypes, dat.describe()

# 2. Missingness: dat.isna().sum() — why is concentration the only column
#    with NaNs? (think: what does a missing well mean?)

# 3. EDA histogram: concentration, raw scale. Choose bins ON PURPOSE.
#    Save: plt.savefig("projects/p2/fig1_hist_raw.png", dpi=300,
#                      bbox_inches="tight")

# 4. Normality: st.shapiro on raw; then np.log10-transform, histogram +
#    shapiro again. Save the log-scale histogram as fig2_hist_log.png

# 5. Q1 — treatment effect:
#    a) Welch t-test on RAW concentration (equal_var=False!)
#    b) Welch t-test on log10 concentration
#    c) median concentration by treatment (groupby)
#    d) broom-style table from the log-scale ols() fit:
#       params / bse / tvalues / pvalues in one DataFrame

# 6. Q2 — two-way ANOVA: ols("log10(conc) ~ C(treatment) * C(peptide)")
#    (drop NAs first), anova_lm(fit, typ=2)

# 7. Q3 — missingness vs batch: pd.crosstab(dat.batch, dat.conc.isna()),
#    then st.chi2_contingency — check the `expected` counts (>= 5?)

# 8. Q4 — power simulation: 500 simulations of n = 50/group drawn from
#    rng.lognormal(log(25), 0.6) (control) and rng.lognormal(log(25) + 0.4, 0.6)
#    (treated = +50%); Welch t-test each; report the fraction p < 0.05

# 9. Figures 3-4: boxplot of log10 conc by treatment (fig3), interaction
#    plot — mean log10 conc by peptide, one line per treatment (fig4)

# 10. Receipt: pd.show_versions(), datetime.date.today(), and
#     hashlib.md5(open("data/peptides_assay_clean.csv", "rb").read()).hexdigest()

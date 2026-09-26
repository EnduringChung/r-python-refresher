# Project 2 — Reference solution (Python)
# Run from the repo root:  python3 projects/p2/solution.py

import datetime
import hashlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as st
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

rng = np.random.default_rng(42)

dat = pd.read_csv("data/peptides_assay_clean.csv")

# 1. Structure ---------------------------------------------------------------
print(f"dimensions: {dat.shape[0]} rows x {dat.shape[1]} cols\n")
dat.info()
print(dat.describe())

# 2. Missingness -------------------------------------------------------------
print("\nNaNs per column:")
print(dat.isna().sum())
print("\nOnly concentration has NaNs — a missing well means the measurement "
      "failed (instrument/pipetting), not that the sample was absent.")

# 3-4. EDA: the skew, and the log fix ----------------------------------------
d = dat.dropna(subset=["concentration_ng_ul"]).copy()
d["log_conc"] = np.log10(d["concentration_ng_ul"])

# 3. raw histogram
fig, ax = plt.subplots(figsize=(5, 4))
ax.hist(dat["concentration_ng_ul"].dropna(),
        bins=np.arange(0, 160, 10), edgecolor="white")
ax.set(xlabel="concentration (ng/uL)", ylabel="wells",
       title=f"Right-skewed — classic lognormal\n"
             f"(Shapiro p={st.shapiro(d['concentration_ng_ul']).pvalue:.2g})")
plt.savefig("projects/p2/fig1_hist_raw.png", dpi=300, bbox_inches="tight")
plt.close()

# 4. log10 histogram
fig, ax = plt.subplots(figsize=(5, 4))
ax.hist(d["log_conc"], bins=np.arange(0.4, 2.2, 0.2), edgecolor="white")
ax.set(xlabel="log10 concentration", ylabel="wells",
       title=f"Log scale: approximately normal\n"
             f"(Shapiro p={st.shapiro(d['log_conc']).pvalue:.2f})")
plt.savefig("projects/p2/fig2_hist_log.png", dpi=300, bbox_inches="tight")
plt.close()

# 5. Q1 — treatment effect ---------------------------------------------------
print("\n== Q1: treatment effect ==")
t_raw = st.ttest_ind(d[d.treatment == "treated"].concentration_ng_ul,
                     d[d.treatment == "control"].concentration_ng_ul,
                     equal_var=False)
t_log = st.ttest_ind(d[d.treatment == "treated"].log_conc,
                     d[d.treatment == "control"].log_conc,
                     equal_var=False)
print(f"raw-scale Welch:  t={t_raw.statistic:.3f}  df={t_raw.df:.1f}  "
      f"p={t_raw.pvalue:.3g}")
print(f"log-scale Welch:  t={t_log.statistic:.3f}  df={t_log.df:.1f}  "
      f"p={t_log.pvalue:.3g}\n")
print(d.groupby("treatment").concentration_ng_ul.agg(["median", "count"]))

fit = ols("log_conc ~ C(treatment)", data=d).fit()
print(pd.DataFrame({"estimate": fit.params, "std.error": fit.bse,
                    "statistic": fit.tvalues, "p.value": fit.pvalues}))

# 6. Q2 — two-way ANOVA (log scale) ------------------------------------------
print("\n== Q2: two-way ANOVA, log10(conc) ~ treatment * peptide ==")
fit2 = ols("log_conc ~ C(treatment) * C(peptide)", data=d).fit()
print(anova_lm(fit2, typ=2))

# 7. Q3 — missingness vs batch ------------------------------------------------
print("\n== Q3: is missingness associated with batch? ==")
tab = pd.crosstab(dat["batch"], dat["concentration_ng_ul"].isna())
tab.columns = ["measured", "missing"]
print(tab)
chi2, p, dof, expected = st.chi2_contingency(
    pd.crosstab(dat["batch"], dat["concentration_ng_ul"].isna()))
print(f"chi2={chi2:.2f}  dof={dof}  p={p:.3g}")
print("expected counts (should be >= 5):")
print(pd.DataFrame(expected, index=tab.index, columns=tab.columns).round(1))

# 8. Q4 — power simulation ----------------------------------------------------
print("\n== Q4: power to detect a +50% effect (n = 50/group) ==")
n = 50
pvals = [st.ttest_ind(rng.lognormal(np.log(25) + 0.4, 0.6, n),
                      rng.lognormal(np.log(25), 0.6, n),
                      equal_var=False).pvalue
         for _ in range(500)]
print(f"empirical power: {np.mean(np.array(pvals) < 0.05)}")

# 9. Figures 3-4 ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5, 4))
groups = [d[d.treatment == g].log_conc for g in ["control", "treated"]]
ax.boxplot(groups)
ax.set_xticks([1, 2])
ax.set_xticklabels(["control", "treated"])
ax.set(ylabel="log10 concentration",
       title=f"Treatment effect, log scale (n = {len(d)} wells)")
plt.savefig("projects/p2/fig3_box_treatment.png", dpi=300, bbox_inches="tight")
plt.close()

means = d.groupby(["peptide", "treatment"]).log_conc.mean().unstack()
fig, ax = plt.subplots(figsize=(6, 4))
for trt in means.columns:
    ax.plot(means.index, means[trt], "o-", label=trt)
ax.set(ylabel="mean log10 concentration",
       title="Interaction check: parallel lines = no interaction")
ax.legend()
plt.savefig("projects/p2/fig4_interaction.png", dpi=300, bbox_inches="tight")
plt.close()

# 10. Receipt ------------------------------------------------------------------
print("\n== Receipt ==")
print("run at:", datetime.datetime.now(datetime.timezone.utc).isoformat())
print("data md5:", hashlib.md5(
    open("data/peptides_assay_clean.csv", "rb").read()).hexdigest())
pd.show_versions()

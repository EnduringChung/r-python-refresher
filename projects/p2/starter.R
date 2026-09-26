# Project 2 — Peptide assay EDA + statistics (R starter)
# Run from the repo root:  Rscript projects/p2/starter.R
# Data: data/peptides_assay_clean.csv (your Project 1 output, committed)

library(tidyverse)
library(broom)

set.seed(42)   # every random step below is reproducible because of this line

dat <- read_csv("data/peptides_assay_clean.csv", show_col_types = FALSE)

# 1. Structure: dimensions + glimpse()

# 2. Missingness: NAs per column; reflect — why does concentration have
#    NAs but nothing else does? (think: what does a missing well mean?)

# 3. EDA histogram: concentration, raw scale. Choose binwidth ON PURPOSE.
#    Save: ggsave("projects/p2/fig1_hist_raw.png", width = 6, height = 4)

# 4. Normality: shapiro.test on raw; then log10-transform, histogram +
#    shapiro again. Save the log-scale histogram as fig2_hist_log.png

# 5. Q1 — treatment effect:
#    a) Welch t-test on RAW concentration (t.test(conc ~ treatment, data))
#    b) Welch t-test on log10 concentration
#    c) median concentration by treatment (aggregate)
#    d) tidy() the log-scale test into a coefficient-style table

# 6. Q2 — two-way ANOVA: log10(conc) ~ treatment * peptide
#    (drop NAs first), tidy() the result

# 7. Q3 — missingness vs batch: build the 3x2 table of batch x is.na(conc),
#    run chisq.test (read the warning if you get one — why is it OK here?)

# 8. Q4 — power simulation: 500 simulations of n = 50/group drawn from
#    rlnorm(log(25), 0.6) (control) and rlnorm(log(25) + 0.4, 0.6)
#    (treated = +50%); Welch t-test each; report the fraction p < 0.05

# 9. Figures 3-4: boxplot of log10 conc by treatment (fig3), interaction
#    plot — mean log10 conc by peptide, one line per treatment (fig4)

# 10. Receipt: print sessionInfo(), Sys.time(), and the md5sum of the data
#     file (tools::md5sum("data/peptides_assay_clean.csv"))

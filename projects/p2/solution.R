# Project 2 — Reference solution (R)
# Run from the repo root:  Rscript projects/p2/solution.R

library(tidyverse)
library(broom)

set.seed(42)

dat <- read_csv("data/peptides_assay_clean.csv", show_col_types = FALSE)

## 1. Structure ----
cat("dimensions:", nrow(dat), "rows x", ncol(dat), "cols\n\n")
glimpse(dat)

## 2. Missingness ----
cat("\nNAs per column:\n")
print(colSums(is.na(dat)))
cat("\nOnly concentration has NAs — a missing well means the measurement
failed (instrument/pipetting), not that the sample was absent.\n\n")

## 3-4. EDA: the skew, and the log fix ----
p1 <- ggplot(dat, aes(concentration_ng_ul)) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "white") +
  labs(x = "concentration (ng/uL)", y = "wells",
       title = "Right-skewed — classic lognormal",
       subtitle = paste0("Shapiro-Wilk p = ",
                         format(shapiro.test(dat$concentration_ng_ul)$p.value,
                                digits = 2)))

dat_log <- dat |>
  mutate(log_conc = log10(concentration_ng_ul))

p2 <- ggplot(dat_log, aes(log_conc)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "white") +
  labs(x = "log10 concentration", y = "wells",
       title = "Log scale: approximately normal",
       subtitle = paste0("Shapiro-Wilk p = ",
                         format(shapiro.test(na.omit(dat_log$log_conc))$p.value,
                                digits = 2)))

ggsave("projects/p2/fig1_hist_raw.png", p1, width = 6, height = 4, dpi = 300)
ggsave("projects/p2/fig2_hist_log.png", p2, width = 6, height = 4, dpi = 300)

## 5. Q1 — treatment effect ----
cat("\n== Q1: treatment effect ==\n")
d <- dat_log |> filter(!is.na(log_conc))

tt_raw <- t.test(concentration_ng_ul ~ treatment, data = d)
tt_log <- t.test(log_conc ~ treatment, data = d)
cat("raw-scale Welch: t =", round(tt_raw$statistic, 3),
    ", df =", round(tt_raw$parameter, 1),
    ", p =", format(tt_raw$p.value, digits = 3), "\n")
cat("log-scale Welch: t =", round(tt_log$statistic, 3),
    ", df =", round(tt_log$parameter, 1),
    ", p =", format(tt_log$p.value, digits = 3), "\n\n")

print(dat |> group_by(treatment) |>
        summarise(median_conc = median(concentration_ng_ul, na.rm = TRUE),
                  n = sum(!is.na(concentration_ng_ul)), .groups = "drop"))

print(tidy(tt_log))

## 6. Q2 — two-way ANOVA (log scale) ----
cat("\n== Q2: two-way ANOVA, log10(conc) ~ treatment * peptide ==\n")
print(tidy(aov(log_conc ~ treatment * peptide, data = d)))

## 7. Q3 — missingness vs batch ----
cat("\n== Q3: is missingness associated with batch? ==\n")
tab <- table(dat$batch, is.na(dat$concentration_ng_ul),
             dnn = list("batch", "missing"))
print(tab)
ct <- chisq.test(tab)
cat("chi2 =", round(ct$statistic, 2), ", df =", ct$parameter,
    ", p =", format(ct$p.value, digits = 3), "\n")
print(round(ct$expected, 1))   # check expected >= 5 yourself

## 8. Q4 — power simulation ----
cat("\n== Q4: power to detect a +50% effect (n = 50/group) ==\n")
n <- 50
pvals <- replicate(500, {
  t.test(rlnorm(n, log(25) + 0.4, 0.6),
         rlnorm(n, log(25),      0.6))$p.value
})
cat("empirical power:", mean(pvals < 0.05), "\n")

## 9. Figures 3-4 ----
p3 <- ggplot(d, aes(treatment, log_conc, fill = treatment)) +
  geom_boxplot(show.legend = FALSE, outlier.alpha = 0.4) +
  labs(x = NULL, y = "log10 concentration",
       title = "Treatment effect, log scale",
       subtitle = paste0("n = ", nrow(d), " wells"))

p4 <- d |>
  group_by(treatment, peptide) |>
  summarise(m = mean(log_conc), .groups = "drop") |>
  ggplot(aes(peptide, m, color = treatment, group = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(x = NULL, y = "mean log10 concentration",
       title = "Interaction check: parallel lines = no interaction")

ggsave("projects/p2/fig3_box_treatment.png", p3, width = 6, height = 4, dpi = 300)
ggsave("projects/p2/fig4_interaction.png", p4, width = 6, height = 4, dpi = 300)

## 10. Receipt ----
cat("\n== Receipt ==\n")
cat("run at:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
cat("data md5:", unname(tools::md5sum("data/peptides_assay_clean.csv")), "\n")
sessionInfo()

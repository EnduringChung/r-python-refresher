# Project 1 — Reference solution (R)
library(tidyverse)
library(lubridate)

raw <- read_csv("data/peptides_assay.csv", show_col_types = FALSE)
n_raw <- nrow(raw)

clean <- raw |>
  rename(                                                     # 1
    sample_id = `Sample ID`,
    treatment = Treatment,
    peptide = Peptide,
    concentration_ng_ul = `Concentration (ng/uL)`,
    replicate = Replicate,
    batch = Batch,
    collection_date = `Collection Date`
  ) |>
  mutate(
    sample_id = str_to_lower(str_trim(sample_id)),            # 2
    treatment = case_match(                                   # 3
      str_to_lower(str_trim(treatment)),
      c("control", "ctrl") ~ "control",
      "treated" ~ "treated"
    ),
    concentration_ng_ul = na_if(as.character(concentration_ng_ul), "n/a"),
    concentration_ng_ul = na_if(concentration_ng_ul, ""),     # 4
    concentration_ng_ul = na_if(concentration_ng_ul, "-999"),
    concentration_ng_ul = as.numeric(concentration_ng_ul),
    collection_date = as_date(                                # 5
      parse_date_time(
        collection_date,
        orders = c("ymd", "mdy", "Bd, Y", "B d, Y")
      )
    )
  ) |>
  distinct() |>                                               # 6
  mutate(                                                     # 7
    across(c(treatment, peptide, batch), as.factor),
    treatment = factor(treatment, levels = c("control", "treated"))
  )

stopifnot(
  "sample_id should have no whitespace" = !any(str_detect(clean$sample_id, "^\\s|\\s$")),
  "treatment has no NAs" = !any(is.na(clean$treatment)),
  "treatment levels fixed" = setequal(levels(clean$treatment), c("control", "treated")),
  "no sentinel -999 left" = !any(clean$concentration_ng_ul == -999, na.rm = TRUE),
  "dates all parsed" = !any(is.na(clean$collection_date))
)

write_csv(clean, "data/peptides_assay_clean.csv")

cat("rows before:", n_raw, "| after:", nrow(clean),
    "| duplicates removed:", n_raw - nrow(clean), "\n\n")
cat("NAs per column:\n")
print(colSums(is.na(clean)))

cat("\nmean concentration (ng/uL) by treatment x peptide:\n")
clean |>
  group_by(treatment, peptide) |>
  summarise(mean_conc = mean(concentration_ng_ul, na.rm = TRUE),
            n = sum(!is.na(concentration_ng_ul)),
            .groups = "drop") |>
  pivot_wider(names_from = treatment, values_from = c(mean_conc, n))

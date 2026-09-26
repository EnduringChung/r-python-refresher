# Project 1 — Clean the peptide assay (R starter)
# Goal: turn data/peptides_assay.csv into a tidy, typed, deduplicated tibble.
# Run from the repo root:  Rscript projects/p1/starter.R

library(tidyverse)

raw <- read_csv("data/peptides_assay.csv", show_col_types = FALSE)

# 1. Clean column names (they have spaces, caps, parentheses)

# 2. Tidy `sample_id`: trim whitespace, lowercase  -> "sample_038"

# 3. Tidy `treatment`: collapse control/CTRL/Control -> "control",
#    treated/treated/Treated -> "treated"

# 4. `concentration_ng_ul`: convert "n/a" and "" to NA, treat -999 as NA,
#    then make it numeric

# 5. `collection_date`: THREE formats mixed together:
#    "2024-03-01", "03/10/2024", "March 11, 2024"  -> parse to Date
#    (hint: readr::parse_date_time + lubridate::as_date)

# 6. Remove exact duplicate rows; report how many you removed

# 7. Convert Treatment, Peptide, Batch to factors; set treatment levels
#    to c("control", "treated")

# 8. Save to data/peptides_assay_clean.csv
# 9. Print: rows before/after, NAs per column, mean concentration by
#    treatment x peptide (watch: some combos may be all-NA!)

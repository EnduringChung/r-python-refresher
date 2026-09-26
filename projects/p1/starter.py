# Project 1 — Clean the peptide assay (Python starter)
# Goal: turn data/peptides_assay.csv into a tidy, typed, deduplicated DataFrame.
# Run from the repo root:  python3 projects/p1/starter.py

import pandas as pd

raw = pd.read_csv("data/peptides_assay.csv")

# 1. Normalize column names -> snake_case, e.g. "Concentration (ng/uL)"
#    becomes concentration_ng_ul  (hint: str.strip + str.lower +
#    regex replace non-alphanumerics)

# 2. sample_id: strip whitespace, lowercase -> "sample_038"

# 3. treatment: lowercase, then collapse {"control", "ctrl"} -> "control"

# 4. concentration_ng_ul: replace "n/a" and "" and -999 with NaN, cast float

# 5. collection_date: parse THREE mixed formats to datetime
#    ("2024-03-01", "03/10/2024", "March 11, 2024")
#    (hint: pd.to_datetime(..., format="mixed") needs pandas >= 2.0)

# 6. Drop exact duplicate rows; report how many you removed

# 7. Convert treatment, peptide, batch to pandas Categorical;
#    treatment ordered ["control", "treated"]

# 8. Save to data/peptides_assay_clean_py.csv (no index)
# 9. Print: rows before/after, NaNs per column, mean concentration by
#    treatment x peptide (pivot_table)

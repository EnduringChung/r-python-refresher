# Project 1 — Reference solution (Python)
import re

import numpy as np
import pandas as pd

raw = pd.read_csv("data/peptides_assay.csv")
n_raw = len(raw)


def snake_case(col: str) -> str:
    col = col.strip().lower()
    return re.sub(r"[^a-z0-9]+", "_", col).strip("_")


clean = raw.rename(columns=snake_case)                          # 1

clean["sample_id"] = clean["sample_id"].str.strip().str.lower()  # 2
clean["treatment"] = clean["treatment"].str.strip().str.lower()
clean["treatment"] = clean["treatment"].replace({"ctrl": "control"})  # 3
clean["concentration_ng_ul"] = (                                # 4
    clean["concentration_ng_ul"]
    .replace(["n/a", "", -999], np.nan)
    .astype(float)
)
clean["collection_date"] = pd.to_datetime(                      # 5
    clean["collection_date"], format="mixed"
)

clean = clean.drop_duplicates().reset_index(drop=True)          # 6

clean["treatment"] = pd.Categorical(                            # 7
    clean["treatment"], categories=["control", "treated"]
)
for col in ["peptide", "batch"]:
    clean[col] = clean[col].astype("category")

# ---- guards (same spirit as R stopifnot) ----
assert not clean["sample_id"].str.match(r"^\s|\s$").any()
assert not clean["treatment"].isna().any()
assert set(clean["treatment"].cat.categories) == {"control", "treated"}
assert not (clean["concentration_ng_ul"] == -999).any()
assert not clean["collection_date"].isna().any()

clean.to_csv("data/peptides_assay_clean_py.csv", index=False)   # 8

print(f"rows before: {n_raw} | after: {len(clean)} | "
      f"duplicates removed: {n_raw - len(clean)}\n")
print("NaNs per column:")
print(clean.isna().sum())

print("\nmean concentration (ng/uL) by treatment x peptide:")
print(clean.pivot_table(
    index="peptide",
    columns="treatment",
    values="concentration_ng_ul",
    aggfunc="mean",
    observed=True,
))

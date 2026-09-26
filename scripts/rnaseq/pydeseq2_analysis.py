# PyDESeq2 differential expression — Python twin of deseq2_analysis.R.
# Requires: pip install pydeseq2  (does NOT run in the browser — local only)
# Run from repo root:  python3 scripts/rnaseq/pydeseq2_analysis.py

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts = pd.read_csv("data/rnaseq/counts.csv", index_col=0).T          # samples × genes
coldata = pd.read_csv("data/rnaseq/coldata.csv", index_col=0)

dds = DeseqDataSet(counts=counts, metadata=coldata, design_factors="condition", seed=2024)
dds.deseq2()

stats = DeseqStats(dds, contrast=("condition", "treated", "control"))
stats.summary()

res = stats.results_df.reset_index().rename(columns={"index": "gene"})
res.to_csv("scripts/rnaseq/output/pydeseq2_results.csv", index=False)

sig = res[res["padj"] < 0.05]
print(f"PyDESeq2: {len(sig)} genes with padj < 0.05")

planted = [f"GENE{i:03d}" for i in range(1, 41)]
print("planted DE genes recovered:", sum(g in set(sig["gene"]) for g in planted), "/ 40")

print(sig.sort_values("padj").head(10)[["gene", "log2FoldChange", "padj"]].to_string(index=False))

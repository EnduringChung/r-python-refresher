# RNA-seq frozen-output pipeline (Week 4)

Bioconductor packages (DESeq2, edgeR) cannot run in the browser — WebAssembly
builds don't exist. This folder is the source of truth for the Week 4
"🧊 frozen output" lessons:

```
learner sees in lesson:   real code + REAL pre-computed output (frozen)
this folder:              the scripts that produced it (runnable locally)
```

## Files

| File | What it does |
|---|---|
| `make_synthetic_counts.R` | generates `../../data/rnaseq/counts.csv` + `coldata.csv` (seed-fixed, planted DE genes) |
| `deseq2_analysis.R` | full DESeq2 run → results CSV + MA/volcano/PCA PNGs in `output/` |
| `edger_analysis.R` | same analysis in edgeR (TMM + quasi-likelihood) |
| `pydeseq2_analysis.py` | same analysis in PyDESeq2 (needs `pip install pydeseq2`) |

**Status:** DESeq2 + edgeR outputs are generated and committed in `output/`
(40 significant genes, 39/40 planted DE genes recovered by both tools — the
cross-tool consistency is a teaching point). PyDESeq2 script is authored and
syntax-checked but its output is NOT yet generated; run it locally after
`pip install pydeseq2` and commit `output/pydeseq2_results.csv`.

## To regenerate outputs

```bash
Rscript scripts/rnaseq/make_synthetic_counts.R
Rscript scripts/rnaseq/deseq2_analysis.R
Rscript scripts/rnaseq/edger_analysis.R
python3 scripts/rnaseq/pydeseq2_analysis.py   # optional, needs pydeseq2
```

Outputs land in `scripts/rnaseq/output/`. The Week 4 author embeds those
artifacts (numbers + PNGs) into the lesson `.qmd` as static output, labeled
🧊, and links the matching script from `scripts/` so learners can run the
real thing on their Mac.

# R + Python Refresher — Master Plan

A 4-week (+1 optional capstone week), 1–2 hrs/day refresher course, authored as Quarto
Live notebooks: read on an iPad in Safari, type code in the browser, run it, see output.
Free, static hosting only, works on Mac/Windows too.

---

## 1. Who this is for (learner profile)

- **R**: knows basics + tidyverse/base R package names, but most written R is
  AI-generated. Main pain: *"R syntax is too flexible — I can't recognize the patterns
  for good code."*
- **Python**: comfortable coding, but weak on classes/`__init__`, OOP architecture.
- **Algorithms**: LeetCode-level beginner.
- **Domains**: statistics, ML/AI, bioinformatics (RNA-seq, single-cell, sequences,
  trajectories, pseudotime, HMMs, prediction, multi-omics).
- **Goal**: step-by-step habits for writing **robust, reproducible** code in both languages.

## 2. Platform — Quarto Live (investigated & verified)

**Extension**: [r-wasm/quarto-live](https://github.com/r-wasm/quarto-live) — Quarto
format `live-html`. Code cells run **in the browser** via **webR** (R on WebAssembly)
and **Pyodide** (Python on WebAssembly). Hosting is any static service — **GitHub
Pages** (free). No Jupyter server, no paid service, no native app. Safari on iPad is a
supported target; work persists in localStorage between sessions.

**What I verified against the actual package indexes (Sep 2026):**

| Runs in-browser (interactive cells) | Does NOT run in-browser |
|---|---|
| R: full tidyverse (dplyr, ggplot2, tidyr, readr, purrr, stringr, forcats, lubridate, broom, patchwork, gt), **entire tidymodels stack**, data.table, testthat, MASS, survival, glmnet, randomForest, e1071, cluster — 23,866 wasm CRAN packages total | R: **all Bioconductor** — DESeq2, edgeR, limma, SingleCellExperiment, scater/scran. (Seurat 5.5 is listed in the wasm repo but untested & heavy — assume local-only) |
| Python: numpy, pandas, matplotlib, scipy, scikit-learn, statsmodels, **biopython**, **pysam** (BAM/VCF!), screed (FASTQ), networkx, pytest, polars, duckdb, h5py, lightgbm, xgboost | Python: **scanpy/anndata workflows with numba**, seaborn is not bundled (but installs via micropip), torch/tensorflow, hmmlearn |

**Consequence — three lesson execution tiers:**

| Tier | How it works | Used for |
|---|---|---|
| **T1 — Live** | `{webr}` / `{pyodide}` cells, fully editable & runnable on iPad | Weeks 1–3 (everything), W4 ML, sequence bioinformatics |
| **T2 — Frozen output** | Real code + **pre-computed outputs embedded**; provided as runnable local scripts | DESeq2/edgeR lessons (Bioconductor) |
| **T3 — Colab/local** | Lesson text + live demo cells where possible + one-click **Google Colab** badge (free) for scanpy/single-cell; local RStudio/Positron scripts for R Bioconductor | Single-cell, trajectories, multi-omics, capstone data work |

This satisfies the fallback requirement: Quarto Live covers ~80% of the course
interactively; heavy bioinformatics uses frozen outputs + Colab/local scripts as the
simplest free alternative (no server to maintain).

**iPad reading loop**: open GitHub Pages site in Safari → read → edit live cell → run →
output appears inline. Solutions are collapsible (`<details>`), exercises have hints.

## 3. Course map

**Format per day (~1–2 h):** every lesson opens with a **cheat sheet** section, then
**drills** (10–15 one-liners to rebuild reflexes), then the **side-by-side lesson**
(R left, Python right, same task), then **exercises with solutions**. Section ends =
**project**. Friday lessons include a **"read & refactor AI code" drill** — messy
AI-generated R/Python refactored into the canonical pattern (directly targets how you
actually work).

### Week 1 — Language core, side-by-side
| Day | Topic | Notes |
|---|---|---|
| 1 | **Setup + The 7 Patterns of R** | Install toolchain; why R feels chaotic; the 7 canonical patterns (vectorize → pipe → df-in/fun-out → formula → `{{}}` → purrr/functionals → S3 dispatch). "Pick one dialect (tidyverse-first), read base fluently." Python: the zen equivalents (comprehensions, EAFP, underscore conventions) |
| 2 | Data structures | vectors/lists vs lists/dicts/tuples; attributes vs attributes; indexing & slicing (`[`, `[[`, `$` vs `[]`, `.loc`, `.iloc`); NA vs None vs NaN |
| 3 | Control flow + functions | if/else/switch, loops when allowed; function anatomy, defaults, `...` vs `*args/**kwargs`; scope & closures; pure functions |
| 4 | **Algorithm Gym #1** | Big-O intuition; arrays & hash maps; 6 starter problems with solutions (two-sum, contains-duplicate, etc.) |
| 5 | Strings, factors, dates + **Project 1** | stringr vs str methods; factor/forcats vs categorical dtype; lubridate vs datetime. **P1**: parse a messy CSV → clean tidy dataset, same result in both languages |

### Week 2 — Data wrangling, visualization, statistics
| Day | Topic | Notes |
|---|---|---|
| 1 | Import/export + tidy data | readr/readxl vs pandas IO; tidy rules; pivot longer/wider; joins (dplyr vs merge/join) |
| 2 | **The grammar of wrangling** | dplyr's 6 verbs ↔ pandas chain (`.query .assign .loc .groupby .agg`); row-wise ops; the "one obvious way" patterns each language wants |
| 3 | Visualization | ggplot2 grammar ↔ matplotlib/seaborn (+ plotnine as the ggplot bridge); the 5 charts that matter; when in-browser plotting works |
| 4 | **Algorithm Gym #2** | Two pointers, sliding window, prefix sums; string problems |
| 5 | Statistics + **Project 2** | Distributions, t-tests, ANOVA, chi², correlation; linear/logistic regression: `lm()/glm()` formula interface ↔ statsmodels formula API; **broom::tidy() ↔ .summary tables**. **P2**: full EDA + statistical report rendered as Quarto HTML (reproducibility thread starts: set.seed, session info) |

### Week 3 — Software engineering + reproducibility
| Day | Topic | Notes |
|---|---|---|
| 1 | **Python classes & OOP** (your weak spot) | `class`, `__init__`, self, instance vs class attrs; `__repr__`, dunders worth knowing; dataclasses; properties. Mirror: R's S3 in 30 min (because you'll read it), R6 as the "real class" |
| 2 | Architecture & robust code | Modules & packages (both languages); virtual envs: **uv** (Python) + **renv** (R); error handling try/except vs tryCatch; logging; type hints; arguments validation (argcheck/stopifnot) |
| 3 | **Reproducibility full stack** | The analysis project template (data/ code/ output/ renv+pyproject); git & GitHub in 90 min; seeds everywhere; **tests: testthat ↔ pytest**; "hand-off checklist" for every analysis |
| 4 | **Algorithm Gym #3** | Stacks, queues, recursion, binary search; intro to dynamic programming (memoization) — with bio flavor (Fibonacci → Fibonacci-like sequence counting) |
| 5 | **Project 3** | Build `bioseqkit`: a tiny R package (testthat + renv) AND a tiny Python package (pytest + uv) that do the same thing (GC content, k-mer counts, reverse complement) — with CI running tests on GitHub |

### Week 4 — Machine learning + bioinformatics I
| Day | Topic | Notes |
|---|---|---|
| 1 | ML foundations side-by-side | Train/test, cross-validation; **tidymodels workflows ↔ sklearn Pipelines**; preprocessing as code (recipes ↔ ColumnTransformer); seeds & the reproducible split |
| 2 | Classifiers & evaluation | Logistic regression, random forest, xgboost/lightgbm; ROC/PR, confusion matrices, calibration; class imbalance; feature importance. Live in-browser (sklearn + tidymodels both wasm ✓) |
| 3 | Sequences & genomic data | FASTA/FASTQ with Biopython + screed (live); k-mers, GC, ORFs; **pysam for BAM/VCF basics (live!)**; GenomicRanges mental model + local R script |
| 4 | **Algorithm Gym #4** | Graphs: BFS/DFS adjacency (networkx live), topological sort; HMM part 1: the Viterbi idea on CpG islands implemented from scratch in numpy (live) |
| 5 | Bulk RNA-seq + **Project 4 kick-off** | Counts → normalization → NB GLM → FDR: **DESeq2/edgeR as frozen-output lessons + local scripts**; PyDESeq2 via Colab badge; variance-stabilized QC plots (ggplot live). **P4**: DE analysis on a public dataset (frozen R) + classifier on the results (live Python) |

### Week 5 (optional capstone) — Bioinformatics II + multi-omics
| Day | Topic | Notes |
|---|---|---|
| 1 | Single-cell core | Scanpy workflow (Colab/local): QC filtering → normalization → HVG → PCA → neighbors → Leiden → UMAP; Seurat comparison script (local); wasm-Seurat flagged as experimental |
| 2 | Trajectory & pseudotime | Diffusion pseudotime & PAGA in scanpy; monocle3/slingshot concepts; interpreting vs over-interpreting trajectories |
| 3 | HMMs + sequence prediction | Full CpG-island HMM (train with Baum-Welch conceptually, Viterbi from scratch — you wrote it in W4); k-mer classifiers for promoter/pathogen sequence prediction (live sklearn) |
| 4 | Multi-omics integration | Concatenation, MOFA/DIABLO concepts (frozen), a live sklearn CCA/PA-geometry mini-example; batch effect awareness |
| 5 | **Capstone** | End-to-end mini-study: public multi-omics-ish dataset → QC → DE → classifier → report, using the full reproducibility stack, published as a Quarto site |

## 4. Standing threads (appear in every week)

1. **Pattern boxes** 🧩 — each R topic reduced to one canonical pattern + "recognize
   it in the wild" (reading messy/AI code) + refactor drill.
2. **Side-by-side tables** — every operation shown in R and Python.
3. **Reproducibility checklist** — seeds, envs, tests, session info; introduced W2,
   enforced in every project, mastered W3D3.
4. **Algorithm Gym** — one dedicated day per week, progressive (complexity →
   two-pointers → recursion/DP → graphs/HMM), all with bio-flavored examples where natural.
5. **Colab/local badges** — any cell that can't run in-browser is labeled with how to run it.

## 5. Repo structure

```
R and Python Refresher/
├── PLAN.md                  ← this file
├── _quarto.yml              ← site config (live-html, navbar by week)
├── _extensions/             ← quarto-live (quarto add r-wasm/quarto-live)
├── index.qmd                ← course home + how-to-use-on-iPad
├── setup.qmd                ← toolchain install (Mac/Windows) + smoke test page
├── cheatsheets/             ← master index of all cheat sheets (printable)
├── w1/ ... w5/              ← one .qmd per day (lesson+drills+exercises+project)
├── projects/                ← project briefs + starter files + solutions
├── solutions/               ← drill & exercise solutions (also collapsible in-page)
├── scripts/                 ← runnable local .R / .py for T2/T3 lessons
├── data/                    ← small datasets (checked in) + download scripts for big ones
└── .github/workflows/       ← render + deploy to GitHub Pages on push; test CI for P3
```

## 6. Build phases

| Phase | Deliverable | Verify |
|---|---|---|
| 0 | Scaffold: `_quarto.yml`, extension install, `setup.qmd` with live R + Python smoke-test cells | `quarto preview` locally; **you test one page on iPad Safari** before mass content build |
| 1 | GitHub repo + Actions → GitHub Pages auto-deploy | Site loads on iPad, live cells run, localStorage persists |
| 2 | Weeks 1–2 (10 lessons, 2 projects, cheat sheets) | You complete W1; feedback loop |
| 3 | Week 3 (OOP, architecture, repro stack, P3 with CI) | P3 CI green |
| 4 | Weeks 4–5 (ML, bioinformatics, frozen-output lessons, Colab badges, capstone) | Colab badges work; frozen outputs match scripts |

## 7. Decisions (resolved)

1. **Hosting**: GitHub Pages, repo `EnduringChung/r-python-refresher`, deploy via Actions on push. Live at https://enduringchung.github.io/r-python-refresher/
2. **Week 5 capstone**: included
3. **Datasets**: mix — toy data checked in + real public datasets for projects
4. **Render engine**: default jupyter (Python on Mac); knitr available for frozen R lessons

---

## 8. Authoring standards (locked in after Week 1 — follow for Weeks 2–5)

### Per-day content checklist (the "depth standard")

- [ ] Cheat sheet table at top (60–90 sec scan)
- [ ] Prose explains the **why**, not just syntax (each pattern gets its motivation)
- [ ] **Drill pad**: ~12 tasks, easy → spicy, one live cell per language
- [ ] **5–6 graded exercises** (at least one Python), each with `.hint` + `.solution` divs linked via `exercise="key"`
- [ ] **Common pitfalls** section — real bugs only, harvested from verification runs
- [ ] **Self-quiz**: 5 predict-the-output questions + collapsed answers
- [ ] Wrap-up checklist with checkbox items

### Language parity rules

- Every concept shown in **both** dialects; Python blocks are **live `{pyodide}` cells, never static ` ```python `** (static only for: conceptual comments, solution listings, ❌ bad-code exhibits)
- Where feasible, R and Python examples produce **identical numbers** on the same data (Project 1 pattern: both solutions print 120 rows / 13 NAs / identical means)
- Pyodide cannot fetch URLs (`pd.read_csv(url)` fails) → inline data or checked-in files; keep in-browser datasets small (Safari memory ceiling)
- Python cells auto-print only the **last** bare expression → `print()` everything

### Code verification workflow (mandatory before any deploy)

1. Extract and run **all** R snippets via `Rscript`, all Python via `python3`, with assertions on claimed outputs
2. Answer keys are guilty until proven: 6 of Week 1's keys were wrong until executed. Never "verify by reasoning"
3. R `identical()` is type-strict (`c(2,3)` vs `c(2L,3L)` fails) — use `all.equal()` for numeric comparisons in tests
4. After every recode/factor relevel: count NAs (the `CTRL` → `"ctrl"` → NA class of bug)

### Exercise self-sufficiency rule (added after the W4/W5 NameError reports)

Every exercise cell MUST be runnable on a freshly-loaded page, without
running any lesson cell first. Quarto-live executes exercise code in its own
environment — lesson state does not carry over. For each exercise add a
**setup cell** immediately before it:

    ```{pyodide}
    #| setup: true
    #| exercise: <key>
    <reload data, rebuild splits/models the exercise needs>
    ```

- Setup cells re-run before every learner evaluation (fast: data loads in
  <1 s, small model fits a few seconds — keep them lean)
- Drill pads and quizzes: prefer comment prompts or inline definitions;
  never reference lesson variables in executable lines
- Verify by extracting and running each setup standalone (see §verification)

### Quarto Live syntax facts (learned the hard way)

- `.hint`/`.solution` divs MUST carry `exercise="key"`; solution code inside is a **static** fenced block, not a cell
- Fill-in blanks need **6+ underscores** (`______`)
- `::: {.callout-collapse}` doesn't exist — use `::: {.callout-tip collapse="true"}`
- `.quartoignore` unreliable for excluding files → use `project.render` list in `_quarto.yml`
- `{webr}`/`{pyodide}` cells share one session per page per language — `library()`/`import` once
- First pyodide package load (statsmodels etc.) can take ~30s — warn learners in a comment
- R facts that bit us: `wday(label=TRUE)` → "Thu" (3 letters); `table()["missing"]` → NA (Python `Counter` → 0); `read_csv` trims header spaces, `pd.read_csv` doesn't; `Date - POSIXct` = garbage, convert with `as_datetime()` first

### Deploy ritual (every content batch)

1. `quarto render` → grep for errors
2. Commit + push → `gh run watch` until success
3. `curl` the changed pages for HTTP 200 (watch URL typos — `.qmd.html` happened)
4. Sanity-grep deployed HTML for cell wiring (`webr-`, `pyodide-data`)

## 9. Remaining build queue

- [ ] **Week 2**: wrangling grammar, ggplot2↔matplotlib, statistics + Project 2 (Quarto report) — at §8 standard *(assigned to parallel builder)*
- [x] Master cheat-sheet index page (`cheatsheets/index.qmd`)
- [x] Print CSS (`styles.css` @media print)
- [x] Project 3 scaffolding: `projects/p3/bioseqkitr` (R pkg, installs, gc tests green) + `projects/p3/bioseqkit_py` (pytest green, stubs skip-marked) + `.github/workflows/bioseqkit.yml` CI (runs on `projects/p3/**` changes)
- [x] Week 4 RNA-seq plumbing: `scripts/rnaseq/` — synthetic counts (seed 2024), DESeq2 + edgeR scripts **run and outputs committed** (40 sig genes, 39/40 planted recovered by both), PyDESeq2 script authored (needs `pip install pydeseq2` locally to generate its output)
- [x] **Week 3 content**: OOP/architecture lessons, repro stack, Algorithm Gym #3, W3D5 lesson wiring up the p3 scaffold (unskip tests → CI green)
- [x] **Week 4 content**: ML lessons, sequences, frozen DESeq2/edgeR lessons embedding `scripts/rnaseq/output/` artifacts, Algorithm Gym #4, Project 4 — P4 dataset recalibrated (10-gene signature, depth jitter; logistic CV ~0.76, boost needs `min_child_samples=5` on n=45)
- [ ] **Week 5 content**: single-cell (Colab badges), trajectories/pseudotime, HMMs from scratch, multi-omics, capstone *(assigned to parallel builder — note W4D4 built the from-scratch Viterbi in R+Python already; extend, don't duplicate)*
- [ ] Consider Posit Cloud path for R-Bioconductor lessons on iPad (open decision)

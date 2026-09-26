# Capstone data generator — synthetic two-view mini-study.
# Run from the repo root:  python3 scripts/capstone/generate_data.py
# Seed 2026; writes three checked-in files under data/.
#
# Planted truth (for the grading audit — learners: read this AFTER your
# own analysis):
#   - 40 RNA genes are up in treated (3x) -> discoverable at FDR < 0.05
#   - 12 proteins follow the same condition contrast
#   - cell_type (A/B/C) marker profiles are the LOUDEST structure in RNA
#     (PC1 = cell type, not condition — the EDA surprise)
#   - batch B2 shifts 15 RNA genes (a smaller, nastier batch effect)
#   - 4 samples have broken library sizes (the QC catch)

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 2026
N = 180
GENES = 1000
PROTEINS = 100


def build(seed: int = SEED):
    """Generate the capstone dataset; returns data + the planted truth."""
    rng = np.random.default_rng(seed)

    sample_id = [f"S{i:03d}" for i in range(1, N + 1)]
    condition = np.array(["control", "treated"] * (N // 2))
    rng.shuffle(condition)
    cell_type = rng.choice(["A", "B", "C"], N, p=[0.4, 0.35, 0.25])
    batch = np.array(["B1", "B2"] * (N // 2))
    rng.shuffle(batch)

    meta = pd.DataFrame({
        "sample_id": sample_id,
        "condition": condition,
        "cell_type": cell_type,
        "batch": batch,
    })

    # ---------------------------------------------------------------- RNA
    base = rng.lognormal(0.0, 1.0, GENES)                # gene abundance
    depth = 1e4 * rng.lognormal(0.0, 0.35, N)            # library sizes

    # break 4 samples' library sizes (QC catch: 3 crash, 1 doubles)
    broken = rng.choice(N, 4, replace=False)
    depth[broken[:3]] *= 0.02
    depth[broken[3:]] *= 2.0

    lam = depth[:, None] * base[None, :] * np.exp(
        (0.35 * np.where(condition == "treated", 1.0, 0.0))[:, None] / 2)

    # cell types = marker profiles: 150 genes shift per type (gene-selective —
    # a global shift on all genes would be normalized away as "depth")
    ct_genes = rng.choice(GENES, 150, replace=False)
    mult = np.ones((3, GENES))
    for t, ct in enumerate(["A", "B", "C"]):
        lo, hi = [(1.5, 3.0), (0.9, 1.4), (0.3, 0.9)][t]
        mult[t, ct_genes] = rng.uniform(lo, hi, 150)
    ct_multiplier = np.array([mult[["A", "B", "C"].index(c)] for c in cell_type])
    lam = lam * ct_multiplier

    # planted DE: 40 genes up in treated (3x — a strong, findable effect)
    de_genes = rng.choice(GENES, 40, replace=False)
    is_treated = (condition == "treated")[:, None]
    lam[:, de_genes] = lam[:, de_genes] * np.where(is_treated, 3.0, 1.0)
    # planted batch effect: 15 genes shifted in B2
    batch_genes = rng.choice(np.setdiff1d(np.arange(GENES), de_genes),
                             15, replace=False)
    lam[:, batch_genes] = lam[:, batch_genes] * np.where(
        (batch == "B2")[:, None], 1.3, 1.0)

    counts = rng.poisson(lam / 50)                       # scaling for realism
    # dropout: low-count entries zero out more often
    drop = rng.random(counts.shape) < 0.08 * (1 - counts / counts.max())
    counts[drop] = 0

    rna = pd.DataFrame(counts, columns=[f"GENE{i:04d}" for i in range(GENES)])
    rna.insert(0, "sample_id", sample_id)

    # ------------------------------------------------------------ protein
    # protein background is condition-free; only the planted 12 shift
    prot = rng.normal(10, 1.5, (N, PROTEINS))
    de_prots = rng.choice(PROTEINS, 12, replace=False)
    prot[:, de_prots] += np.where(condition == "treated", 1.0, 0.0)[:, None]
    prot = pd.DataFrame(prot, columns=[f"PROT{i:03d}" for i in range(PROTEINS)])
    prot.insert(0, "sample_id", sample_id)

    truth = {
        "de_genes": sorted(int(g) for g in de_genes),
        "de_proteins": sorted(int(p) for p in de_prots),
        "batch_genes": sorted(int(g) for g in batch_genes),
        "broken_samples": [sample_id[b] for b in broken],
    }
    return meta, rna, prot, truth


if __name__ == "__main__":
    meta, rna, prot, truth = build(SEED)

    out = Path("data")
    out.mkdir(exist_ok=True)
    meta.to_csv(out / "capstone_meta.csv", index=False)
    rna.to_csv(out / "capstone_rna_counts.csv", index=False)
    prot.to_csv(out / "capstone_protein.csv", index=False)

    counts = rna.iloc[:, 1:].to_numpy()
    print(f"meta:     {meta.shape}  ({meta.condition.value_counts().to_dict()})")
    print(f"rna:      {rna.shape}, totals {counts.sum(1).min():.0f}-"
          f"{counts.sum(1).max():.0f} (broken samples present)")
    print(f"protein:  {prot.shape}")
    print("\nSHA-256 receipts:")
    for f in ["capstone_meta.csv", "capstone_rna_counts.csv",
              "capstone_protein.csv"]:
        h = hashlib.sha256((out / f).read_bytes()).hexdigest()
        print(f"  {f}: {h[:32]}...")
    print("\nplanted: 40 DE genes, 12 DE proteins, 15 batch-shifted genes,",
          "4 broken libraries")

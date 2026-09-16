#!/usr/bin/env python3
"""Family-wise error across the whole six-stage screen, by permutation.

CLADE's stages are individually conservative and jointly convergent by design,
but that is an argument, not a number. This measures the thing actually
claimed: under no true candidate-phenotype relationship, how often does ANY
candidate in the family reach convergence?

The shuffle preserves lineage. Phenotype is permuted within each sequence type,
so every lineage keeps its exact resistant/susceptible mix and only the link to
the candidate is broken. A whole-cohort shuffle would destroy the population
structure that CLADE exists to handle, making the null trivially easy to pass
and understating the real false-positive risk.

Stages are evaluated sequentially, exactly as the framework applies them: a
candidate that fails Stage 1 or Stage 3 is never carried into the expensive
Stage 4 and Stage 5 evaluations. That is not an optimisation bolted on for
speed -- it is what the framework does, so the null has to do it too, or it
would be calibrating a different procedure.

Stage 1 uses an independent likelihood-ratio test, so no firthlogist is
required. Stage 5 uses conditional logistic regression stratified by control
genome, which handles the reuse that matched-with-replacement pairing creates.

Usage:
    fwer_permutation.py --genotypes g.csv --phenotype p.tsv --lineages l.tsv \
        --tree t.nwk --distances d.tsv --n-permutations 20 --outdir out/
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter("ignore")


def stratified_shuffle(y, strata, rng):
    out = y.copy()
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        out[idx] = rng.permutation(y[idx])
    return out


def stage1_lrt(X0, y, g):
    """Independent LRT for one candidate, plus separation status."""
    import statsmodels.api as sm
    from scipy.stats import chi2

    n11 = int(((g == 1) & (y == 1)).sum()); n10 = int(((g == 1) & (y == 0)).sum())
    n01 = int(((g == 0) & (y == 1)).sum()); n00 = int(((g == 0) & (y == 0)).sum())
    if min(n11, n10, n01, n00) == 0:
        return np.nan, np.nan, "complete"
    try:
        full = sm.Logit(y, sm.add_constant(np.hstack([X0, g.reshape(-1, 1)]),
                                           has_constant="add")).fit(disp=0, maxiter=200)
        red = sm.Logit(y, sm.add_constant(X0, has_constant="add")).fit(disp=0, maxiter=200)
        stat = 2.0 * (float(full.llf) - float(red.llf))
        if not np.isfinite(stat) or stat < 0:
            return np.nan, np.nan, "failed"
        return float(chi2.sf(stat, 1)), float(full.params[-1]), "none"
    except Exception:  # noqa: BLE001
        return np.nan, np.nan, "failed"


def run_screen(geno, y, strata, tree, dist, X0, alpha, baseline, margin, min_origins):
    """One full pass of the screen. Returns the convergent candidate count."""
    from statsmodels.stats.multitest import fdrcorrection

    # --- Stage 1 -----------------------------------------------------------
    names, pvals, coefs = [], [], []
    for c in geno.columns:
        g = geno[c].to_numpy(dtype=float)
        p, b, sep = stage1_lrt(X0, y, g)
        if sep != "none" or not np.isfinite(p):
            continue                      # unidentified: no p-value, cannot pass
        names.append(c); pvals.append(p); coefs.append(b)
    if not names:
        return 0, {}
    _, q = fdrcorrection(np.array(pvals), alpha=alpha)

    # --- Stage 3 (direction) applied to Stage 1 survivors ------------------
    survivors = [n for n, qq, b in zip(names, q, coefs) if qq < alpha and b > 0]
    if not survivors:
        return 0, {}

    # --- Stage 4 (timing) --------------------------------------------------
    from clade.validation.stage4_temporal_order import temporal_ordering

    pmap = dict(zip(geno.index, y))
    cmaps = {c: dict(zip(geno.index, geno[c])) for c in survivors}
    s4 = temporal_ordering(tree, pmap, cmaps, ambiguity_resolution=0).set_index("Candidate")
    survivors = [c for c in survivors
                 if s4.loc[c, "Total_gains"] >= min_origins
                 and s4.loc[c, "Pct_post_resistance"] >= baseline + margin]
    if not survivors:
        return 0, {}

    # --- Stage 5 (matched pairs, conditional) ------------------------------
    from clade.validation.stage5_matched_neighbors import nearest_negative_neighbor
    from lifelines import CoxPHFitter

    pos = [s for s in geno.index if pmap[s] == 1]
    neg = [s for s in geno.index if pmap[s] == 0]
    if not pos or not neg:
        return 0, {}
    m = nearest_negative_neighbor(dist.loc[geno.index, geno.index], pos, neg)
    m = m.reset_index()
    # The index name varies with how the genotype frame was loaded
    # (Sample_ID, positive, or unnamed), so take the first column positionally
    # rather than guessing at a name.
    pairs = m.rename(columns={m.columns[0]: "case", "nearest_neighbor": "control"})[
        ["case", "control"]
    ]
    passed = {}
    for c in survivors:
        gser = geno[c]
        rows = []
        for control, grp in pairs.groupby("control"):
            for s in list(grp["case"]) + [control]:
                v = gser.get(s)
                if pd.isna(v):
                    continue
                rows.append({"stratum": control, "is_case": int(s in set(grp["case"])),
                             "carrier": float(v)})
        df = pd.DataFrame(rows)
        if df.empty:
            continue
        keep = [s for s, gg in df.groupby("stratum")
                if gg.carrier.nunique() > 1 and gg.is_case.nunique() > 1]
        df = df[df.stratum.isin(keep)]
        if df.empty or df.carrier.nunique() < 2:
            continue
        df = df.assign(T=1.0, E=df.is_case.astype(int))
        try:
            cph = CoxPHFitter()
            cph.fit(df[["T", "E", "carrier", "stratum"]], duration_col="T",
                    event_col="E", strata=["stratum"], robust=True)
            p = float(cph.summary.loc["carrier", "p"])
            b = float(cph.summary.loc["carrier", "coef"])
            if p < alpha and b > 0:
                passed[c] = p
        except Exception:  # noqa: BLE001
            continue
    return len(passed), passed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--lineages", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--distances", required=True)
    ap.add_argument("--n-permutations", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--temporal-margin", type=float, default=0.0)
    ap.add_argument("--min-origins", type=int, default=2)
    ap.add_argument("--top-lineages", type=int, default=20)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    import dendropy

    geno = pd.read_csv(args.genotypes, index_col=0)
    pheno = pd.read_csv(args.phenotype, sep="\t", index_col=0).iloc[:, 0].astype(int)
    lin = pd.read_csv(args.lineages, sep="\t", index_col=0).iloc[:, 0].astype(str)
    tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)
    tips = {t.label for t in tree.taxon_namespace}
    dist = pd.read_csv(args.distances, sep="\t", index_col=0)
    dist.index = dist.index.astype(str); dist.columns = dist.columns.astype(str)

    common = [s for s in geno.index.intersection(pheno.index).intersection(lin.index)
              if s in tips and s in dist.index]
    geno, pheno, lin = geno.loc[common], pheno.loc[common], lin.loc[common]
    dist = dist.loc[common, common]
    geno = geno.apply(pd.to_numeric, errors="coerce")

    top = lin.value_counts().head(args.top_lineages).index.tolist()
    X0 = pd.get_dummies(lin.where(lin.isin(top), "POOLED"), prefix="L",
                        drop_first=True).values.astype(float)
    baseline = float(pheno.mean())

    print(f"samples={len(common)}  candidates={geno.shape[1]}  baseline={baseline:.4f}")
    print(f"permutations={args.n_permutations}  seed={args.seed}  alpha={args.alpha}")

    t0 = time.time()
    obs, obs_detail = run_screen(geno, pheno.values.astype(int), lin.values, tree, dist,
                                 X0, args.alpha, baseline, args.temporal_margin,
                                 args.min_origins)
    print(f"observed convergent: {obs}  {list(obs_detail)}  ({time.time()-t0:.0f}s)")

    rng = np.random.default_rng(args.seed)
    counts = []
    for i in range(args.n_permutations):
        ysh = stratified_shuffle(pheno.values.astype(int), lin.values, rng)
        k, _ = run_screen(geno, ysh, lin.values, tree, dist, X0, args.alpha,
                          float(ysh.mean()), args.temporal_margin, args.min_origins)
        counts.append(k)
        print(f"  permutation {i+1}/{args.n_permutations}: {k} convergent", flush=True)

    counts = np.array(counts)
    any_rate = float((counts >= 1).mean())
    emp_p = (int((counts >= obs).sum()) + 1) / (args.n_permutations + 1)
    lines = [
        "FAMILY-WISE ERROR ACROSS THE SIX-STAGE SCREEN",
        "=" * 64,
        f"Permutations              : {args.n_permutations}  (seed {args.seed})",
        f"Shuffle                   : within lineage",
        f"Stages evaluated          : 1 (independent LRT) -> 3 -> 4 -> 5 (conditional)",
        "",
        f"Observed convergent       : {obs}",
        f"Null mean convergent      : {counts.mean():.3f}",
        f"Null distribution         : {[int(c) for c in counts]}",
        "",
        f"FAMILY-WISE ERROR RATE    : {any_rate:.3f}",
        f"   (fraction of null runs in which ANY candidate reached convergence)",
        f"Empirical p for observed  : {emp_p:.4f}",
    ]
    text = "\n".join(lines) + "\n"
    print("\n" + text)
    with open(f"{args.outdir}/fwer_permutation.txt", "w") as fh:
        fh.write(text)
    pd.DataFrame({"permutation": range(len(counts)), "n_convergent": counts}).to_csv(
        f"{args.outdir}/fwer_permutation.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())

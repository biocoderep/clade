#!/usr/bin/env python3
"""Sweep the number of lineage covariates and measure Stage 1's null calibration
at each, using the independent likelihood-ratio test now in stage1_association.

Earlier probes (Wald statistic, plain logistic regression) showed false
positives under a lineage-preserving null falling sharply as more lineages get
their own covariate rather than being pooled: ~9 of 37 at 20 dummies, ~1.3 at
60. This repeats that sweep with the LRT-based Stage 1 that replaced the Wald
path, to find how many covariates the case-study cohort actually needs.

Usage:
    stage1_calibration_sweep.py --genotypes g.csv --phenotype p.tsv \
        --lineages l.tsv --dimensions 20 40 60 80 --n-permutations 15
"""
from __future__ import annotations

import argparse
import sys
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


def lrt_significant_count(geno, y, X0, alpha, min_carriers, max_freq):
    import statsmodels.api as sm
    from scipy.stats import chi2
    from statsmodels.stats.multitest import fdrcorrection

    pvals, coefs, names = [], [], []
    for c in geno.columns:
        g = geno[c].to_numpy(dtype=float)
        n_carriers = int(np.nansum(g == 1))
        if n_carriers < min_carriers or n_carriers / len(g) > max_freq:
            continue
        n11 = int(((g == 1) & (y == 1)).sum()); n10 = int(((g == 1) & (y == 0)).sum())
        n01 = int(((g == 0) & (y == 1)).sum()); n00 = int(((g == 0) & (y == 0)).sum())
        if min(n11, n10, n01, n00) == 0:
            continue  # separated: no p-value, cannot count as a false positive
        try:
            full = sm.Logit(y, sm.add_constant(np.hstack([X0, g.reshape(-1, 1)]),
                                               has_constant="add")).fit(disp=0, maxiter=200)
            red = sm.Logit(y, sm.add_constant(X0, has_constant="add")).fit(disp=0, maxiter=200)
            stat = 2.0 * (float(full.llf) - float(red.llf))
            if not np.isfinite(stat) or stat < 0:
                continue
            pvals.append(float(chi2.sf(stat, 1)))
            coefs.append(float(full.params[-1]))
            names.append(c)
        except Exception:  # noqa: BLE001, S112 - a non-converging refit for one
            # candidate under one permutation is an expected outcome at this volume,
            # not something to log individually.
            continue
    if not pvals:
        return 0, 0
    _, q = fdrcorrection(np.array(pvals), alpha=alpha)
    return int((q < alpha).sum()), len(pvals)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--lineages", required=True)
    ap.add_argument("--dimensions", type=int, nargs="+", default=[20, 40, 60, 80, 100])
    ap.add_argument("--n-permutations", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-carriers", type=int, default=5)
    ap.add_argument("--max-freq", type=float, default=0.85)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    geno = pd.read_csv(args.genotypes, index_col=0).apply(pd.to_numeric, errors="coerce")
    pheno = pd.read_csv(args.phenotype, sep="\t", index_col=0).iloc[:, 0].astype(int)
    lin = pd.read_csv(args.lineages, sep="\t", index_col=0).iloc[:, 0].astype(str)
    common = geno.index.intersection(pheno.index).intersection(lin.index)
    geno, pheno, lin = geno.loc[common], pheno.loc[common], lin.loc[common]
    y0 = pheno.to_numpy(dtype=float)
    strata = lin.to_numpy()
    n_lineages = lin.nunique()

    print(f"samples={len(common)}  candidates={geno.shape[1]}  distinct lineages={n_lineages}")
    print(f"dimensions to test: {args.dimensions}  permutations each: {args.n_permutations}\n")

    results = []
    for k in args.dimensions:
        top = lin.value_counts().head(k).index.tolist()
        X0 = pd.get_dummies(lin.where(lin.isin(top), "POOLED"), prefix="L",
                            drop_first=True).to_numpy(dtype=float)
        pooled_n = int((~lin.isin(top)).sum())

        rng = np.random.default_rng(args.seed)
        counts, n_tested = [], None
        for i in range(args.n_permutations):
            ysh = stratified_shuffle(y0, strata, rng)
            k_sig, n_t = lrt_significant_count(geno, ysh, X0, args.alpha,
                                               args.min_carriers, args.max_freq)
            counts.append(k_sig)
            n_tested = n_t
        counts = np.array(counts)
        # Under BH with a complete null, expected discoveries are ~0, not
        # alpha*n_tested -- the naive nominal figure would be misleading, so it is
        # not computed. Report the raw mean and the fraction of permutations with
        # ANY hit instead.
        any_hit = float((counts >= 1).mean())
        print(f"k={k:>3d} covariates ({pooled_n} pooled)  tested={n_tested}  "
              f"mean_significant={counts.mean():.3f}  "
              f"fraction_with_any_hit={any_hit:.3f}  counts={list(counts)}")
        results.append({"k_covariates": k, "n_pooled": pooled_n, "n_tested": n_tested,
                        "mean_significant": counts.mean(), "fraction_with_any_hit": any_hit,
                        "counts": list(counts)})

    df = pd.DataFrame(results)
    df.to_csv(f"{args.outdir}/stage1_calibration_sweep.tsv", sep="\t", index=False)

    best = df.loc[df.mean_significant.idxmin()]
    print(f"\nBest calibration: k={int(best.k_covariates)} covariates, "
          f"mean_significant={best.mean_significant:.3f}")
    print(f"Wrote {args.outdir}/stage1_calibration_sweep.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())

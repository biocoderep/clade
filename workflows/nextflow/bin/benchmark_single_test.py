#!/usr/bin/env python3
"""Benchmark CLADE against a conventional single-test association pipeline.

The question a reader of a framework paper asks first: what does the extra
machinery buy? This answers it by running the standard approach on the same
cohort and comparing dispositions candidate by candidate.

The comparator is pyseer's fixed-effects, distance-based population-structure
correction: classical multidimensional scaling of the patristic distance
matrix, the leading components used as covariates in logistic regression, then
Benjamini-Hochberg across candidates. That is what `pyseer --distances
--max-dimensions k` computes.

It is a reimplementation, not pyseer itself -- pyseer is distributed through
bioconda and was not installable in this environment -- and it is labelled that
way wherever the result appears. The model is specified in pyseer's
documentation and is simple enough to reproduce exactly; what cannot be claimed
is bit-identical agreement with that tool.

Output: one row per candidate with the single-test verdict, CLADE's verdict,
and which stage CLADE used to reject it.
"""
from __future__ import annotations

import argparse
import sys
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter("ignore")


def classical_mds(dist: pd.DataFrame, k: int) -> np.ndarray:
    """Classical MDS (Torgerson): double-centre the squared distances and take
    the leading eigenvectors. This is the projection pyseer uses to summarise
    population structure from a distance matrix."""
    d2 = dist.to_numpy(dtype=float) ** 2
    n = d2.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1][:k]
    vals, vecs = vals[order], vecs[:, order]
    vals = np.clip(vals, 0, None)
    return vecs * np.sqrt(vals)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--distances", required=True)
    ap.add_argument("--clade-results", help="TSV of CLADE dispositions to compare against")
    ap.add_argument("--max-dimensions", type=int, default=10)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    import statsmodels.api as sm
    from scipy.stats import chi2
    from statsmodels.stats.multitest import fdrcorrection

    geno = pd.read_csv(args.genotypes, index_col=0)
    pheno = pd.read_csv(args.phenotype, sep="\t", index_col=0).iloc[:, 0].astype(int)
    dist = pd.read_csv(args.distances, sep="\t", index_col=0)
    dist.index = dist.index.astype(str); dist.columns = dist.columns.astype(str)
    common = [s for s in geno.index.intersection(pheno.index) if s in dist.index]
    geno, pheno, dist = geno.loc[common], pheno.loc[common], dist.loc[common, common]
    geno = geno.apply(pd.to_numeric, errors="coerce")

    print(f"samples={len(common)}  candidates={geno.shape[1]}")
    print(f"computing {args.max_dimensions}-dimensional MDS of the distance matrix...")
    mds = classical_mds(dist, args.max_dimensions)
    y = pheno.to_numpy(dtype=float)

    rows = []
    for c in geno.columns:
        g = geno[c].fillna(0).to_numpy(dtype=float)
        n11 = int(((g == 1) & (y == 1)).sum()); n10 = int(((g == 1) & (y == 0)).sum())
        n01 = int(((g == 0) & (y == 1)).sum()); n00 = int(((g == 0) & (y == 0)).sum())
        sep = min(n11, n10, n01, n00) == 0
        try:
            full = sm.Logit(y, sm.add_constant(np.hstack([mds, g.reshape(-1, 1)]),
                                               has_constant="add")).fit(disp=0, maxiter=200)
            red = sm.Logit(y, sm.add_constant(mds, has_constant="add")).fit(disp=0, maxiter=200)
            stat = 2.0 * (float(full.llf) - float(red.llf))
            p = float(chi2.sf(stat, 1)) if np.isfinite(stat) and stat >= 0 else np.nan
            b = float(full.params[-1])
        except Exception:  # noqa: BLE001
            p, b = np.nan, np.nan
        rows.append({"Candidate": c, "singletest_coef": b, "singletest_p": p,
                     "separation": sep})

    df = pd.DataFrame(rows)
    ok = df.singletest_p.notna()
    q = np.full(len(df), np.nan)
    if ok.any():
        _, qq = fdrcorrection(df.loc[ok, "singletest_p"].to_numpy(), alpha=args.alpha)
        q[ok.to_numpy()] = qq
    df["singletest_q"] = q
    df["singletest_hit"] = (df.singletest_q < args.alpha) & (df.singletest_coef > 0)

    print(f"\nSINGLE-TEST PIPELINE (pyseer-equivalent, {args.max_dimensions} MDS dimensions)")
    print(f"  significant at q<{args.alpha}          : {int((df.singletest_q < args.alpha).sum())}")
    print(f"  ...and positively signed (a 'hit'): {int(df.singletest_hit.sum())}")
    print(f"  of those hits, separated          : {int((df.singletest_hit & df.separation).sum())}")

    if args.clade_results:
        cl = pd.read_csv(args.clade_results, sep="\t")
        col = "Disposition" if "Disposition" in cl.columns else cl.columns[1]
        df = df.merge(cl[["Candidate", col]].rename(columns={col: "clade_disposition"}),
                      on="Candidate", how="left")
        hits = df[df.singletest_hit]
        print("\n  What CLADE does with the single-test hits:")
        if len(hits):
            print(hits.clade_disposition.value_counts().to_string())
            print(f"\n  Hits CLADE does NOT call convergent: "
                  f"{int((hits.clade_disposition != 'convergent').sum())} of {len(hits)}")
    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

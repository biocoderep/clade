#!/usr/bin/env python3
"""Effective independent test count from candidate correlation structure.

Benjamini-Hochberg assumes the tests it corrects are (near) independent. In a
clonal population they frequently are not: paralogous loci, genes in linkage,
and outright duplicate columns all produce candidate vectors that carry the
same information. Counting each as a separate test inflates the denominator and
distorts every q-value in the family.

This reports the redundancy explicitly: the pairwise correlation matrix, the
clusters that collapse at a threshold, and the resulting effective number of
independent genetic units. It does NOT silently drop candidates -- which ones
to keep is a judgement about biology, not a computation, so the decision stays
with the analyst and the evidence for it is written down.

Outputs
    candidate_correlation.tsv     full pairwise |r| matrix
    effective_tests.tsv           one row per candidate: cluster id, representative
    correlation_summary.txt       human-readable summary
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd


def greedy_cluster(corr: pd.DataFrame, threshold: float) -> dict[str, str]:
    """Collapse candidates into clusters at |r| >= threshold.

    Greedy and order-dependent by construction, so candidates are visited in
    sorted name order to make the result deterministic: the same input must
    give the same clusters on every run and every machine.
    """
    assigned: dict[str, str] = {}
    for name in sorted(corr.columns):
        if name in assigned:
            continue
        assigned[name] = name  # this candidate represents its own cluster
        for other in sorted(corr.columns):
            if other == name or other in assigned:
                continue
            if abs(corr.loc[name, other]) >= threshold:
                assigned[other] = name
    return assigned


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--threshold", type=float, default=0.9,
                    help="absolute correlation at or above which two candidates "
                         "are treated as one effective test (default: 0.9)")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    if not 0 < args.threshold <= 1:
        sys.exit(f"--threshold must be in (0,1]; got {args.threshold}")

    geno = pd.read_csv(args.genotypes, index_col=0)
    geno = geno.select_dtypes(include=[np.number])
    if geno.shape[1] < 2:
        sys.exit("need at least two candidates to assess correlation structure")

    # Constant columns have undefined correlation; report them rather than
    # letting them appear as silent NaN in the matrix.
    constant = [c for c in geno.columns if geno[c].nunique(dropna=True) <= 1]
    variable = [c for c in geno.columns if c not in constant]
    if len(variable) < 2:
        sys.exit("fewer than two non-constant candidates; nothing to cluster")

    corr = geno[variable].corr().fillna(0.0)
    corr.to_csv(f"{args.outdir}/candidate_correlation.tsv", sep="\t",
                float_format="%.6f")

    clusters = greedy_cluster(corr, args.threshold)
    rows = [{"Candidate": c, "Cluster": clusters[c], "Is_representative": clusters[c] == c}
            for c in sorted(variable)]
    rows += [{"Candidate": c, "Cluster": "constant", "Is_representative": False}
             for c in sorted(constant)]
    table = pd.DataFrame(rows)
    table.to_csv(f"{args.outdir}/effective_tests.tsv", sep="\t", index=False)

    n_nominal = len(geno.columns)
    n_effective = len({clusters[c] for c in variable})

    # Pairs worth naming in the write-up: near-identical candidates.
    pairs = []
    cols = list(corr.columns)
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            r = corr.loc[a, b]
            if abs(r) >= args.threshold:
                pairs.append((a, b, r))
    pairs.sort(key=lambda t: -abs(t[2]))

    lines = [
        "Candidate correlation structure",
        "=" * 60,
        f"Nominal candidates          : {n_nominal}",
        f"Constant (no variation)     : {len(constant)}",
        f"Effective independent units : {n_effective}  (|r| >= {args.threshold})",
        "",
        "Interpretation: BH-FDR across the nominal count treats correlated",
        "candidates as independent tests, inflating the denominator. The",
        "effective count is the honest number of independent genetic units.",
        "",
        f"Correlated pairs at |r| >= {args.threshold}: {len(pairs)}",
    ]
    for a, b, r in pairs[:40]:
        flag = "  <-- IDENTICAL" if abs(r) > 0.9999 else ""
        lines.append(f"    {a}  vs  {b}   r = {r:+.4f}{flag}")
    if len(pairs) > 40:
        lines.append(f"    ... and {len(pairs) - 40} more (see candidate_correlation.tsv)")
    if constant:
        lines += ["", "Constant candidates (excluded from clustering):"]
        lines += [f"    {c}" for c in constant]

    text = "\n".join(lines) + "\n"
    with open(f"{args.outdir}/correlation_summary.txt", "w") as fh:
        fh.write(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())

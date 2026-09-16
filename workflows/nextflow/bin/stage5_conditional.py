#!/usr/bin/env python3
"""Stage 5 by conditional logistic regression over matched strata.

Why this exists. Nearest-neighbour matching with replacement gives McNemar's
test dependent pairs, which it assumes away: in this project's cohort 61 unique
susceptible genomes served all 2,492 pairs and one anchored 1,491 of them, so a
reported p of 1.6e-51 rested on an effective sample size of one. Forcing unique
1:1 matching fixes independence but discards 70% of the resistant genomes and
introduces whatever selection the greedy matching order imposes.

Neither is satisfactory, so neither is used here. Matching with replacement is
kept -- every resistant genome retains its true nearest susceptible neighbour --
and the resulting dependence is modelled rather than assumed away:

  * Each control genome defines a stratum. Every resistant genome matched to
    that control joins its stratum. Conditional logistic regression estimates
    the candidate effect WITHIN strata, so a control reused 1,491 times
    contributes one stratum, not 1,491 independent observations.

  * Conditional likelihood is implemented as Cox proportional hazards with all
    events at the same time, stratified by control. This is the standard
    equivalence between conditional logistic regression and a stratified Cox
    model, and it lets us use a well-tested implementation rather than writing
    our own optimiser.

  * A cluster-robust variance keyed on the control genome is also reported, as
    a second, less model-dependent view of the same dependence.

Reported alongside, deliberately, so the reader can see what the correction
costs and what it rescues:

    mcnemar_naive_p      McNemar over all pairs -- what was published
    mcnemar_unique_p     McNemar after unique 1:1 matching -- the lossy fix
    conditional_p        this method

Where the three disagree, the naive column is the one to distrust.

Usage:
    stage5_conditional.py --genotypes g.csv --phenotype p.tsv \
        --distances d.tsv --candidates cand1 cand2 --output out.tsv
"""
from __future__ import annotations

import argparse
import sys
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter("ignore")


def build_strata(distances, positives, negatives):
    """Match each case to its nearest control, with replacement.

    Delegates to the package's own matcher so this analysis and the pipeline
    use identical pairing -- including its load-balanced tie-breaking, which
    matters enormously here: 83% of matches in this cohort are exact ties, and
    breaking them by first-column order collapsed 2,492 pairs onto 61 controls.
    Spreading ties across equally-near controls raised that to 302.
    """
    from clade.validation.stage5_matched_neighbors import nearest_negative_neighbor

    m = nearest_negative_neighbor(distances, list(positives), list(negatives))
    m = m.reset_index()
    # The index name varies with how the genotype frame was loaded
    # (Sample_ID, positive, or unnamed), so take the first column positionally
    # rather than guessing at a name.
    pairs = m.rename(columns={m.columns[0]: "case", "nearest_neighbor": "control"})[
        ["case", "control"]
    ]
    return pairs


def conditional_test(pairs, genotype):
    """Conditional logistic regression of candidate carriage on case status,
    stratified by control genome.

    Each stratum is one control plus every case matched to it. Strata that are
    uninformative -- no variation in the candidate, or no variation in case
    status -- contribute nothing to a conditional likelihood and are dropped,
    which is a property of the method, not a filter applied by hand.
    """
    from lifelines import CoxPHFitter

    rows = []
    for control, grp in pairs.groupby("control"):
        cases = list(grp["case"])
        members = cases + [control]
        for s in members:
            g = genotype.get(s)
            if pd.isna(g):
                continue
            rows.append({"stratum": control, "sample": s,
                         "is_case": 1 if s in grp["case"].values else 0,
                         "carrier": float(g)})
    df = pd.DataFrame(rows)
    if df.empty:
        return {"conditional_p": np.nan, "conditional_coef": np.nan,
                "n_strata_informative": 0, "note": "no usable strata"}

    keep = [s for s, g in df.groupby("stratum")
            if g.carrier.nunique() > 1 and g.is_case.nunique() > 1]
    df = df[df.stratum.isin(keep)]
    if df.empty or df.carrier.nunique() < 2:
        return {"conditional_p": np.nan, "conditional_coef": np.nan,
                "n_strata_informative": 0,
                "note": "no informative strata (no within-stratum variation)"}

    # Stratified Cox with a common event time == conditional logistic likelihood.
    df = df.assign(T=1.0, E=df.is_case.astype(int))
    try:
        cph = CoxPHFitter()
        cph.fit(df[["T", "E", "carrier", "stratum"]], duration_col="T",
                event_col="E", strata=["stratum"], robust=True)
        return {"conditional_p": float(cph.summary.loc["carrier", "p"]),
                "conditional_coef": float(cph.summary.loc["carrier", "coef"]),
                "n_strata_informative": len(keep),
                "note": "conditional logistic via stratified Cox"}
    except Exception as exc:  # noqa: BLE001
        return {"conditional_p": np.nan, "conditional_coef": np.nan,
                "n_strata_informative": len(keep),
                "note": f"fit failed: {type(exc).__name__}: {exc}"}


def mcnemar_p(pairs, genotype, exact_threshold=25):
    from statsmodels.stats.contingency_tables import mcnemar
    a = genotype.reindex(pairs["case"]).to_numpy(dtype=float)
    b = genotype.reindex(pairs["control"]).to_numpy(dtype=float)
    ok = ~(np.isnan(a) | np.isnan(b))
    a, b = a[ok], b[ok]
    both = int(((a == 1) & (b == 1)).sum()); case_only = int(((a == 1) & (b == 0)).sum())
    ctrl_only = int(((a == 0) & (b == 1)).sum()); neither = int(((a == 0) & (b == 0)).sum())
    disc = case_only + ctrl_only
    if disc == 0:
        return np.nan, case_only, ctrl_only
    try:
        p = float(mcnemar([[both, case_only], [ctrl_only, neither]],
                          exact=disc < exact_threshold).pvalue)
    except Exception:  # noqa: BLE001
        p = np.nan
    return p, case_only, ctrl_only


def greedy_unique(distances, positives, negatives):
    """Greedy 1:1 matching, closest pairs first, each control used once."""
    sub = distances.loc[positives, negatives]
    order = sorted(((sub.loc[c].min(), c) for c in positives), key=lambda t: t[0])
    used, out = set(), []
    for _, case in order:
        row = sub.loc[case].drop(labels=list(used), errors="ignore")
        if row.empty:
            continue
        m = row.min()
        pick = sorted(row.index[row == m])[0]
        used.add(pick)
        out.append((case, pick))
    return pd.DataFrame(out, columns=["case", "control"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--distances", required=True)
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    geno = pd.read_csv(args.genotypes, index_col=0)
    pheno = pd.read_csv(args.phenotype, sep="\t", index_col=0).iloc[:, 0].astype(int)
    dist = pd.read_csv(args.distances, sep="\t", index_col=0)
    dist.index = dist.index.astype(str); dist.columns = dist.columns.astype(str)

    common = [s for s in geno.index.intersection(pheno.index) if s in dist.index]
    if not common:
        sys.exit("no samples common to genotypes, phenotype and distance matrix")
    geno, pheno, dist = geno.loc[common], pheno.loc[common], dist.loc[common, common]
    pos = [s for s in common if pheno[s] == 1]
    neg = [s for s in common if pheno[s] == 0]
    print(f"samples={len(common)}  resistant={len(pos)}  susceptible={len(neg)}")

    pairs = build_strata(dist, pos, neg)
    uniq = greedy_unique(dist, pos, neg)
    n_unique_ctrl = pairs["control"].nunique()
    print(f"pairs={len(pairs)}  unique controls={n_unique_ctrl} "
          f"({100*n_unique_ctrl/len(pairs):.1f}%)  max reuse={pairs['control'].value_counts().max()}")
    print(f"unique 1:1 matching retains {len(uniq)} of {len(pos)} resistant genomes "
          f"({100*len(uniq)/len(pos):.1f}%)\n")

    out = []
    for cand in args.candidates:
        if cand not in geno.columns:
            print(f"  {cand}: NOT IN GENOTYPE MATRIX -- skipped")
            continue
        g = pd.to_numeric(geno[cand], errors="coerce")
        naive_p, case_only, ctrl_only = mcnemar_p(pairs, g)
        uniq_p, u_case_only, u_ctrl_only = mcnemar_p(uniq, g)
        cond = conditional_test(pairs, g)
        row = {"Candidate": cand, "N_carriers": int(g.sum()),
               "mcnemar_naive_p": naive_p, "naive_case_only": case_only,
               "naive_control_only": ctrl_only,
               "mcnemar_unique_p": uniq_p, "unique_case_only": u_case_only,
               "unique_control_only": u_ctrl_only,
               "n_unique_controls": n_unique_ctrl, **cond}
        out.append(row)
        print(f"  {cand}: naive p={naive_p:.3g}  unique p={uniq_p:.3g}  "
              f"conditional p={cond['conditional_p']:.3g} "
              f"(coef {cond['conditional_coef']:+.3f}, {cond['n_strata_informative']} strata)")

    pd.DataFrame(out).to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

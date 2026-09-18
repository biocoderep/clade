#!/usr/bin/env python3
"""Record both available Firth p-value estimators, side by side, on real data.

Firth bias-reduced logistic regression admits two p-values, and they are not
interchangeable:

  * the penalised likelihood-ratio p-value the library computes, which refits a
    null model per coefficient; and
  * a Wald p-value built from the fitted coefficient and its standard error.

Standard guidance for Firth regression favours the penalised LRT, because Wald
inference degrades under separation. But an implementation can fail to
converge, and a Wald computation can lose all precision, and neither announces
itself in the output -- both just return a number.

So this records both for every candidate and reports where they disagree,
rather than picking one and hoping. Two specific failure signatures are
detected and flagged, because both were observed on real data in this project:

  * a penalised-LRT value repeated identically across unrelated candidates,
    which is a convergence floor rather than a statistic; and
  * a Wald p-value underflowing to exactly zero.

The Wald computation uses the normal survival function. Computing it as
1 - Phi(|z|) suffers catastrophic cancellation and floors at machine epsilon
(2.22e-16), a value indistinguishable from a genuine p-value in output.

Outputs
    estimator_comparison.tsv    per-candidate: both p-values, both q-values
    estimator_summary.txt       counts, disagreements, failure signatures
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import norm


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--lineages", required=True)
    ap.add_argument("--min-carriers", type=int, default=5)
    ap.add_argument("--max-freq", type=float, default=0.85)
    ap.add_argument("--top-lineages", type=int, default=20,
                    help="lineages given individual dummy covariates; the rest "
                         "are pooled (default: 20)")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    try:
        from firthlogist import FirthLogisticRegression
    except ImportError:
        sys.exit("estimator_diagnostic.py requires firthlogist "
                 "(Python <3.11, scikit-learn <1.6). Install the 'firth' extra.")
    from statsmodels.stats.multitest import fdrcorrection

    geno = pd.read_csv(args.genotypes, index_col=0)
    pheno = pd.read_csv(args.phenotype, sep="\t", index_col=0).iloc[:, 0].astype(int)
    lineages = pd.read_csv(args.lineages, sep="\t", index_col=0).iloc[:, 0].astype(str)

    common = geno.index.intersection(pheno.index).intersection(lineages.index)
    if len(common) == 0:
        sys.exit("no samples in common across genotypes, phenotype and lineages")
    geno, pheno, lineages = geno.loc[common], pheno.loc[common], lineages.loc[common]

    top = lineages.value_counts().head(args.top_lineages).index.tolist()
    covariates = pd.get_dummies(lineages, prefix="L").reindex(
        columns=[f"L_{s}" for s in top], fill_value=0).values.astype(float)
    y = pheno.values.astype(int)

    rows = []
    for name in geno.columns:
        g = pd.to_numeric(geno[name], errors="coerce").values
        n_carriers = int(np.nansum(g == 1))
        prevalence = n_carriers / len(g) if len(g) else 0.0

        if n_carriers < args.min_carriers:
            rows.append({"Candidate": name, "N_carriers": n_carriers,
                         "Status": f"skipped_too_rare(n={n_carriers})"})
            continue
        if prevalence > args.max_freq:
            rows.append({"Candidate": name, "N_carriers": n_carriers,
                         "Status": f"skipped_near_fixed(freq={prevalence:.3f})"})
            continue

        X = np.hstack([covariates, np.nan_to_num(g).reshape(-1, 1)])
        try:
            # skip_pvals=False: we explicitly WANT the library's own estimator,
            # even though it costs a refit per coefficient. Comparing them is
            # the entire purpose of this step.
            model = FirthLogisticRegression(max_iter=200, skip_pvals=False)
            model.fit(X, y)
            coef = float(model.coef_[-1])
            bse = float(model.bse_[-1])
            try:
                lib_p = float(model.pvals_[-1])
            except Exception:  # noqa: BLE001 - a missing/invalid pvals_ just means 'no library estimate'
                lib_p = np.nan
            if bse == 0 or np.isnan(bse):
                wald_z, wald_p = np.nan, np.nan
            else:
                wald_z = coef / bse
                wald_p = 2 * norm.sf(abs(wald_z))   # survival function, not 1-cdf
            rows.append({"Candidate": name, "N_carriers": n_carriers,
                         "Coef": coef, "BSE": bse, "Wald_z": wald_z,
                         "Wald_p": wald_p, "Library_p": lib_p, "Status": "tested"})
        except Exception as exc:  # noqa: BLE001 - a per-candidate failure is data, not a crash
            rows.append({"Candidate": name, "N_carriers": n_carriers,
                         "Status": f"failed: {type(exc).__name__}: {exc}"})

    df = pd.DataFrame(rows)
    ok = df["Status"] == "tested"

    if ok.sum() == 0:
        sys.exit("no candidate could be tested; check --min-carriers / --max-freq")

    # FDR within each estimator separately, so the two are comparable like-for-like.
    _, q_wald = fdrcorrection(df.loc[ok, "Wald_p"].fillna(1.0).values, alpha=args.alpha)
    df.loc[ok, "q_wald"] = q_wald
    df.loc[ok, "Significant_wald"] = q_wald < args.alpha
    if df.loc[ok, "Library_p"].notna().all():
        _, q_lib = fdrcorrection(df.loc[ok, "Library_p"].values, alpha=args.alpha)
        df.loc[ok, "q_library"] = q_lib
        df.loc[ok, "Significant_library"] = q_lib < args.alpha
        df.loc[ok, "Disagreement"] = df.loc[ok, "Significant_wald"] != df.loc[ok, "Significant_library"]
        df.loc[ok, "log10_ratio"] = (np.log10(df.loc[ok, "Library_p"].clip(lower=1e-320))
                                     - np.log10(df.loc[ok, "Wald_p"].clip(lower=1e-320)))

    df.to_csv(f"{args.outdir}/estimator_comparison.tsv", sep="\t", index=False,
              float_format="%.6g")

    # --- failure signatures -------------------------------------------------
    sub = df[ok]
    underflow = sub[sub["Wald_p"] == 0.0]["Candidate"].tolist()
    floors = []
    if "Library_p" in sub:
        counts = Counter(np.round(sub["Library_p"].dropna(), 12))
        for value, n in counts.items():
            if n >= 3:   # three unrelated candidates sharing a p-value exactly
                members = sub[np.isclose(sub["Library_p"], value)]["Candidate"].tolist()
                floors.append((value, members))

    lines = [
        "Firth p-value estimator comparison",
        "=" * 64,
        f"Candidates tested                     : {int(ok.sum())} of {len(df)}",
    ]
    if "Significant_library" in df:
        lines += [
            f"Significant at q<{args.alpha} (Wald)          : {int(sub['Significant_wald'].sum())}",
            f"Significant at q<{args.alpha} (penalised LRT) : {int(sub['Significant_library'].sum())}",
            f"Candidates where the two disagree     : {int(sub['Disagreement'].sum())}",
            f"Median log10(library_p / wald_p)      : {sub['log10_ratio'].median():+.2f}",
            "    negative => library systematically SMALLER (anticonservative)",
            "    positive => library systematically LARGER  (conservative)",
        ]
    lines += ["", "Numerical failure signatures:"]
    if underflow:
        lines.append(f"  Wald p underflowed to exactly 0 for {len(underflow)}: " + ", ".join(underflow[:6]))
    if floors:
        for value, members in floors:
            lines.append(f"  penalised-LRT returned {value:.6g} identically for "
                         f"{len(members)} unrelated candidates -- a convergence "
                         f"floor, not a statistic: " + ", ".join(members[:6]))
    if not underflow and not floors:
        lines.append("  none detected")
    lines += ["",
              "Neither estimator should be trusted where a failure signature",
              "appears. Candidates that also pass Stage 4 (tree topology) and",
              "Stage 5 (pairwise distance) do not depend on this p-value."]

    text = "\n".join(lines) + "\n"
    with open(f"{args.outdir}/estimator_summary.txt", "w") as fh:
        fh.write(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())

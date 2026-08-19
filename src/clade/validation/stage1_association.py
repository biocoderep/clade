"""
Stage 1 — Structure-corrected statistical association.

CLADE's preferred Stage 1 method is kinship-matrix-corrected regression
(`pyseer`, https://github.com/mgalardini/pyseer) — an external, compiled
tool CLADE does not reimplement. Where kinship correction fails to converge
(a real and, in our experience, not-infrequent occurrence in strongly
clonal datasets — see docs/reproducibility/original_discovery_pipeline.md
and the CLADE manuscript Section 3.2), this module provides the fallback
path that CLADE actually ships as Python: Firth bias-reduced logistic
regression with lineage-indicator covariates, and a Fisher's-exact +
BH-FDR + Lasso screen for the full-cohort case.

These are direct refactors of the original firth_final.py and
fisher_fdr_lasso.py scripts into testable functions — the statistical
logic is unchanged, only the file I/O and print statements were removed.
"""
from __future__ import annotations

import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from statsmodels.stats.multitest import multipletests


def firth_association(df: pd.DataFrame, candidate: str, phenotype_col: str, st_dummies: pd.DataFrame) -> dict:
    """Firth bias-reduced logistic regression for one candidate, with ST
    lineage covariates. The candidate is always the last column, so its
    coefficient/p-value are at index `st_dummies.shape[1]`.

    Requires the `firthlogist` package (Python <3.11, scikit-learn <1.6 as
    of this writing — see docs/reproducibility/open_items.md for the exact
    dependency conflict this project hit and how it was resolved).
    """
    from firthlogist import FirthLogisticRegression

    n_st_cols = st_dummies.shape[1]
    X = pd.concat([st_dummies, df[[candidate]].astype(float)], axis=1).values
    y = df[phenotype_col].astype(int).values
    cand_idx = n_st_cols

    n_carriers = int(df[candidate].sum())
    try:
        fl = FirthLogisticRegression(test_vars=cand_idx, skip_ci=True)
        fl.fit(X, y)
        return {
            "Candidate": candidate,
            "N_carriers": n_carriers,
            "Coef": float(fl.coef_[cand_idx]),
            "Firth_p": float(fl.pvals_[cand_idx]),
            "Status": "OK",
        }
    except Exception as e:  # noqa: BLE001 - real convergence failures are the point
        return {
            "Candidate": candidate,
            "N_carriers": n_carriers,
            "Coef": None,
            "Firth_p": None,
            "Status": f"FAILED: {e}",
        }


def fisher_fdr_screen(
    df: pd.DataFrame,
    candidate_cols: list[str],
    phenotype_col: str,
    min_carriers: int = 10,
    max_freq: float = 0.80,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Fisher's exact test per candidate + Benjamini-Hochberg FDR correction,
    after a prevalence pre-filter (drops near-fixed/degenerate candidates).

    This is a pre-filter only — a candidate passing this screen is not yet
    evidence on its own; it identifies which candidates are worth carrying
    into Stage 1's kinship/Firth regression or Stage 1's Lasso follow-up.
    """
    n = len(df)
    valid = [f for f in candidate_cols if min_carriers <= df[f].sum() and df[f].sum() / n <= max_freq]

    p_values = []
    for f in valid:
        ct = pd.crosstab(df[f], df[phenotype_col]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        _, p = stats.fisher_exact(ct)
        p_values.append(p)

    if not valid:
        return pd.DataFrame(columns=["Candidate", "Fisher_p", "FDR_q", "FDR_significant"])

    reject, qvals, _, _ = multipletests(p_values, alpha=alpha, method="fdr_bh")
    result = pd.DataFrame({"Candidate": valid, "Fisher_p": p_values, "FDR_q": qvals, "FDR_significant": reject})
    return result.sort_values("Fisher_p").reset_index(drop=True)


def lasso_st_adjusted(df: pd.DataFrame, sig_features: list[str], st_dummies: pd.DataFrame, phenotype_col: str, C: float = 0.1) -> pd.Series:
    """L1-penalized logistic regression on FDR-significant candidates plus
    ST covariates. A coefficient of exactly 0 means the candidate's signal
    was fully absorbed by lineage structure (consistent with a clonal
    artifact); non-zero means it survives structure control.
    """
    X = pd.concat([df[sig_features], st_dummies], axis=1).astype(float)
    y = df[phenotype_col].astype(int)

    X_std = (X - X.mean()) / X.std(ddof=0)
    X_std = X_std.fillna(0)

    lasso = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=2000)
    lasso.fit(X_std, y)

    coefs = pd.Series(lasso.coef_[0], index=X.columns)
    return coefs[sig_features].sort_values(key=abs, ascending=False)

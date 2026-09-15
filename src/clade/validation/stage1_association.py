"""
Stage 1 — Structure-corrected statistical association.

CLADE's preferred Stage 1 method is kinship-matrix-corrected regression
(`pyseer`, https://github.com/mgalardini/pyseer) — an external tool CLADE does
not reimplement. Where kinship correction fails to converge (a real and, in this
project's experience, not-infrequent occurrence in strongly clonal datasets),
this module provides the fallback CLADE actually ships as Python: Firth
bias-reduced logistic regression with lineage-indicator covariates, plus a
Fisher's-exact + BH-FDR + Lasso screen for the full-cohort case.

The statistical logic is unchanged from the original firth_final.py and
fisher_fdr_lasso.py. What has been added:

* `fisher_fdr_screen` now returns **every** candidate with an explicit
  `Screen_status`, instead of silently dropping the ones that fail the
  prevalence pre-filter. Downstream, "excluded by the prevalence filter" and
  "tested and not significant" and "never tested" are three different things,
  and conflating them is how a screened-out candidate ended up reported as
  "not individually re-tested" — CLADE's phrase for a coverage limitation
  rather than a result.
* `firth_association` checks that the covariate matrix and the data frame are
  aligned on the same index before `pd.concat`, which would otherwise align on
  index and silently introduce NaN rows.
* `firth_association` no longer reads `firthlogist`'s `pvals_` attribute (the
  library's built-in profile-likelihood-ratio p-value). A documented forensic
  investigation (`CLADE_TECHNICAL_LOG_PART2.md` in this project's working
  history, sections 12-15) found that output numerically unreliable once the
  covariate design matrix reaches the dimensionality this module actually
  uses -- 20 lineage-dummy columns plus the tested variant -- returning
  p-values roughly 100-300x too small on real data, on pure noise, and under
  permutation alike; a null-permutation calibration run that should show ~5%
  of candidates passing at alpha=0.05 instead showed 91.7% (33/36) with the
  library p-value, against ~28% (10/36) with the fix below -- still elevated,
  attributed there to candidate paralog/linkage correlation and coarse ST
  pooling, not to this defect. The fix computes a manual Wald p-value from
  the library's independently-checked `coef_`/`bse_` outputs instead:
  z = coef/SE, p = 2*(1-Phi(|z|)). This has not been independently
  re-verified against a live `firthlogist` run from this codebase (the
  package requires Python <3.11, unavailable in the environment this fix was
  written in) -- it is applied on the strength of that documented diagnostic,
  not re-derived from scratch. Anyone re-verifying it: check `bse_` is
  actually populated when `test_vars=` is passed (the current construction
  below) -- the log's own verified fix used `skip_pvals=False` without
  `test_vars`, so this is worth confirming they compose the same way.
"""
from __future__ import annotations

import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from statsmodels.stats.multitest import multipletests

from clade.io.validation import CladeInputError, is_missing, validate_binary


def firth_association(
    df: pd.DataFrame,
    candidate: str,
    phenotype_col: str,
    st_dummies: pd.DataFrame,
) -> dict:
    """Firth bias-reduced logistic regression for one candidate, with ST
    lineage covariates.

    The candidate is appended as the last column, so its coefficient and
    p-value sit at index `st_dummies.shape[1]`.

    Requires the `firthlogist` package (Python <3.11, scikit-learn <1.6 as of
    this writing — see docs/reproducibility/open_items.md for the dependency
    conflict this project hit and how it was resolved).
    """
    from firthlogist import FirthLogisticRegression

    if candidate not in df.columns:
        raise CladeInputError(f"candidate '{candidate}' is not a column in the input frame.")
    if not df.index.equals(st_dummies.index):
        raise CladeInputError(
            "the covariate matrix and the data frame are not aligned on the same index; "
            "concatenating them would silently introduce missing rows. Build the ST dummies "
            "from the same frame you pass here."
        )

    cand = validate_binary(df[candidate], f"genotype column '{candidate}'")
    pheno = validate_binary(df[phenotype_col], "phenotype", allow_missing=False)

    complete = cand.notna()
    n_dropped = int((~complete).sum())

    n_st_cols = st_dummies.shape[1]
    X = pd.concat([st_dummies.loc[complete], cand[complete].rename(candidate)], axis=1).to_numpy(
        dtype=float
    )
    y = pheno[complete].to_numpy(dtype=int)
    cand_idx = n_st_cols

    n_carriers = int(cand.sum())
    if n_carriers == 0:
        return {
            "Candidate": candidate,
            "N_carriers": 0,
            "N_dropped_missing": n_dropped,
            "Coef": None,
            "Firth_p": None,
            "Status": "SKIPPED: no carriers",
        }

    try:
        fl = FirthLogisticRegression(test_vars=cand_idx)
        fl.fit(X, y)
        coef = fl.coef_[cand_idx]
        # Deliberately NOT fl.pvals_[cand_idx] -- see the module docstring.
        # firthlogist's profile-likelihood-ratio p-value is numerically
        # unreliable at this covariate dimensionality; a manual Wald p-value
        # from coef_/bse_ is used instead.
        bse = fl.bse_[cand_idx] if hasattr(fl, "bse_") else float("nan")
        if is_missing(coef) or is_missing(bse) or bse == 0:
            pval = float("nan")
        else:
            from scipy.stats import norm
            wald_z = coef / bse
            pval = float(2 * (1 - norm.cdf(abs(wald_z))))
        if is_missing(coef) or is_missing(pval):
            return {
                "Candidate": candidate,
                "N_carriers": n_carriers,
                "N_dropped_missing": n_dropped,
                "Coef": None,
                "Firth_p": None,
                "Status": "FAILED: non-finite estimate",
            }
        return {
            "Candidate": candidate,
            "N_carriers": n_carriers,
            "N_dropped_missing": n_dropped,
            "Coef": float(coef),
            "Firth_p": float(pval),
            "Status": "OK",
        }
    except Exception as e:  # noqa: BLE001 - real convergence failures are the point
        return {
            "Candidate": candidate,
            "N_carriers": n_carriers,
            "N_dropped_missing": n_dropped,
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
    after a prevalence pre-filter that drops near-fixed and too-rare candidates.

    This is a pre-filter only — passing it is not evidence on its own; it
    identifies which candidates are worth carrying into Stage 1's kinship/Firth
    regression.

    Returns one row per *input* candidate. `Screen_status` is one of:

        tested                      Fisher's exact ran; see Fisher_p / FDR_q
        excluded_too_rare           fewer than `min_carriers` carriers
        excluded_near_fixed         carrier frequency above `max_freq`
        excluded_no_variation       candidate is constant in this dataset

    Only `tested` rows take part in the FDR correction, and only they carry a
    q-value. An excluded row has `FDR_significant = None` — unknown, not False.
    """
    pheno = validate_binary(df[phenotype_col], "phenotype", allow_missing=False)
    n = len(df)
    if n == 0:
        raise CladeInputError("no samples available for the Stage 1 screen.")

    statuses: dict[str, str] = {}
    carriers: dict[str, int] = {}
    tested: list[str] = []

    for f in candidate_cols:
        if f not in df.columns:
            raise CladeInputError(f"candidate '{f}' is not a column in the input frame.")
        col = validate_binary(df[f], f"genotype column '{f}'")
        n_car = int(col.sum())
        carriers[f] = n_car
        n_obs = int(col.notna().sum())
        freq = n_car / n_obs if n_obs else 0.0

        if col.nunique(dropna=True) < 2:
            statuses[f] = "excluded_no_variation"
        elif n_car < min_carriers:
            statuses[f] = "excluded_too_rare"
        elif freq > max_freq:
            statuses[f] = "excluded_near_fixed"
        else:
            statuses[f] = "tested"
            tested.append(f)

    p_values = []
    for f in tested:
        col = validate_binary(df[f], f"genotype column '{f}'")
        complete = col.notna()
        ct = (
            pd.crosstab(col[complete], pheno[complete])
            .reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        )
        _, p = stats.fisher_exact(ct.to_numpy())
        p_values.append(float(p))

    qvals: dict[str, float] = {}
    rejects: dict[str, bool] = {}
    if tested:
        reject, q, _, _ = multipletests(p_values, alpha=alpha, method="fdr_bh")
        qvals = dict(zip(tested, q))
        rejects = dict(zip(tested, reject))

    p_by_cand = dict(zip(tested, p_values))
    rows = [
        {
            "Candidate": f,
            "N_carriers": carriers[f],
            "Screen_status": statuses[f],
            "Fisher_p": p_by_cand.get(f),
            "FDR_q": qvals.get(f),
            "FDR_significant": bool(rejects[f]) if f in rejects else None,
        }
        for f in candidate_cols
    ]
    result = pd.DataFrame(rows)
    return result.sort_values(
        ["Screen_status", "Fisher_p"], na_position="last"
    ).reset_index(drop=True)


def lasso_st_adjusted(
    df: pd.DataFrame,
    sig_features: list[str],
    st_dummies: pd.DataFrame,
    phenotype_col: str,
    C: float = 0.1,
    random_state: int = 0,
) -> pd.Series:
    """L1-penalised logistic regression on FDR-significant candidates plus ST
    covariates. A coefficient of exactly 0 means the candidate's signal was
    fully absorbed by lineage structure (consistent with a clonal artifact);
    non-zero means it survives this (coarse) structure control.

    Note this yields shrunk point estimates, not p-values: a non-zero
    coefficient is weaker evidence than a kinship-corrected test and should not
    be reported as though it were one.

    Zero-variance columns are dropped rather than standardised into NaN and
    then filled with 0, which silently retained them as all-zero covariates.
    """
    if not sig_features:
        return pd.Series(dtype=float)
    if not df.index.equals(st_dummies.index):
        raise CladeInputError(
            "the covariate matrix and the data frame are not aligned on the same index."
        )

    X = pd.concat([df[sig_features], st_dummies], axis=1).astype(float)
    y = validate_binary(df[phenotype_col], "phenotype", allow_missing=False).astype(int)

    sd = X.std(ddof=0)
    constant = sd[sd == 0].index.tolist()
    if constant:
        X = X.drop(columns=constant)
        sd = sd.drop(index=constant)

    X_std = (X - X.mean()) / sd
    if X_std.isna().any().any():
        X_std = X_std.fillna(0.0)

    lasso = LogisticRegression(
        penalty="l1", solver="liblinear", C=C, max_iter=2000, random_state=random_state
    )
    lasso.fit(X_std, y)

    coefs = pd.Series(lasso.coef_[0], index=X.columns)
    retained = [f for f in sig_features if f in coefs.index]
    return coefs[retained].sort_values(key=abs, ascending=False)

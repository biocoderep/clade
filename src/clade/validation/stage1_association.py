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
* `firth_association` uses neither of `firthlogist`'s inference outputs.

  An earlier revision of this module replaced the library's profile-likelihood
  p-value (`pvals_`) with a manual Wald p-value computed from `coef_`/`bse_`.
  That was wrong, and the note it left asking a future reader to "check `bse_`
  is actually populated" has now been acted on. Measured against an
  independent, unpenalised likelihood-ratio test on 37 real candidates
  (n=3,249, 20 lineage covariates):

      firthlogist `bse_` across all 37 candidates : 0.1200 - 0.1559  (1.3x)
      true standard error over the same candidates: 0.144  - 23531   (163000x)
      correlation between them (separation cases excluded): r = -0.18

  `bse_` is effectively a constant. It does not measure the precision of the
  estimate, so a Wald statistic built from it reduces to |coef| rescaled and
  carries no information about identifiability. Three unrelated candidates
  shared `bse_ = 0.119952` to six decimal places. A candidate with complete
  separation -- coefficient unidentified, true SE ~23,000 -- was assigned
  p = 3.7e-182.

  A 100-permutation lineage-preserving null confirmed the consequence: with
  the Wald path, 20.5 of 37 candidates were called significant at q<0.05 when
  no true association existed. Under Benjamini-Hochberg the expected count
  under a complete null is approximately zero, not the 5% figure sometimes
  quoted. The same arbitration found the library's own `pvals_` tracked the
  independent LRT more closely than the Wald did (median |log10| gap 2.42 vs
  3.79) -- so the earlier "fix" was worse than what it replaced.

  This module therefore computes its own likelihood-ratio test by refitting
  the reduced model, and reports the Firth coefficient for direction only.
  `check_separation` runs first, because where a candidate perfectly predicts
  the phenotype no method yields a meaningful p-value and the honest output
  is a flag, not a number.
"""
from __future__ import annotations

import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from statsmodels.stats.multitest import multipletests

from clade.io.validation import CladeInputError, is_missing, validate_binary


def _independent_fit(X, y, cand_idx: int) -> dict:
    """The base Stage 1 computation: an unpenalised logistic fit, its
    likelihood-ratio p-value, and its coefficient -- all from `statsmodels`,
    none from `firthlogist`.

    This exists as the *required* computation, not a fallback, because it is
    the only one of the two available fits whose behaviour at this covariate
    dimensionality has actually been measured (module docstring): `firthlogist`
    has no role in it at all, so a machine without that package installed
    (Python >=3.11, where it cannot even be built -- see §5.14) still gets a
    real Stage 1 result rather than a hard failure. When `firthlogist` *is*
    importable, `firth_association` additionally fits it and prefers its
    coefficient (Firth's penalty gives a more stable point estimate near
    separation); the p-value below is used either way.

    Returns `coef`/`se`/`lrt_p` as NaN, not raising, when the refit itself
    fails to converge -- a per-candidate fit failure is data, not a crash.
    """
    import numpy as np
    import statsmodels.api as sm
    from scipy.stats import chi2

    reduced_cols = [i for i in range(X.shape[1]) if i != cand_idx]
    Xf = sm.add_constant(X, has_constant="add")
    Xr = sm.add_constant(X[:, reduced_cols], has_constant="add")
    try:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            full = sm.Logit(y, Xf).fit(disp=0, maxiter=200)
            red = sm.Logit(y, Xr).fit(disp=0, maxiter=200)
        coef = float(full.params[-1])
        se = float(full.bse[-1])
        stat = 2.0 * (float(full.llf) - float(red.llf))
        lrt_p = float(chi2.sf(stat, 1)) if np.isfinite(stat) and stat >= 0 else float("nan")
        return {"coef": coef, "se": se, "lrt_p": lrt_p}
    except Exception:  # noqa: BLE001 - a non-converging refit is a real outcome here
        return {"coef": float("nan"), "se": float("nan"), "lrt_p": float("nan")}


def _likelihood_ratio_p(X, y, cand_idx: int) -> float:
    """Backwards-compatible wrapper around `_independent_fit` returning only
    the p-value. Prefer `_independent_fit` in new code -- it does the same
    refit once and also returns the coefficient, instead of fitting twice.
    """
    import numpy as np
    import statsmodels.api as sm
    from scipy.stats import chi2

    reduced_cols = [i for i in range(X.shape[1]) if i != cand_idx]
    Xf = sm.add_constant(X, has_constant="add")
    Xr = sm.add_constant(X[:, reduced_cols], has_constant="add")
    try:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            full = sm.Logit(y, Xf).fit(disp=0, maxiter=200)
            red = sm.Logit(y, Xr).fit(disp=0, maxiter=200)
        stat = 2.0 * (float(full.llf) - float(red.llf))
        if not np.isfinite(stat) or stat < 0:
            return float("nan")
        return float(chi2.sf(stat, 1))
    except Exception:  # noqa: BLE001 - a non-converging refit is a real outcome here
        return float("nan")


def check_separation(genotype: pd.Series, phenotype: pd.Series) -> dict:
    """Detect complete or quasi-complete separation from the 2x2 table.

    Separation is when a candidate perfectly (or almost perfectly) predicts the
    phenotype. The maximum-likelihood coefficient then diverges: it is
    unbounded, its standard error explodes, and every inference method produces
    an impressive-looking number that means nothing. Firth's penalty keeps the
    estimate finite, which is what it is for -- but finite is not the same as
    identified, and a finite estimate from a separated design still supports no
    claim about effect size.

    This is checked *before* fitting because every spurious Stage 1 result in
    this project's case study traced back to it. One candidate with zero
    susceptible carriers was reported at p = 1.6e-51; its true standard error
    was 23,425.

    Returns the cell counts, the minimum cell, and a status:

        none        every cell >= 10; standard inference is reasonable
        quasi       some cell in 1..9; estimates are unstable
        complete    some cell is 0; the coefficient is not identified
    """
    g = validate_binary(genotype, "genotype")
    y = validate_binary(phenotype, "phenotype", allow_missing=False)
    both = g.notna() & y.notna()
    g, y = g[both], y[both]

    n11 = int(((g == 1) & (y == 1)).sum())
    n10 = int(((g == 1) & (y == 0)).sum())
    n01 = int(((g == 0) & (y == 1)).sum())
    n00 = int(((g == 0) & (y == 0)).sum())
    min_cell = min(n11, n10, n01, n00)

    if min_cell == 0:
        status = "complete"
    elif min_cell < 10:
        status = "quasi"
    else:
        status = "none"

    return {
        "n_carrier_positive": n11,
        "n_carrier_negative": n10,
        "n_noncarrier_positive": n01,
        "n_noncarrier_negative": n00,
        "min_cell": min_cell,
        "separation": status,
        "identifiable": status == "none",
    }


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

    Does NOT require `firthlogist`. The base fit and the p-value are always
    computed by `_independent_fit` (plain `statsmodels`, no Firth penalty) --
    see that function's docstring for why this is the required computation,
    not a fallback. Where `firthlogist` (Python <3.11, scikit-learn <1.6 --
    see §5.14) *is* importable, it is fit as well and its coefficient is
    preferred, since Firth's penalty gives a more stable point estimate close
    to separation; the independent p-value is used either way. Its own
    p-value/standard-error outputs are never read (module docstring: neither
    was found reliable at this covariate dimensionality).
    """
    try:
        from firthlogist import FirthLogisticRegression
    except ImportError:
        FirthLogisticRegression = None

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

    # Separation is checked before fitting: where a candidate perfectly
    # predicts the phenotype the coefficient is not identified, and every
    # method returns an impressive number that supports no claim.
    sep = check_separation(cand, pheno)

    try:
        # Base computation: always run, never dependent on firthlogist.
        # Supplies both the p-value (used unconditionally) and the
        # coefficient (used unless a working Firth fit overrides it below).
        fit = _independent_fit(X, y, cand_idx)
        coef = fit["coef"]
        pval = fit["lrt_p"]
        engine = "statsmodels (unpenalised, independent)"

        # Optional refinement: prefer Firth's coefficient when available and
        # finite. Its p-value/bse_ are still never read (see module docstring).
        if FirthLogisticRegression is not None:
            try:
                fl = FirthLogisticRegression(test_vars=cand_idx)
                fl.fit(X, y)
                firth_coef = fl.coef_[cand_idx]
                if not is_missing(firth_coef):
                    coef = firth_coef
                    engine = "firthlogist coefficient + independent LRT p-value"
            except Exception:  # noqa: BLE001, S110 - Firth failing to converge is
                # not fatal; the independent fit above already stands, and there
                # is nothing candidate-specific worth logging at this volume.
                pass

        if is_missing(coef):
            return {
                "Candidate": candidate,
                "N_carriers": n_carriers,
                "N_dropped_missing": n_dropped,
                "Coef": None,
                "Firth_p": None,
                "Status": "FAILED: non-finite estimate",
                **sep,
            }

        status = "OK"
        if sep["separation"] == "complete":
            status = "SEPARATION: coefficient not identified; p-value not reportable"
            pval = None
        elif sep["separation"] == "quasi":
            status = f"QUASI_SEPARATION: min cell {sep['min_cell']}; estimate unstable"

        return {
            "Candidate": candidate,
            "N_carriers": n_carriers,
            "N_dropped_missing": n_dropped,
            "Coef": float(coef),
            "Firth_p": None if pval is None or is_missing(pval) else float(pval),
            "Status": status,
            # Which fit produced Coef, for audit -- distinct from Status, which
            # the CLI matches on exactly ("OK") to decide Stage 1's verdict and
            # must therefore stay a fixed vocabulary, not carry free text.
            "Engine": engine,
            **sep,
        }
    except Exception as e:  # noqa: BLE001 - real convergence failures are the point
        return {
            "Candidate": candidate,
            "N_carriers": n_carriers,
            "N_dropped_missing": n_dropped,
            "Coef": None,
            "Firth_p": None,
            "Status": f"FAILED: {e}",
            **sep,
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

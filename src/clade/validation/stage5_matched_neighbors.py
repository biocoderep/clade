"""
Stage 5 — Matched-neighbour comparison.

Each phenotype-positive genome is paired with its nearest phenotype-negative
genome by real patristic distance, and candidate carriage is compared within
pairs by McNemar's test. This operates at a fundamentally different scale from
Stages 1 and 4 — immediate local relatedness rather than population-wide or
whole-tree pattern — which is why a Stage 5 contradiction is reported as an
open finding rather than resolved by preferring another stage.

Three properties of this design are now measured and returned rather than left
implicit, because each can materially affect how the p-value should be read:

**Neighbour reuse (pseudoreplication).** Nothing stops one negative genome from
being the nearest neighbour of many positives. McNemar's test assumes the pairs
are independent; when a single negative anchors hundreds of pairs they are not,
and the p-value is anticonservative. The original implementation did this
silently. `matched_neighbor_test` now reports `Max_neighbor_reuse` and
`N_unique_neighbors` so the degree of non-independence is visible, and
`nearest_negative_neighbor(unique=True)` offers a one-to-one matching as a
sensitivity check. **The default remains reuse-permitted**, matching the
published analysis — this is instrumentation, not a silent change of method.

**Distance ties.** Exact ties are common in patristic distance matrices
(identical or near-identical sequences). `idxmin` resolves them by column
order, so the chosen neighbour depended on the order of the caller's sample
list. Ties are now broken by sorted sample ID, making the result reproducible
regardless of input ordering, and `n_ties` is reported.

**Missing genotypes.** An uncalled genotype was previously counted as absence,
because `NaN == 1` is False — so missing data silently became negative
evidence. Pairs where either member has a missing call for the candidate are
now excluded from that candidate's test and counted in `N_pairs_dropped_missing`.
"""
from __future__ import annotations

import pandas as pd
from statsmodels.stats.contingency_tables import mcnemar

from clade.io.validation import CladeInputError, validate_binary


def nearest_negative_neighbor(
    distances: pd.DataFrame,
    positive_samples: list[str],
    negative_samples: list[str],
    unique: bool = False,
) -> pd.DataFrame:
    """For each phenotype-positive sample, find its nearest phenotype-negative
    sample by real patristic distance.

    `distances`: a square, sample x sample distance matrix of real tree
    distance — not Euclidean distance in some other feature space.

    Ties are broken deterministically by sorted sample ID, so the pairing does
    not depend on the order of `negative_samples`.

    `unique=True` performs greedy one-to-one matching (closest pairs first,
    each negative used at most once), which removes the pseudoreplication in
    the default behaviour at the cost of dropping positives once negatives run
    out. Intended as a sensitivity analysis alongside the default, not as a
    replacement for it.

    Returns a DataFrame indexed by positive sample with `nearest_neighbor`,
    `distance`, and `tied` (whether the minimum distance was not unique).
    """
    if not positive_samples:
        raise CladeInputError("no phenotype-positive samples supplied to Stage 5.")
    if not negative_samples:
        raise CladeInputError(
            "no phenotype-negative samples available to match against; Stage 5 cannot "
            "form any pairs."
        )

    ordered_negatives = sorted(negative_samples)
    sub = distances.loc[positive_samples, ordered_negatives]

    min_dist = sub.min(axis=1)
    tied = sub.eq(min_dist, axis=0).sum(axis=1) > 1

    if not unique:
        return pd.DataFrame(
            {"nearest_neighbor": sub.idxmin(axis=1), "distance": min_dist, "tied": tied}
        )

    # Greedy one-to-one: consider candidate pairs in ascending distance order.
    long = (
        sub.stack()
        .rename("distance")
        .reset_index()
        .rename(columns={"level_0": "positive", "level_1": "negative"})
        .sort_values(["distance", "positive", "negative"], kind="mergesort")
    )
    used_pos, used_neg, rows = set(), set(), []
    for row in long.itertuples(index=False):
        if row.positive in used_pos or row.negative in used_neg:
            continue
        used_pos.add(row.positive)
        used_neg.add(row.negative)
        rows.append(
            {
                "positive": row.positive,
                "nearest_neighbor": row.negative,
                "distance": row.distance,
                "tied": bool(tied.get(row.positive, False)),
            }
        )
    return pd.DataFrame(rows).set_index("positive")


def matched_neighbor_test(
    genotypes: pd.DataFrame,
    positive_samples: list[str],
    matches: pd.DataFrame,
    candidate_cols: list[str],
    exact_threshold: int = 25,
) -> pd.DataFrame:
    """Paired McNemar's test of candidate carriage between each
    phenotype-positive sample and its matched nearest phenotype-negative
    neighbour (from `nearest_negative_neighbor`).

    Uses the exact McNemar variant when the number of discordant pairs is below
    `exact_threshold`, matching the original script.

    Pairs with a missing genotype call on either side are dropped for that
    candidate and counted, rather than being silently scored as absence.
    """
    paired = [s for s in positive_samples if s in matches.index]
    if not paired:
        raise CladeInputError(
            "none of the phenotype-positive samples appear in the matched-pair table."
        )

    neighbors = matches.loc[paired, "nearest_neighbor"]
    reuse = neighbors.value_counts()
    max_reuse = int(reuse.iloc[0]) if len(reuse) else 0
    n_unique = int(neighbors.nunique())
    n_tied = int(matches.loc[paired, "tied"].sum()) if "tied" in matches.columns else 0

    missing_geno = [s for s in list(paired) + list(neighbors) if s not in genotypes.index]
    if missing_geno:
        raise CladeInputError(
            f"{len(set(missing_geno))} sample(s) in the matched pairs are absent from the "
            f"genotype matrix, e.g. {sorted(set(missing_geno))[:5]}."
        )

    rows = []
    for cand in candidate_cols:
        if cand not in genotypes.columns:
            raise CladeInputError(f"candidate '{cand}' is not a column in the genotype matrix.")

        col = validate_binary(genotypes[cand], f"genotype column '{cand}'")
        pos_has = col.loc[paired].to_numpy()
        neighbor_has = col.loc[neighbors.to_numpy()].to_numpy()

        complete = pd.notna(pos_has) & pd.notna(neighbor_has)
        n_dropped = int((~complete).sum())
        p_ok = pos_has[complete]
        n_ok = neighbor_has[complete]

        both = int(((p_ok == 1) & (n_ok == 1)).sum())
        pos_only = int(((p_ok == 1) & (n_ok == 0)).sum())
        neighbor_only = int(((p_ok == 0) & (n_ok == 1)).sum())
        neither = int(((p_ok == 0) & (n_ok == 0)).sum())

        discordant = pos_only + neighbor_only
        table = [[both, pos_only], [neighbor_only, neither]]

        if discordant == 0:
            # McNemar is undefined with no discordant pairs; NaN, not 1.0.
            pval = float("nan")
            note = "no discordant pairs; McNemar undefined"
        else:
            try:
                pval = float(mcnemar(table, exact=discordant < exact_threshold).pvalue)
                note = "exact" if discordant < exact_threshold else "asymptotic"
            except Exception as exc:  # noqa: BLE001 - surfaced, not swallowed
                pval = float("nan")
                note = f"McNemar failed: {type(exc).__name__}: {exc}"

        rows.append(
            {
                "Candidate": cand,
                "N_pairs": int(complete.sum()),
                "Resistant_carries": pos_only + both,
                "Neighbor_carries": neighbor_only + both,
                "Resistant_only": pos_only,
                "Neighbor_only": neighbor_only,
                "Discordant_pairs": discordant,
                # Direction: enriched in the resistant member of the pair.
                # Plain Python bool: numpy's np.bool_ fails `is True` identity checks
                # in the disposition engine.
                "Enriched_in_resistant": (bool(pos_only > neighbor_only) if discordant else None),
                "McNemar_p": pval,
                "McNemar_note": note,
                "N_pairs_dropped_missing": n_dropped,
                "N_unique_neighbors": n_unique,
                "Max_neighbor_reuse": max_reuse,
                "N_tied_matches": n_tied,
            }
        )

    return pd.DataFrame(rows)

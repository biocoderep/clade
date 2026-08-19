"""
Stage 5 — Matched-neighbor comparison.

Direct refactor of matched_neighbor.py into testable functions. This stage
operates at a fundamentally different scale than Stages 1 and 4 (immediate
local relatedness, via real patristic distance, rather than population-wide
or whole-tree patterns) — in the CLADE case study this was the one stage
that produced genuine, unresolved disagreement with the others for three
candidates. Report a Stage 5 contradiction as an open finding, not
something to resolve by preferring another stage's verdict by convention.
"""
from __future__ import annotations

import pandas as pd
from statsmodels.stats.contingency_tables import mcnemar


def nearest_negative_neighbor(distances: pd.DataFrame, positive_samples: list[str], negative_samples: list[str]) -> pd.DataFrame:
    """For each phenotype-positive sample, find its single nearest
    phenotype-negative sample by real patristic distance.

    `distances`: a square, sample x sample distance matrix (real tree
    distance, e.g. from a precomputed pairwise distance matrix — not
    Euclidean distance in some other feature space).
    Returns a DataFrame indexed by positive sample, with columns
    `nearest_neighbor` and `distance`.
    """
    sub = distances.loc[positive_samples, negative_samples]
    return pd.DataFrame({
        "nearest_neighbor": sub.idxmin(axis=1),
        "distance": sub.min(axis=1),
    })


def matched_neighbor_test(
    genotypes: pd.DataFrame,
    positive_samples: list[str],
    matches: pd.DataFrame,
    candidate_cols: list[str],
) -> pd.DataFrame:
    """Paired McNemar's test of candidate carriage between each
    phenotype-positive sample and its matched nearest phenotype-negative
    neighbor (from `nearest_negative_neighbor`).

    Uses the exact McNemar variant when the number of discordant pairs is
    small (<25), matching the original script's threshold.
    """
    rows = []
    for cand in candidate_cols:
        pos_has = genotypes.loc[positive_samples, cand].values
        neighbor_has = genotypes.loc[matches["nearest_neighbor"].values, cand].values

        both = int(((pos_has == 1) & (neighbor_has == 1)).sum())
        pos_only = int(((pos_has == 1) & (neighbor_has == 0)).sum())
        neighbor_only = int(((pos_has == 0) & (neighbor_has == 1)).sum())
        neither = int(((pos_has == 0) & (neighbor_has == 0)).sum())

        table = [[both, pos_only], [neighbor_only, neither]]
        try:
            result = mcnemar(table, exact=(pos_only + neighbor_only) < 25)
            pval = result.pvalue
        except Exception:  # noqa: BLE001
            pval = None

        rows.append({
            "Candidate": cand,
            "Resistant_carries": pos_only + both,
            "Neighbor_carries": neighbor_only + both,
            "Resistant_only": pos_only,
            "Neighbor_only": neighbor_only,
            "McNemar_p": pval,
        })

    return pd.DataFrame(rows)

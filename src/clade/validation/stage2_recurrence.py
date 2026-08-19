"""
Stage 2 — Clonal-recurrence check.

Direct refactor of st_crosstab_remaining.py into a testable function. Does
not reject candidates by itself — a single-lineage candidate remains
eligible for later stages, but its Stage-2 result is a required caveat on
how any Stage-1 significance should be read (see CLADE_Framework_Specification.md).
"""
from __future__ import annotations

import pandas as pd


def lineage_recurrence(genotypes: pd.DataFrame, lineages: pd.Series, candidate_cols: list[str] | None = None) -> pd.DataFrame:
    """For each candidate column, count carriers, distinct lineages (STs)
    among carriers, and what fraction sit in the single largest lineage.

    `genotypes`: samples x candidates (0/1), indexed by sample ID.
    `lineages`: sample ID -> ST, same index space.
    """
    if candidate_cols is None:
        candidate_cols = list(genotypes.columns)

    df = genotypes.join(lineages.rename("ST"), how="inner")

    rows = []
    for col in candidate_cols:
        carriers = df[df[col] == 1]
        n = len(carriers)
        if n == 0:
            rows.append({"Candidate": col, "Carriers": 0, "Distinct_STs": 0, "Pct_dominant_ST": None, "Dominant_ST": None})
            continue
        vc = carriers["ST"].value_counts()
        rows.append({
            "Candidate": col,
            "Carriers": n,
            "Distinct_STs": carriers["ST"].nunique(),
            "Pct_dominant_ST": round(vc.iloc[0] / n, 3),
            "Dominant_ST": vc.index[0],
        })

    return pd.DataFrame(rows).sort_values("Distinct_STs").reset_index(drop=True)

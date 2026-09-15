"""
Stage 2 — Clonal-recurrence check.

Asks whether a candidate appears across many lineages or is confined to one.
A variant restricted to a single sequence type may simply have been inherited
by every descendant of one successful clone; one recurring across many lineages
is harder to explain that way.

Stage 2 does **not** reject candidates by itself, and the disposition engine
treats it as informative rather than as a gate: a single-lineage candidate can
still be genuine if it arose independently several times *within* that lineage
(the SecA case, 100% ST2 but five independent origins by Stage 4). Read
`Distinct_STs` together with Stage 4's `Total_gains`, never alone.

The previous implementation inner-joined genotypes to lineages and silently
dropped every sample missing a lineage assignment, so `Carriers` could be far
lower than the true carrier count with no indication. That loss is now
reported.
"""
from __future__ import annotations

import pandas as pd

from clade.io.validation import validate_genotype_frame


def lineage_recurrence(
    genotypes: pd.DataFrame,
    lineages: pd.Series,
    candidate_cols: list[str] | None = None,
    warn_on_loss: bool = True,
) -> pd.DataFrame:
    """For each candidate column, count carriers, distinct lineages (STs)
    among carriers, and what fraction sit in the single largest lineage.

    `genotypes`: samples x candidates (0/1), indexed by sample ID.
    `lineages`: sample ID -> ST, over the same index space.

    Returned columns include `Carriers_without_ST`, the number of carriers that
    had no lineage assignment and so could not contribute to the ST counts.
    `Pct_dominant_ST` is computed over carriers *with* a known ST, so its
    denominator is stated explicitly rather than mixed.
    """
    validated = validate_genotype_frame(genotypes, candidate_cols)
    cols = list(validated.columns)

    st = lineages.rename("ST")
    df = validated.join(st, how="left")

    rows = []
    for col in cols:
        carriers = df[df[col] == 1]
        n_total = len(carriers)
        with_st = carriers[carriers["ST"].notna()]
        n_with_st = len(with_st)

        if n_total == 0:
            rows.append(
                {
                    "Candidate": col,
                    "Carriers": 0,
                    "Carriers_with_ST": 0,
                    "Carriers_without_ST": 0,
                    "Distinct_STs": 0,
                    "Pct_dominant_ST": None,
                    "Dominant_ST": None,
                }
            )
            continue

        if n_with_st == 0:
            rows.append(
                {
                    "Candidate": col,
                    "Carriers": n_total,
                    "Carriers_with_ST": 0,
                    "Carriers_without_ST": n_total,
                    "Distinct_STs": 0,
                    "Pct_dominant_ST": None,
                    "Dominant_ST": None,
                }
            )
            continue

        vc = with_st["ST"].value_counts()
        rows.append(
            {
                "Candidate": col,
                "Carriers": n_total,
                "Carriers_with_ST": n_with_st,
                "Carriers_without_ST": n_total - n_with_st,
                "Distinct_STs": int(with_st["ST"].nunique()),
                "Pct_dominant_ST": round(vc.iloc[0] / n_with_st, 3),
                "Dominant_ST": vc.index[0],
            }
        )

    result = pd.DataFrame(rows).sort_values("Distinct_STs").reset_index(drop=True)

    if warn_on_loss:
        n_missing = int(df["ST"].isna().sum())
        if n_missing:
            result.attrs["lineage_warning"] = (
                f"{n_missing} of {len(df)} samples have no lineage assignment and were "
                "excluded from the ST counts (but are still counted as carriers)."
            )

    return result

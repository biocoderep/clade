"""
Build the final, per-candidate evidence table — the single artifact CLADE
is meant to produce. Every candidate that entered the framework must appear
here, whatever its disposition; a table that only lists survivors cannot be
audited (see CLADE_Framework_Specification.md, Section 4: "Reporting rule").
"""
from __future__ import annotations

import pandas as pd

from clade.classification.disposition import Disposition, StageResults, classify_candidate


def build_evidence_table(candidates: dict[str, StageResults]) -> pd.DataFrame:
    """`candidates`: candidate name -> StageResults. Returns one row per
    candidate with every stage's raw values plus the final Disposition —
    nothing filtered out.
    """
    rows = []
    for name, r in candidates.items():
        disposition = classify_candidate(r)
        rows.append({
            "Candidate": name,
            "Stage1_tested": r.stage1_tested,
            "Stage1_significant": r.stage1_significant,
            "Stage2_distinct_STs": r.stage2_distinct_sts,
            "Stage2_pct_dominant_ST": r.stage2_pct_dominant_st,
            "Stage3_correct_direction": r.stage3_correct_direction,
            "Stage4_pct_post_resistance": r.stage4_pct_post_resistance,
            "Stage4_baseline": r.stage4_baseline,
            "Stage5_significant": r.stage5_significant,
            "Stage5_concordant": r.stage5_concordant,
            "Cohort_frequency": r.cohort_frequency,
            "Disposition": disposition.value,
        })
    return pd.DataFrame(rows)


def evidence_table_to_markdown(table: pd.DataFrame) -> str:
    """Render the evidence table as Markdown, grouped by Disposition tier —
    same tiered structure as 14_FINAL_evidence_table.md, generated instead
    of hand-written."""
    lines = ["# CLADE Evidence Table\n"]
    tier_order = [d.value for d in Disposition]
    for tier in tier_order:
        subset = table[table["Disposition"] == tier]
        if subset.empty:
            continue
        lines.append(f"## {tier} ({len(subset)} candidates)\n")
        lines.append(subset.drop(columns=["Disposition"]).to_markdown(index=False))
        lines.append("")
    return "\n".join(lines)

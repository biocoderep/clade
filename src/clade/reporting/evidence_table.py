"""
Build the final, per-candidate evidence table — the single artifact CLADE is
meant to produce.

Every candidate that entered the framework must appear here, whatever its
disposition; a table that only lists survivors cannot be audited (see
CLADE_Framework_Specification.md, Section 4: "Reporting rule").

Two columns exist specifically so a reader can tell *why* a verdict was
reached without re-running anything: `Stages_missing` names the stages that
produced no usable result, and `Reason` carries the disposition engine's own
justification. A candidate is never promoted by silence, and the table now
shows where the silence was.
"""
from __future__ import annotations

import pandas as pd

from clade.classification.disposition import (
    Disposition,
    DispositionResult,
    StageResults,
    classify_candidate,
)

#: Tri-state values are rendered as text so that "not run" is visually distinct
#: from "False" in the Markdown output, where an empty cell is ambiguous.
_TRISTATE = {True: "pass", False: "FAIL", None: "not run"}


def _fmt_tristate(v) -> str:
    return _TRISTATE[v] if v in _TRISTATE else str(v)


def build_evidence_table(
    candidates: dict[str, StageResults],
    temporal_margin: float = 0.10,
    **classify_kwargs,
) -> pd.DataFrame:
    """`candidates`: candidate name -> StageResults.

    Returns one row per candidate with every stage's raw values, the derived
    per-stage pass/fail/not-run status, the final disposition and the reason
    for it. Nothing is filtered out.
    """
    rows = []
    for name, r in candidates.items():
        verdict: DispositionResult = classify_candidate(
            r, temporal_margin=temporal_margin, **classify_kwargs
        )
        status = verdict.stage_status
        rows.append(
            {
                "Candidate": name,
                "Disposition": verdict.disposition.value,
                "Stage1_status": _fmt_tristate(status.get("stage1")),
                "Stage3_status": _fmt_tristate(status.get("stage3")),
                "Stage4_status": _fmt_tristate(status.get("stage4")),
                "Stage5_status": _fmt_tristate(status.get("stage5")),
                "Stage1_tested": r.stage1_tested,
                "Stage1_significant": r.stage1_significant,
                "Stage2_distinct_STs": r.stage2_distinct_sts,
                "Stage2_pct_dominant_ST": r.stage2_pct_dominant_st,
                "Stage3_correct_direction": r.stage3_correct_direction,
                "Stage4_pct_post_resistance": r.stage4_pct_post_resistance,
                "Stage4_baseline": r.stage4_baseline,
                "Stage4_total_gains": r.stage4_total_gains,
                "Stage5_significant": r.stage5_significant,
                "Stage5_enriched_in_resistant": r.stage5_enriched_in_resistant,
                "Stage6_mechanism_found": r.stage6_mechanism_found,
                "Cohort_frequency": r.cohort_frequency,
                "Stages_missing": ", ".join(verdict.stages_missing) or "-",
                "Reason": "; ".join(verdict.reasons),
            }
        )
    return pd.DataFrame(rows)


#: Narrative shown under each disposition heading, so the table is readable
#: without the specification to hand.
_TIER_NOTES = {
    Disposition.CONVERGENT.value: (
        "Explicit positive result at every stage required for convergence. "
        "This is convergent *computational* evidence, not experimental validation."
    ),
    Disposition.UNRESOLVED.value: (
        "Stages genuinely disagree — local matched-pair evidence points one way, "
        "population-wide or whole-tree evidence the other. Reported as an open "
        "question, not resolved by preferring one stage."
    ),
    Disposition.REJECTED.value: "At least one stage returned an explicit negative result.",
    Disposition.CLONAL_ARTIFACT.value: (
        "Signal is concentrated in a single lineage with no evidence of independent "
        "recurrence, and at least one stage contradicts."
    ),
    Disposition.UNINFORMATIVE.value: (
        "Near-fixed across the cohort; no association test on it is informative."
    ),
    Disposition.INSUFFICIENT_EVIDENCE.value: (
        "No stage contradicts, but one or more required stages produced no usable "
        "result. Absence of contradiction is not support — these candidates are "
        "under-tested, not validated."
    ),
    Disposition.NOT_RETESTED.value: (
        "Stage 1 was never run for these candidates. A coverage limitation, not a "
        "scientific verdict on them."
    ),
}

#: Columns shown in the Markdown summary. The full frame keeps everything.
_SUMMARY_COLUMNS = [
    "Candidate",
    "Stage1_status",
    "Stage3_status",
    "Stage4_status",
    "Stage5_status",
    "Stage2_distinct_STs",
    "Stage2_pct_dominant_ST",
    "Cohort_frequency",
    "Stages_missing",
]


def evidence_table_to_markdown(table: pd.DataFrame, include_reasons: bool = True) -> str:
    """Render the evidence table as Markdown, grouped by disposition tier."""
    lines = [
        "# CLADE Evidence Table",
        "",
        "Every candidate that entered the framework appears below, whatever the outcome.",
        "",
        ("Stage status: `pass` = explicit positive, `FAIL` = explicit negative, "
         "`not run` = no usable result for that stage."),
        "",
    ]

    counts = table["Disposition"].value_counts().to_dict()
    lines.append("| Disposition | Candidates |")
    lines.append("|---|---:|")
    for tier in [d.value for d in Disposition]:
        if tier in counts:
            lines.append(f"| {tier} | {counts[tier]} |")
    lines.append(f"| **total** | **{len(table)}** |")
    lines.append("")

    for tier in [d.value for d in Disposition]:
        subset = table[table["Disposition"] == tier]
        if subset.empty:
            continue
        lines.append(f"## {tier} ({len(subset)} candidates)")
        lines.append("")
        note = _TIER_NOTES.get(tier)
        if note:
            lines.append(f"*{note}*")
            lines.append("")
        cols = [c for c in _SUMMARY_COLUMNS if c in subset.columns]
        lines.append(subset[cols].to_markdown(index=False))
        lines.append("")
        if include_reasons and "Reason" in subset.columns:
            for _, row in subset.iterrows():
                lines.append(f"- **{row['Candidate']}** — {row['Reason']}")
            lines.append("")

    return "\n".join(lines)

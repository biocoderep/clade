"""
Stage 3 — Effect-direction check.

In the original project this was never a separate script — it was a
manual read of Stage 1's regression coefficient sign, applied for the
first time only partway through the analysis (see
docs/reproducibility/original_discovery_pipeline.md: this single check
disqualified the majority of candidates that had passed Stage 1). CLADE
codifies it as an explicit, required step rather than an optional glance,
which is the actual point of naming it as its own stage.
"""
from __future__ import annotations


def direction_check(coefficient: float | None) -> dict:
    """A genuine compensatory mutation should be enriched (positive
    coefficient), not depleted, among resistant genomes.

    Returns a dict with the raw coefficient, a boolean `correct_direction`,
    and a human-readable `verdict`. `coefficient=None` (e.g. a failed
    Stage-1 fit) yields `correct_direction=None` and a verdict explaining why.
    """
    if coefficient is None:
        return {"Coefficient": None, "correct_direction": None, "verdict": "No Stage 1 coefficient available"}
    correct = coefficient > 0
    verdict = "correct (positive)" if correct else "wrong (negative) — anti-correlated with resistance"
    return {"Coefficient": coefficient, "correct_direction": correct, "verdict": verdict}

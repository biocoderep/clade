"""
Final candidate disposition — combining all six stages into one verdict.

Honesty note: in the CLADE case study, this classification was originally
done by hand while writing the evidence table, and involved some
case-by-case judgment for borderline candidates (e.g. Phosphodiesterase was
called "weak" rather than cleanly convergent or rejected, because it passed
direction but failed the naive full-cohort screen and was near-single-clone).
This module codifies the *general* rule that covers the clear-cut majority
of cases — it will not perfectly reproduce every nuanced judgment call from
the original manual table, and callers working through a genuinely
borderline candidate should read the underlying stage results themselves,
not just trust the returned enum.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Disposition(Enum):
    CONVERGENT = "convergent"
    UNRESOLVED = "unresolved"
    REJECTED = "rejected"
    CLONAL_ARTIFACT = "clonal_artifact"
    UNINFORMATIVE = "uninformative"
    NOT_RETESTED = "not_individually_retested"


@dataclass
class StageResults:
    stage1_tested: bool
    stage1_significant: bool | None = None
    stage2_distinct_sts: int | None = None
    stage2_pct_dominant_st: float | None = None
    stage3_correct_direction: bool | None = None
    stage4_pct_post_resistance: float | None = None
    stage4_baseline: float | None = None
    stage5_significant: bool | None = None
    stage5_concordant: bool | None = None
    cohort_frequency: float | None = None


def classify_candidate(
    r: StageResults,
    uninformative_freq_threshold: float = 0.85,
    single_lineage_threshold: float = 0.95,
    temporal_margin: float = 0.10,
) -> Disposition:
    """Apply CLADE's disposition rule to one candidate's stage results.

    A stage that was not applicable/not run (None) does not block
    CONVERGENT — only an explicit negative result does. This matches the
    framework's own definition: "passes every stage for which a valid test
    could be run."
    """
    if r.cohort_frequency is not None and r.cohort_frequency >= uninformative_freq_threshold:
        return Disposition.UNINFORMATIVE

    if not r.stage1_tested:
        return Disposition.NOT_RETESTED

    if r.stage1_significant is False:
        return Disposition.REJECTED

    passes_temporal = None
    if r.stage4_pct_post_resistance is not None and r.stage4_baseline is not None:
        passes_temporal = r.stage4_pct_post_resistance >= (r.stage4_baseline + temporal_margin)

    any_explicit_failure = (
        r.stage3_correct_direction is False
        or passes_temporal is False
        or r.stage5_concordant is False
    )

    if not any_explicit_failure:
        return Disposition.CONVERGENT

    stage5_contradicts = (
        r.stage5_significant is True
        and r.stage5_concordant is False
        and (r.stage3_correct_direction is False or passes_temporal is False)
    )
    if stage5_contradicts:
        return Disposition.UNRESOLVED

    if r.stage2_pct_dominant_st is not None and r.stage2_pct_dominant_st >= single_lineage_threshold:
        return Disposition.CLONAL_ARTIFACT

    return Disposition.REJECTED

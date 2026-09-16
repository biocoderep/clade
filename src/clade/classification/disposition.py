"""
Final candidate disposition — combining the six stages into one verdict.

Design rule (the thing this module exists to enforce): **missing evidence is
never treated as evidence.** Three states are kept distinct throughout:

    True   an explicit positive result for that stage
    False  an explicit negative result for that stage
    None   the stage was not run, or could not produce a valid result

A candidate reaches CONVERGENT only when the stages that actually constitute
convergence produced explicit positive results. A candidate with Stage 1
significance and nothing else is *not* convergent — it is under-tested, and is
reported as such (INSUFFICIENT_EVIDENCE) rather than being promoted by the
absence of contradicting data.

Honesty note (carried forward from the original implementation): in the CLADE
case study this classification was originally done by hand while writing the
evidence table, and involved case-by-case judgment for borderline candidates
(e.g. Phosphodiesterase was called "weak" rather than cleanly convergent or
rejected, because it passed direction but failed the naive full-cohort screen
and was near-single-clone). This module codifies the *general* rule that covers
the clear-cut majority of cases. It will not reproduce every nuanced judgment
call from the original manual table, and callers working through a genuinely
borderline candidate should read `DispositionResult.reasons` and the underlying
stage results, not just the returned enum.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from clade.io.validation import is_missing


class Disposition(Enum):
    CONVERGENT = "convergent"
    UNRESOLVED = "unresolved"
    REJECTED = "rejected"
    CLONAL_ARTIFACT = "clonal_artifact"
    UNINFORMATIVE = "uninformative"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"
    NOT_RETESTED = "not_individually_retested"


#: Stages that must have produced an explicit positive result before a
#: candidate can be called CONVERGENT. Stage 2 is deliberately excluded: a
#: single-lineage candidate can still be convergent if it recurs independently
#: within that lineage (the SecA case), so Stage 2 informs rather than gates.
#: Stage 6 is excluded because it is a manual process that may legitimately be
#: outstanding; it is reported separately via `stage6_mechanism_found`.
CONVERGENCE_REQUIRED_STAGES = ("stage1", "stage3", "stage4", "stage5")


@dataclass
class StageResults:
    """Per-candidate stage evidence.

    Every stage field is tri-state. Do not coerce a not-run stage to False
    before constructing this object — that is precisely the bug this type is
    shaped to prevent.
    """

    # Stage 1 — structure-corrected association
    stage1_tested: bool = False
    stage1_significant: bool | None = None

    # Stage 2 — clonal recurrence (informative, not a gate)
    stage2_distinct_sts: int | None = None
    stage2_pct_dominant_st: float | None = None
    stage2_independent_origins: int | None = None

    # Stage 3 — effect direction
    stage3_correct_direction: bool | None = None

    # Stage 4 — temporal ordering
    stage4_pct_post_resistance: float | None = None
    stage4_baseline: float | None = None
    stage4_total_gains: int | None = None

    # Stage 5 — matched neighbours
    stage5_significant: bool | None = None
    stage5_enriched_in_resistant: bool | None = None

    # Stage 6 — external corroboration (manual)
    stage6_mechanism_found: bool | None = None

    cohort_frequency: float | None = None

    # ---- derived, tri-state ----

    def stage4_supports(
        self, temporal_margin: float = 0.10, min_origins: int = 2
    ) -> bool | None:
        """True/False if Stage 4 ran and is interpretable, else None.

        Requires at least `min_origins` independently reconstructed gain
        events. One origin is not evidence of temporal ordering: a single
        event necessarily yields 0% or 100% post-resistance, and 100% from one
        event is indistinguishable from coincidence. Three candidates in this
        project's case study were reported at "100% post-resistance" on the
        strength of exactly one reconstructed origin.

        Zero or one origin therefore returns None -- not evaluable -- never
        False. Treating "too few events to judge" as a failure is the same
        missing-as-negative error this class exists to prevent, in the other
        direction.
        """
        if self.stage4_pct_post_resistance is None or self.stage4_baseline is None:
            return None
        if self.stage4_total_gains is not None and self.stage4_total_gains < min_origins:
            return None
        pct = self.stage4_pct_post_resistance
        if is_missing(pct):
            return None
        return pct >= (self.stage4_baseline + temporal_margin)

    def stage5_supports(self) -> bool | None:
        """Stage 5 supports the candidate only if it is both significant and
        enriched in the resistant member of each pair. Significance alone is
        direction-free and cannot support a compensatory claim."""
        if self.stage5_significant is None:
            return None
        if not self.stage5_significant:
            return False
        if self.stage5_enriched_in_resistant is None:
            return None
        return bool(self.stage5_enriched_in_resistant)

    def stage_status(
        self, temporal_margin: float = 0.10, min_origins: int = 2
    ) -> dict[str, bool | None]:
        return {
            "stage1": self.stage1_significant if self.stage1_tested else None,
            "stage3": self.stage3_correct_direction,
            "stage4": self.stage4_supports(temporal_margin, min_origins),
            "stage5": self.stage5_supports(),
        }


@dataclass
class DispositionResult:
    """An auditable verdict: the category, why it was reached, and which
    stages contributed."""

    disposition: Disposition
    reasons: list[str] = field(default_factory=list)
    stage_status: dict[str, bool | None] = field(default_factory=dict)
    stages_supporting: tuple[str, ...] = ()
    stages_contradicting: tuple[str, ...] = ()
    stages_missing: tuple[str, ...] = ()

    @property
    def value(self) -> str:
        return self.disposition.value


def classify_candidate(
    r: StageResults,
    uninformative_freq_threshold: float = 0.85,
    single_lineage_threshold: float = 0.95,
    temporal_margin: float = 0.10,
    min_origins: int = 2,
) -> DispositionResult:
    """Apply CLADE's disposition rule to one candidate's stage results.

    Returns a `DispositionResult`. For backwards compatibility the result
    exposes `.value`, and `classify_candidate_enum` returns the bare enum.

    Order of adjudication:

    1. Near-fixed across the cohort      -> UNINFORMATIVE (no test is informative)
    2. Stage 1 never run                 -> NOT_RETESTED  (coverage limit, not a verdict)
    3. Stage 1 ran and was negative      -> REJECTED
    4. Stage 1 ran, result indeterminate -> INSUFFICIENT_EVIDENCE
    5. Explicit contradictions present   -> UNRESOLVED / CLONAL_ARTIFACT / REJECTED
    6. All required stages positive      -> CONVERGENT
    7. Otherwise (gaps, no conflict)     -> INSUFFICIENT_EVIDENCE
    """
    reasons: list[str] = []
    status = r.stage_status(temporal_margin, min_origins)

    # 1. Near-fixed candidates are uninformative regardless of any statistic.
    if r.cohort_frequency is not None and r.cohort_frequency >= uninformative_freq_threshold:
        reasons.append(
            f"cohort frequency {r.cohort_frequency:.3f} >= {uninformative_freq_threshold}: "
            "near-fixed, no test on it is informative"
        )
        return _result(Disposition.UNINFORMATIVE, reasons, status)

    # 2. Never tested is a coverage limitation, not a scientific verdict.
    if not r.stage1_tested:
        reasons.append("Stage 1 was not run for this candidate (coverage limitation, not a verdict)")
        return _result(Disposition.NOT_RETESTED, reasons, status)

    # 3. Explicit Stage 1 negative.
    if r.stage1_significant is False:
        reasons.append("Stage 1 ran and returned a non-significant association")
        return _result(Disposition.REJECTED, reasons, status)

    # 4. Stage 1 ran but produced no usable answer (e.g. non-convergence).
    if r.stage1_significant is None:
        reasons.append(
            "Stage 1 was attempted but produced no valid result (e.g. model non-convergence); "
            "this is missing evidence, not negative evidence"
        )
        return _result(Disposition.INSUFFICIENT_EVIDENCE, reasons, status)

    supporting = tuple(s for s in CONVERGENCE_REQUIRED_STAGES if status.get(s) is True)
    contradicting = tuple(s for s in CONVERGENCE_REQUIRED_STAGES if status.get(s) is False)
    missing = tuple(s for s in CONVERGENCE_REQUIRED_STAGES if status.get(s) is None)

    # 5. Contradictions.
    if contradicting:
        # Stage 5 (local, matched-pair) disagreeing with the population-wide or
        # whole-tree stages is a genuine open finding, not something to resolve
        # by preferring one stage by convention.
        population_wide_fail = any(s in contradicting for s in ("stage3", "stage4"))
        if r.stage5_supports() is True and population_wide_fail:
            reasons.append(
                "Stage 5 (local matched-pair) supports the candidate while "
                f"{'/'.join(s for s in contradicting if s in ('stage3', 'stage4'))} "
                "(population-wide / whole-tree) contradict it — reported as an open "
                "disagreement rather than resolved by convention"
            )
            return _result(Disposition.UNRESOLVED, reasons, status, supporting, contradicting, missing)

        if (
            r.stage2_pct_dominant_st is not None
            and r.stage2_pct_dominant_st >= single_lineage_threshold
            and (r.stage2_independent_origins is None or r.stage2_independent_origins <= 1)
        ):
            reasons.append(
                f"{r.stage2_pct_dominant_st:.1%} of carriers sit in a single lineage with no "
                "evidence of independent recurrence, and at least one stage contradicts: "
                "signal is explained by clonal expansion"
            )
            return _result(Disposition.CLONAL_ARTIFACT, reasons, status, supporting, contradicting, missing)

        reasons.append(f"explicit negative result at: {', '.join(contradicting)}")
        return _result(Disposition.REJECTED, reasons, status, supporting, contradicting, missing)

    # 6. Convergence requires explicit positives, not merely an absence of failures.
    if not missing:
        reasons.append(
            "explicit positive result at every stage required for convergence: "
            + ", ".join(CONVERGENCE_REQUIRED_STAGES)
        )
        if r.stage6_mechanism_found is True:
            reasons.append("Stage 6 found mechanistic/literature support")
        elif r.stage6_mechanism_found is False:
            reasons.append(
                "Stage 6 search performed and found no mechanistic support — convergent on "
                "computational evidence only"
            )
        else:
            reasons.append("Stage 6 (manual corroboration) not recorded for this candidate")
        return _result(Disposition.CONVERGENT, reasons, status, supporting, contradicting, missing)

    # 7. No contradiction, but the evidence base is incomplete.
    reasons.append(
        "no stage contradicts, but these stages produced no usable result: "
        + ", ".join(missing)
        + " — insufficient evidence for convergence (absence of contradiction is not support)"
    )
    return _result(Disposition.INSUFFICIENT_EVIDENCE, reasons, status, supporting, contradicting, missing)


def _result(disposition, reasons, status, supporting=(), contradicting=(), missing=()):
    return DispositionResult(
        disposition=disposition,
        reasons=reasons,
        stage_status=status,
        stages_supporting=supporting,
        stages_contradicting=contradicting,
        stages_missing=missing,
    )


def classify_candidate_enum(r: StageResults, **kwargs) -> Disposition:
    """The bare enum, for callers that do not need the audit trail."""
    return classify_candidate(r, **kwargs).disposition

"""
Scientific control tests: small synthetic datasets whose correct answer is
known by construction.

These are the regression tests that matter. A refactor that keeps every unit
test green but flips one of these has broken the science, not the code.

Nothing here is biological evidence — these are software validation datasets
with hand-derived expected answers.
"""
from __future__ import annotations

import dendropy
import numpy as np
import pandas as pd
import pytest

from clade.classification.disposition import (
    Disposition,
    StageResults,
    classify_candidate,
)
from clade.io.validation import CladeInputError
from clade.phylogeny.fitch import fitch_reconstruct
from clade.validation.stage2_recurrence import lineage_recurrence
from clade.validation.stage3_direction import direction_check
from clade.validation.stage4_temporal_order import (
    cohort_baseline_prevalence,
    temporal_ordering,
    temporal_ordering_sensitivity,
)
from clade.validation.stage5_matched_neighbors import (
    matched_neighbor_test,
    nearest_negative_neighbor,
)

# --------------------------------------------------------------------------
# Stage 4 / Fitch — known evolutionary histories
#
#            /-A          Resistance: A, B      (susceptible: C, D)
#         /-|
#        |   \-B
#      --|
#        |   /-C
#         \-|
#            \-D
# --------------------------------------------------------------------------

FOUR_TAXON = "((A,B),(C,D));"
PHENO = {"A": 1, "B": 1, "C": 0, "D": 0}


def _tree(newick: str = FOUR_TAXON) -> dendropy.Tree:
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def _gains(cand: dict, pheno: dict = PHENO, newick: str = FOUR_TAXON, **kw) -> pd.Series:
    out = temporal_ordering(_tree(newick), pheno, {"c": cand}, **kw)
    return out.iloc[0]


def test_clonal_inheritance_is_one_origin():
    """Present in both members of one clade -> a single ancestral gain.

    This is the confound CLADE exists to catch: perfect association with
    resistance, one evolutionary event."""
    r = _gains({"A": 1, "B": 1, "C": 0, "D": 0})
    assert r["Total_gains"] == 1


def test_two_independent_origins_are_counted_separately():
    """Same mutation in two unrelated clades -> two independent gains."""
    r = _gains({"A": 1, "B": 0, "C": 1, "D": 0})
    assert r["Total_gains"] == 2


def test_gain_on_resistant_background_is_post_resistance():
    """A gain restricted to the resistant clade arises on a resistant parent."""
    r = _gains({"A": 1, "B": 0, "C": 0, "D": 0})
    assert r["Total_gains"] == 1
    assert r["Gains_on_resistant_bg"] == 1
    assert r["Pct_post_resistance"] == 1.0


def test_gain_on_susceptible_background_is_pre_resistance():
    """A gain restricted to the susceptible clade cannot be compensatory."""
    r = _gains({"A": 0, "B": 0, "C": 1, "D": 0})
    assert r["Total_gains"] == 1
    assert r["Gains_on_susceptible_bg"] == 1
    assert r["Pct_post_resistance"] == 0.0


def test_ancestral_mutation_is_not_called_post_resistance():
    """A mutation present in every genome predates everything: zero gains, and
    an undefined (NaN) percentage that must not be read as 0% or 100%."""
    r = _gains({"A": 1, "B": 1, "C": 1, "D": 1})
    assert r["Total_gains"] == 0
    assert np.isnan(r["Pct_post_resistance"])


def test_absent_candidate_yields_no_gains():
    r = _gains({"A": 0, "B": 0, "C": 0, "D": 0})
    assert r["Total_gains"] == 0


def test_mixed_origins_split_by_background():
    """One gain in the resistant clade, one in the susceptible clade."""
    r = _gains({"A": 1, "B": 0, "C": 1, "D": 0})
    assert r["Gains_on_resistant_bg"] == 1
    assert r["Gains_on_susceptible_bg"] == 1
    assert r["Pct_post_resistance"] == 0.5


def test_cohort_baseline_ignores_missing():
    assert cohort_baseline_prevalence({"A": 1, "B": 1, "C": 0, "D": None}) == pytest.approx(2 / 3)
    assert cohort_baseline_prevalence({"A": 1, "B": 0}) == 0.5


# --------------------------------------------------------------------------
# Stage 4 / Fitch — input integrity
# --------------------------------------------------------------------------


def test_nan_genotype_is_rejected_not_treated_as_a_third_state():
    with pytest.raises(CladeInputError, match="non-binary|not 0 or 1"):
        fitch_reconstruct(_tree(), {"A": 1, "B": 2, "C": 0, "D": 0})


def test_missing_leaf_is_unknown_not_absent():
    """None/NaN is genuinely unknown and must not be silently scored as 0."""
    t = fitch_reconstruct(_tree(), {"A": 1, "B": None, "C": 0, "D": 0})
    leaves = {n.taxon.label: n.fitch_set for n in t.leaf_node_iter()}
    assert leaves["B"] == {0, 1}


def test_tree_sample_mismatch_raises_instead_of_returning_zeros():
    """Label-convention mismatch used to yield a clean table of zero gains."""
    with pytest.raises(CladeInputError, match="leaves|present"):
        temporal_ordering(_tree("((X1,X2),(X3,X4));"), PHENO, {"c": PHENO})


def test_ambiguity_resolution_is_reported_as_a_sensitivity():
    """A candidate whose verdict depends on the tie-breaking policy must be
    identifiable rather than reported at face value."""
    out = temporal_ordering_sensitivity(
        _tree(), PHENO, {"c": {"A": 1, "B": 0, "C": 1, "D": 0}}, baseline=0.5
    )
    assert "verdict_stable" in out.columns
    assert len(out) == 1


# --------------------------------------------------------------------------
# Stage 2 — recurrence
# --------------------------------------------------------------------------


def test_single_lineage_candidate_is_flagged():
    g = pd.DataFrame({"c": [1, 1, 1, 0]}, index=list("abcd"))
    lin = pd.Series(["ST1", "ST1", "ST1", "ST2"], index=list("abcd"))
    r = lineage_recurrence(g, lin).iloc[0]
    assert r["Distinct_STs"] == 1
    assert r["Pct_dominant_ST"] == 1.0


def test_multi_lineage_candidate_is_distinguished():
    g = pd.DataFrame({"c": [1, 1, 1, 0]}, index=list("abcd"))
    lin = pd.Series(["ST1", "ST2", "ST3", "ST4"], index=list("abcd"))
    r = lineage_recurrence(g, lin).iloc[0]
    assert r["Distinct_STs"] == 3


def test_carriers_without_lineage_are_counted_not_dropped_silently():
    g = pd.DataFrame({"c": [1, 1, 0]}, index=["s1", "s2", "s3"])
    lin = pd.Series(["ST1"], index=["s1"])
    r = lineage_recurrence(g, lin).iloc[0]
    assert r["Carriers"] == 2
    assert r["Carriers_without_ST"] == 1


def test_non_binary_genotype_is_rejected():
    g = pd.DataFrame({"c": [1, 2, 0]}, index=list("abc"))
    lin = pd.Series(["ST1", "ST1", "ST2"], index=list("abc"))
    with pytest.raises(CladeInputError, match="outside"):
        lineage_recurrence(g, lin)


# --------------------------------------------------------------------------
# Stage 3 — direction
# --------------------------------------------------------------------------


def test_direction_signs():
    assert direction_check(2.5)["correct_direction"] is True
    assert direction_check(-1.59)["correct_direction"] is False


def test_zero_and_missing_coefficients_are_not_negative_evidence():
    assert direction_check(0.0)["correct_direction"] is None
    assert direction_check(None)["correct_direction"] is None
    assert direction_check(float("nan"))["correct_direction"] is None


# --------------------------------------------------------------------------
# Stage 5 — matched neighbours
# --------------------------------------------------------------------------


def _dist(samples, pairs):
    d = pd.DataFrame(10.0, index=samples, columns=samples, dtype=float)
    for s in samples:
        d.loc[s, s] = 0.0
    for (a, b), v in pairs.items():
        d.loc[a, b] = d.loc[b, a] = v
    return d


def test_tie_breaking_is_independent_of_input_order():
    samples = ["p1", "n1", "n2"]
    d = _dist(samples, {("p1", "n1"): 0.5, ("p1", "n2"): 0.5})
    a = nearest_negative_neighbor(d, ["p1"], ["n1", "n2"])
    b = nearest_negative_neighbor(d, ["p1"], ["n2", "n1"])
    assert a.loc["p1", "nearest_neighbor"] == b.loc["p1", "nearest_neighbor"]
    assert bool(a.loc["p1", "tied"]) is True


def test_neighbor_reuse_is_reported():
    samples = ["p1", "p2", "p3", "n1"]
    d = _dist(samples, {})
    nn = nearest_negative_neighbor(d, ["p1", "p2", "p3"], ["n1"])
    g = pd.DataFrame({"c": [1, 1, 1, 0]}, index=samples)
    r = matched_neighbor_test(g, ["p1", "p2", "p3"], nn, ["c"]).iloc[0]
    assert r["Max_neighbor_reuse"] == 3
    assert r["N_unique_neighbors"] == 1


def test_unique_matching_removes_reuse():
    samples = ["p1", "p2", "n1", "n2"]
    d = _dist(samples, {("p1", "n1"): 0.1, ("p2", "n2"): 0.2})
    nn = nearest_negative_neighbor(d, ["p1", "p2"], ["n1", "n2"], unique=True)
    assert nn["nearest_neighbor"].nunique() == 2


def test_missing_genotype_excludes_the_pair_rather_than_scoring_absence():
    samples = ["p1", "p2", "n1", "n2"]
    d = _dist(samples, {("p1", "n1"): 0.1, ("p2", "n2"): 0.2})
    nn = nearest_negative_neighbor(d, ["p1", "p2"], ["n1", "n2"], unique=True)
    g = pd.DataFrame({"c": [1.0, np.nan, 0.0, 0.0]}, index=samples)
    r = matched_neighbor_test(g, ["p1", "p2"], nn, ["c"]).iloc[0]
    assert r["N_pairs_dropped_missing"] == 1
    assert r["N_pairs"] == 1


def test_direction_of_stage5_effect_is_reported():
    samples = ["p1", "p2", "n1", "n2"]
    d = _dist(samples, {("p1", "n1"): 0.1, ("p2", "n2"): 0.2})
    nn = nearest_negative_neighbor(d, ["p1", "p2"], ["n1", "n2"], unique=True)
    g = pd.DataFrame({"c": [1, 1, 0, 0]}, index=samples)
    r = matched_neighbor_test(g, ["p1", "p2"], nn, ["c"]).iloc[0]
    # pandas returns np.bool_ on row access; compare by value.
    assert bool(r["Enriched_in_resistant"]) is True

    g2 = pd.DataFrame({"c": [0, 0, 1, 1]}, index=samples)
    r2 = matched_neighbor_test(g2, ["p1", "p2"], nn, ["c"]).iloc[0]
    assert bool(r2["Enriched_in_resistant"]) is False


def test_no_discordant_pairs_gives_nan_not_a_significant_p():
    samples = ["p1", "n1"]
    d = _dist(samples, {("p1", "n1"): 0.1})
    nn = nearest_negative_neighbor(d, ["p1"], ["n1"])
    g = pd.DataFrame({"c": [1, 1]}, index=samples)
    r = matched_neighbor_test(g, ["p1"], nn, ["c"]).iloc[0]
    assert np.isnan(r["McNemar_p"])


def test_stage5_without_any_negative_samples_raises():
    d = _dist(["p1", "p2"], {})
    with pytest.raises(CladeInputError, match="negative"):
        nearest_negative_neighbor(d, ["p1", "p2"], [])


# --------------------------------------------------------------------------
# Verdict engine — missing evidence must never become support
# --------------------------------------------------------------------------


def _full_pass(**over):
    base = {
        "stage1_tested": True,
        "stage1_significant": True,
        "stage3_correct_direction": True,
        "stage4_pct_post_resistance": 0.95,
        "stage4_baseline": 0.50,
        "stage4_total_gains": 5,
        "stage5_significant": True,
        "stage5_enriched_in_resistant": True,
    }
    base.update(over)
    return StageResults(**base)


def test_all_stages_positive_is_convergent():
    assert classify_candidate(_full_pass()).disposition is Disposition.CONVERGENT


def test_stage1_only_is_not_convergent():
    """The central regression: one passing test is not six-stage convergence."""
    v = classify_candidate(StageResults(stage1_tested=True, stage1_significant=True))
    assert v.disposition is Disposition.INSUFFICIENT_EVIDENCE
    assert set(v.stages_missing) == {"stage3", "stage4", "stage5"}


@pytest.mark.parametrize("drop", ["stage3", "stage4", "stage5"])
def test_any_single_missing_stage_blocks_convergence(drop):
    over = {
        "stage3": {"stage3_correct_direction": None},
        "stage4": {"stage4_pct_post_resistance": None, "stage4_baseline": None},
        "stage5": {"stage5_significant": None, "stage5_enriched_in_resistant": None},
    }[drop]
    v = classify_candidate(_full_pass(**over))
    assert v.disposition is Disposition.INSUFFICIENT_EVIDENCE
    assert drop in v.stages_missing


def test_indeterminate_stage1_is_not_convergent():
    v = classify_candidate(StageResults(stage1_tested=True, stage1_significant=None))
    assert v.disposition is Disposition.INSUFFICIENT_EVIDENCE


def test_untested_is_distinct_from_rejected():
    assert classify_candidate(StageResults(stage1_tested=False)).disposition is Disposition.NOT_RETESTED
    assert (
        classify_candidate(StageResults(stage1_tested=True, stage1_significant=False)).disposition
        is Disposition.REJECTED
    )


def test_zero_gain_stage4_is_missing_not_failing():
    """A candidate with no reconstructed origin has an undefined percentage;
    that must not be scored as a temporal failure."""
    r = _full_pass(stage4_pct_post_resistance=float("nan"), stage4_total_gains=0)
    v = classify_candidate(r)
    assert v.stage_status["stage4"] is None
    assert v.disposition is Disposition.INSUFFICIENT_EVIDENCE


def test_dnaa_profile_is_unresolved():
    """Local matched-pair evidence supports; direction and timing contradict."""
    v = classify_candidate(
        _full_pass(
            stage3_correct_direction=False,
            stage4_pct_post_resistance=0.23,
            stage4_baseline=0.767,
            stage4_total_gains=30,
            stage5_significant=True,
            stage5_enriched_in_resistant=True,
        )
    )
    assert v.disposition is Disposition.UNRESOLVED


def test_wrong_direction_alone_is_rejected():
    v = classify_candidate(_full_pass(stage3_correct_direction=False,
                                      stage5_significant=False,
                                      stage5_enriched_in_resistant=False))
    assert v.disposition is Disposition.REJECTED


def test_near_fixed_is_uninformative():
    v = classify_candidate(_full_pass(cohort_frequency=0.97))
    assert v.disposition is Disposition.UNINFORMATIVE


def test_single_lineage_with_contradiction_is_clonal_artifact():
    v = classify_candidate(
        _full_pass(
            stage2_pct_dominant_st=1.0,
            stage2_independent_origins=1,
            stage3_correct_direction=False,
            stage5_significant=False,
            stage5_enriched_in_resistant=False,
        )
    )
    assert v.disposition is Disposition.CLONAL_ARTIFACT


def test_seca_profile_single_lineage_but_recurrent_stays_convergent():
    """SecA is 100% ST2 yet convergent, because it arose independently five
    times within that lineage. Stage 2 must inform, not gate."""
    v = classify_candidate(
        _full_pass(stage2_pct_dominant_st=1.0, stage2_independent_origins=5)
    )
    assert v.disposition is Disposition.CONVERGENT


def test_stage5_significant_but_wrong_direction_does_not_support():
    r = _full_pass(stage5_significant=True, stage5_enriched_in_resistant=False)
    assert r.stage5_supports() is False
    assert classify_candidate(r).disposition is Disposition.REJECTED


def test_verdict_carries_an_audit_trail():
    v = classify_candidate(_full_pass())
    assert v.reasons
    assert set(v.stage_status) == {"stage1", "stage3", "stage4", "stage5"}

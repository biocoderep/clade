from clade.classification.disposition import Disposition, StageResults, classify_candidate


def test_convergent_when_all_applicable_stages_agree():
    r = StageResults(
        stage1_tested=True, stage1_significant=True,
        stage2_distinct_sts=5, stage2_pct_dominant_st=0.4,
        stage3_correct_direction=True,
        stage4_pct_post_resistance=1.0, stage4_baseline=0.767,
        stage5_significant=True, stage5_concordant=True,
        cohort_frequency=0.09,
    )
    assert classify_candidate(r) == Disposition.CONVERGENT


def test_convergent_even_when_single_lineage_if_everything_else_passes():
    """Mirrors SecA in the CLADE case study: 100% one lineage, but recurs
    independently within it and passes every other stage."""
    r = StageResults(
        stage1_tested=True, stage1_significant=True,
        stage2_distinct_sts=1, stage2_pct_dominant_st=1.0,
        stage3_correct_direction=True,
        stage4_pct_post_resistance=1.0, stage4_baseline=0.767,
        stage5_significant=True, stage5_concordant=True,
        cohort_frequency=0.09,
    )
    assert classify_candidate(r) == Disposition.CONVERGENT


def test_unresolved_when_stage5_contradicts_stage3_and_4():
    """Mirrors DnaA: wrong direction + mostly pre-resistance timing, but
    Stage 5 shows a significant positive local signal."""
    r = StageResults(
        stage1_tested=True, stage1_significant=True,
        stage2_distinct_sts=47, stage2_pct_dominant_st=0.3,
        stage3_correct_direction=False,
        stage4_pct_post_resistance=0.233, stage4_baseline=0.767,
        stage5_significant=True, stage5_concordant=False,
        cohort_frequency=0.05,
    )
    assert classify_candidate(r) == Disposition.UNRESOLVED


def test_rejected_when_stage1_not_significant():
    r = StageResults(stage1_tested=True, stage1_significant=False, cohort_frequency=0.1)
    assert classify_candidate(r) == Disposition.REJECTED


def test_rejected_when_wrong_direction_and_no_stage5_contradiction():
    r = StageResults(
        stage1_tested=True, stage1_significant=True,
        stage2_distinct_sts=48, stage2_pct_dominant_st=0.17,
        stage3_correct_direction=False,
        stage4_pct_post_resistance=0.143, stage4_baseline=0.767,
        stage5_significant=False, stage5_concordant=None,
        cohort_frequency=0.06,
    )
    assert classify_candidate(r) == Disposition.REJECTED


def test_clonal_artifact_when_single_lineage_and_fails():
    r = StageResults(
        stage1_tested=True, stage1_significant=True,
        stage2_distinct_sts=1, stage2_pct_dominant_st=0.986,
        stage3_correct_direction=False,
        cohort_frequency=0.02,
    )
    assert classify_candidate(r) == Disposition.CLONAL_ARTIFACT


def test_uninformative_when_near_fixed():
    r = StageResults(stage1_tested=True, stage1_significant=True, cohort_frequency=0.95)
    assert classify_candidate(r) == Disposition.UNINFORMATIVE


def test_not_retested_when_stage1_never_ran():
    r = StageResults(stage1_tested=False, cohort_frequency=0.05)
    assert classify_candidate(r) == Disposition.NOT_RETESTED

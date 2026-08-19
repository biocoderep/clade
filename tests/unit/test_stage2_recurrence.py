import pandas as pd
import pytest

from clade.validation.stage2_recurrence import lineage_recurrence


def test_single_lineage_candidate_flagged():
    genotypes = pd.DataFrame(
        {"geneA": [1, 1, 1, 0, 0]},
        index=["s1", "s2", "s3", "s4", "s5"],
    )
    lineages = pd.Series(["ST2", "ST2", "ST2", "ST9", "ST9"], index=genotypes.index)

    result = lineage_recurrence(genotypes, lineages)
    row = result.set_index("Candidate").loc["geneA"]

    assert row["Carriers"] == 3
    assert row["Distinct_STs"] == 1
    assert row["Dominant_ST"] == "ST2"
    assert row["Pct_dominant_ST"] == 1.0


def test_multi_lineage_candidate_not_flagged():
    genotypes = pd.DataFrame(
        {"geneB": [1, 1, 0, 1, 0]},
        index=["s1", "s2", "s3", "s4", "s5"],
    )
    lineages = pd.Series(["ST2", "ST9", "ST9", "ST41", "ST41"], index=genotypes.index)

    result = lineage_recurrence(genotypes, lineages)
    row = result.set_index("Candidate").loc["geneB"]

    assert row["Carriers"] == 3
    assert row["Distinct_STs"] == 3
    assert row["Pct_dominant_ST"] == pytest.approx(1 / 3, abs=0.01)


def test_zero_carriers_handled():
    genotypes = pd.DataFrame({"geneC": [0, 0, 0]}, index=["s1", "s2", "s3"])
    lineages = pd.Series(["ST2", "ST2", "ST2"], index=genotypes.index)

    result = lineage_recurrence(genotypes, lineages)
    row = result.set_index("Candidate").loc["geneC"]

    assert row["Carriers"] == 0
    assert row["Distinct_STs"] == 0
    assert row["Dominant_ST"] is None

import pandas as pd

from clade.validation.stage5_matched_neighbors import matched_neighbor_test, nearest_negative_neighbor


def test_nearest_negative_neighbor_picks_closest():
    distances = pd.DataFrame(
        {"n1": [1.0, 5.0], "n2": [5.0, 1.0]},
        index=["p1", "p2"],
    )
    result = nearest_negative_neighbor(distances, ["p1", "p2"], ["n1", "n2"])

    assert result.loc["p1", "nearest_neighbor"] == "n1"
    assert result.loc["p2", "nearest_neighbor"] == "n2"
    assert result.loc["p1", "distance"] == 1.0


def test_matched_neighbor_carriage_counts():
    genotypes = pd.DataFrame(
        {"geneX": [1, 1, 1, 0]},
        index=["p1", "p2", "n1", "n2"],
    )
    matches = pd.DataFrame(
        {"nearest_neighbor": ["n1", "n2"]},
        index=["p1", "p2"],
    )

    result = matched_neighbor_test(genotypes, ["p1", "p2"], matches, ["geneX"])
    row = result.set_index("Candidate").loc["geneX"]

    # p1 (carrier) matched to n1 (carrier) -> concordant
    # p2 (carrier) matched to n2 (non-carrier) -> discordant, resistant-only
    assert row["Resistant_carries"] == 2
    assert row["Neighbor_carries"] == 1
    assert row["Resistant_only"] == 1
    assert row["Neighbor_only"] == 0
    assert row["McNemar_p"] is not None
    assert 0.0 <= row["McNemar_p"] <= 1.0

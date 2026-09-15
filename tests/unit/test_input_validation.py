"""
Input-validation tests.

Each case here previously passed silently and produced a plausible-looking
number from corrupt input. The requirement is not that CLADE tolerates bad
data — it is that CLADE refuses it, loudly, with a message naming the problem.
"""
from __future__ import annotations

import dendropy
import numpy as np
import pandas as pd
import pytest

from clade.io.validation import (
    CladeInputError,
    require_sample_overlap,
    require_unique_index,
    validate_binary,
    validate_distance_matrix,
    validate_genotype_frame,
    validate_phenotype,
    validate_tree_against_samples,
)


def test_duplicate_sample_ids_are_rejected():
    """Duplicates misalign every join and mispair matched neighbours."""
    df = pd.DataFrame({"c": [1, 0, 1]}, index=["a", "a", "b"])
    with pytest.raises(CladeInputError, match="duplicated sample IDs"):
        require_unique_index(df, "genotypes")


def test_non_binary_genotype_is_rejected():
    with pytest.raises(CladeInputError, match="outside"):
        validate_binary(pd.Series([0, 1, 2]), "genotype")


def test_non_numeric_genotype_is_rejected():
    with pytest.raises(CladeInputError, match="non-numeric"):
        validate_binary(pd.Series([0, 1, "present"]), "genotype")


def test_missing_genotype_is_allowed_and_preserved_as_missing():
    out = validate_binary(pd.Series([0, 1, np.nan]), "genotype")
    assert out.isna().sum() == 1
    assert out.notna().sum() == 2


def test_float_and_bool_codings_are_accepted():
    assert validate_binary(pd.Series([0.0, 1.0]), "g").tolist() == [0.0, 1.0]
    assert validate_binary(pd.Series([True, False]), "g").tolist() == [1.0, 0.0]


def test_empty_dataset_is_rejected():
    with pytest.raises(CladeInputError, match="empty"):
        validate_genotype_frame(pd.DataFrame())


def test_candidate_absent_from_genotypes_is_named_in_the_error():
    g = pd.DataFrame({"present": [0, 1]}, index=["a", "b"])
    with pytest.raises(CladeInputError, match="ghost"):
        validate_genotype_frame(g, ["ghost"])


def test_phenotype_with_one_class_is_rejected():
    """Association testing is undefined when every sample shares a phenotype."""
    with pytest.raises(CladeInputError, match="two distinct values"):
        validate_phenotype(pd.Series([1, 1, 1], index=list("abc")))


def test_phenotype_missing_values_are_rejected_by_default():
    """Coercing an unknown phenotype to 0 would inflate the susceptible group."""
    with pytest.raises(CladeInputError, match="missing"):
        validate_phenotype(pd.Series([1, 0, np.nan], index=list("abc")))


def test_no_sample_overlap_is_rejected_with_examples():
    with pytest.raises(CladeInputError, match="share 0 sample"):
        require_sample_overlap(["a", "b"], ["x", "y"], "phenotype", "genotypes")


def test_sample_overlap_returns_shared_ids():
    assert require_sample_overlap(["a", "b", "c"], ["b", "c", "d"], "l", "r") == ["b", "c"]


# ---- tree ----------------------------------------------------------------


def _tree(newick="((A,B),(C,D));"):
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)


def test_tree_with_no_matching_labels_is_rejected():
    with pytest.raises(CladeInputError, match="are present"):
        validate_tree_against_samples(_tree(), ["S1", "S2", "S3", "S4"])


def test_tree_with_partial_overlap_below_threshold_is_rejected():
    with pytest.raises(CladeInputError, match="below the required"):
        validate_tree_against_samples(_tree(), ["A", "S2", "S3", "S4"])


def test_tree_with_sufficient_overlap_reports_the_difference():
    info = validate_tree_against_samples(_tree(), ["A", "B", "C", "X"])
    assert info["n_shared"] == 3
    assert info["samples_not_in_tree"] == ["X"]
    assert info["leaves_not_in_samples"] == ["D"]


def test_unlabelled_tree_is_rejected():
    with pytest.raises(CladeInputError, match="no labelled leaves"):
        validate_tree_against_samples(_tree("((,),(,));"), ["A"])


# ---- distance matrix -----------------------------------------------------


def test_non_square_distance_matrix_is_rejected():
    d = pd.DataFrame([[0.0, 1.0], [1.0, 0.0]], index=["a", "b"], columns=["a", "c"])
    with pytest.raises(CladeInputError, match="not square"):
        validate_distance_matrix(d)


def test_duplicate_labels_in_distance_matrix_are_rejected():
    d = pd.DataFrame(0.0, index=["a", "a"], columns=["a", "a"])
    with pytest.raises(CladeInputError, match="duplicated"):
        validate_distance_matrix(d)


def test_column_order_is_normalised_to_row_order():
    d = pd.DataFrame(
        [[0.0, 1.0], [1.0, 0.0]], index=["a", "b"], columns=["b", "a"]
    )
    out = validate_distance_matrix(d)
    assert list(out.columns) == list(out.index)


# ---- loaders -------------------------------------------------------------


def test_duplicate_candidate_columns_across_files_are_rejected(tmp_path):
    from clade.io.genotypes import load_genotypes

    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    a.write_text("Sample_ID,geneX\ns1,1\ns2,0\n")
    b.write_text("Sample_ID,geneX\ns1,0\ns2,1\n")
    with pytest.raises(CladeInputError, match="redefines"):
        load_genotypes(str(a), str(b))


def test_genotype_file_without_index_column_is_rejected(tmp_path):
    from clade.io.genotypes import load_genotypes

    p = tmp_path / "g.csv"
    p.write_text("wrong_id,geneX\ns1,1\n")
    with pytest.raises(CladeInputError, match="Sample_ID"):
        load_genotypes(str(p))


def test_phenotype_with_multiple_columns_is_rejected(tmp_path):
    from clade.io.genotypes import load_phenotype

    p = tmp_path / "p.tsv"
    p.write_text("sample\ta\tb\ns1\t1\t0\n")
    with pytest.raises(CladeInputError, match="single phenotype column"):
        load_phenotype(str(p))


def test_lineage_file_missing_required_column_is_rejected(tmp_path):
    from clade.io.genotypes import load_lineages

    p = tmp_path / "l.tsv"
    p.write_text("sample\twrong\ns1\tST1\n")
    with pytest.raises(CladeInputError, match="no 'ST' column"):
        load_lineages(str(p))

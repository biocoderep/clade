import os

import pytest

from clade.provenance.phenotype import derive_phenotype_from_amr_matrix, verify_phenotype_provenance

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "..", "fixtures")


@pytest.fixture
def tiny_rtab(tmp_path):
    path = tmp_path / "feature_matrix.Rtab"
    path.write_text(
        "Gene\ts1\ts2\ts3\ts4\n"
        "amr_blaOXA-23\t1\t1\t0\t1\n"
        "amr_blaOXA-66\t0\t0\t1\t0\n"
    )
    return str(path)


def test_derive_phenotype_reads_correct_row(tiny_rtab):
    pheno = derive_phenotype_from_amr_matrix(tiny_rtab, "amr_blaOXA-23")
    assert pheno.to_dict() == {"s1": 1, "s2": 1, "s3": 0, "s4": 1}


def test_derive_phenotype_missing_row_raises(tiny_rtab):
    with pytest.raises(KeyError):
        derive_phenotype_from_amr_matrix(tiny_rtab, "amr_doesnotexist")


def test_verify_provenance_detects_identity_mismatch(tiny_rtab):
    """Reproduces the class of error this module exists to catch: a
    'claimed' phenotype that actually encodes a different gene entirely."""
    claimed = derive_phenotype_from_amr_matrix(tiny_rtab, "amr_blaOXA-66")  # wrong gene, on purpose
    derived = derive_phenotype_from_amr_matrix(tiny_rtab, "amr_blaOXA-23")  # correct source

    report = verify_phenotype_provenance(claimed, derived)
    assert report["n_compared"] == 4
    assert report["n_mismatched"] > 0
    assert report["mismatch_rate"] > 0.5


def test_verify_provenance_passes_when_source_matches(tiny_rtab):
    derived = derive_phenotype_from_amr_matrix(tiny_rtab, "amr_blaOXA-23")
    report = verify_phenotype_provenance(derived, derived)
    assert report["n_mismatched"] == 0
    assert report["mismatch_rate"] == 0.0

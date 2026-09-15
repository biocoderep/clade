"""
Tests for the Nextflow pipeline's bridging scripts.

These sit outside the installable package (they are Nextflow `bin/` scripts),
so they are loaded by path. They are still worth testing: `build_matrix.py`
produces the genotype matrix every CLADE stage then consumes, so an encoding
mistake here propagates through the entire framework.
"""
from __future__ import annotations

import importlib.util
import pathlib

import pandas as pd
import pytest

BIN = pathlib.Path(__file__).resolve().parents[2] / "workflows" / "nextflow" / "bin"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


build_matrix = pytest.importorskip("pandas") and _load("build_matrix")


VCF = """\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t1\t0\t.\t./.
chr1\t200\t.\tG\tC\t.\tPASS\t.\tGT\t0/0\t1/1\t0\t1
"""


def _write(tmp_path, text=VCF):
    p = tmp_path / "core.vcf"
    p.write_text(text)
    return str(p)


def test_alt_call_is_one_and_ref_call_is_zero(tmp_path):
    df = build_matrix.extract_genotypes(
        _write(tmp_path), [("candA", "chr1", 100)]
    ).set_index("Sample_ID")
    assert df.loc["s1", "candA"] == 1.0
    assert df.loc["s2", "candA"] == 0.0


def test_no_call_is_missing_not_absent(tmp_path):
    """`.` and `./.` mean the position could not be genotyped. Recording them
    as 0 would assert the sample does not carry the mutation."""
    df = build_matrix.extract_genotypes(
        _write(tmp_path), [("candA", "chr1", 100)]
    ).set_index("Sample_ID")
    assert pd.isna(df.loc["s3", "candA"])
    assert pd.isna(df.loc["s4", "candA"])


def test_phased_genotypes_are_understood(tmp_path):
    vcf = VCF.replace("\t0/0\t1/1\t0\t1", "\t0|0\t1|1\t.|.\t1")
    df = build_matrix.extract_genotypes(
        _write(tmp_path, vcf), [("candB", "chr1", 200)]
    ).set_index("Sample_ID")
    assert df.loc["s1", "candB"] == 0.0
    assert df.loc["s2", "candB"] == 1.0
    assert pd.isna(df.loc["s3", "candB"])


def test_position_absent_from_vcf_is_all_missing_not_all_zero(tmp_path):
    """A candidate position that never appears in the VCF tells us nothing
    about any sample — it does not tell us nobody carries it."""
    df = build_matrix.extract_genotypes(
        _write(tmp_path), [("ghost", "chr1", 999)]
    ).set_index("Sample_ID")
    assert df["ghost"].isna().all()


def test_multiple_candidates_are_extracted_independently(tmp_path):
    df = build_matrix.extract_genotypes(
        _write(tmp_path), [("candA", "chr1", 100), ("candB", "chr1", 200)]
    ).set_index("Sample_ID")
    assert list(df.columns) == ["candA", "candB"]
    assert df.loc["s2", "candA"] == 0.0
    assert df.loc["s2", "candB"] == 1.0


def test_output_is_consumable_by_clade_validation(tmp_path):
    """The matrix this produces must pass CLADE's own input validation."""
    from clade.io.validation import validate_genotype_frame

    df = build_matrix.extract_genotypes(
        _write(tmp_path), [("candA", "chr1", 100), ("candB", "chr1", 200)]
    ).set_index("Sample_ID")
    out = validate_genotype_frame(df)
    assert set(out.columns) == {"candA", "candB"}

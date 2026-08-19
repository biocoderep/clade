"""
End-to-end smoke test for `clade validate`, using only the tiny synthetic
fixtures in tests/fixtures/ — never the real 3,261-genome dataset (CI must
not depend on it; see docs/reproducibility/open_items.md).

This dataset is deliberately too small to pass the default Fisher/FDR
prevalence filter (min_carriers=10), so it exercises the "nothing reaches
Stage 1 significance -> NOT_RETESTED" path rather than a convergent one.
That is itself a real, useful thing to test: the pipeline must not crash
when nothing survives Stage 1, and every candidate must still appear in
the output table (CLADE's reporting rule).
"""
import os
import subprocess
import sys

FIXTURES = os.path.join(os.path.dirname(__file__), "..", "fixtures")


def test_version():
    result = subprocess.run([sys.executable, "-m", "clade.cli", "--version"], capture_output=True, text=True, check=False)
    assert result.returncode == 0
    assert "clade" in result.stdout.lower()


def test_help():
    result = subprocess.run([sys.executable, "-m", "clade.cli", "--help"], capture_output=True, text=True, check=False)
    assert result.returncode == 0
    assert "validate" in result.stdout


def test_validate_end_to_end_on_tiny_fixtures(tmp_path):
    output = tmp_path / "evidence_table.md"
    cmd = [
        sys.executable, "-m", "clade.cli", "validate",
        "--genotypes", os.path.join(FIXTURES, "example_genotypes.csv"),
        "--phenotype", os.path.join(FIXTURES, "example_phenotype.tsv"),
        "--lineages", os.path.join(FIXTURES, "example_lineages.tsv"),
        "--candidates", "geneA", "geneB",
        "--tree", os.path.join(FIXTURES, "tiny_tree.nwk"),
        "--output", str(output),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)

    assert result.returncode == 0, result.stderr
    assert output.exists()

    content = output.read_text()
    assert "geneA" in content
    assert "geneB" in content
    assert "not_individually_retested" in content  # expected disposition given min_carriers filter

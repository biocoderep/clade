"""
End-to-end tests for `clade validate`, using only the synthetic fixtures in
tests/fixtures/ — never the real case-study dataset (CI must not depend on it;
see docs/reproducibility/open_items.md).

The substantive test here is `test_planted_truths_are_recovered`: the synthetic
cohort has four candidates whose correct dispositions are known by
construction, and the pipeline must recover all four. That is a scientific
regression test, not a smoke test — if a refactor flips one of those verdicts,
the framework is broken regardless of whether it still runs.

Note for anyone reading the history: an earlier version of this file asserted
that both candidates in the tiny fixture came back
`not_individually_retested`, and rationalised it as "exercising the nothing
survives Stage 1 path". It was in fact locking in a defect — a candidate
excluded by a hard-coded prevalence filter was being reported with CLADE's
phrase for *a coverage limitation*, which is not what had happened.
"""
import os
import subprocess
import sys

import pandas as pd
import pytest

FIXTURES = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "fixtures"))
COHORT = os.path.join(FIXTURES, "synthetic_cohort")


def _run(args, **kw):
    return subprocess.run(
        [sys.executable, "-m", "clade.cli", *args], capture_output=True, text=True, check=False, **kw
    )


def test_version():
    r = _run(["--version"])
    assert r.returncode == 0
    assert "clade" in r.stdout.lower()


def test_help_lists_every_stage_input():
    r = _run(["validate", "--help"])
    assert r.returncode == 0
    for flag in ("--genotypes", "--phenotype", "--lineages", "--tree", "--distances", "--stage6"):
        assert flag in r.stdout


def test_no_subcommand_prints_help_and_fails():
    r = _run([])
    assert r.returncode == 1
    assert "validate" in r.stdout


def _validate(tmp_path, extra=(), candidates=None):
    out = tmp_path / "evidence_table.md"
    tsv = tmp_path / "evidence_table.tsv"
    cands = candidates or [
        "true_compensatory", "clonal_marker", "wrong_direction", "near_fixed"
    ]
    r = _run(
        [
            "validate",
            "--genotypes", os.path.join(COHORT, "genotypes.csv"),
            "--phenotype", os.path.join(COHORT, "phenotype.tsv"),
            "--lineages", os.path.join(COHORT, "lineages.tsv"),
            "--candidates", *cands,
            "--output", str(out),
            "--output-tsv", str(tsv),
            *extra,
        ]
    )
    return r, out, tsv


FULL = (
    "--tree", os.path.join(COHORT, "tree.nwk"),
    "--distances", os.path.join(COHORT, "distances.tsv"),
    "--stage1-results", os.path.join(COHORT, "stage1_results.tsv"),
    "--stage6", os.path.join(COHORT, "stage6.json"),
)


def test_planted_truths_are_recovered(tmp_path):
    """The scientific regression test. Each candidate's correct disposition is
    fixed by how the fixture was constructed."""
    r, _out, tsv = _validate(tmp_path, FULL)
    assert r.returncode == 0, r.stderr
    table = pd.read_csv(tsv, sep="\t").set_index("Candidate")

    # Planted: independent origins in several lineages, all post-resistance,
    # enriched in resistant genomes, with recorded Stage 6 support.
    assert table.loc["true_compensatory", "Disposition"] == "convergent"
    assert table.loc["true_compensatory", "Stage2_distinct_STs"] > 1
    assert table.loc["true_compensatory", "Stage4_total_gains"] > 1

    # Planted: one ancestral origin inherited by an entire lineage. Its
    # apparent association is absorbed by the lineage covariate.
    assert table.loc["clonal_marker", "Disposition"] != "convergent"
    assert table.loc["clonal_marker", "Stage2_distinct_STs"] == 1
    assert table.loc["clonal_marker", "Stage4_total_gains"] == 1

    # Planted: depleted among resistant genomes.
    assert table.loc["wrong_direction", "Disposition"] == "rejected"
    assert table.loc["wrong_direction", "Stage3_correct_direction"] is False or (
        table.loc["wrong_direction", "Stage3_status"] == "FAIL"
    )

    # Planted: present in almost every genome.
    assert table.loc["near_fixed", "Disposition"] == "uninformative"


def test_every_candidate_appears_whatever_the_verdict(tmp_path):
    """CLADE's reporting rule: a table listing only survivors cannot be audited."""
    _r, out, tsv = _validate(tmp_path, FULL)
    content = out.read_text()
    for cand in ("true_compensatory", "clonal_marker", "wrong_direction", "near_fixed"):
        assert cand in content
    assert len(pd.read_csv(tsv, sep="\t")) == 4


def test_missing_stages_do_not_produce_convergence(tmp_path):
    """Without a tree or distances, Stages 4 and 5 cannot run — and no
    candidate may be promoted to convergent on the strength of their silence."""
    r, _out, tsv = _validate(
        tmp_path, ("--stage1-results", os.path.join(COHORT, "stage1_results.tsv"))
    )
    assert r.returncode == 0, r.stderr
    table = pd.read_csv(tsv, sep="\t")
    assert "convergent" not in set(table["Disposition"])
    row = table.set_index("Candidate").loc["true_compensatory"]
    assert row["Disposition"] == "insufficient_evidence"
    assert "stage4" in row["Stages_missing"] and "stage5" in row["Stages_missing"]


def test_skipped_stage_is_reported_as_not_run_not_as_failure(tmp_path):
    _r, _out, tsv = _validate(
        tmp_path, ("--stage1-results", os.path.join(COHORT, "stage1_results.tsv"))
    )
    table = pd.read_csv(tsv, sep="\t").set_index("Candidate")
    assert table.loc["true_compensatory", "Stage4_status"] == "not run"
    assert table.loc["true_compensatory", "Stage5_status"] == "not run"


def test_unknown_candidate_fails_cleanly(tmp_path):
    r, _out, _tsv = _validate(tmp_path, FULL, candidates=["no_such_gene"])
    assert r.returncode == 2
    assert "not present" in r.stderr


def test_missing_file_fails_cleanly(tmp_path):
    out = tmp_path / "o.md"
    r = _run(
        [
            "validate",
            "--genotypes", os.path.join(COHORT, "does_not_exist.csv"),
            "--phenotype", os.path.join(COHORT, "phenotype.tsv"),
            "--lineages", os.path.join(COHORT, "lineages.tsv"),
            "--candidates", "true_compensatory",
            "--output", str(out),
        ]
    )
    assert r.returncode == 2
    assert "not found" in r.stderr.lower() or "no such file" in r.stderr.lower()


def test_mismatched_tree_labels_fail_cleanly(tmp_path):
    bad = tmp_path / "bad_tree.nwk"
    bad.write_text("((X1:0.1,X2:0.1):0.1,(X3:0.1,X4:0.1):0.1);\n")
    r, _out, _tsv = _validate(
        tmp_path,
        ("--tree", str(bad), "--stage1-results", os.path.join(COHORT, "stage1_results.tsv")),
    )
    assert r.returncode == 2
    assert "leaves" in r.stderr or "tree" in r.stderr.lower()


def test_output_is_deterministic(tmp_path):
    _a, out_a, _ta = _validate(tmp_path / "a", FULL)
    _b, out_b, _tb = _validate(tmp_path / "b", FULL)
    assert out_a.read_text() == out_b.read_text()


def test_unique_neighbors_sensitivity_runs(tmp_path):
    r, _out, tsv = _validate(tmp_path, (*FULL, "--unique-neighbors"))
    assert r.returncode == 0, r.stderr
    assert len(pd.read_csv(tsv, sep="\t")) == 4


@pytest.fixture(autouse=True)
def _tmpdirs(tmp_path):
    (tmp_path / "a").mkdir(exist_ok=True)
    (tmp_path / "b").mkdir(exist_ok=True)
    yield

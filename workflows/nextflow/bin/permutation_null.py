#!/usr/bin/env python3
"""Empirical false-discovery calibration by lineage-preserving permutation.

A numerically correct p-value is not the same thing as a well-calibrated one.
This answers the question that actually matters: if there were no true
candidate-phenotype relationship at all, how often would this pipeline still
call a candidate convergent?

The shuffle is done WITHIN each lineage, not across the cohort. That choice is
the whole point. A naive whole-cohort shuffle destroys population structure
entirely, which makes the null far too easy to pass and understates the real
false-positive risk -- exactly the confounding CLADE exists to handle would be
permuted away. Shuffling within lineage keeps each lineage's resistant /
susceptible mix intact and breaks only the link to the candidate.

Two scopes:

  --scope stage1   Permute and re-run the association test only. Fast, and
                   diagnoses the association estimator's calibration.
  --scope full     Permute and re-run the whole convergence criterion
                   (Stages 1, 3, 4 and 5). Slower, and this is the number that
                   belongs next to a claim of "k convergent candidates",
                   because it is the null for the same statistic.

--scope full is the default. A Stage-1-only null answers a narrower question
than the convergent count it is usually quoted beside, and quoting the narrow
number for the broad claim overstates the calibration evidence.

Outputs
    permutation_raw.tsv        every permutation x candidate result
    permutation_summary.txt    observed vs null, empirical p, verdict
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def stratified_shuffle(y: np.ndarray, strata: np.ndarray,
                       rng: np.random.Generator) -> np.ndarray:
    """Permute y within each stratum, preserving each lineage's case count."""
    out = y.copy()
    for s in np.unique(strata):
        idx = np.where(strata == s)[0]
        out[idx] = rng.permutation(y[idx])
    return out


def count_convergent(genotypes, phenotype, lineages, workdir, tag, cli_opts,
                     candidates_file, tree, distances):
    """Run the framework once, via the same CLI the pipeline runs, and return
    the number of convergent candidates.

    This shells out to `clade validate` rather than reimplementing the stages.
    That is the point: a permutation null is only meaningful if it exercises
    the identical estimator as the observed statistic. A reimplementation that
    drifts -- a different Stage 1 backend, a different disposition rule --
    produces a null for a statistic nobody reported. An earlier version of this
    script did exactly that and silently returned 0 observed convergent for a
    fixture the CLI calls convergent.
    """
    import subprocess

    pheno_path = workdir / f"pheno_{tag}.tsv"
    out_md = workdir / f"evidence_{tag}.md"
    out_tsv = workdir / f"evidence_{tag}.tsv"

    pd.DataFrame({"sample": phenotype.index, "phenotype": phenotype.values}) \
      .to_csv(pheno_path, sep="\t", index=False)

    # sys.executable -m clade.cli, not a bare "clade": the console script is
    # only on PATH if the active environment's bin directory happens to be, and
    # under a pipeline or a venv-qualified interpreter it frequently is not.
    # Invoking through the running interpreter guarantees the same environment
    # that imported this script is the one that runs the validation.
    cmd = [
        sys.executable, "-m", "clade.cli", "validate",
        "--genotypes", str(cli_opts["genotypes_path"]),
        "--phenotype", str(pheno_path),
        "--lineages", str(cli_opts["lineages_path"]),
        "--candidates", *candidates_file,
        "--output", str(out_md),
        "--output-tsv", str(out_tsv),
        "--alpha", str(cli_opts["alpha"]),
        "--min-carriers", str(cli_opts["min_carriers"]),
        "--max-freq", str(cli_opts["max_freq"]),
        "--temporal-margin", str(cli_opts["temporal_margin"]),
        "--ambiguity-resolution", str(cli_opts["ambiguity_resolution"]),
    ]
    if cli_opts["scope"] == "full":
        cmd += ["--tree", str(tree), "--distances", str(distances)]
    if cli_opts.get("stage1_results"):
        cmd += ["--stage1-results", str(cli_opts["stage1_results"])]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise RuntimeError(f"could not launch the CLADE CLI: {exc}") from exc
    if not out_tsv.exists():
        raise RuntimeError(
            "clade validate produced no evidence table"
            f" (exit {proc.returncode}).\nstderr:\n{proc.stderr[-2000:]}"
        )

    table = pd.read_csv(out_tsv, sep="\t")
    col = "Disposition" if "Disposition" in table.columns else table.columns[1]
    dispositions = table[col].astype(str).str.lower()
    n_convergent = int((dispositions == "convergent").sum())
    detail = dict(zip(table.iloc[:, 0].astype(str), dispositions))

    # Whether Stage 1 actually produced a verdict for anything. A run in which
    # Stage 1 never ran still yields a tidy table of zeros, which is
    # indistinguishable from a real null unless this is checked.
    stage1_ran = False
    if "Stage1_status" in table.columns:
        stage1_ran = bool((table["Stage1_status"].astype(str).str.lower() != "not run").any())
    elif "Stage1_significant" in table.columns:
        stage1_ran = bool(table["Stage1_significant"].notna().any())
    return n_convergent, detail, stage1_ran


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--lineages", required=True)
    ap.add_argument("--tree")
    ap.add_argument("--distances")
    ap.add_argument("--n-permutations", type=int, default=100)
    ap.add_argument("--scope", choices=["stage1", "full"], default="full")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-carriers", type=int, default=5)
    ap.add_argument("--max-freq", type=float, default=0.85)
    ap.add_argument("--temporal-margin", type=float, default=0.0)
    ap.add_argument("--ambiguity-resolution", type=int, choices=[0, 1], default=0)
    ap.add_argument("--stage1-results",
                    help="external structure-corrected Stage 1 results; "
                         "note these are FIXED across permutations and so "
                         "cannot form a valid null -- only use when "
                         "recomputing Stage 1 per permutation is intended")
    ap.add_argument("--keep-work", action="store_true",
                    help="keep per-permutation evidence tables for inspection")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    if args.n_permutations < 1:
        sys.exit("--n-permutations must be >= 1")
    if args.scope == "full" and not (args.tree and args.distances):
        sys.exit("--scope full needs --tree and --distances; use --scope stage1 without them")

    genotypes = pd.read_csv(args.genotypes, index_col=0)
    phenotype = pd.read_csv(args.phenotype, sep="\t", index_col=0).iloc[:, 0].astype(int)
    lineages = pd.read_csv(args.lineages, sep="\t", index_col=0).iloc[:, 0].astype(str)

    common = genotypes.index.intersection(phenotype.index).intersection(lineages.index)
    if len(common) == 0:
        sys.exit("no samples common to genotypes, phenotype and lineages")
    genotypes = genotypes.loc[common]
    phenotype = phenotype.loc[common]
    lineages = lineages.loc[common]
    candidates = list(genotypes.columns)

    # The seed is printed, not just used: a calibration figure that cannot be
    # regenerated exactly is not evidence.
    print("Lineage-preserving permutation null")
    print(f"  samples={len(common)}  candidates={len(candidates)}  "
          f"lineages={lineages.nunique()}")
    print(f"  permutations={args.n_permutations}  scope={args.scope}  seed={args.seed}")

    workdir = Path(args.outdir) / "_perm_work"
    workdir.mkdir(parents=True, exist_ok=True)
    cli_opts = {
        "genotypes_path": Path(args.genotypes).resolve(),
        "lineages_path": Path(args.lineages).resolve(),
        "alpha": args.alpha,
        "min_carriers": args.min_carriers,
        "max_freq": args.max_freq,
        "temporal_margin": args.temporal_margin,
        "ambiguity_resolution": args.ambiguity_resolution,
        "scope": args.scope,
        "stage1_results": args.stage1_results,
    }

    observed, _observed_detail, stage1_ran = count_convergent(
        genotypes, phenotype, lineages, workdir, "observed", cli_opts,
        candidates, args.tree, args.distances)
    print(f"  observed convergent (real phenotype): {observed}")

    # Refuse to produce a null that only looks like one.
    #
    # If Stage 1 cannot run -- firthlogist absent, no kinship backend -- every
    # candidate comes back insufficient_evidence, the observed count is 0, and
    # every permutation is 0 too. That yields a clean-looking table, an
    # empirical p of 1.0, and no indication anywhere that the association test
    # never executed. Failing loudly here is the whole point.
    if not stage1_ran:
        sys.exit(
            "ABORT: Stage 1 did not run for any candidate, so there is no\n"
            "statistic to calibrate and any null computed here would be\n"
            "meaningless (every permutation would return zero convergent\n"
            "candidates for the same reason the observed run did).\n\n"
            "Stage 1 needs either:\n"
            "  * firthlogist installed (requires Python <3.11 and\n"
            "    scikit-learn <1.6) so the built-in Firth fit is available; or\n"
            "  * a kinship-corrected backend such as pyseer.\n\n"
            "Note that --stage1-results does NOT solve this: externally\n"
            "computed results are fixed across permutations, so Stage 1 would\n"
            "not be recomputed under the shuffled phenotype and the resulting\n"
            "'null' would not be a null at all."
        )

    # A fixed Stage 1 cannot produce a valid null for the convergent count:
    # the association result would not be recomputed under the shuffled
    # phenotype, so Stage 1 contributes the same verdict every permutation and
    # only Stages 3-5 vary. That is a narrower, conditional null. It is allowed
    # -- it answers a real question -- but it is labelled everywhere it appears
    # so it cannot be quoted as the calibration of the full statistic.
    conditional = bool(args.stage1_results)
    if conditional:
        print()
        print("  WARNING: --stage1-results is fixed across permutations.")
        print("  Stage 1 will NOT be recomputed under the shuffled phenotype, so")
        print("  this is a CONDITIONAL null over Stages 3-5 only, not a")
        print("  calibration of the convergent count. Labelled as such in the")
        print("  summary. For the full null, install firthlogist or use pyseer.")
        print()

    rng = np.random.default_rng(args.seed)
    strata = lineages.values
    y = phenotype.values
    null_counts, records = [], []

    for i in range(args.n_permutations):
        shuffled = pd.Series(stratified_shuffle(y, strata, rng), index=phenotype.index)
        k, detail, _ = count_convergent(
            genotypes, shuffled, lineages, workdir, f"perm{i}", cli_opts,
            candidates, args.tree, args.distances)
        null_counts.append(k)
        for name, disp in detail.items():
            records.append({"permutation": i, "Candidate": name, "Disposition": disp})
        if (i + 1) % 10 == 0 or i == 0:
            print(f"    permutation {i+1}/{args.n_permutations}: {k} convergent")

    pd.DataFrame(records).to_csv(f"{args.outdir}/permutation_raw.tsv",
                                 sep="\t", index=False)

    null = np.array(null_counts)
    null_mean = float(null.mean())
    # Add-one empirical p: the observed run is itself one draw, so a p of
    # exactly zero is never claimable from a finite number of permutations.
    emp_p = (int((null >= observed).sum()) + 1) / (args.n_permutations + 1)

    lines = [
        "Empirical false-discovery calibration",
        "=" * 64,
        f"Scope                      : {args.scope}"
        + ("  (Stages 1+3+4+5)" if args.scope == "full" else "  (association only)")
        + ("   [CONDITIONAL: Stage 1 fixed, not recomputed]" if conditional else ""),
        f"Permutations               : {args.n_permutations}   seed={args.seed}",
        "Shuffle                    : within lineage (population structure preserved)",
        "",
        f"Observed convergent        : {observed}",
        f"Null mean convergent       : {null_mean:.3f}",
        f"Null max                   : {int(null.max())}",
        f"Permutations >= observed   : {int((null >= observed).sum())} of {args.n_permutations}",
        f"Empirical p                : {emp_p:.4f}",
        "",
        f"Null distribution          : {[int(x) for x in null]}",
        "",
    ]
    if conditional:
        lines += [
            "CONDITIONAL NULL -- READ BEFORE QUOTING",
            "Stage 1 results were supplied externally and held fixed across",
            "every permutation, so the association test was not recomputed",
            "under the shuffled phenotype. Only Stages 3-5 varied. This is NOT",
            "a calibration of the convergent count and must not be reported as",
            "one. Recompute Stage 1 per permutation for that number.",
            "",
        ]
    if args.scope == "stage1":
        lines += [
            "NOTE: this null covers the association stage only. If a convergent",
            "count from all stages is reported, this is not its null -- re-run",
            "with --scope full so the null matches the statistic being claimed.",
            "",
        ]
    if emp_p <= 0.01:
        verdict = "observed count is well outside the null; the convergent result is not explained by population structure alone"
    elif emp_p <= 0.05:
        verdict = "observed count exceeds the null at conventional significance"
    else:
        verdict = "observed count is NOT distinguishable from the lineage-preserving null; do not report these candidates as convergent on this evidence"
    lines.append(f"Verdict: {verdict}")

    text = "\n".join(lines) + "\n"
    with open(f"{args.outdir}/permutation_summary.txt", "w") as fh:
        fh.write(text)
    print()
    print(text)

    if not args.keep_work:
        shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

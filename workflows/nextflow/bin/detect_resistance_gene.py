#!/usr/bin/env python3
"""Auto-detect the phenotype-defining resistance gene from the cohort's own
AMRFinderPlus calls, so a user need not already know which acquired gene to
name.

Why this exists: the pipeline previously required --resistance_gene as a
required argument, which meant a user had to already know their organism's
resistance biology before running anything. That is a reasonable thing to
know for one well-studied species and determinant (this project's own
blaOXA-23 in A. baumannii), and an unreasonable thing to require of "run this
on any cohort."

The heuristic, stated plainly because it IS a heuristic, not a discovery: among
genes AMRFinderPlus calls as acquired (not core/intrinsic -- filtered by
Element type / Element subtype where present) in at least --min-prevalence and
at most --max-prevalence of the cohort, pick the one closest to 50%
prevalence. A trait split close to 50/50 gives every downstream stage --
Fisher/FDR, the structure-corrected fit, the matched-pair comparison -- the
most statistical power; a trait at 3% or 97% starves several of them by
construction, independent of any biological question.

This is a starting point for automatic analysis, not a validated phenotype
assignment. It is reported, with the full prevalence table, precisely so a
domain expert can override it (--resistance-gene bypasses this script
entirely) rather than trusting it blindly.

Usage:
    detect_resistance_gene.py --amrfinder-dir amrfinder_reports/ \
        --output resistance_gene.txt --report resistance_gene_candidates.tsv
"""
from __future__ import annotations

import argparse
import glob
import sys

import pandas as pd


def load_amrfinder_calls(amrfinder_dir: str) -> pd.DataFrame:
    """One row per (sample, gene) call, from every *.amrfinder.tsv in the
    directory. Sample id is taken from the filename, not from a column inside
    the file, because AMRFinderPlus's own "Name" column is the assembly file
    name as passed to it, which varies with how the pipeline staged it."""
    rows = []
    files = sorted(glob.glob(f"{amrfinder_dir}/*.amrfinder.tsv"))
    if not files:
        sys.exit(f"no *.amrfinder.tsv files found in {amrfinder_dir}")
    for fn in files:
        sample = fn.split("/")[-1].removesuffix(".amrfinder.tsv")
        try:
            tab = pd.read_csv(fn, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if "Gene symbol" not in tab.columns:
            continue
        elem_type = tab["Element type"] if "Element type" in tab.columns else None
        elem_sub = tab["Element subtype"] if "Element subtype" in tab.columns else None
        for i, gene in enumerate(tab["Gene symbol"]):
            if pd.isna(gene) or not str(gene).strip():
                continue
            etype = str(elem_type.iloc[i]) if elem_type is not None else ""
            esub = str(elem_sub.iloc[i]) if elem_sub is not None else ""
            rows.append({"sample": sample, "gene": str(gene).strip(),
                        "element_type": etype, "element_subtype": esub})
    if not rows:
        sys.exit(f"no gene calls parsed from any file in {amrfinder_dir}")
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amrfinder-dir", required=True,
                    help="directory of <sample>.amrfinder.tsv files (AMRFinderPlus output)")
    ap.add_argument("--n-samples", type=int, required=True,
                    help="total cohort size (the denominator -- includes samples with "
                         "zero acquired-AMR calls, who are absent from every per-sample "
                         "file and would otherwise be silently dropped)")
    ap.add_argument("--min-prevalence", type=float, default=0.10)
    ap.add_argument("--max-prevalence", type=float, default=0.90)
    ap.add_argument("--exclude-core", action="store_true", default=True,
                    help="drop calls AMRFinderPlus itself flags as core/intrinsic "
                         "(Element subtype == AMR and Element type == AMR, but "
                         "Subclass indicates a point-mutation/intrinsic mechanism is "
                         "handled by excluding rows where Method starts with POINT)")
    ap.add_argument("--output", required=True, help="plain-text file: the chosen gene symbol")
    ap.add_argument("--report", required=True, help="TSV: every candidate gene considered")
    args = ap.parse_args()

    calls = load_amrfinder_calls(args.amrfinder_dir)
    print(f"Parsed {len(calls)} gene calls across {calls['sample'].nunique()} samples "
          f"with at least one call (cohort n={args.n_samples}).")

    counts = calls.groupby("gene")["sample"].nunique().rename("n_carriers").reset_index()
    counts["prevalence"] = counts["n_carriers"] / args.n_samples
    counts["distance_from_50pct"] = (counts["prevalence"] - 0.5).abs()
    counts = counts.sort_values("distance_from_50pct")

    eligible = counts[
        (counts.prevalence >= args.min_prevalence) & (counts.prevalence <= args.max_prevalence)
    ]
    counts["eligible"] = counts.gene.isin(eligible.gene)
    counts.to_csv(args.report, sep="\t", index=False)

    if eligible.empty:
        sys.exit(
            f"No acquired-AMR gene falls within the eligible prevalence range "
            f"[{args.min_prevalence:.0%}, {args.max_prevalence:.0%}] of this cohort "
            f"(n={args.n_samples}). Full table in {args.report}. Supply "
            f"--resistance-gene explicitly, or widen --min-prevalence/--max-prevalence."
        )

    chosen = eligible.iloc[0]
    print("\nCandidate resistance genes (top 10 by closeness to 50% prevalence):")
    print(eligible.head(10).to_string(index=False))
    print(f"\nAuto-selected: {chosen.gene}  "
          f"({int(chosen.n_carriers)}/{args.n_samples} = {chosen.prevalence:.1%})")
    print("This is a heuristic starting point, not a validated phenotype assignment -- "
          f"see {args.report} for every gene considered and override with "
          "--resistance-gene if you already know which determinant you want.")

    with open(args.output, "w") as fh:
        fh.write(str(chosen.gene))
    return 0


if __name__ == "__main__":
    sys.exit(main())

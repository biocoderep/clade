#!/usr/bin/env python3
"""Invoke the CLADE six-stage validation CLI from the pipeline.

A thin adapter, deliberately. It reads candidate names out of the same
name,position CSV that build_matrix.py consumes, and forwards every stage
option through to `clade validate` unchanged.

It stays thin on purpose: the pipeline and the permutation null must both run
the *identical* estimator, and the surest way to guarantee that is for neither
to reimplement it.

Exit codes are passed through, not swallowed:
    0 at least one convergent candidate
    1 ran cleanly, nothing convergent
    2 insufficient evidence
    3 error
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys


def read_candidate_names(path: str) -> list[str]:
    names: list[str] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames and "name" in reader.fieldnames:
            for row in reader:
                value = (row.get("name") or "").strip()
                if value:
                    names.append(value)
        else:
            fh.seek(0)
            for row in csv.reader(fh):
                if row and row[0].strip() and row[0].strip().lower() != "name":
                    names.append(row[0].strip())
    return names


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genotypes", required=True)
    p.add_argument("--phenotype", required=True)
    p.add_argument("--lineages", required=True)
    p.add_argument("--candidates-file", required=True)
    p.add_argument("--tree", required=True)
    p.add_argument("--distances", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--output-tsv")
    p.add_argument("--stage1-results")
    p.add_argument("--stage6")
    p.add_argument("--alpha", type=float)
    p.add_argument("--min-carriers", type=int)
    p.add_argument("--max-freq", type=float)
    p.add_argument("--temporal-margin", type=float)
    p.add_argument("--ambiguity-resolution", type=int, choices=[0, 1])
    p.add_argument("--unique-neighbors", action="store_true")
    args = p.parse_args()

    names = read_candidate_names(args.candidates_file)
    if not names:
        print(f"ERROR: no candidate names found in {args.candidates_file}", file=sys.stderr)
        return 3

    cmd = [
        "clade", "validate",
        "--genotypes", args.genotypes,
        "--phenotype", args.phenotype,
        "--lineages", args.lineages,
        "--candidates", *names,
        "--tree", args.tree,
        "--distances", args.distances,
        "--output", args.output,
    ]
    # Forward only what was actually supplied, so the CLI's own defaults stand
    # rather than being silently overwritten with this script's.
    for flag, value in (
        ("--output-tsv", args.output_tsv),
        ("--stage1-results", args.stage1_results),
        ("--stage6", args.stage6),
        ("--alpha", args.alpha),
        ("--min-carriers", args.min_carriers),
        ("--max-freq", args.max_freq),
        ("--temporal-margin", args.temporal_margin),
        ("--ambiguity-resolution", args.ambiguity_resolution),
    ):
        if value is not None:
            cmd += [flag, str(value)]
    if args.unique_neighbors:
        cmd.append("--unique-neighbors")

    print("Running:", " ".join(cmd), file=sys.stderr)
    return subprocess.run(cmd, check=False).returncode


if __name__ == "__main__":
    sys.exit(main())

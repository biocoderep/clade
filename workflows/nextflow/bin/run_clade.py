#!/usr/bin/env python3
"""
Thin wrapper around `clade validate`.

Exists for one real reason (caught during plan review, before this was built
wrong): the actual `clade validate --candidates` CLI takes candidate names as
separate arguments (`nargs="+"`), not a file path. This reads the same
`name,position` candidates CSV that `build_matrix.py` uses, extracts just the
`name` column, and expands it correctly.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys


def read_candidate_names(path: str) -> list[str]:
    names = []
    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].strip().lower() == "name":
                continue
            names.append(row[0].strip())
    return names


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genotypes", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--lineages", required=True)
    parser.add_argument("--candidates-file", required=True, help="Same name,position CSV build_matrix.py uses")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    names = read_candidate_names(args.candidates_file)
    if not names:
        print(f"ERROR: no candidate names found in {args.candidates_file}", file=sys.stderr)
        return 1

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
    print("Running:", " ".join(cmd), file=sys.stderr)
    result = subprocess.run(cmd, check=False)
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())

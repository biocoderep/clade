#!/usr/bin/env python3
"""
Compute a pairwise patristic distance matrix from a Newick tree.

This is a genuinely new script (not adapted from an existing one) — nothing
in this project previously computed a distance matrix reusable outside the
original analysis server. It exists specifically so CLADE's Stage 5
(matched-neighbor comparison) has a `--distances` input to work with; without
it, Stage 5 silently never runs (caught during the plan-review recheck for
this pipeline, before any code was written).

Uses dendropy, already a CLADE dependency (clade.phylogeny.fitch), not a new
library for this project.
"""
from __future__ import annotations

import argparse

import dendropy


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, help="Newick phylogeny")
    parser.add_argument("--output", required=True, help="Output TSV path (square matrix, tab-separated)")
    args = parser.parse_args()

    tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)
    pdm = tree.phylogenetic_distance_matrix()
    taxa = [t.label for t in tree.taxon_namespace]

    with open(args.output, "w") as f:
        f.write("\t" + "\t".join(taxa) + "\n")
        for t1 in tree.taxon_namespace:
            row = [t1.label]
            for t2 in tree.taxon_namespace:
                dist = 0.0 if t1 is t2 else pdm.patristic_distance(t1, t2)
                row.append(f"{dist:.6f}")
            f.write("\t".join(row) + "\n")

    print(f"Wrote {len(taxa)}x{len(taxa)} distance matrix to {args.output}")


if __name__ == "__main__":
    main()

"""
Phenotype provenance.

This module exists because of a real, caught error in the CLADE case
study: a hand-maintained phenotype file (`phenotypes_final.tsv`) was
used through an earlier round of analysis before being discovered — by
direct per-sample value comparison against an independently-derived
phenotype, not by any automated check — to represent a different AMR
determinant than its filename and documentation claimed. See
docs/reproducibility/phenotype_correction.md for the full account.

The fix encoded here: derive phenotype vectors directly and reproducibly
from their generating genotype source (an AMR gene presence/absence
matrix row) rather than trusting a separately-maintained file, and provide
a verification function that would have caught the original error
automatically had it existed at the time.
"""
from __future__ import annotations

import pandas as pd


def derive_phenotype_from_amr_matrix(rtab_path: str, gene_row: str) -> pd.Series:
    """Derive a binary phenotype directly from one row of an AMR
    presence/absence matrix (Rtab format: first column is the gene/feature
    name, remaining columns are sample IDs with 0/1 values).

    This is the source-of-truth path CLADE recommends: e.g.
    `derive_phenotype_from_amr_matrix("feature_matrix_final.Rtab", "amr_blaOXA-23")`
    reproduces exactly the `pheno.tsv` used throughout the CLADE case study.
    """
    with open(rtab_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == gene_row:
                values = pd.Series([int(x) for x in parts[1:]], index=samples, name=gene_row)
                return values
    raise KeyError(f"Row {gene_row!r} not found in {rtab_path}")


def verify_phenotype_provenance(claimed_phenotype: pd.Series, derived_phenotype: pd.Series) -> dict:
    """Compare a claimed (hand-maintained, or inherited) phenotype file
    against one freshly derived from its stated source, per-sample.

    Returns a dict with `n_compared`, `n_mismatched`, `mismatch_rate`, and
    `mismatched_samples` (capped at 20 for readability). A non-trivial
    mismatch rate is exactly the signal that caught the original
    phenotypes_final.tsv / blaOXA-66-vs-blaOXA-23 error — treat any
    mismatch rate above noise level as a provenance failure requiring
    investigation before the claimed file is used for anything.
    """
    common = claimed_phenotype.index.intersection(derived_phenotype.index)
    claimed = claimed_phenotype.loc[common]
    derived = derived_phenotype.loc[common]

    mismatched = common[claimed.values != derived.values]
    n = len(common)
    n_mismatch = len(mismatched)

    return {
        "n_compared": n,
        "n_mismatched": n_mismatch,
        "mismatch_rate": (n_mismatch / n) if n else float("nan"),
        "mismatched_samples": list(mismatched[:20]),
    }

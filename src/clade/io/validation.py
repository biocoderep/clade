"""
Input validation for CLADE.

Every CLADE stage previously accepted whatever it was handed and failed
silently on malformed input: non-binary genotype codes were counted as absent,
NaN became its own ancestral state, duplicated sample IDs misaligned matched
pairs, and a tree whose leaf labels did not match any sample ID produced a
clean-looking table of zeros. None of those raised.

The functions here fail loudly instead. They are deliberately strict: a
scientific pipeline that silently reinterprets its input is worse than one that
refuses to run.

`CladeInputError` is raised for all of it, so the CLI can turn any validation
failure into a clear message and a non-zero exit status rather than a traceback.
"""
from __future__ import annotations

import math

import pandas as pd


class CladeInputError(ValueError):
    """Raised when input data is malformed, inconsistent, or unusable."""


def is_missing(value) -> bool:
    """True for None and NaN.

    CLADE distinguishes missing from zero everywhere, so this predicate is
    used rather than an inline `x != x` NaN idiom, which reads as a typo.
    """
    if value is None:
        return True
    try:
        return bool(math.isnan(value))
    except (TypeError, ValueError):
        return False


def _fmt(items, limit: int = 5) -> str:
    items = list(items)
    shown = ", ".join(str(i) for i in items[:limit])
    if len(items) > limit:
        shown += f", ... ({len(items)} total)"
    return shown


def require_unique_index(obj, name: str) -> None:
    """Duplicated sample IDs silently misalign every downstream join and,
    in Stage 5, silently mispair resistant genomes with neighbours."""
    idx = obj.index
    if idx.duplicated().any():
        dups = idx[idx.duplicated()].unique()
        raise CladeInputError(
            f"{name} contains duplicated sample IDs: {_fmt(dups)}. "
            "Deduplicate before running CLADE — duplicates misalign joins and matched pairs."
        )


def require_non_empty(obj, name: str) -> None:
    if len(obj) == 0:
        raise CladeInputError(f"{name} is empty; nothing to analyse.")


def validate_binary(
    values: pd.Series,
    name: str,
    allow_missing: bool = True,
) -> pd.Series:
    """Check a 0/1 coded series.

    Missing values are permitted (they are meaningful — an uncalled genotype is
    not an absent genotype) but any *other* value is an error. Booleans and
    floats that are exactly 0.0/1.0 are accepted and normalised.
    """
    s = values
    if s.dtype == bool:
        return s.astype("float64")

    numeric = pd.to_numeric(s, errors="coerce")
    became_nan = numeric.isna() & s.notna()
    if became_nan.any():
        bad = s[became_nan].unique()
        raise CladeInputError(
            f"{name} contains non-numeric values: {_fmt(bad)}. Expected 0, 1, or missing."
        )

    present = numeric.dropna()
    invalid = present[~present.isin([0, 1])]
    if len(invalid) > 0:
        raise CladeInputError(
            f"{name} contains values outside {{0, 1}}: {_fmt(sorted(invalid.unique()))}. "
            "CLADE treats candidate carriage and phenotype as binary; recode before running."
        )

    if not allow_missing and numeric.isna().any():
        n = int(numeric.isna().sum())
        raise CladeInputError(f"{name} contains {n} missing value(s), which are not permitted here.")

    return numeric.astype("float64")


def validate_genotype_frame(
    genotypes: pd.DataFrame,
    candidate_cols: list[str] | None = None,
    allow_missing: bool = True,
) -> pd.DataFrame:
    """Validate and normalise a samples x candidates genotype matrix."""
    require_non_empty(genotypes, "genotype matrix")
    require_unique_index(genotypes, "genotype matrix")

    if candidate_cols is not None:
        missing = [c for c in candidate_cols if c not in genotypes.columns]
        if missing:
            raise CladeInputError(
                f"candidate column(s) not present in the genotype matrix: {_fmt(missing)}. "
                f"Available columns: {_fmt(genotypes.columns, limit=10)}"
            )
        cols = candidate_cols
    else:
        cols = list(genotypes.columns)

    dup_cols = [c for c in pd.Index(cols)[pd.Index(cols).duplicated()]]
    if dup_cols:
        raise CladeInputError(f"duplicated candidate column name(s): {_fmt(set(dup_cols))}")

    out = {}
    for c in cols:
        col = genotypes[c]
        if isinstance(col, pd.DataFrame):
            raise CladeInputError(
                f"candidate '{c}' appears more than once in the genotype matrix; "
                "column names must be unique."
            )
        out[c] = validate_binary(col, f"genotype column '{c}'", allow_missing=allow_missing)
    return pd.DataFrame(out, index=genotypes.index)


def validate_phenotype(phenotype: pd.Series, allow_missing: bool = False) -> pd.Series:
    """Phenotype must be binary. Missing phenotype is disallowed by default:
    a sample with no known resistance status cannot contribute to any stage,
    and silently coercing it to 0 inflates the susceptible group."""
    require_non_empty(phenotype, "phenotype")
    require_unique_index(phenotype, "phenotype")
    validated = validate_binary(phenotype, "phenotype", allow_missing=allow_missing)
    if validated.nunique(dropna=True) < 2:
        raise CladeInputError(
            "phenotype has fewer than two distinct values; association testing is "
            "undefined when every sample shares the same phenotype."
        )
    return validated


def require_sample_overlap(
    left_index,
    right_index,
    left_name: str,
    right_name: str,
    min_overlap: int = 1,
) -> list[str]:
    """Return the shared sample IDs, raising if the overlap is unusably small.

    Catches the common failure where two files use different ID conventions
    (accession vs. isolate name, or one with a suffix) and the inner join
    silently yields an empty or tiny dataset.
    """
    shared = [s for s in left_index if s in set(right_index)]
    if len(shared) < min_overlap:
        raise CladeInputError(
            f"{left_name} and {right_name} share {len(shared)} sample ID(s), "
            f"fewer than the required {min_overlap}. "
            f"{left_name} examples: {_fmt(list(left_index))}; "
            f"{right_name} examples: {_fmt(list(right_index))}. "
            "Check that both files use the same sample identifier convention."
        )
    return shared


def report_join_loss(before: int, after: int, name: str) -> str | None:
    """Describe how many samples a join dropped, or None if it dropped none."""
    if before == 0 or after == before:
        return None
    lost = before - after
    frac = lost / before
    return f"{name}: {lost} of {before} samples ({frac:.1%}) dropped by the join"


def validate_tree_against_samples(
    tree,
    sample_ids,
    min_overlap_fraction: float = 0.5,
) -> dict:
    """Check that a phylogeny's leaf labels actually correspond to the samples.

    Without this, a label-convention mismatch makes every leaf an unknown state,
    Fitch resolves the whole tree to the root default, and Stage 4 returns a
    confident-looking table of zero gains.
    """
    leaf_labels = {
        leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon is not None
    }
    samples = set(sample_ids)
    shared = leaf_labels & samples

    if not leaf_labels:
        raise CladeInputError("the supplied tree has no labelled leaves.")

    frac_of_samples = len(shared) / len(samples) if samples else 0.0
    if frac_of_samples < min_overlap_fraction:
        raise CladeInputError(
            f"only {len(shared)} of {len(samples)} samples ({frac_of_samples:.1%}) are present "
            f"as leaves in the tree, below the required {min_overlap_fraction:.0%}. "
            f"Tree leaf examples: {_fmt(sorted(leaf_labels))}; "
            f"sample examples: {_fmt(sorted(samples))}. "
            "Stage 4 would otherwise reconstruct ancestral states from almost no data."
        )

    return {
        "n_tree_leaves": len(leaf_labels),
        "n_samples": len(samples),
        "n_shared": len(shared),
        "samples_not_in_tree": sorted(samples - leaf_labels),
        "leaves_not_in_samples": sorted(leaf_labels - samples),
    }


def validate_distance_matrix(distances: pd.DataFrame) -> pd.DataFrame:
    """Check a pairwise distance matrix is square, labelled consistently, and
    free of duplicate labels."""
    require_non_empty(distances, "distance matrix")
    require_unique_index(distances, "distance matrix (rows)")
    if distances.columns.duplicated().any():
        dups = distances.columns[distances.columns.duplicated()].unique()
        raise CladeInputError(f"distance matrix has duplicated column labels: {_fmt(dups)}")
    if list(distances.index) != list(distances.columns):
        only_rows = set(distances.index) - set(distances.columns)
        only_cols = set(distances.columns) - set(distances.index)
        if only_rows or only_cols:
            raise CladeInputError(
                "distance matrix is not square/consistently labelled. "
                f"Only in rows: {_fmt(only_rows)}; only in columns: {_fmt(only_cols)}"
            )
        distances = distances.loc[:, list(distances.index)]
    return distances

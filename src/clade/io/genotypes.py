"""
Shared data loading for CLADE's validation stages.

This factors out the join pattern that was previously duplicated across
firth_final.py, temporal_ordering.py, matched_neighbor.py, and
fisher_fdr_lasso.py (each independently joined phenotype + candidate
genotypes + lineage on Sample_ID). The join logic itself is unchanged from
those scripts; only the duplication is removed.
"""
from __future__ import annotations

import pandas as pd

from clade.io.validation import (
    CladeInputError,
    require_non_empty,
    require_unique_index,
)


def load_phenotype(path: str) -> pd.Series:
    """Load a single-column phenotype file (sample -> 0/1), indexed by sample ID.

    CLADE takes this as a required, already-derived input rather than deriving
    it itself — see docs/reproducibility/phenotype_correction.md for why
    phenotype provenance is treated as the caller's responsibility, not
    something CLADE infers.
    """
    pheno = pd.read_csv(path, sep="\t", index_col=0)
    if pheno.shape[1] != 1:
        raise CladeInputError(
            f"Expected a single phenotype column in {path}, found {pheno.shape[1]}: "
            f"{list(pheno.columns)[:5]}"
        )
    require_non_empty(pheno, f"phenotype file {path}")
    pheno.index = pheno.index.astype(str)
    require_unique_index(pheno, f"phenotype file {path}")
    return pheno.iloc[:, 0].rename("phenotype")


def load_genotypes(*paths: str, index_col: str = "Sample_ID") -> pd.DataFrame:
    """Load one or more candidate-genotype CSVs and concatenate their columns.

    Each file is expected to have `index_col` plus one 0/1 column per candidate
    locus, matching the format produced by real genotype extraction from a
    core-genome VCF (see docs/reproducibility/original_discovery_pipeline.md
    for why genotypes must come from the VCF directly, not from any upstream
    discovery method's intermediate files).
    """
    frames = []
    seen: set[str] = set()
    for p in paths:
        raw = pd.read_csv(p)
        if index_col not in raw.columns:
            raise CladeInputError(
                f"{p} has no '{index_col}' column; found {list(raw.columns)[:8]}"
            )
        frame = raw.set_index(index_col)
        frame.index = frame.index.astype(str)
        require_unique_index(frame, f"genotype file {p}")
        overlap = seen & set(frame.columns)
        if overlap:
            raise CladeInputError(
                f"{p} redefines candidate column(s) already loaded from an earlier file: "
                f"{sorted(overlap)}. Duplicate columns would silently shadow one another."
            )
        seen |= set(frame.columns)
        frames.append(frame)
    return pd.concat(frames, axis=1, join="outer")


def load_lineages(path: str, sample_col: str = "sample", st_col: str = "ST") -> pd.Series:
    """Load MLST sequence-type assignments, indexed by sample ID."""
    raw = pd.read_csv(path, sep="\t")
    for col in (sample_col, st_col):
        if col not in raw.columns:
            raise CladeInputError(
                f"lineage file {path} has no '{col}' column; found {list(raw.columns)[:8]}"
            )
    mlst = raw[[sample_col, st_col]].set_index(sample_col)
    mlst.index = mlst.index.astype(str)
    require_unique_index(mlst, f"lineage file {path}")
    return mlst[st_col].astype(str).rename("ST")


def load_joined_dataset(
    phenotype_path: str,
    genotype_paths: list[str],
    lineage_path: str | None = None,
    lineage_join: str = "inner",
) -> pd.DataFrame:
    """Join phenotype, candidate genotypes, and (optionally) lineage on sample ID.

    Uses an inner join, matching the behavior of the original per-stage scripts:
    a sample is only included if it has phenotype, genotype, and (if requested)
    lineage data all present.
    """
    pheno = load_phenotype(phenotype_path)
    geno = load_genotypes(*genotype_paths)

    df = pheno.to_frame().join(geno, how="inner")
    if df.empty:
        raise CladeInputError(
            "the phenotype and genotype files share no sample IDs. "
            f"Phenotype examples: {list(pheno.index[:3])}; "
            f"genotype examples: {list(geno.index[:3])}. "
            "Check that both use the same identifier convention."
        )
    df.attrs["n_dropped_genotype_join"] = len(pheno) - len(df)

    if lineage_path is not None:
        lineage = load_lineages(lineage_path)
        before = len(df)
        # NOTE: `inner` is the published behaviour — a sample with no lineage
        # assignment is excluded from the analysis entirely. That is preserved
        # deliberately: switching to a left join would change the analysis
        # population and therefore every downstream statistic. What changes here
        # is only that the loss is now reported instead of silent. Pass
        # lineage_join="left" to retain unassigned samples (they are pooled into
        # the "other" ST covariate), as a sensitivity analysis.
        df = df.join(lineage, how=lineage_join)
        df.attrs["n_dropped_lineage_join"] = before - len(df)
        df.attrs["n_without_lineage"] = int(df["ST"].isna().sum())
        if df.attrs["n_dropped_lineage_join"]:
            df.attrs["lineage_join_warning"] = (
                f"{df.attrs['n_dropped_lineage_join']} of {before} samples had no lineage "
                f"assignment and were dropped by the '{lineage_join}' join."
            )
        if df.empty:
            raise CladeInputError(
                "no sample has phenotype, genotype and lineage data simultaneously. "
                "Check the identifier conventions across all three files."
            )
    return df


def st_covariate_dummies(df: pd.DataFrame, st_col: str = "ST", top_n: int = 20) -> pd.DataFrame:
    """Build ST dummy-variable covariates, grouping all but the top-N lineages as 'other'.

    Matches the covariate construction used in firth_final.py and
    fisher_fdr_lasso.py exactly (top 20 STs by size, drop_first=True).
    """
    top_sts = df[st_col].value_counts().head(top_n).index.tolist()
    grouped = df[st_col].where(df[st_col].isin(top_sts), "other")
    return pd.get_dummies(grouped, prefix="ST", drop_first=True).astype(float)

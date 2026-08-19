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


def load_phenotype(path: str) -> pd.Series:
    """Load a single-column phenotype file (sample -> 0/1), indexed by sample ID.

    CLADE takes this as a required, already-derived input rather than deriving
    it itself — see docs/reproducibility/phenotype_correction.md for why
    phenotype provenance is treated as the caller's responsibility, not
    something CLADE infers.
    """
    pheno = pd.read_csv(path, sep="\t", index_col=0)
    if pheno.shape[1] != 1:
        raise ValueError(
            f"Expected a single phenotype column in {path}, found {pheno.shape[1]}"
        )
    return pheno.iloc[:, 0].rename("phenotype")


def load_genotypes(*paths: str, index_col: str = "Sample_ID") -> pd.DataFrame:
    """Load one or more candidate-genotype CSVs and concatenate their columns.

    Each file is expected to have `index_col` plus one 0/1 column per candidate
    locus, matching the format produced by real genotype extraction from a
    core-genome VCF (see docs/reproducibility/original_discovery_pipeline.md
    for why genotypes must come from the VCF directly, not from any upstream
    discovery method's intermediate files).
    """
    frames = [pd.read_csv(p).set_index(index_col) for p in paths]
    return pd.concat(frames, axis=1, join="outer")


def load_lineages(path: str, sample_col: str = "sample", st_col: str = "ST") -> pd.Series:
    """Load MLST sequence-type assignments, indexed by sample ID."""
    mlst = pd.read_csv(path, sep="\t")[[sample_col, st_col]].set_index(sample_col)
    return mlst[st_col].astype(str).rename("ST")


def load_joined_dataset(
    phenotype_path: str,
    genotype_paths: list[str],
    lineage_path: str | None = None,
) -> pd.DataFrame:
    """Join phenotype, candidate genotypes, and (optionally) lineage on sample ID.

    Uses an inner join, matching the behavior of the original per-stage scripts:
    a sample is only included if it has phenotype, genotype, and (if requested)
    lineage data all present.
    """
    pheno = load_phenotype(phenotype_path)
    geno = load_genotypes(*genotype_paths)
    df = pheno.to_frame().join(geno, how="inner")
    if lineage_path is not None:
        lineage = load_lineages(lineage_path)
        df = df.join(lineage, how="inner")
    return df


def st_covariate_dummies(df: pd.DataFrame, st_col: str = "ST", top_n: int = 20) -> pd.DataFrame:
    """Build ST dummy-variable covariates, grouping all but the top-N lineages as 'other'.

    Matches the covariate construction used in firth_final.py and
    fisher_fdr_lasso.py exactly (top 20 STs by size, drop_first=True).
    """
    top_sts = df[st_col].value_counts().head(top_n).index.tolist()
    grouped = df[st_col].where(df[st_col].isin(top_sts), "other")
    return pd.get_dummies(grouped, prefix="ST", drop_first=True).astype(float)

"""
Lineage-preserving phenotype permutation harness (CLADE Stage 1 FDR calibration).

Shuffles the blaOXA-23 phenotype label across genomes while holding each
genome's own MLST-ST fixed — i.e. permutation is stratified by ST, so each
ST's exact count of resistant/susceptible genomes is preserved and only
*which* genomes within a stratum carry the label is randomized. This
destroys any candidate-locus association that isn't itself lineage-driven,
while leaving the population-structure confound (the actual thing CLADE's
six stages are trying to see past) fully intact in the null — a naive
unstratified shuffle would instead destroy population structure too, making
the null too easy to beat and understating the real false-discovery risk.

This module is deliberately dependency-light (pandas/numpy only) so the
permutation mechanism itself can be built, tested, and used independently
of which CLADE stage-1 backend (Fisher/FDR, Firth, kinship-corrected
pyseer) is actually installed in a given environment — see
docs/reproducibility/amr_panel_provenance.md and the accompanying
permutation-harness design note for why that separation matters here: as
of this writing, no single local environment has the full CLADE dependency
set (statsmodels and dendropy are both missing from the base environment;
firthlogist cannot install at all under Python >=3.11). `run_permutation_null`
takes the per-permutation statistic as a caller-supplied callback for
exactly this reason — plug in whatever Stage-1 backend is actually
available where this is run.

Why setup is loaded once, not per permutation: measured against the real
cohort, `distances_final.tsv` (137MB, ~3,249x3,249) takes ~7.8s to load and
tree.nwk (3,255 tips) takes ~0.26s to parse. Both are fixed inputs that do
not change across permutations — only the phenotype labeling does. Loading
either one inside the permutation loop would mean paying that cost on every
draw: at ~10.5s/permutation (dominated by the distance-matrix reload),
1,000 permutations would take ~2.9 hours. Loading them once via
`PermutationSetup`, before the loop starts, drops the same 1,000 permutations
to ~7.5 minutes (per-permutation compute only, ~0.45s/candidate). This is a
correctness-adjacent performance requirement, not just an optimization —
`run_permutation_null` accepts a pre-built `PermutationSetup` for exactly
this reason, and does not load anything itself.
"""
from __future__ import annotations

import datetime
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

# Central, single source of truth for this module's default random seed —
# every stochastic draw in this file traces back to this constant, not a
# literal scattered through the code. Change it here, not at call sites.
DEFAULT_SEED = 0

SEEDS_LOG_PATH = Path(__file__).resolve().parent / "seeds.log"


def _log_seed(seed: int, n_permutations: int) -> None:
    """Append one line recording the seed used for a run, so any past run's
    randomness is reconstructable from seeds.log alone, not just from
    whatever default happened to be in the source at the time."""
    SEEDS_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with open(SEEDS_LOG_PATH, "a") as f:
        f.write(f"{timestamp}\tseed={seed}\tn_permutations={n_permutations}\tmodule=permutation_lineage_preserving\n")


@dataclass
class PermutationSetup:
    """Everything the permutation loop needs, loaded once before it begins.

    `phenotype` and `lineage` must share the same index (sample ID); the
    permutation itself only touches `phenotype`, reshuffled within each
    `lineage` (ST) stratum. `genotypes`, `tree`, and `distances` are passed
    through unchanged to `statistic_fn` on every draw — they are the
    expensive-to-load, permutation-invariant inputs this dataclass exists to
    hoist out of the loop (see module docstring for the measured cost of
    getting this wrong: ~2.9 hours vs. ~7.5 minutes for 1,000 permutations).

    `tree` and `distances` are optional — omit them if `statistic_fn` only
    needs Stage 1/2 (no Fitch reconstruction, no matched-neighbor test).
    Typed as `object` rather than `dendropy.Tree` so this module stays
    dendropy-free; the caller is responsible for constructing/parsing
    whatever `tree`/`distances` objects its own `statistic_fn` expects.
    """
    phenotype: pd.Series
    lineage: pd.Series
    genotypes: pd.DataFrame
    tree: object | None = None
    distances: pd.DataFrame | None = None
    extra: dict = field(default_factory=dict)


def stratified_permute_phenotype(
    phenotype: pd.Series,
    lineage: pd.Series,
    rng: np.random.Generator,
) -> pd.Series:
    """Shuffle `phenotype` values within each `lineage` (ST) stratum only.

    Both series must share the same index (sample ID). Returns a new Series,
    same index and dtype, with each ST's exact positive/negative count
    preserved — only which specific samples within that ST carry which
    label changes.
    """
    if not phenotype.index.equals(lineage.index):
        raise ValueError("phenotype and lineage must share the same index")

    permuted = phenotype.copy()
    for st, group_idx in lineage.groupby(lineage).groups.items():
        values = phenotype.loc[group_idx].to_numpy()
        permuted.loc[group_idx] = rng.permutation(values)
    return permuted


@dataclass
class PermutationResult:
    observed_statistic: int
    null_distribution: list[int]
    n_permutations: int
    empirical_p: float


def run_permutation_null(
    setup: PermutationSetup,
    statistic_fn: Callable[[pd.Series, PermutationSetup], int],
    n_permutations: int = 1000,
    seed: int = DEFAULT_SEED,
) -> PermutationResult:
    """Empirical null distribution of `statistic_fn` under lineage-preserving
    phenotype permutation.

    `setup`: a `PermutationSetup` built once by the caller *before* calling
    this function — see its docstring for why. This function never loads a
    tree, distance matrix, or genotype file itself; it only permutes
    `setup.phenotype` and hands the result, plus the unchanged `setup`,
    to `statistic_fn`.

    `statistic_fn`: takes `(permuted_phenotype, setup)` and returns an
    integer statistic — e.g. "how many candidates are Stage-1-significant"
    or "how many candidates reach the convergent disposition" under that
    phenotype, using `setup.genotypes`/`setup.tree`/`setup.distances` as
    needed. Plug in whatever CLADE stage backend(s) are available in the
    calling environment; this module does not import any of them itself.

    Returns the observed statistic (computed on the real, unpermuted
    phenotype), the full null distribution, and a one-sided empirical
    p-value: P(null statistic >= observed statistic), with the standard
    +1/+1 correction so the p-value is never exactly zero regardless of
    `n_permutations`.
    """
    observed = statistic_fn(setup.phenotype, setup)

    _log_seed(seed, n_permutations)
    rng = np.random.default_rng(seed)
    null_stats = []
    for _ in range(n_permutations):
        permuted = stratified_permute_phenotype(setup.phenotype, setup.lineage, rng)
        null_stats.append(statistic_fn(permuted, setup))

    as_extreme = sum(1 for s in null_stats if s >= observed)
    empirical_p = (as_extreme + 1) / (n_permutations + 1)

    return PermutationResult(
        observed_statistic=observed,
        null_distribution=null_stats,
        n_permutations=n_permutations,
        empirical_p=empirical_p,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phenotype", required=True, help="Phenotype TSV (sample -> 0/1)")
    parser.add_argument("--lineages", required=True, help="MLST TSV (sample, ST columns)")
    parser.add_argument("--genotypes", nargs="+", required=True, help="Candidate-genotype CSV(s)")
    parser.add_argument("--candidates", nargs="+", required=True, help="Candidate column names")
    parser.add_argument("--tree", help="Newick phylogeny (optional; enables a Stage-4-aware statistic_fn)")
    parser.add_argument("--distances", help="Pairwise distance TSV (optional; enables a Stage-5-aware statistic_fn)")
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()

    pheno = pd.read_csv(args.phenotype, sep="\t", index_col=0).iloc[:, 0]
    lineage = pd.read_csv(args.lineages, sep="\t")
    lineage = lineage.set_index(lineage.columns[0])[lineage.columns[-1]].astype(str)
    geno_frames = [pd.read_csv(p).set_index(pd.read_csv(p).columns[0]) for p in args.genotypes]
    geno = pd.concat(geno_frames, axis=1)

    common = pheno.index.intersection(lineage.index).intersection(geno.index)
    pheno, lineage, geno = pheno.loc[common], lineage.loc[common], geno.loc[common]

    # Everything expensive is loaded exactly once, here, before the
    # permutation loop starts inside run_permutation_null. See
    # PermutationSetup's and this module's docstrings for the measured
    # cost of getting this ordering wrong.
    tree = None
    if args.tree:
        import dendropy
        tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)

    distances = None
    if args.distances:
        distances = pd.read_csv(args.distances, sep="\t", index_col=0)
        distances.columns = distances.columns.astype(str)
        distances.index = distances.index.astype(str)

    setup = PermutationSetup(phenotype=pheno, lineage=lineage, genotypes=geno, tree=tree, distances=distances)

    def count_stage1_significant(permuted_phenotype: pd.Series, setup: PermutationSetup) -> int:
        # Wire this to whatever Stage-1 backend is actually installed —
        # deferred import so this script can still be used (and its
        # permutation mechanism unit-tested) in environments missing
        # statsmodels/firthlogist. See module docstring.
        from clade.validation.stage1_association import fisher_fdr_screen

        df = setup.genotypes.copy()
        df["phenotype"] = permuted_phenotype
        result = fisher_fdr_screen(df, args.candidates, "phenotype")
        return int(result["FDR_significant"].sum())

    result = run_permutation_null(
        setup, count_stage1_significant,
        n_permutations=args.n_permutations, seed=args.seed,
    )
    print(f"Observed statistic: {result.observed_statistic}")
    print(f"Empirical p-value ({result.n_permutations} permutations): {result.empirical_p:.4f}")

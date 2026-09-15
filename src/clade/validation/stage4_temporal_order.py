"""
Stage 4 — Phylogenetic temporal ordering.

Refactor of temporal_ordering.py built on clade.phylogeny.fitch (the same
two-pass Fitch parsimony algorithm, unchanged for well-formed input).

Reads as: for each candidate, reconstruct where on the tree it was gained, and
ask what fraction of those independent gains occurred on a branch whose parent
was already resistant. A compensatory mutation should arise *after* resistance,
so its percentage should exceed the cohort's background resistance prevalence.

Two cautions are built into the API rather than left to the reader:

* `Pct_post_resistance` is NaN when no gain event was reconstructed. That is
  undefined, not zero, and `StageResults.stage4_supports()` returns None for it
  rather than treating it as a failure.
* The result depends on tree topology *and* on how ambiguous ancestral states
  are resolved. `temporal_ordering_sensitivity` re-runs under both resolution
  policies so a candidate whose verdict flips can be identified rather than
  reported at face value. See
  docs/limitations/phylogenetic_alignment_sensitivity.md.
"""
from __future__ import annotations

import dendropy
import pandas as pd

from clade.io.validation import is_missing, validate_tree_against_samples
from clade.phylogeny.fitch import count_gain_events, fitch_reconstruct


def temporal_ordering(
    tree: dendropy.Tree,
    phenotype_map: dict,
    candidate_maps: dict[str, dict],
    ambiguity_resolution: int = 0,
    validate: bool = True,
    min_overlap_fraction: float = 0.5,
) -> pd.DataFrame:
    """Reconstruct ancestral resistance state, then for each candidate
    reconstruct its own ancestral state and count 0->1 gain events by whether
    resistance was already present at the parent node.

    `phenotype_map` / each value in `candidate_maps`: taxon label -> 0/1/None.

    Returns one row per candidate with Total_gains, Gains_on_resistant_bg,
    Gains_on_susceptible_bg and Pct_post_resistance — the same schema as the
    original temporal_ordering_results.csv.
    """
    if validate:
        validate_tree_against_samples(
            tree, phenotype_map.keys(), min_overlap_fraction=min_overlap_fraction
        )

    tree = fitch_reconstruct(tree, phenotype_map, ambiguity_resolution=ambiguity_resolution)
    # Freeze the phenotype reconstruction onto every node before the candidate
    # passes overwrite .fitch_final. This attribute is what count_gain_events
    # reads as the resistance background.
    for node in tree.preorder_node_iter():
        node.oxa_final = node.fitch_final

    rows = []
    for candidate, cand_map in candidate_maps.items():
        reconstructed = fitch_reconstruct(
            tree, cand_map, ambiguity_resolution=ambiguity_resolution, require_overlap=False
        )
        stats = count_gain_events(reconstructed, background_attr="oxa_final")
        rows.append({"Candidate": candidate, **stats})

    return pd.DataFrame(rows)


def temporal_ordering_sensitivity(
    tree: dendropy.Tree,
    phenotype_map: dict,
    candidate_maps: dict[str, dict],
    baseline: float | None = None,
    temporal_margin: float = 0.10,
    **kwargs,
) -> pd.DataFrame:
    """Run Stage 4 under both ambiguity-resolution policies and report whether
    each candidate's verdict is stable.

    `resolution=0` (default, DELTRAN-like) dates gains as late as possible;
    `resolution=1` (ACCTRAN-like) as early as possible. A candidate whose
    Pct_post_resistance crosses the decision threshold between the two is
    topology/policy-sensitive and should not be reported from Stage 4 alone.
    """
    frames = {}
    for resolution in (0, 1):
        # dendropy trees are mutated in place; re-parse so the two runs are independent.
        t = dendropy.Tree.get(
            data=tree.as_string(schema="newick"), schema="newick", preserve_underscores=True
        )
        frames[resolution] = temporal_ordering(
            t, phenotype_map, candidate_maps, ambiguity_resolution=resolution, **kwargs
        ).set_index("Candidate")

    out = pd.DataFrame(
        {
            "Pct_post_resistance_delayed": frames[0]["Pct_post_resistance"],
            "Pct_post_resistance_accelerated": frames[1]["Pct_post_resistance"],
            "Total_gains_delayed": frames[0]["Total_gains"],
            "Total_gains_accelerated": frames[1]["Total_gains"],
        }
    )

    if baseline is not None:
        threshold = baseline + temporal_margin
        passes_delayed = out["Pct_post_resistance_delayed"] >= threshold
        passes_accel = out["Pct_post_resistance_accelerated"] >= threshold
        out["verdict_stable"] = passes_delayed == passes_accel

    return out.reset_index()


def cohort_baseline_prevalence(phenotype_map: dict) -> float:
    """The cohort-wide resistance prevalence, used as the reference baseline
    against which each candidate's Pct_post_resistance should be compared."""
    values = [
        v
        for v in phenotype_map.values()
        if not is_missing(v)
    ]
    return sum(values) / len(values) if values else float("nan")

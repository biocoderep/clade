"""
Stage 4 — Phylogenetic temporal ordering.

Direct refactor of temporal_ordering.py into a testable function, built on
clade.phylogeny.fitch (the same two-pass Fitch parsimony algorithm,
unchanged). See docs/limitations/phylogenetic_alignment_sensitivity.md —
this stage's output depends entirely on tree topology, which CLADE's own
150-genome comparison found is not fully robust to core-genome alignment
choice. Report Stage 4 results alongside at least one other stage, not
in isolation.
"""
from __future__ import annotations

import dendropy
import pandas as pd

from clade.phylogeny.fitch import count_gain_events, fitch_reconstruct


def temporal_ordering(tree: dendropy.Tree, phenotype_map: dict, candidate_maps: dict[str, dict]) -> pd.DataFrame:
    """Reconstruct ancestral resistance state, then for each candidate,
    reconstruct its own ancestral state and count 0->1 gain events by
    whether resistance was already present at the parent node.

    `phenotype_map` / each value in `candidate_maps`: taxon label -> 0/1.
    Returns one row per candidate with Total_gains, Gains_on_resistant_bg,
    Gains_on_susceptible_bg, Pct_post_resistance — same schema as the
    original temporal_ordering_results.csv.
    """
    tree = fitch_reconstruct(tree, phenotype_map)
    for node in tree.preorder_node_iter():
        node.oxa_final = node.fitch_final

    rows = []
    for candidate, cand_map in candidate_maps.items():
        tree2 = fitch_reconstruct(tree, cand_map)
        stats = count_gain_events(tree2, background_attr="oxa_final")
        rows.append({"Candidate": candidate, **stats})

    return pd.DataFrame(rows)


def cohort_baseline_prevalence(phenotype_map: dict) -> float:
    """The cohort-wide resistance prevalence, used as the reference baseline
    against which each candidate's Pct_post_resistance should be compared."""
    values = [v for v in phenotype_map.values() if v is not None]
    return sum(values) / len(values) if values else float("nan")

"""
Fitch-parsimony ancestral-state reconstruction for a binary character.

This is the exact algorithm originally written directly inside
temporal_ordering.py (two-pass: postorder state-set intersection/union,
preorder final-state resolution), factored out here unchanged so it can be
unit-tested independently and reused for both the phenotype and each
candidate's reconstruction (Stage 4).

Known limitation: this operates on tree *topology* only, and CLADE's own
150-genome comparison found topology is not as robust to core-genome
alignment choice (SNP-only vs. full-alignment) as initially assumed — see
docs/limitations/phylogenetic_alignment_sensitivity.md before treating
Stage 4 results as bias-proof in isolation.
"""
from __future__ import annotations

import dendropy


def fitch_reconstruct(tree: dendropy.Tree, state_map: dict) -> dendropy.Tree:
    """Two-pass Fitch parsimony for a binary (0/1) character.

    `state_map`: taxon label -> 0 or 1. Taxa missing from `state_map` are
    treated as unknown ({0, 1}) at the leaf, matching the original script's
    behavior. Mutates and returns `tree`, setting `.fitch_set` and
    `.fitch_final` on every node.
    """
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label if leaf.taxon else None
        val = state_map.get(label)
        leaf.fitch_set = {val} if val is not None else {0, 1}

    for node in tree.postorder_internal_node_iter():
        child_sets = [c.fitch_set for c in node.child_nodes()]
        inter = set.intersection(*child_sets)
        node.fitch_set = inter if inter else set.union(*child_sets)

    root = tree.seed_node
    root.fitch_final = min(root.fitch_set)
    for node in tree.preorder_node_iter():
        if node is root:
            continue
        parent_state = node.parent_node.fitch_final
        node.fitch_final = parent_state if parent_state in node.fitch_set else min(node.fitch_set)

    return tree


def count_gain_events(tree: dendropy.Tree, background_attr: str = "oxa_final") -> dict:
    """Count 0->1 transitions on `tree` (after `fitch_reconstruct` has set
    `.fitch_final`) and classify each by the resistance background at the
    parent node (read from `background_attr`, set by a prior Fitch pass on
    the phenotype).

    Returns a dict with Total_gains, Gains_on_resistant_bg,
    Gains_on_susceptible_bg, Pct_post_resistance — same fields and
    definitions as the original temporal_ordering.py output.
    """
    gains_on_resistant_bg = 0
    gains_on_susceptible_bg = 0
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        parent_state = node.parent_node.fitch_final
        child_state = node.fitch_final
        if parent_state == 0 and child_state == 1:
            bg = getattr(node.parent_node, background_attr)
            if bg == 1:
                gains_on_resistant_bg += 1
            else:
                gains_on_susceptible_bg += 1

    total_gains = gains_on_resistant_bg + gains_on_susceptible_bg
    pct = gains_on_resistant_bg / total_gains if total_gains else float("nan")
    return {
        "Total_gains": total_gains,
        "Gains_on_resistant_bg": gains_on_resistant_bg,
        "Gains_on_susceptible_bg": gains_on_susceptible_bg,
        "Pct_post_resistance": pct,
    }

"""
Fitch-parsimony ancestral-state reconstruction for a binary character.

This is the two-pass algorithm originally written inside temporal_ordering.py
(postorder state-set intersection/union, preorder final-state resolution),
factored out so it can be unit-tested independently and reused for both the
phenotype and each candidate's reconstruction (Stage 4).

**The reconstruction itself is unchanged** — for well-formed binary input this
module returns exactly what the original script returned. What has been added
is refusal to proceed on input the original silently mis-handled:

* a non-binary or NaN leaf state used to become its own Fitch state, so a
  missing genotype call was reconstructed as a third character state and could
  manufacture or suppress gain events. It now raises.
* a `state_map` sharing no labels with the tree used to return a clean table of
  zero gains. `require_overlap` now catches that.

Two conventions matter for interpretation and are named explicitly rather than
left implicit in a `min()` call:

**Root resolution.** Where the root's Fitch set is ambiguous ({0, 1}), the root
is resolved to `ambiguity_resolution` (default 0, the original behaviour). For a
presence/absence character this assumes the ancestor lacked the mutation, which
*maximises* the number of inferred 0->1 gains. This is the conventional choice
for derived-state counting, but it is a choice, and a candidate whose gain count
is sensitive to it should be treated with caution.

**Internal ties.** Where a node's state is ambiguous and its parent's state is
not in its Fitch set, the node takes `ambiguity_resolution`. With the default of
0 this is a DELTRAN-like policy: changes are pushed toward the tips, so gains
are dated as late as possible. ACCTRAN (early gains) would be `1`. The two can
differ in how many gains land on a resistant background, which is exactly what
Stage 4 measures — see `docs/limitations/phylogenetic_alignment_sensitivity.md`.

Known limitation: this operates on tree *topology* only, and CLADE's own
150-genome comparison found topology is not as robust to core-genome alignment
choice (SNP-only vs. full-alignment) as initially assumed. Do not treat Stage 4
as bias-proof in isolation.
"""
from __future__ import annotations

import dendropy

from clade.io.validation import CladeInputError, is_missing

VALID_STATES = (0, 1)


def _coerce_state(value, label: str):
    """Return 0, 1, or None (unknown). Anything else is an error."""
    if value is None:
        return None
    if is_missing(value):
        return None
    if isinstance(value, bool):
        return int(value)
    try:
        as_int = int(value)
    except (TypeError, ValueError):
        raise CladeInputError(
            f"leaf '{label}' has non-binary state {value!r}; Fitch reconstruction "
            "requires 0, 1, or a missing value."
        ) from None
    if as_int != value or as_int not in VALID_STATES:
        raise CladeInputError(
            f"leaf '{label}' has state {value!r}, which is not 0 or 1. A non-binary state "
            "would be reconstructed as a third character state and corrupt the gain count."
        )
    return as_int


def fitch_reconstruct(
    tree: dendropy.Tree,
    state_map: dict,
    ambiguity_resolution: int = 0,
    require_overlap: bool = True,
    min_known_leaves: int = 1,
) -> dendropy.Tree:
    """Two-pass Fitch parsimony for a binary (0/1) character.

    `state_map`: taxon label -> 0, 1, or None. Taxa missing from `state_map`,
    or mapped to None/NaN, are treated as unknown ({0, 1}) at the leaf —
    the original behaviour, and the correct one for genuinely uncalled data.

    Any *other* leaf value raises `CladeInputError` rather than becoming a
    third state.

    Mutates and returns `tree`, setting `.fitch_set` and `.fitch_final` on
    every node, plus `.fitch_n_known_leaves` on the tree.
    """
    if ambiguity_resolution not in VALID_STATES:
        raise CladeInputError(f"ambiguity_resolution must be 0 or 1, got {ambiguity_resolution!r}")

    n_known = 0
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label if leaf.taxon else None
        val = _coerce_state(state_map.get(label), str(label))
        if val is None:
            leaf.fitch_set = set(VALID_STATES)
        else:
            leaf.fitch_set = {val}
            n_known += 1

    if require_overlap and n_known < min_known_leaves:
        raise CladeInputError(
            f"none of the tree's leaf labels are present in the supplied state map "
            f"(0 of {len(list(tree.leaf_node_iter()))} leaves have a known state). "
            "Ancestral reconstruction would be vacuous. Check that the tree's leaf "
            "labels use the same identifiers as the genotype/phenotype files."
        )

    tree.fitch_n_known_leaves = n_known

    for node in tree.postorder_internal_node_iter():
        child_sets = [c.fitch_set for c in node.child_nodes()]
        inter = set.intersection(*child_sets)
        node.fitch_set = inter if inter else set.union(*child_sets)

    root = tree.seed_node
    root.fitch_final = (
        next(iter(root.fitch_set)) if len(root.fitch_set) == 1 else ambiguity_resolution
    )
    for node in tree.preorder_node_iter():
        if node is root:
            continue
        parent_state = node.parent_node.fitch_final
        if parent_state in node.fitch_set:
            node.fitch_final = parent_state
        elif len(node.fitch_set) == 1:
            node.fitch_final = next(iter(node.fitch_set))
        else:
            node.fitch_final = ambiguity_resolution

    return tree


def count_gain_events(tree: dendropy.Tree, background_attr: str = "oxa_final") -> dict:
    """Count 0->1 transitions on `tree` (after `fitch_reconstruct` has set
    `.fitch_final`) and classify each by the resistance background at the
    parent node (read from `background_attr`, set by a prior Fitch pass on
    the phenotype).

    Returns Total_gains, Gains_on_resistant_bg, Gains_on_susceptible_bg and
    Pct_post_resistance — the same fields and definitions as the original
    temporal_ordering.py output. `Pct_post_resistance` is NaN when no gain
    event was reconstructed; that is *undefined*, not zero, and callers must
    not read it as a failed stage.
    """
    gains_on_resistant_bg = 0
    gains_on_susceptible_bg = 0
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        parent_state = node.parent_node.fitch_final
        child_state = node.fitch_final
        if parent_state == 0 and child_state == 1:
            if not hasattr(node.parent_node, background_attr):
                raise CladeInputError(
                    f"node is missing the background attribute '{background_attr}'; "
                    "reconstruct the phenotype first and copy it onto every node "
                    "before counting candidate gains."
                )
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

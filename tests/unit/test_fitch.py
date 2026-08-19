"""
Hand-verified Fitch-parsimony test on a 5-tip synthetic tree.

Tree: ((A,B),(C,(D,E)))
Phenotype (resistance): A=1, B=1, C=0, D=0, E=0
Candidate: A=1, B=0, C=0, D=0, E=1

Worked by hand (see PR/commit description for the full derivation):
- Candidate gains at leaf A (parent clade (A,B), phenotype-positive background) -> resistant-background gain
- Candidate gains at leaf E (parent clade (D,E), phenotype-negative background) -> susceptible-background gain
Expected: Total_gains=2, Gains_on_resistant_bg=1, Gains_on_susceptible_bg=1, Pct_post_resistance=0.5
"""
import os

import dendropy
import pytest

from clade.phylogeny.fitch import count_gain_events, fitch_reconstruct

FIXTURE = os.path.join(os.path.dirname(__file__), "..", "fixtures", "tiny_tree.nwk")


@pytest.fixture
def tree():
    return dendropy.Tree.get(path=FIXTURE, schema="newick", preserve_underscores=True)


def test_phenotype_reconstruction(tree):
    pheno_map = {"A": 1, "B": 1, "C": 0, "D": 0, "E": 0}
    tree = fitch_reconstruct(tree, pheno_map)
    finals = {leaf.taxon.label: leaf.fitch_final for leaf in tree.leaf_node_iter()}
    assert finals == pheno_map  # leaves with known values must be preserved exactly


def test_candidate_gain_events_against_resistance_background(tree):
    pheno_map = {"A": 1, "B": 1, "C": 0, "D": 0, "E": 0}
    candidate_map = {"A": 1, "B": 0, "C": 0, "D": 0, "E": 1}

    tree = fitch_reconstruct(tree, pheno_map)
    for node in tree.preorder_node_iter():
        node.oxa_final = node.fitch_final

    tree2 = fitch_reconstruct(tree, candidate_map)
    stats = count_gain_events(tree2, background_attr="oxa_final")

    assert stats["Total_gains"] == 2
    assert stats["Gains_on_resistant_bg"] == 1
    assert stats["Gains_on_susceptible_bg"] == 1
    assert stats["Pct_post_resistance"] == pytest.approx(0.5)


def test_no_gains_when_candidate_is_ancestral(tree):
    """A candidate present at every tip has zero 0->1 transitions."""
    pheno_map = {"A": 1, "B": 1, "C": 0, "D": 0, "E": 0}
    candidate_map = {"A": 1, "B": 1, "C": 1, "D": 1, "E": 1}

    tree = fitch_reconstruct(tree, pheno_map)
    for node in tree.preorder_node_iter():
        node.oxa_final = node.fitch_final

    tree2 = fitch_reconstruct(tree, candidate_map)
    stats = count_gain_events(tree2, background_attr="oxa_final")
    assert stats["Total_gains"] == 0

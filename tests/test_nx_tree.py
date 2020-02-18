import pytest
from splitp import *


@pytest.fixture(scope='class')
def get_trees():
    matrix_1 = [[0.95, 0.05 / 3, 0.05 / 3, 0.05 / 3],
                [0.05 / 3, 0.95, 0.05 / 3, 0.05 / 3],
                [0.05 / 3, 0.05 / 3, 0.95, 0.05 / 3],
                [0.05 / 3, 0.05 / 3, 0.05 / 3, 0.95]]
    trees = []
    t1 = NXTree("(((A,B),C),(D,(E,F)));")
    t1.reassign_all_transition_matrices(np.array(matrix_1))
    trees.append(t1)
    return trees


def test_trivial_parsimony(get_trees):
    """ Testing that all trivial splits have parsimony of 1 """
    for tree in get_trees:
        splits = generate_all_splits(tree.get_num_taxa(), True, True)
        scores = [tree.parsimony_score(s) for s in splits]
        assert all([s == 1 for s in scores])


def test_subflats_equal(get_trees):
    for tree in get_trees:
        splits = generate_all_splits(tree.get_num_taxa(), trivial=False)
        pattern_probs = tree.get_pattern_probabilities()
        for sp in splits[::6]:
            flat = tree.flattening(sp, pattern_probs)
            sf1 = tree.subflattening(tree.transformed_flattening(flat))
            sf2 = tree.subflattening_alt(flat)
            assert np.all(np.isclose(sf1, sf2))

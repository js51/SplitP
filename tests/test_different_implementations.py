import pytest
from splitp import *

def test_both_implementations_give_same_probabilities():
    """ Code the same tree under both implementations, and check that the same
    site-pattern probability distribution is produced"""
    # 6 Taxa Tree (for old tree implementation)
    initial_dist = [0.25, 0.25, 0.25, 0.25]
    matrix = [[0.95, 0.05 / 3, 0.05 / 3, 0.05 / 3],
              [0.05 / 3, 0.95, 0.05 / 3, 0.05 / 3],
              [0.05 / 3, 0.05 / 3, 0.95, 0.05 / 3],
              [0.05 / 3, 0.05 / 3, 0.05 / 3, 0.95]]
    T_old = old_tree(11, 4)
    T_old.initDist = initial_dist
    T_old.addNode(node(None, 0))
    T_old.addNode(node(matrix, 1), 0)
    T_old.addNode(node(matrix, 2), 0)
    T_old.addNode(node(matrix, 3), 1)
    T_old.addNode(node(matrix, 7), 1)
    T_old.addNode(node(matrix, 8), 2)
    T_old.addNode(node(matrix, 4), 2)
    T_old.addNode(node(matrix, 5), 3)
    T_old.addNode(node(matrix, 6), 3)
    T_old.addNode(node(matrix, 9), 4)
    T_old.addNode(node(matrix, 10), 4)
    #################################

    T_new = NXTree("(((A,B),C),(D,(E,F)));")
    T_new.reassign_all_transition_matrices(np.array(matrix))

    pattern_probs_old = T_old.getLikelihoods()
    pattern_probs_new = T_new.get_pattern_probabilities()

    assert pattern_probs_old.equals(pattern_probs_new)

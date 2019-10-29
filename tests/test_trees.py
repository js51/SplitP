import pytest
from SplitP import *

@pytest.fixture(scope='class')
def get_trees():
    matrix = [[0.95, 0.05],[0.05, 0.95]]
    initial_dist = [0.50, 0.50]
    ### Tree 1 (6 Taxon binary) ###
    T6_2 = tree(11, 2)
    T6_2.initDist = initial_dist
    T6_2.addNode(node(None, 0))
    T6_2.addNode(node(matrix, 1),  0)
    T6_2.addNode(node([[1,0],[0,1]], 2),  0)
    T6_2.addNode(node(matrix, 3),  1)
    T6_2.addNode(node(matrix, 7),  1)
    T6_2.addNode(node(matrix, 8),  2)
    T6_2.addNode(node(matrix, 4),  2)
    T6_2.addNode(node(matrix, 5),  3)
    T6_2.addNode(node(matrix, 6),  3)
    T6_2.addNode(node(matrix, 9),  4)
    T6_2.addNode(node(matrix, 10), 4)
    return [T6_2]

def test_trivial_parsimony(get_trees):
    """ Testing that all trivial splits have parsimony of 1 """
    for tree in get_trees:
        splits = generateAllSplits(tree.getNumTaxa(), True, True)
        scores = [tree.getParsimony(s) for s in splits]
        assert all([s == 1 for s in scores])


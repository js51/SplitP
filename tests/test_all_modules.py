from splitp.nx_tree import NXTree
from splitp.tree_helper_functions import all_splits


import pytest
from splitp import *

@pytest.fixture(scope='class')
def test_cases():
    cases = [
        { 
            'newick_string'  : "(0:0.1,1:0.4);",
            'number_of_taxa' : 2,
            'true_splits'    : [], 
            'model'          : None, 
        },
        { 
            'newick_string'  : "((0:0.1,1:0.3),2:0.03);",
            'number_of_taxa' : 3,
            'true_splits'    : [], 
            'model'          : None, 
        },
        { 
            'newick_string'  : "((0,2),(1,3));",
            'number_of_taxa' : 4,
            'true_splits'    : ['02|13'], 
            'model'          : ('JC', 0.05), 
        }
    ]
    return cases


def test_trees(test_cases):
    """ Test a bunch of things about these trees """
    for case in test_cases:
        tree = NXTree(case['newick_string'])
        true_splits = list(tree.true_splits())
        assert true_splits == case['true_splits']
        false_splits = list(tree.false_splits())
        assert tree.get_num_taxa() == case['number_of_taxa']
        splits = list(all_splits(tree.get_num_taxa()))
        assert set(false_splits + true_splits) == set(splits)
import pytest
from splitp.parsers import newick, fasta
from splitp.phylogeny import Phylogeny
import splitp.phylogenetics

def test_true_splits(test_cases_trees):
    for case in test_cases_trees:
        if "newick_string" in case and "true_splits" in case:
            newick_string = case['newick_string']
            tree = Phylogeny(newick_string)
            assert set(tree.splits(as_strings=True)) == set(case['true_splits'])
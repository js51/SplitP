import pytest
from splitp.parsers import newick, fasta
from splitp.phylogeny import Phylogeny
import splitp.phylogenetics

def test_newick_string_from_splits(test_cases_trees):
    for test_case in test_cases_trees:
        if 'newick_string' in test_case and test_case['number_of_taxa'] > 3:
            true_newick_string = test_case['newick_string']
            phylogeny = Phylogeny(true_newick_string)
            splits = phylogeny.splits()
            newick_string = splitp.phylogenetics.newick_string_from_splits(splits)
            assert newick_string == true_newick_string
import pytest
from splitp.parsers import newick, fasta
import networkx as nx
from splitp.phylogeny import Phylogeny
import splitp.phylogenetics

def test_newick_string_from_splits(test_cases_trees):
    for test_case in test_cases_trees:
        if 'newick_string' in test_case and test_case['number_of_taxa'] > 3:
            true_newick_string = test_case['newick_string']
            phylogeny = Phylogeny(true_newick_string)
            splits = list(phylogeny.splits())
            newick_string = splitp.phylogenetics.newick_string_from_splits(splits)
            # Now check if the two Phylogenies are isomorphic
            ph1 = Phylogeny(newick.strip_newick(newick_string)).unrooted_networkx_graph()
            ph2 = Phylogeny(newick.strip_newick(true_newick_string)).unrooted_networkx_graph()
            print(ph1.edges())
            print(ph2.edges())
            assert nx.is_isomorphic(ph1, ph2)
""" A collection of tree reconstruction methods. 
    Each one of these methods accepts a sequence alignment in the form of a dictionary, 
    and returns a collection of splits which have been chosen as `most likely' to be true.
    Some methods will return additional information. Not all methods guarantee that the splits returned will be compatible.
"""

import splitp as sp
from splitp import tree_helper_functions as hf

def erickson_SVD():
    pass

def all_split_scores():
    pass

def split_tree_parsimony(alignment):
    alignment_dict = {}
    for table_pattern, value in alignment.itertuples(index=False, name=None):
        alignment_dict[table_pattern] = value
    num_taxa = len(list(alignment_dict.keys())[0]) # Length of first pattern
    all_splits = list(hf.all_splits(num_taxa))
    scores = {split : 0 for split in all_splits}
    for split in all_splits:
        newick_string = []
        for part in split.split('|'):
            newick_string.append(f'({",".join(c for c in part)})')
        newick_string = f"({newick_string[0]},{newick_string[1]});"
        split_tree = sp.NXTree(newick_string, taxa_ordering='sorted')
        for pattern, value in alignment_dict.items():
            scores[split] += value * split_tree.hartigan_algorithm(pattern)
    return scores
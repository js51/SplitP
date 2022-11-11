""" A collection of tree reconstruction methods. 
    Each one of these methods accepts a sequence alignment in the form of a dictionary, 
    and returns a collection of splits which have been chosen as `most likely' to be true.
    Some methods will return additional information. Not all methods guarantee that the splits returned will be compatible.
"""

from build.lib.splitp.nx_tree import NXTree
import splitp as sp
from splitp import tree_helper_functions as hf
import numpy as np

def erickson_SVD(alignment, method=sp.Method.flattenings):
    num_taxa = len(list(alignment.keys())[0]) # Length of first pattern
    taxa = [str(i) for i in range(num_taxa)]
    tree = NXTree.dummy_tree(taxa=taxa)
    twoSplits = tree.all_splits(size=2)
    minimum_score = 1
    chosen_split = None
    for split in twoSplits:
        subflat = tree.signed_sum_subflattening(split, alignment)
        score = tree.split_score(subflat)
        if score < minimum_score:
            minimum_score = score
            chosen_split = split
    chosen_pair = min(chosen_split)
    print(chosen_pair)
    

def all_split_scores():
    pass

def split_tree_parsimony(alignment, splits=None):
    if type(alignment) is dict:
        alignment_dict = alignment
    else:
        alignment_dict = {}
        for table_pattern, value in alignment.itertuples(index=False, name=None):
            alignment_dict[table_pattern] = value
    num_taxa = len(list(alignment_dict.keys())[0]) # Length of first pattern
    all_splits = list(hf.all_splits(num_taxa)) if splits==None else splits
    scores = {split : 0 for split in all_splits}
    for split in all_splits:
        newick_string = []
        for part in split.split('|'):
            newick_string.append(f'({",".join(c for c in part)})')
        newick_string = f"({newick_string[0]},{newick_string[1]});"
        split_tree = sp.NXTree(newick_string, taxa_ordering='sorted')
        for pattern, value in alignment_dict.items():
            scores[split] += value * (split_tree.hartigan_algorithm(pattern)/(num_taxa-1))
    return scores

def euclidean_split_distance(alignment, splits):
    print("assuming sorted taxa")
    states = ('A', 'G', 'C', 'T')
    alignment_dict = {}
    for table_pattern, value in alignment.itertuples(index=False, name=None):
        alignment_dict[table_pattern] = value
    num_taxa = len(list(alignment_dict.keys())[0]) # Length of first pattern
    all_splits = list(hf.all_splits(num_taxa)) if splits==None else splits
    scores = {split : 0 for split in all_splits}
    for split in all_splits:
        split_list = split.split('|')
        for pattern, value in alignment_dict.items():
            part_a = "".join(pattern[int(s, base=num_taxa+1)] for s in split_list[0])
            part_b = "".join(pattern[int(s, base=num_taxa+1)] for s in split_list[1])
            vec_a = np.array([ part_a.count(state) for state in states ])
            vec_b = np.array([ part_b.count(state) for state in states ])
            vec_a = vec_a / np.linalg.norm(vec_a)
            vec_b = vec_b / np.linalg.norm(vec_b)
            scores[split] += value * (2-np.linalg.norm(vec_a - vec_b))/2
    return scores
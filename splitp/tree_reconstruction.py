""" A collection of tree reconstruction methods. 
    Each one of these methods accepts a sequence alignment in the form of a dictionary, 
    and returns a collection of splits which have been chosen as `most likely' to be true.
    Some methods will return additional information. Not all methods guarantee that the splits returned will be compatible.
"""

import splitp as sp
from splitp import tree_helper_functions as hf
import numpy as np
from itertools import combinations

def erickson_SVD(alignment, taxa=None, method=sp.Method.flattenings):
    all_scores = {}
    bundle = {}
    labels = {}

    def _erickstep(all_taxa, alignment):
        scores = {}
        for pair in combinations(taxa, 2):
            flat_pair = tuple(sorted(element for tup in pair for element in (tup if not isinstance(tup, str) else (tup,))))
            other = tuple(sorted(element for tup in all_taxa for element in (tup if not isinstance(tup, str) else (tup,)) if element not in flat_pair))
            split = (flat_pair, other)
            try: 
                score = all_scores[split]
            except KeyError:
                if method == sp.Method.flattenings:
                    xflat = tree.reduced_sparse_flattening(split, alignment)
                elif method == sp.Method.subflattenings:
                    xflat = tree.fast_signed_sum_subflattening(split, alignment, bundle=bundle, labels=labels)
                score = tree.split_score(xflat)
                all_scores[split] = score
            scores[pair] = (pair, split, score)
        best_pair, best_split, best_score = min(scores.values(), key=lambda x: x[2])
        return best_pair, best_split, best_score

    num_taxa = len(list(alignment.keys())[0]) # Length of first pattern
    if taxa is None: 
        taxa = [str(np.base_repr(i, base=max(i+1,2))) if num_taxa <= 36 else f't{str(i)}' for i in range(num_taxa)]
    tree = sp.NXTree.dummy_tree(taxa=taxa)
    true_splits = []
    while len(true_splits) < num_taxa - 2:
        best_pair, best_split, best_score = _erickstep(taxa, alignment)
        true_splits.append(tuple(sorted(best_split)))
        taxa = tuple([element for element in taxa if (element not in best_split[0] and not set(element).issubset(best_split[0]))] + [best_split[0]])
    return true_splits

def newick_string_from_splits(splits):
    def _consolidate(tup, smaller_halves):
        if len(tup) == 2: return tup
        if len(tup) == 1: return tup[0]
        if isinstance(tup, str): return tup
        for smaller_half in smaller_halves:
            if set(smaller_half).issubset(tup) and len(smaller_half) < len(tup):
                # then the consolidation is made up of the smaller half and what is left over
                left_over = tuple(set(tup).difference(set(smaller_half)))
                try: 
                    smaller_halves.remove(smaller_half)
                    smaller_halves.remove(left_over)
                except ValueError:
                    pass
                return tuple((_consolidate(left_over, smaller_halves), _consolidate(smaller_half, smaller_halves)))

    splits = iter(sorted(splits, key=lambda x: min(len(x[0]), len(x[1])), reverse=True))
    first_split = next(splits)
    smaller_halves = [min(split, key=len) for split in splits]
    consolidated_split = tuple(_consolidate(half, smaller_halves) for half in first_split)
    return str(consolidated_split).replace("'", "").replace(" ", "") + ";"

def tree_from_splits(splits):
    return sp.NXTree(newick_string_from_splits(splits))

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

def neighbor_joining(distance_matrix):
    D = distance_matrix
    # Create the Q matrix
    Q = np.zeros((len(D), len(D)))
    for i in range(len(D)):
        for j in range(len(D)):
            if i != j:
                Q[i][j] = (len(D)-2)*D[i][j] - sum(D[i]) - sum(D[j])
    # Find the minimum value in the Q matrix
    min_i, min_j = np.unravel_index(np.argmin(Q, axis=None), Q.shape)

def JC_corrected_distance_matrix(alignment):
    num_taxa = len(list(alignment.keys())[0]) # Length of first pattern
    taxa = [str(i) for i in range(num_taxa)]
    distance_matrix = [ [0 for _ in range(num_taxa)] for _ in range(num_taxa)]
    for i in range(num_taxa):
        for j in range(i+1, num_taxa):
            for pattern, value in alignment.items():
                if pattern[i] != pattern[j]:
                    distance_matrix[i][j] += value
                    distance_matrix[j][i] += value
    for i in range(num_taxa):
        for j in range(i+1, num_taxa):
            distance_matrix[i][j] = -3.0/4.0 * np.log(1 - 4.0/3.0 * distance_matrix[i][j])
            distance_matrix[j][i] = distance_matrix[i][j]
    return distance_matrix

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
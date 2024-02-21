import numpy as np
from itertools import combinations
from numpy.random import shuffle
from math import floor


def split_balance(s, asTuple=False):
    """Returns a string formatted 'X|X' which describes the balance of a given split string"""
    s = s.split("|")
    if not asTuple:
        return str(len(s[0])) + "|" + str(len(s[1]))
    else:
        return (len(s[0]), len(s[1]))

def format_split(tree, split):
    if isinstance(split, str):
        return split  # If already a string, just send it back
    if len(split[0]) + len(split[1]) > 35:
        raise ValueError(
            "Cannot produce string format for split with more than 35 taxa."
        )
    if not all(len(taxon) == 1 for taxon in tree.get_taxa()):
        raise ValueError("Cannot produce string format for split with taxa name of length > 1.")
    else:
        return f'{"".join(split[0])}|{"".join(split[1])}'

def all_splits(tree, trivial=False, size=None, randomise=False, string_format=False):
    taxa = tree.taxa
    if string_format and len(taxa) > 35:
        raise ValueError(
            "Cannot generate splits for more than 35 taxa in string format. Use string_format=False."
        )
    if size is not None:
        sizes_to_do = [size]
    else:
        sizes_to_do = list(
            range(1 if trivial else 2, floor(len(taxa) / 2) + 1)
        )
    for bal in sizes_to_do:
        even_split = bal == len(taxa) / 2
        if not even_split:
            combos = combinations(taxa, bal)
        else: # We don't want to double up by selecting 012 and then 345
            combos = combinations(taxa[1:], bal - 1)  # Don't select 0
        if randomise:
            combos = list(combos)
            shuffle(combos)
        for left_taxa in combos:
            if even_split:
                left_taxa = (taxa[0],) + left_taxa
            right_taxa = sorted(tuple(set(taxa) - set(left_taxa)), key=taxa.index)
            left_taxa = sorted(left_taxa, key=taxa.index)
            left, right = tuple(left_taxa), tuple(right_taxa) 
            if (taxa[0] in right):  # Put taxa[0] on left side of split
                left, right = right, left
            if string_format:
                yield format_split(tree, (left, right))
            else:
                yield (left, right)
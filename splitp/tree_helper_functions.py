import numpy as np
import pandas as pd
import collections
from scipy.sparse import issparse

def balanced_newick_tree(num_taxa):
    if num_taxa%2 != 0:
        raise ValueError("There is no balanced tree on {num_taxa} taxa. Please specify an even number.")
    from math import floor
    def _balanced_newick_subtree(nt, left=False):
        if nt == 2:
            return "(_,_)"
        elif nt==3:
            return "((_,_),_)" if left else "(_,(_,_))"
        else:
            if nt%2==0:
                return f"({_balanced_newick_subtree(nt/2, True)},{_balanced_newick_subtree(nt/2)})"
            else: 
                return f"({_balanced_newick_subtree(floor(nt/2) + int(left), True)},{_balanced_newick_subtree(floor(nt/2) + int(not left))})"
    newick_string = f"({_balanced_newick_subtree(num_taxa/2, True)},{_balanced_newick_subtree(num_taxa/2)});"
    for i in range(0, num_taxa):
        newick_string = newick_string.replace('_', str(np.base_repr(i, base=i+3)), 1)
    return newick_string

def get_balance(s, asTuple=False):
    """Returns a string formatted 'X|X' which describes the balance of a given split string"""
    s = s.split("|")
    if not asTuple:
        return str(len(s[0])) + '|' + str(len(s[1]))
    else:
        return (len(s[0]), len(s[1]))

def is_sparse(matrix):
    return issparse(matrix)

def frob_norm(matrix, data_table=None):
    """Calculates the Frobenius Norm for a given matrix"""
    if data_table is not None:
        return sum(val**2 for _, val in data_table.itertuples(index=False))**(1/2)
    if is_sparse(matrix):
        return sum(matrix[i,j]**2 for i, j in zip(*matrix.nonzero()))**(1/2)
    else:
        return np.sqrt(sum(val**2 for val in np.nditer(matrix)))
        
def make_substitution_matrix(subs_prob, k):
    matrix = []
    for i in range(k):
        matrix.append([1-subs_prob if j==i else subs_prob/(k-1) for j in range(k)])
    return matrix

def all_splits(num_taxa, trivial=False, only_balance=None, randomise=False):
    """Generates all splits as string-representations

    Args:
        num_taxa: The number of taxa on the tree.
        trivial: Whether or not to calculate trivial splits, default True.
        only_trivial: Whether to ONLY create trivial splits, default False.

    Returns:
        A list of string-representations of splits (using '|'-notation)
    """
    k = only_balance
    n = num_taxa
    taxa_string = "".join(np.base_repr(i, base=n) for i in range(n))
    r = 0 if trivial else 1
    loop_over = range(r, 2**(n-1) - r)
    if randomise:
        import random
        loop_over = [i for i in loop_over]
        random.shuffle(loop_over)
    for i in loop_over:
        template = format(i, f'0{n}b')
        if not only_balance or sum(int(b) for b in template) in [only_balance, n-only_balance]:
            if r < sum(int(b) for b in template) < n-r:
                left  = ""
                right = ""
                for t, b in zip(taxa_string, template):
                    if b == '0':
                        left += t
                    else:
                        right += t
                if '0' in right:
                    left, right = right, left
                yield f'{left}|{right}'
###
# Functions for reading in a sequence alignment

def read_alignment_from_file(pathToFile, check=True):
    file = open(pathToFile, 'r')
    alignment = collections.OrderedDict()  # Keeps insertion order
    currentSeq = ""
    for line in file:
        if '>' in line:
            currentSeq = line.replace('>', '').replace('\n', '')
        else:
            if currentSeq in alignment:
                alignment[currentSeq] += line.replace('\n', '')
            else:
                alignment[currentSeq] = line.replace('\n', '')
    if check:
        return alignment if __valid_alignment(alignment) else None
    else:
        return alignment


def __valid_alignment(alignmentDict):
    ordered = isinstance(alignmentDict, collections.OrderedDict)  # Must be ordered to avoid incorrect patterns
    sameLength = len(set(len(value) for key, value in alignmentDict.items())) == 1  # Must be same length
    return sameLength and ordered


def get_pattern_counts(alignment, asNumpyArray=False):
    patterns = {}
    sequences = list(alignment.values())
    sequenceLength = len(sequences[0])
    usableSequenceLength = 0
    for i in range(sequenceLength):
        patternAtSitei = ''.join(s[i] for s in sequences).upper()
        if all(c in ('A', 'C', 'G', 'T') for c in patternAtSitei):  # Need to handle U (let U = T)
            usableSequenceLength += 1
            if patternAtSitei in patterns:
                patterns[patternAtSitei] += 1
            else:
                patterns[patternAtSitei] = 1

    if asNumpyArray:
        patterns = np.array([[key, val] for key, val in patterns.items()])

    return patterns, usableSequenceLength


def pattern_counts_to_probs(patterns, seqLen):
    newCounts = [float(float(count) / seqLen) for count in patterns[:, 1]]
    patternList = patterns[:, 0]
    return np.array([patternList, newCounts]).transpose()


def pattern_probs_from_alignment(pathToFile, check=True):
    alignment = read_alignment_from_file(pathToFile, check)
    counts, sequenceLength = get_pattern_counts(alignment, True)
    probs = pattern_counts_to_probs(counts, sequenceLength)
    probs = pd.DataFrame(probs, index=probs[:, 0])
    probs[[1]] = probs[[1]].astype("float")
    return probs, len(alignment), sequenceLength


def scores_to_weights(scores):
    scores = [s if s != 0 else 10 ** (-100) for s in scores]
    inverses = [1 / s for s in scores]
    weights = [i / sum(inverses) for i in inverses]
    return weights


def split_equivalence_classes(splits, group):
    """Breaks up list of splits into equivalence classes under the canonical action of the specified group."""
    def __action(perm, split):
        index_of_sep = split.index('|')
        split = split.replace('|', '')
        split = split = ''.join(str(perm(int(s))) for s in split)
        split = sorted(["".join(sorted(split[:index_of_sep])), "".join(sorted(split[index_of_sep:]))])
        split = f'{split[0]}|{split[1]}'
        return split
    def __orbit(split, group):
        orbit = set()
        for g in group.elements:
            orbit.add(__action(g, split))
        return frozenset(orbit)
    def __orbits(splits, group):
        orbits = set()
        for split in splits:
            orbits.add(__orbit(split, group))
        return sorted(list(set(orbit) for orbit in orbits), key=len)
    return __orbits(splits, group)

def normalised(scores):
    scores = scores.copy()
    if type(scores) is dict:
        minimum = min(scores.values())
        maximum = max(scores.values())
        for key, value in scores.items():
            scores[key] = (value - minimum)/(maximum - minimum)
    elif type(scores) is list:
        minimum = min(scores)
        maximum = max(scores)
        scores = [(score-minimum)/(maximum-minimum) for score in scores]
    return scores

def check_splits_representatives(split_reps, orbits):
    represented_orbits = set()
    for orbit in orbits:
        for split in split_reps:
            if split in orbit:
                represented_orbits.add(frozenset(orbit))
    return len(represented_orbits) == len(orbits)

# Function for scaling the S part of the hard-coded H matrix
def scaled_h_matrix(_lambda):
    if _lambda == None:
        return None, None
    if _lambda == "identity":
        return np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 1, 1, 1]]), "identity"
    H = np.matrix([[1, -1], [1, 1]])
    H = np.kron(H, H)
    S = np.empty((0, 4))
    for row in H:
        if len(set(np.asarray(row)[0])) != 1:
            S = np.append(S, row, axis=0)
    S = _lambda * S
    H = np.append(S, [[1, 1, 1, 1]], axis=0)
    return (H, _lambda)

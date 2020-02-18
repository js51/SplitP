import numpy as np
import copy as cp
import pandas as pd
import itertools
import collections
from itertools import permutations

def get_balance(s, asTuple=False):
    """Returns a string formatted 'X|X' which describes the balance of a given split string"""
    s = s.split("|")
    if not asTuple:
        return str(len(s[0])) + '|' + str(len(s[1]))
    else:
        return (len(s[0]), len(s[1]))

def frob_norm(m):
    """Calculates the Frobenius Norm for a given numpy array"""
    norm = 0
    for i in np.nditer(m):
        norm += i ** 2
    norm = np.sqrt(norm)
    return norm

def make_substitution_matrix(subs_prob, k):
    matrix = []
    for i in range(k):
        matrix.append([1-subs_prob if j==i else subs_prob/(k-1) for j in range(k)])
    return matrix

def generate_all_splits(num_taxa, trivial=True, onlyTrivial=False, onlyBalance=-1):
    """Generates all splits as string-representations

    Args:
        num_taxa: The number of taxa on the tree.
        trivial: Whether or not to calculate trivial splits, default True.
        onlyTrivial: Whether to ONLY create trivial splits, default False.

    Returns:
        A list of string-representations of splits (using '|'-notation)
    """
    s = [str(i) for i in range(num_taxa)]
    s.append('|')
    perms = [''.join(p) for p in permutations(s)]
    perms = set(perms)
    aList = []
    complete_set = set()
    for p in perms:
        halves = p.split('|')
        limit = 0
        if not trivial: limit = 1
        if onlyBalance == -1:
            criterion = (len(halves[0]) > limit) and (len(halves[1]) > limit)
        else:
            criterion = (len(halves[0]) == onlyBalance) or (len(halves[1]) == onlyBalance)
        if criterion:
            a = ''.join(sorted(halves[0]))
            b = ''.join(sorted(halves[1]))
            if a not in aList and b not in aList:
                if '0' in a:
                    complete_set.add(a + '|' + b)
                else:
                    complete_set.add(b + '|' + a)
                aList.append(a)
                aList.append(b)

    if onlyTrivial:
        new_set = complete_set.copy()
        for s in complete_set:
            halves = s.split('|')
            if len(halves[0]) != 1 and len(halves[1]) != 1:
                new_set.remove(s)
        return list(new_set)

    return sorted(list(complete_set))


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

        file = open("output/patternCounts", 'w')
    for p in patterns:
        file.write(str(p))
    file.close()

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


####

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

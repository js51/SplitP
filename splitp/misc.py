import numpy as np
import pandas as pd
import collections
from math import floor
from scipy.sparse import issparse
from splitp import enums
from warnings import warn
from math import exp


def make_substitution_matrix(subs_prob, k):
    matrix = []
    for i in range(k):
        matrix.append(
            [1 - subs_prob if j == i else subs_prob /
                (k - 1) for j in range(k)]
        )
    return matrix


def scores_to_weights(scores):
    scores = [s if s != 0 else 10 ** (-100) for s in scores]
    inverses = [1 / s for s in scores]
    weights = [i / sum(inverses) for i in inverses]
    return weights


def normalised(scores):
    scores = scores.copy()
    if type(scores) is dict:
        minimum = min(scores.values())
        maximum = max(scores.values())
        for key, value in scores.items():
            scores[key] = (value - minimum) / (maximum - minimum)
    elif type(scores) is list:
        minimum = min(scores)
        maximum = max(scores)
        scores = [(score - minimum) / (maximum - minimum) for score in scores]
    return scores


def build_JC_matrix(tree, l):
        matrix = [[0 for _ in range(tree.num_bases)] for n in range(tree.num_bases)]
        for r, row in enumerate(matrix):
            for c, _ in enumerate(row):
                matrix[r][c] = (1 / 4) + (3 / 4) * exp((-4 * l) / 3) if r == c else (1 / 4) - (1 / 4) * exp(
                    (-4 * l) / 3)
        return np.array(matrix).T
    

def build_K2ST_matrix(tree, transition, transversion):
    if tree.num_bases != 4:
        warn(f"K2ST matrices are 4x4 but your model has {tree.num_bases} states!" )
    if transversion > transition: 
        warn(f"transitions are known to be more likely than transversions!")
    purines = ('A', 'G')
    pyrimidines = ('C', 'T')
    matrix = [[0 for i in range(tree.num_bases)] for n in range(tree.num_bases)]
    for r, row in enumerate(matrix):
        from_state = tree.state_space[r]
        for c, _ in enumerate(row):
            to_state = tree.state_space[c]
            if from_state == to_state:
                # No change
                matrix[r][c] = 1-(transition+2*transversion)
            elif from_state in purines and to_state in purines:
                matrix[r][c] = transition 
            elif from_state in pyrimidines and to_state in pyrimidines:
                matrix[r][c] = transition
            else:
                matrix[r][c] = transversion
    return np.array(matrix).T


def rate_matrix(tree, model):
    def _JC_rate_matrix(mutation_rate=None):
        if mutation_rate:
            a = mutation_rate
            return [[-3*a, a, a, a],
                    [a, -3*a, a, a],
                    [a, a, -3*a, a],
                    [a, a, a, -3*a]]
    def _K2ST_rate_matrix(rate_transition=None, rate_transversion=None, ratio=None):
        if tree.num_bases != 4:
            warn(f"K2ST matrices are 4x4 but your model has {tree.num_bases} states!" )
        purines = ('A', 'G')
        pyrimidines = ('C', 'T')
        matrix = [[0 for i in range(tree.num_bases)] for n in range(tree.num_bases)]
        if rate_transition and rate_transversion:
            transition = rate_transition
            transversion = rate_transversion
            if transversion > transition: 
                warn(f"transitions are known to be more likely than transversions!")
        elif ratio:
            transition = ratio
            transversion = 1
        for r, row in enumerate(matrix):
            from_state = tree.state_space[r]
            for c, _ in enumerate(row):
                to_state = tree.state_space[c]
                if from_state == to_state:
                    # No change
                    matrix[r][c] = -(transition+2*transversion)
                elif from_state in purines and to_state in purines:
                    matrix[r][c] = transition 
                elif from_state in pyrimidines and to_state in pyrimidines:
                    matrix[r][c] = transition
                else:
                    matrix[r][c] = transversion
        return np.array(matrix).T
    if   model is enums.Model.JC:   return _JC_rate_matrix
    elif model is enums.Model.K2ST: return _K2ST_rate_matrix


def scale_TR_rate_matrix(tree, Q, return_scale_factor=False):
    scale_factor = 1/(-sum((1/4)*Q[i][i] for i in range(4)))
    if return_scale_factor:
        return scale_factor
    else:
        return scale_factor*np.array(Q)

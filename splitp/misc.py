import numpy as np
import pandas as pd
import collections
from math import floor
from scipy.sparse import issparse


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

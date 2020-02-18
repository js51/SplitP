# imports
from Tree import tree
from Tree import node
import trees
import pandas as pd
import TreeHelperFunctions as hf
import matplotlib.pyplot as plt
import numpy as np

trees = trees.trees

trials = 1000

lambdas = [0,1,2]
matNames = [None, 1, "identity"]
HMats = [hf.scaled_h_matrix(i) for i in matNames]

plt.figure(0)

key = None
value = None
for k, v in trees.items():
    if k == 6:
        key = k
        value = v

# Constants for each tree
splits = hf.generate_all_splits(key, trivial=False)
sitePatternProbs = value[0].get_pattern_probabilities()
numRuns = 1000  # number of lines to plot

scoreDiffsCollection = []
scoreDiffs2Collection = []

for i in range(numRuns):
    print("Run " + str(i) + " of " + str(numRuns))
    print(str(key), "taxon tree")
    scoreDiffs = []
    scoreDiffs2 = []
    empericalProbs = value[0].draw_from_multinomial(sitePatternProbs, trials)

    for H in HMats:
        # Preparation
        splitScores = []
        VecF = None
        L, R = None, None

        for split in splits:
            # print("Split:", split)
            F = value[0].flattening(split, empericalProbs)
            if H[1] != None:
                SF = value[0].subflattening_alt(F, S=H[0], returnLRMats=False)
            else:
                SF = F
            score = value[0].split_score(SF)
            splitScores.append(score)
        trueScores = sorted(splitScores)[0:value[1]]
        falseScores = sorted(splitScores)[value[1]:]
        trueFalseSep = falseScores[0] - trueScores[-1]
        minMaxSep = sorted(splitScores)[-1] / sorted(splitScores)[0]
        scoreDiffs.append(trueFalseSep)
        scoreDiffs2.append(trueFalseSep / minMaxSep)

    scoreDiffs2Collection.append(scoreDiffs2)
    scoreDiffsCollection.append(scoreDiffs)
    plt.plot(lambdas, scoreDiffs)
    plt.plot(lambdas, scoreDiffs2)


_scoreDiffs1 = [[] for n in lambdas]  # a list of scores for each lambda
for _diffs1 in scoreDiffsCollection:  # loop over runs
    for i in range(len(_diffs1)):
        _scoreDiffs1[i].append(_diffs1[i])
_scoreDiffs1Means = [np.mean(lamScores) for lamScores in _scoreDiffs1]
_scoreDiffs1StdevL = [_scoreDiffs1Means[i] - np.std(_scoreDiffs1[i]) for i in range(len(_scoreDiffs1))]
_scoreDiffs1StdevU = [_scoreDiffs1Means[i] + np.std(_scoreDiffs1[i]) for i in range(len(_scoreDiffs1))]

_scoreDiffs2 = [[] for n in lambdas]  # a list of scores for each lambda
for _diffs2 in scoreDiffs2Collection:  # loop over runs
    for i in range(len(_diffs2)):
        _scoreDiffs2[i].append(_diffs2[i])
_scoreDiffs2Means = [np.mean(lamScores) for lamScores in _scoreDiffs2]
_scoreDiffs2StdevL = [_scoreDiffs2Means[i] - np.std(_scoreDiffs2[i]) for i in range(len(_scoreDiffs2))]
_scoreDiffs2StdevU = [_scoreDiffs2Means[i] + np.std(_scoreDiffs2[i]) for i in range(len(_scoreDiffs2))]

plt.figure()

plt.plot(lambdas, _scoreDiffs1Means, color='red')
plt.plot(lambdas, [0.017119617763233518 for i in range(len(lambdas))], color='gray')
plt.plot(lambdas, _scoreDiffs1StdevL, color='blue',alpha = 0.2)
plt.plot(lambdas, _scoreDiffs1StdevU, color='blue',alpha = 0.2)
plt.fill_between(lambdas, _scoreDiffs1StdevL, _scoreDiffs1StdevU, color='blue', alpha = 0.05)

plt.show()

plt.figure()

yerr = [ np.std(_scoreDiffs1[i]) for i in range(len(_scoreDiffs1)) ]
plt.errorbar(lambdas, _scoreDiffs1Means, yerr = yerr, marker='o', capsize=5)
plt.plot(lambdas, [0.017119617763233518 for i in range(len(lambdas))], color='gray')

plt.show()

plt.figure()
plt.scatter(lambdas, _scoreDiffs2Means)
plt.fill_between(lambdas, _scoreDiffs2StdevL, _scoreDiffs2StdevU, color='blue', alpha=0.05)
#

plt.show()

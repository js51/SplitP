# imports
from SplitP import *
import trees
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

trees = trees.trees
treeTuple = trees[6]
tree = treeTuple[0]
trueSplits = treeTuple[2]
numSpecies = tree.getNumTaxa()

sequenceLength = 1000
matNames = [n for n in range(-20, 20+1)]
matNames.remove(0)
matNames.append(-0.01)
matNames.append(0.01)
matNames.sort()
lambdas = [i for i in range(len(matNames))]
HMats = [hf.scaledHMatrix(i) for i in matNames]

all_splits = hf.generateAllSplits(numSpecies, trivial=False)
patternProbs = tree.getLikelihoods()
numRuns = 1000  # number of lines to plot

scoreDiffsCollection = []
scoreDiffs2Collection = []

finalResults = { matName : [[],[]] for matName in matNames }
singleRun = []
# { Matrix : [[qualityMeasure1s], [qualityMeasure2s]] }
for run in range(numRuns):
    print("Run " + str(run) + " of " + str(numRuns))
    empericalProbs = tree.drawFromMultinomialFast(patternProbs, sequenceLength)

    scoreLists = {mat : [] for mat in matNames}
    for split in all_splits:
        F = tree.flattening(split, empericalProbs)
        for i, H in enumerate(HMats):
            SF = tree.subFlatteningAltFast(F, S=H[0])
            score = tree.splitScore(SF)
            scoreLists[matNames[i]].append(score)

    for matName, scoreList in scoreLists.items():
        trueScores = []
        falseScores = []
        for s, split in enumerate(all_splits):
            if split in trueSplits:
                trueScores.append(scoreList[s])
            else:
                falseScores.append(scoreList[s])
        quality_measure_1 = min(falseScores) - max(trueScores)
        quality_measure_2 = quality_measure_1/(max(scoreList) - min(scoreList))
        finalResults[matName][0].append(quality_measure_1)
        finalResults[matName][1].append(quality_measure_2)

finalResultsTable = []
for matName, measurePair in finalResults.items():
    finalResultsTable.append([np.mean(measurePair[0]), np.std(measurePair[0]), np.mean(measurePair[1]), np.std(measurePair[1])])
finalResultsTable = pd.DataFrame(finalResultsTable, index=matNames)


f = open("SMatrixScaling6tax2.csv", "w")
f.write(pd.DataFrame(finalResultsTable).to_csv(index=False))
f.close()

nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        "axes.labelsize": 11,
        "font.size": 11,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
}
import seaborn as sns
import matplotlib as mpl
sns.set()
plt.style.use('default')
mpl.rcParams.update(nice_fonts)

with sns.color_palette("colorblind"):
    plt.figure(figsize=(4,3))
    plt.errorbar(matNames, list(finalResultsTable[0]), yerr = list(finalResultsTable[1]), markersize=4, marker='o', capsize=2.5, elinewidth=0.75)
    plt.xlabel("$\lambda$ for scaled matrix $S=H^{(\lambda)}$")
    plt.ylabel("Measure of performance 1")
    plt.tight_layout()
    plt.savefig('MOP1_LONG_6TAX.pdf')
    plt.show()

    plt.figure(figsize=(4,3))
    plt.errorbar(matNames, list(finalResultsTable[2]), yerr = list(finalResultsTable[3]), markersize=4, marker='o', capsize=2.5, elinewidth=0.75)
    plt.xlabel("$\lambda$ for scaled matrix $S=H^{(\lambda)}$")
    plt.ylabel("Measure of performance 2")
    plt.tight_layout()
    plt.savefig('MOP2_LONG_6TAX.pdf')
    plt.show()

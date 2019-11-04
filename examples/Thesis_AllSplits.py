import SplitP as sp
from SplitP import *
import trees
import matplotlib as mpl
import matplotlib.pyplot as plt
from more_itertools import sort_together
import numpy as np

treeTuple = trees.trees[6]
tree = treeTuple[0]

def sort(string):
    return ''.join(sorted(string))

# Settings
true_splits = treeTuple[2]
scalingFactors = [None, 1, 2, "H/2", "Alt Identity", "Identity"]
numRuns = 1000

def scaledHMatrix(lam):
    if lam == "Alt Identity":
        return (np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 1, 1, 1]]), "Alt Identity")
    elif lam == "Identity":
        return (np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), "Identity")
    else:
        return sp.scaledHMatrix(lam)

numSpecies = tree.getNumTaxa()
all_splits = generateAllSplits(numSpecies, trivial=False)
patternProbs = tree.getLikelihoods()

results = {str(lam) : {split : [] for split in all_splits} for lam in scalingFactors}

for i in range(numRuns):
    print("Run " + str(i))
    DTable = tree.drawFromMultinomial(patternProbs, 1000)
    for split in all_splits:
        print(split)
        F = tree.flattening(split, DTable)
        for lam in scalingFactors:
            if lam == None:
                print('Flattening')
                results[str(lam)][split].append(tree.splitScore(F))
            elif lam == "Identity":
                print('Identity Matrix')
                SF = tree.subFlattening(tree.transformedFlattening(F, S=np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])))
                score = tree.splitScore(SF)
                results[str(lam)][split].append(score)
            elif lam == "H/2":
                print('Full H matrix divided by 2')
                SF = tree.subFlattening(tree.transformedFlattening(F, S=(1/2)*scaledHMatrix(1)[0]))
                score = tree.splitScore(SF)
                results[str(lam)][split].append(score)
            else:
                print('Scaling by:', scaledHMatrix(lam)[1])
                SF = tree.subFlatteningAlt(F, scaledHMatrix(lam)[0])
                score = tree.splitScore(SF)
                results[str(lam)][split].append(score)


# Collect mean and variance
finalResults = {str(lam) : [[], [], []] for lam in scalingFactors}
for lam, splitDict in results.items():
    for split, scoreList in splitDict.items():
        meanScore = np.mean(scoreList)
        stddev = np.std(scoreList)
        finalResults[str(lam)][0].append(split)
        finalResults[str(lam)][1].append(meanScore)
        finalResults[str(lam)][2].append(stddev)

nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        "axes.labelsize": 11,
        "font.size": 11,
        "legend.fontsize": 9,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
}
import seaborn as sns
sns.set()
plt.style.use('default')
mpl.rcParams.update(nice_fonts)

flatScores = finalResults['None'][1]
with sns.color_palette("colorblind"):
    plt.figure(figsize=(6,3.5))
    for lam, triple in finalResults.items():
        triple.append(flatScores)
        triple = sort_together(triple, key_list=(len(triple) - 1,))
        if lam not in ['None', "Alt Identity", "Identity"]:
            plt.plot([x.replace('|', '\\textbar') for x in triple[0]], triple[1], marker='o', markersize=2, label="Flattening" if lam == 'None' else "$\lambda =\ $" + lam)
        else:
            plt.plot([x.replace('|', '\\textbar') for x in triple[0]], triple[1], marker='o', markersize=2, label="Flattening" if lam == 'None' else lam)

plt.legend()
plt.xlabel("Split", labelpad=10)
plt.ylabel("Split Score", labelpad=10)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('tryTreesPP.pdf')
plt.show()

from splitp import *
import trees
import matplotlib as mpl
import matplotlib.pyplot as plt
from more_itertools import sort_together
import numpy as np

treeTuple = trees.trees["T6B"]
tree = treeTuple[0]
subflattenings = ["Flattening", "Subflattening", (1,2), (2,1), (2,2)]


def sort(string):
    return ''.join(sorted(string))

# Settings
true_splits = treeTuple[2]
sequenceLengths = [100*n for n in range(1, 11)]
numRuns = 1000

numSpecies = tree.get_num_taxa()
all_splits = generateAllSplits(numSpecies, trivial=False)
patternProbs = tree.get_pattern_probabilities()

results = {str(subflat) : {split : [] for split in all_splits} for subflat in subflattenings}

for i in range(numRuns):
    print("Run " + str(i))
    DTable = tree.draw_from_multinomial(patternProbs, 1000)
    for split in all_splits:
        print(split)
        F = tree.flattening(split, DTable)
        TF = tree.transformed_flattening(F)
        for subflat in subflattenings:
            score = None
            if subflat == "Flattening":
                print('Flattening')
                score = tree.split_score(F)
            elif subflat == "Subflattening":
                print('Subflat')
                SF = tree.subflattening(TF)
                score = tree.split_score(SF)
            else:
                print('Subflattening: ', str(subflat))
                SF = tree.subflattening(TF, type=subflat)
                score = tree.split_score(SF)
            results[str(subflat)][split].append(score)


# Collect mean and variance
finalResults = {str(subflat) : [[], [], []] for subflat in subflattenings}
for subflat, splitDict in results.items():
    for split, scoreList in splitDict.items():
        meanScore = np.mean(scoreList)
        stddev = np.std(scoreList)
        finalResults[str(subflat)][0].append(split)
        finalResults[str(subflat)][1].append(meanScore)
        finalResults[str(subflat)][2].append(stddev)

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

flatScores = finalResults['Flattening'][1]
with sns.color_palette("colorblind"):
    plt.figure(figsize=(6,3.5))
    for subflat, triple in finalResults.items():
        triple.append(flatScores)
        triple = sort_together(triple, key_list=(len(triple) - 1,))
        plt.plot([x.replace('|', '\\textbar') for x in triple[0]], triple[1], marker='o', markersize=2, label=str(subflat))

plt.legend()
plt.xlabel("Split", labelpad=10)
plt.ylabel("Split Score", labelpad=10)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('GSFALLSPLITS_EP_6B.pdf') # GET THIS ONE
plt.show()

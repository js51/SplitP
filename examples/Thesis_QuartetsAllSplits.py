import SplitP as sp
import trees
import matplotlib as mpl
import matplotlib.pyplot as plt
from more_itertools import sort_together
import numpy as np

treeTuple = trees.trees['T4A']

def sort(string):
    return ''.join(sorted(string))

# Settings
true_splits = treeTuple[2]
scalingFactors = [None, 1, 2, 3, 10]
sequenceLengths = [100*n for n in range(1, 11)]
numRuns = 100
bigResults = []
theQuartets = [trees.trees['T4A'][0], trees.trees['T4B'][0], trees.trees['T4C'][0], trees.trees['T4D'][0]]

for t, tree in enumerate(theQuartets):

    finalResults = {str(lam): [0 for sequenceLength in sequenceLengths] for lam in scalingFactors}
    numSpecies = tree.getNumTaxa()
    all_splits = sp.generateAllSplits(numSpecies, trivial=False)
    patternProbs = tree.getLikelihoods()
    results = {str(lam): {split: [] for split in all_splits} for lam in scalingFactors}

    for i in range(numRuns):
        print("Run " + str(i))
        DTable = tree.drawFromMultinomialFast(patternProbs, 1000)
        for split in all_splits:
            print(split)
            F = tree.flattening(split, DTable)
            for lam in scalingFactors:
                if lam == None:
                    print('Flattening')
                    results[str(lam)][split].append(tree.splitScore(F))
                else:
                    print('Scaling by:', sp.scaledHMatrix(lam)[1])
                    SF = tree.subFlatteningAltFast(F, sp.scaledHMatrix(lam)[0])
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

    bigResults.append(finalResults)

for t, finalResults in enumerate(bigResults):
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
        plt.figure(figsize=(6,3))
        plt.title(str(t))
        for lam, triple in finalResults.items():
            triple.append(flatScores)
            triple = sort_together(triple, key_list=(len(triple) - 1,))
            plt.plot([x.replace('|', '\\textbar') for x in triple[0]], triple[1], marker='o', markersize=2, label="Flattening" if lam == 'None' else "$\lambda =\ $" + lam)


    plt.legend()
    plt.xlabel("Split", labelpad=10)
    plt.ylabel("Split Score", labelpad=10)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('6TaxaMultinomialScalings_1000seq_1000runs.pdf')
    plt.show()

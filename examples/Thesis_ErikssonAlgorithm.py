import pandas as pd
import SplitP as sp
import time
import trees
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt

treeTuple = trees.trees['T6B']
tree = treeTuple[0]

def sort(string):
    return ''.join(sorted(string))

# Settings
true_splits = treeTuple[2]
true_splits = treeTuple[2]
scalingFactors = [None, 1, 2, 3]
sequenceLengths = [100*n for n in range(1, 11)]
numRuns = 1000

finalResults = {str(lam) : [0 for sequenceLength in sequenceLengths] for lam in scalingFactors}
numSpecies = tree.getNumTaxa()
all_splits = sp.generateAllSplits(numSpecies, trivial=False)
patternProbs = tree.getLikelihoods()
for i, sequenceLength in enumerate(sequenceLengths):
    start = time.time()
    for run in range(numRuns):
        print("Run", run, "of", numRuns, end='')
        print(" (length: " + str(sequenceLength) + ")")
        empericalProbs = tree.drawFromMultinomialFast(patternProbs, sequenceLength)
        flattenings = {}
        for lam in scalingFactors:
            print("\t" + str(lam), end='')
            taxa = [str(i) for i in range(numSpecies)]
            chosenSplits = set()
            done = False
            for c in range(numSpecies, 3, -1):
                if not done:
                    results = dict()
                    for pair in itertools.combinations(taxa, 2):
                        split = ("".join(t for t in taxa if t not in pair[0]+pair[1]) , pair[0] + pair[1])
                        split = sort(split[0]) + '|' + sort(split[1]) if '0' in split[0] else sort(split[1]) + '|' + sort(split[0])
                        F = None
                        if split in flattenings:
                            F = flattenings[split]
                        else:
                            F = tree.flattening(split, empericalProbs)
                            flattenings[split] = F
                        if lam != None:
                            F = tree.subFlatteningAltFast(F, sp.scaledHMatrix(lam)[0])
                        results[(pair, split)] = tree.splitScore(F)
                    bestSplit = min(results, key=results.get)
                    if bestSplit[1] not in true_splits:
                        done = True
                    for p in bestSplit[0]:
                        taxa.remove(p)
                    taxa.append("".join(bestSplit[0]))
                    chosenSplits.add(bestSplit[1])
            if chosenSplits == true_splits:
                finalResults[str(lam)][i] += 1
                print("\t\t" + "Success")
            else:
                print("\t\t" + "Failure")
    end = time.time()
    print("Time taken:", end - start)

f = open("ErikksonAlgorithmResultsHARD.csv", "w")
f.write(pd.DataFrame(finalResults).to_csv(index=False))
f.close()

data = pd.read_csv("ErikksonAlgorithmResultsHARD.csv")
finalResults = {}
for column in data:
    finalResults[column] = data[column].to_list()

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
sns.set()
plt.style.use('default')
mpl.rcParams.update(nice_fonts)

with sns.color_palette("colorblind"):
    plt.figure(figsize=(6,3))
    for k in finalResults.keys():
        plt.plot(sequenceLengths, [n*100/numRuns for n in finalResults[k]], marker='o', markersize=3, linewidth=0.75, label="Flattening" if k == 'None' else "$\lambda =\ $" + k)
    plt.xticks([100*n for n in range(1,11)])
    plt.yticks([10*n for n in range(1,11)])
    plt.xlabel('Sequence Length')
    plt.ylabel('\% Trees Correctly Reconstructed')
    plt.legend()

plt.tight_layout()
plt.savefig('Eriksson_Alg_Scaling_1.pdf')
plt.show()

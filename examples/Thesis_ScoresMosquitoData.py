from splitp import *
import matplotlib.pyplot as plt
from more_itertools import sort_together

DTable = hf.patternProbsFromAlignment("input/infile_mosquito.fa")[0]
tree = tree(1,4,"Primate_Data_Empty_Tree")
splits = hf.generateAllSplits(len(DTable[0][0]), trivial=False)

splitScores1 = [splits, [], [], [], [], []]
for sp in splits:
    print(sp)
    F = tree.flattening(sp, DTable)
    TF = tree.transformedFlattening(F)
    SF = tree.subFlattening(TF, specialState='T', type=(1, 1))
    SF_1_2 = tree.subFlattening(TF, specialState='T', type=(1, 2))
    SF_2_1 = tree.subFlattening(TF, specialState='T', type=(2, 1))
    SF_2_2 = tree.subFlattening(TF, specialState='T', type=(2, 2))
    splitScores1[1].append(tree.splitScore(F))
    splitScores1[2].append(tree.splitScore(SF))
    splitScores1[3].append(tree.splitScore(SF_1_2))
    splitScores1[4].append(tree.splitScore(SF_2_1))
    splitScores1[5].append(tree.splitScore(SF_2_2))

splitScores1 = sort_together(splitScores1, key_list=(1,2))
splits1 = splitScores1[0]


splitScores2 = [splits, [], [], [], [], []]
for sp in splits:
    print(sp)
    F = tree.flattening(sp, DTable)
    SF = tree.subFlatteningAlt(F)
    SF_2 = tree.subFlatteningAlt(F, S=scaledHMatrix(0.5)[0])
    SF_3 = tree.subFlatteningAlt(F, S=scaledHMatrix(3)[0])
    SF_10 = tree.subFlatteningAlt(F, S=scaledHMatrix(15)[0])
    splitScores2[1].append(tree.splitScore(F))
    splitScores2[2].append(tree.splitScore(SF))
    splitScores2[3].append(tree.splitScore(SF_2))
    splitScores2[4].append(tree.splitScore(SF_3))
    splitScores2[5].append(tree.splitScore(SF_10))

splitScores2 = sort_together(splitScores2, key_list=(1,2))
splits2 = splitScores2[0]

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

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6,4))

with sns.color_palette("colorblind"):
    splits = [x.replace('|', '\\textbar') for x in splits]
    ax1.plot(splits, splitScores1[1], marker='o', markersize=3, linewidth=0.75, label="Flat")
    ax1.plot(splits, splitScores1[2], marker='o', markersize=3, linewidth=0.75, label="(1,1)-Subflat")
    ax1.plot(splits, splitScores1[3], marker='o', markersize=3, linewidth=0.75, label="(1,2)-Subflat")
    ax1.plot(splits, splitScores1[4], marker='o', markersize=3, linewidth=0.75, label="(2,1)-Subflat")
    #ax1.plot(splits, splitScores1[5], marker='o', markersize=3, linewidth=0.75, label="(2,2)-Subflat")
    ax1.set_ylim(0)

    ax2.plot(splits, splitScores2[1], marker='o', markersize=3, linewidth=0.75, label="Flat")
    ax2.plot(splits, splitScores2[2], marker='o', markersize=3, linewidth=0.75, label="Subflat")
    ax2.plot(splits, splitScores2[3], marker='o', markersize=3, linewidth=0.75, label="$\lambda =\ 0.5$")
    ax2.plot(splits, splitScores2[4], marker='o', markersize=3, linewidth=0.75, label="$\lambda =\ 3$")
    ax2.plot(splits, splitScores2[5], marker='o', markersize=3, linewidth=0.75, label="$\lambda =\ 15$")
    ax2.set_ylim(0)
    ax1.legend()
    ax2.legend()

ax3 = fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

ax3.set_xlabel("Split")
ax3.set_ylabel("Split Score", labelpad= 25)
fig.tight_layout()
fig.savefig('Scores_mosquito.pdf')
plt.show()

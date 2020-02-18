from splitp import *
import matplotlib.pyplot as plt
from more_itertools import sort_together

DTable = hf.pattern_probs_from_alignment("input/infile_primate_sorted.fa")[0]
tree = old_tree(5, 4, "Primate_Data_Empty_Tree")
splits = hf.generate_all_splits(5, trivial=False)

splitScores = [splits, [], [], [], [], []]
for sp in splits:
    print(sp)
    F = tree.flattening(sp, DTable)
    TF = tree.transformedFlattening(F)
    SF = tree.subFlattening(TF, specialState='T', type=(1, 1))
    SF_1_2 = tree.subFlattening(TF, specialState='T', type=(1, 2))
    SF_2_1 = tree.subFlattening(TF, specialState='T', type=(2, 1))
    SF_2_2 = tree.subFlattening(TF, specialState='T', type=(2, 2))
    splitScores[1].append(tree.splitScore(F))
    splitScores[2].append(tree.splitScore(SF))
    splitScores[3].append(tree.splitScore(SF_1_2))
    splitScores[4].append(tree.splitScore(SF_2_1))
    splitScores[5].append(tree.splitScore(SF_2_2))

splitScores = sort_together(splitScores, key_list=(2, 1))
splits = splitScores[0]

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

plt.figure(figsize=(6, 2.8))
with sns.color_palette("colorblind"):
    splits = [x.replace('|', '\\textbar') for x in splits]
    plt.plot(splits, splitScores[1], marker='o', markersize=3, linewidth=0.75, label="Flat")
    plt.plot(splits, splitScores[2], marker='o', markersize=3, linewidth=0.75,label="(1,1)-Subflat")
    plt.plot(splits, splitScores[3], marker='o', markersize=3, linewidth=0.75,label="(1,2)-Subflat")
    plt.plot(splits, splitScores[4], marker='o', markersize=3, linewidth=0.75,label="(2,1)-Subflat")
    plt.plot(splits, splitScores[5], marker='o', markersize=3, linewidth=0.75,label="(2,2)-Subflat")
    plt.legend()
    plt.xlabel("Split")
    plt.ylabel("Split Score")
    plt.tight_layout()
    plt.savefig('SVDAlg_Primate.pdf')
plt.show()

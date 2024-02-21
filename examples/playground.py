# %% Imports
import numpy as np
import splitp as sp
import matplotlib.pyplot as plt

# %%
# Define a 8-taxon rooted binary tree in newick format with random branch lengths
newick = "((((A:0.1,B:0.04):0.06,(C:0.13,D:0.05):0.06):0.1,(E:0.1,F:0.1):0.1):0.15,(G:0.1,H:0.2):0.08);"
#newick = "((A:0.05,B:0.4):0.025,(C:0.05,D:0.4):0.025)"
tree = sp.Phylogeny(newick)
tree.draw()

# %%
# Define the model
model = sp.model.GTR.JukesCantor(0.5)

# %%
# Generate a sequence alignment
alignment = sp.simulation.generate_alignment(tree, model, 10000)

# %%
# Get all the true splits
splits = list(tree.splits())
all_splits = list(sp.splits.all_splits(tree, string_format=False))

# %%
# For each split, compute the flattening and store the result
flattenings = []
for split in all_splits:
    flattenings.append(sp.flattening(split, alignment, sp.FlatFormat.reduced))

# %%
# Compute split scores and KL divergences for each flattening
from splitp import phylogenetics
split_scores = []
KL_scores = []
for flat in flattenings:
    split_scores.append(sp.phylogenetics.split_score(flat))
    KL_scores.append(sp.phylogenetics.flattening_rank_1_approximation_divergence(np.array(flat)))

# %%
# Create new alignments, do an erikson SVD for each method and check if the sets of splits are the true ones. Tally the scores and decide the winner
score_flattening = 0
score_mutual_information = 0
for a in range(100):
    print(f"Alignment {a}")
    alignment = sp.simulation.generate_alignment(tree, model, 200)
    splits_flat = sp.phylogenetics.erickson_SVD(alignment, taxa=tree.get_taxa(), method=sp.Method.flattening)
    splits_KL = sp.phylogenetics.erickson_SVD(alignment, taxa=tree.get_taxa(), method=sp.Method.mutual_information)
    if set(splits_flat) >= set(splits):
        print("Flattening")
        score_flattening += 1
    if set(splits_KL) >= set(splits):
        print("Mutual information")
        score_mutual_information += 1


# %%

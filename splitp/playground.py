# %% Imports
import numpy as np
import splitp as sp

# %%
# Define a 8-taxon rooted binary tree in newick format with random branch lengths
newick = "((((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6):0.4,(E:0.7,F:0.8):0.9):0.3,(G:0.9,H:0.9):0.6);"
tree = sp.Phylogeny(newick)
tree.draw()

# %%
# Define the model
model = sp.model.GTR.JukesCantor()

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
    flattenings.append(sp.flattening(split, alignment))

# %%
# Compute split scores and KL divergences for each flattening
from splitp import phylogenetics
split_scores = []
KL_scores = []
for flat in flattenings:
    split_scores.append(sp.phylogenetics.split_score(flat))
    KL_scores.append(sp.phylogenetics.flattening_rank_1_approximation_divergence(np.array(flat.todense())))

# %%

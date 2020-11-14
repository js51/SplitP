import splitp
import numpy as np

# Mosquito
tree = splitp.NXTree("(A,B);")
DTable, _, _ = splitp.pattern_probs_from_alignment("./infile_mosquito.fa")
splits = list(splitp.all_splits(4, trivial=False))
flattening = tree.flattening(splits[1], DTable)
sparse_flattening = tree.sparse_flattening(splits[1], DTable, num_taxa=4)
subflattening1 = tree.subflattening_alt(flattening)
#subflattening2 = tree.sparse_subflattening(splits[0], DTable)

print(tree.split_score(subflattening1), tree.split_score(subflattening1), tree.split_score(subflattening1, force_frob_norm=True))
print(tree.split_score(flattening), tree.split_score(sparse_flattening), tree.split_score(flattening, force_frob_norm=True))


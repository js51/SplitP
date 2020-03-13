import splitp
import pandas as pd

tree = splitp.NXTree("(A,B);")
DTable, _, _ = splitp.pattern_probs_from_alignment("./20_taxa_infile.fa")

splits = ["01235|46789ABCDEFGHIJKLMNO"]

flattening = tree.sparse_flattening(splits[0], DTable, 25)

from scipy.sparse.linalg import svds, eigs
from scipy.sparse import csc_matrix

print(svds(flattening,4,return_singular_vectors=False))

flattening = flattening.todok()
subflattening = tree.sparse_subflattening(splits[0], DTable)

print(svds(subflattening,4,return_singular_vectors=False))


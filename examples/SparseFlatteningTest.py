import splitp
import pandas as pd
import time
import numpy as np

tree = splitp.NXTree("(A,B);")
DTable, _, _ = splitp.pattern_probs_from_alignment("./20_taxa_infile.fa")
num_taxa = 20
splits = splitp.all_splits(num_taxa, trivial=False)
split = "0123|4E6789ABCD5FGHIJ"

# Construction
start = time.time()
coo_flattening = tree.sparse_flattening(split, DTable, format='coo')
end = time.time()
print("coo flattening constructed in", end-start)

start = time.time()
dok_flattening = tree.sparse_flattening(split, DTable, format='dok')
end = time.time()
print("dok flattening constructed in", end-start)

# Computing scores
start = time.time()
coo_score = tree.split_score(coo_flattening, data_table_for_frob_norm=DTable)
end = time.time()
print("coo score calculated in", end-start)

start = time.time()
dok_score = tree.split_score(dok_flattening, data_table_for_frob_norm=DTable)
end = time.time()
print("dok score calculated in", end-start)

#dictionary = { key : val for key, val in DTable.itertuples(index=False)}
start = time.time()
subflattening = tree.signed_sum_subflattening(split, DTable)
subflattening_score = tree.split_score(subflattening)
end = time.time()
print("subflattening and score computed in", end-start)

print(subflattening_score)
print("Scores", coo_score, dok_score)

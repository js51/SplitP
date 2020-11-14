import splitp
import itertools
import random
n = 20
splits_of_size = {}
results = {}
results['Flat'] = {}
results['Subflat'] = {}

for i in range(2, 11):
	splits_of_size[i] = list(itertools.islice(splitp.all_splits(n, only_balance=i, randomise=True), 190))

for key, value in splits_of_size.items():
	print("Size ", str(key))
	
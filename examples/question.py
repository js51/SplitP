I am attempting to create large, sparse matrices for the purpose of computing their singular values (at least the first four of them) in Python using SciPy sparse matrices.

Fix `n=20`, then these matrices come in the shape `(4^m, 4^(n-m))` and have a small number of non-zero values (as few as 500). I am successfully creating these matrices and obtaining the singular values when `m` is close to `n` (when the matrices are 'close' to being square) very quickly and with barely any memory. 

But when they are short-and-wide or tall-and-skinny, (`m=3` for example), The Python process is being killed due to too much memory usage.

Here is my code:

	def sparse_flattening(self, split, table, n=20, format='dok'):
	    import scipy
	
	    if format == 'coo':
		    from scipy.sparse import coo_matrix
		    split = split.split('|')
		    rows = []
		    cols = []
		    data = []
		    for r in table.itertuples(index=False):
			    pattern = r[0]
			    # Compute the correct location using the pattern and the split.
			    row = self.__index_of(''.join([str(pattern[int(s, n)]) for s in split[0]]))
			    col = self.__index_of(''.join([str(pattern[int(s, n)]) for s in split[1]]))
			    rows.append(row)
			    cols.append(col)
			    data.append(r[1])
		    return coo_matrix((data, (rows, cols)), shape=(4**len(split[0]),4**len(split[1])))
		
	    elif format == 'dok':
		    from scipy.sparse import dok_matrix
		    split = split.split('|')
		    flattening = dok_matrix((4**len(split[0]),4**len(split[1])))
		    for r in table.itertuples(index=False):
			    pattern = r[0]
			    row = self.__index_of(''.join([str(pattern[int(s, n)]) for s in split[0]]))
			    col = self.__index_of(''.join([str(pattern[int(s, n)]) for s in split[1]]))
			    flattening[row, col] = r[1]
		    return flattening


I have tried a few of the different SciPy sparse matrix formats.

My question is, why does the shape of the matrix affect the memory usage so dramatically despite there being the same number of entries in each of these matrices, and how can I store these tall or short matrices efficiently in python?

**EDIT:** On further inspection, it seems like the slow and memory demanding part of my code was computing the frobenius norm of the matrices using `scipy.sparse.linalg.norm(matrix, 'fro')`. I was normalising a sum over the singular values and assumed the singular value computation was the slow part. Not understanding why the norm function could eat up so much memory, I looked at the implementation of `_sparse_frobenius_norm` SciPy and discovered it was faulty. I have written my own code for the norm and no longer have memory troubles.
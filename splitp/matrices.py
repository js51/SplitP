import scipy
import numpy as np

def is_sparse(matrix):
    return scipy.sparse.issparse(matrix)


def frobenius_norm(matrix, data_table=None):
    """Calculates the Frobenius Norm for a given matrix"""
    if data_table is not None:
        return sum(val**2 for _, val in data_table.itertuples(index=False))**(1/2)
    if is_sparse(matrix):
        return sum(matrix[i, j]**2 for i, j in zip(*matrix.nonzero()))**(1/2)
    else:
        return np.sqrt(sum(val**2 for val in np.nditer(matrix)))

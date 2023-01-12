from splitp import Phylogeny
from math import floor
import numpy as np

def balanced_newick_tree(num_taxa):
    if num_taxa % 2 != 0:
        raise ValueError(
            "There is no balanced tree on {num_taxa} taxa. Please specify an even number.")
    if num_taxa == 2:
        return "(0,1);"

    newick_string = f"({__balanced_newick_subtree(num_taxa/2, True)},{__balanced_newick_subtree(num_taxa/2)});"
    for i in range(0, num_taxa):
        leaf_name = str(np.base_repr(i, base=max(i+1, 2))
                        ) if num_taxa <= 36 else f't{str(i)}'
        newick_string = newick_string.replace('_', leaf_name, 1)
    return Phylogeny(newick_string)

def __balanced_newick_subtree(nt, left=False):
    if nt == 2:
        return "(_,_)"
    elif nt == 3:
        return "((_,_),_)" if left else "(_,(_,_))"
    else:
        if nt % 2 == 0:
            return f"({__balanced_newick_subtree(nt/2, True)},{__balanced_newick_subtree(nt/2)})"
        else:
            return f"({__balanced_newick_subtree(floor(nt/2) + int(left), True)},{__balanced_newick_subtree(floor(nt/2) + int(not left))})"
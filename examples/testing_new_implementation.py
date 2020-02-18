from splitp import *
import time

# 6 Taxa Tree
num_nodes 		= 11					         # The total number of nodes in the tree
num_bases 		= 4   				         # Size of the character set
initial_dist 	= [0.25, 0.25, 0.25, 0.25]    # Initial Distribution for root node
matrix 			= [[0.95, 0.05/3, 0.05/3, 0.05/3],		# ACGT
		  		   [0.05/3, 0.95, 0.05/3, 0.05/3],
                    [0.05/3, 0.05/3, 0.95, 0.05/3],
                    [0.05/3, 0.05/3, 0.05/3, 0.95]]
T1 = tree(num_nodes, num_bases)
T1.initDist = initial_dist
# Adding nodes to the tree, each is given a matrix for it's incoming edge and an ID.
# The second parameter is the id of the parent node.
T1.addNode(node(None, 0)) # Root Node
T1.addNode(node(matrix, 1),  0)
T1.addNode(node(matrix, 2),  0)
T1.addNode(node(matrix, 3),  1)
T1.addNode(node(matrix, 7),  1)
T1.addNode(node(matrix, 8),  2)
T1.addNode(node(matrix, 4),  2)
T1.addNode(node(matrix, 5),  3)
T1.addNode(node(matrix, 6),  3)
T1.addNode(node(matrix, 9),  4)
T1.addNode(node(matrix, 10), 4)
###

T = nx_tree("(((A,B),C),(D,(E,F)));")

T.reassign_all_transition_matrices(np.array([[0.95, 0.05/3, 0.05/3, 0.05/3],
                                    [0.05/3, 0.95, 0.05/3, 0.05/3],
                                    [0.05/3, 0.05/3, 0.95, 0.05/3],
                                    [0.05/3, 0.05/3, 0.05/3, 0.95]]))


start1 = time.time()
LTable_old = T1.getLikelihoods()
end1 = time.time()
print(end1 - start1)

start2 = time.time()
LTable = T.get_pattern_probabilities()
end2 = time.time()
print(end2 - start2)

print(LTable_old.equals(LTable))

from splitp import *

def makeSubsMatrix(subs_prob, k):
    matrix = []
    for i in range(k):
        matrix.append([1-subs_prob if j==i else subs_prob/(k-1) for j in range(k)])
    return matrix

trees = dict()

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
T1.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 2),  0)
T1.addNode(node(matrix, 3),  1)
T1.addNode(node(matrix, 7),  1)
T1.addNode(node(matrix, 8),  2)
T1.addNode(node(matrix, 4),  2)
T1.addNode(node(matrix, 5),  3)
T1.addNode(node(matrix, 6),  3)
T1.addNode(node(matrix, 9),  4)
T1.addNode(node(matrix, 10), 4)
###

# 8 Taxa Tree
T2 = tree(15, 4)
T2.initDist = initial_dist
T2.addNode(node(None, 0))
T2.addNode(node(matrix, 1),  0)
T2.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 2),  0)
T2.addNode(node(matrix, 3),  1)
T2.addNode(node(matrix, 4),  1)
T2.addNode(node(matrix, 5),  2)
T2.addNode(node(matrix, 6),  2)
T2.addNode(node(matrix, 7),  3)
T2.addNode(node(matrix, 8),  3)
T2.addNode(node(matrix, 9),  4)
T2.addNode(node(matrix, 10), 4)
T2.addNode(node(matrix, 11), 5)
T2.addNode(node(matrix, 12), 5)
T2.addNode(node(matrix, 13), 6)
T2.addNode(node(matrix, 14), 6)
###

# 6 Taxa Tree (2)
matrixXS = makeSubsMatrix(0.025, 4)
m = np.matrix(matrixXS)
T6B = tree(11, 4)
T6B.initDist = initial_dist
T6B.addNode(node(None, 0))
T6B.addNode(node(matrixXS, 1),  0)
T6B.addNode(node(matrixXS, 2),  0)
T6B.addNode(node((m@m).tolist(), 3),  1)
T6B.addNode(node((m@m@m@m@m).tolist(), 7),  1)
T6B.addNode(node((m@m@m@m@m).tolist(), 8),  2)
T6B.addNode(node((m@m).tolist(), 4),  2)
T6B.addNode(node((m@m@m@m@m@m@m@m).tolist(), 5),  3)
T6B.addNode(node((m@m@m@m@m@m@m@m).tolist(), 6),  3)
T6B.addNode(node((m@m@m@m@m@m@m@m).tolist(), 9),  4)
T6B.addNode(node((m@m@m@m@m@m@m@m).tolist(), 10), 4)
###


matrixLong		 = [[0.7, 0.3/3, 0.3/3, 0.3/3],		# ACGT
		  		    [0.3/3, 0.7, 0.3/3, 0.3/3],
                    [0.3/3, 0.3/3, 0.7, 0.3/3],
                    [0.3/3, 0.3/3, 0.3/3, 0.7]]
matrixShort = [[0.95, 0.05/3, 0.05/3, 0.05/3],		# ACGT
		  		  [0.05/3, 0.95, 0.05/3, 0.05/3],
                    [0.05/3, 0.05/3, 0.95, 0.05/3],
                    [0.05/3, 0.05/3, 0.05/3, 0.95]]
# 4 taxon tree
T3 = tree(7, 4)
T3.initDist = initial_dist
T3.addNode(node(None, 0))
T3.addNode(node(matrix, 1),  0)
T3.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 2),  0)
T3.addNode(node(matrixShort, 3),  1)
T3.addNode(node(matrixLong, 4),  1)
T3.addNode(node(matrixShort, 5),  2)
T3.addNode(node(matrixLong, 6),  2)
###

# 7 Taxa Tree
num_nodes 		= 13
num_bases 		= 4
initial_dist 	= [0.25, 0.25, 0.25, 0.25]
matrix_a	    = makeSubsMatrix(0.01, 4)
matrix_b	    = makeSubsMatrix(0.03, 4)
matrix_c	    = makeSubsMatrix(0.08, 4)
T7 = tree(num_nodes, num_bases, name="7TaxonTreeA")
T7.initDist = initial_dist

T7.addNode(node(None, 0)) # Root Node
T7.addNode(node(matrix_a, 1),  0)
T7.addNode(node(matrix_a, 2),  0)

# Internal Nodes:
T7.addNode(node(matrix_b, 4),  1)
T7.addNode(node(matrix_b, 6),  2)
T7.addNode(node(matrix_b, 10), 6)

# Leaves:
T7.addNode(node(matrix_c, 3),  1)
T7.addNode(node(matrix_c, 5),  4)
T7.addNode(node(matrix_c, 7),  4)
T7.addNode(node(matrix_c, 8),  2)
T7.addNode(node(matrix_c, 9),  6)
T7.addNode(node(matrix_c, 11),  10)
T7.addNode(node(matrix_c, 12), 10)
###

##### Long branch attraction trees #####
matrix = [[0.95, 0.05/3, 0.05/3, 0.05/3],
		  [0.05/3, 0.95, 0.05/3, 0.05/3],
          [0.05/3, 0.05/3, 0.95, 0.05/3],
          [0.05/3, 0.05/3, 0.05/3, 0.95]]
import numpy as np
m = np.matrix(matrix)
matrixLong = (m@m@m@m@m@m@m@m).tolist()
LBATrueSplits = {'01|23'}
# A
T4A = tree(7, 4)
T4A.initDist = initial_dist
T4A.addNode(node(None, 0))
T4A.addNode(node(matrix, 1),  0)
T4A.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 2),  0)
T4A.addNode(node(matrix, 3),  1)
T4A.addNode(node(matrix, 4),  1)
T4A.addNode(node(matrix, 5),  2)
T4A.addNode(node(matrix, 6),  2)
###
# B
T4B = tree(7, 4)
T4B.initDist = initial_dist
T4B.addNode(node(None, 0))
T4B.addNode(node(matrix, 1),  0)
T4B.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 2),  0)
T4B.addNode(node(matrix, 3),  1)
T4B.addNode(node(matrixLong, 4),  1)
T4B.addNode(node(matrix, 5),  2)
T4B.addNode(node(matrixLong, 6),  2)
###
# C
T4C = tree(7, 4)
T4C.initDist = initial_dist
T4C.addNode(node(None, 0))
T4C.addNode(node(matrix, 1),  0)
T4C.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 2),  0)
T4C.addNode(node(matrixLong, 3),  1)
T4C.addNode(node(matrixLong, 4),  1)
T4C.addNode(node(matrix, 5),  2)
T4C.addNode(node(matrix, 6),  2)
###
# D
T4D = tree(7, 4)
T4D.initDist = initial_dist
T4D.addNode(node(None, 0))
T4D.addNode(node(matrix, 1),  0)
T4D.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 2),  0)
T4D.addNode(node(matrixLong, 3),  1)
T4D.addNode(node(matrixLong, 4),  1)
T4D.addNode(node(matrixLong, 5),  2)
T4D.addNode(node(matrixLong, 6),  2)
###
# Star
star = tree(7, 4)
star.initDist = initial_dist
star.addNode(node(None, 0))
star.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 1),  0)
star.addNode(node([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], 2),  0)
star.addNode(node(matrixLong, 3),  1)
star.addNode(node(matrixLong, 4),  1)
star.addNode(node(matrix, 5),  2)
star.addNode(node(matrix, 6),  2)
###



# (tree, num true splits)
trees[6] = (T1, 3, {"01|2345", "012|345", "0123|45"})
trees[4] = (T3, 1, {"01|23"})
trees[8] = (T2, 5, {"01|234567", "0123|4567", "012345|67", "014567|23", "012367|45"})
trees["T6B"] = (T6B, 3, {"01|2345", "012|345", "0123|45"})
trees[7] = (T7, 4, {"012|3456", "03456|12", "0123|456", "01234|56"})
trees['T4A'] = (T4A, 1, LBATrueSplits)
trees['T4B'] = (T4B, 1, LBATrueSplits)
trees['T4C'] = (T4C, 1, LBATrueSplits)
trees['T4D'] = (T4D, 1, LBATrueSplits)
trees['star'] = (star, 3, generateAllSplits(4,trivial=False))
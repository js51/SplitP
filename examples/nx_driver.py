from splitp import *

newick_tree_string = "(((A,B),C),(D,(E,F)));"
T = NXTree(newick_tree_string)
T.reassign_all_transition_matrices(np.array([[0.95, 0.05/3, 0.05/3, 0.05/3],
                                    [0.05/3, 0.95, 0.05/3, 0.05/3],
                                    [0.05/3, 0.05/3, 0.95, 0.05/3],
                                    [0.05/3, 0.05/3, 0.05/3, 0.95]]))
#T.addNode("Hello", 100)
LTable = T.get_pattern_probabilities()

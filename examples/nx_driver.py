from SplitP import *

newick_tree_string = "((A,B),(C,D));"
json_tree_string = newick_to_json(newick_tree_string, generate_names=True)
print(json_to_newick(json_tree_string))

print(json_tree_string)

T = nxtree(newick_tree_string)
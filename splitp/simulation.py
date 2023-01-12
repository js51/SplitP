import itertools
import numpy as np
from warnings import warn
from numpy.random import choices
from networkx import dfs_tree


def evolve_pattern(tree, model):
    def __evolve_on_subtree(subtree, state):
        root_node = [n for n, d in subtree.in_degree() if d == 0][0]
        children = list(subtree.successors(root_node))
        probs = list(
            tree.networkx_graph.nodes[root_node]["transition_matrix"][
                :, tree.state_space.index(state)
            ]
        )
        new_state = choices(tree.state_space, probs)[0]
        if len(children) == 0:
            return f"{str(root_node)}:{new_state}"
        else:
            subtrees = [dfs_tree(tree.networkx_graph, child) for child in children]
            return ",".join(
                __evolve_on_subtree(subtree, new_state) for subtree in subtrees
            )

    root_state = choices(tree.state_space, tree.initial_dist)[0]
    root_node = [n for n, d in tree.networkx_graph.in_degree() if d == 0][0]
    children = list(tree.networkx_graph.successors(root_node))
    subtrees = [dfs_tree(tree.networkx_graph, child) for child in children]
    result_string = ",".join(
        __evolve_on_subtree(subtree, root_state) for subtree in subtrees
    )
    result = {
        pair.split(":")[0]: pair.split(":")[1] for pair in result_string.split(",")
    }
    return "".join(result[k] for k in sorted(result.keys(), key=tree.taxa.index))


def generate_alignment(tree, model, sequence_length):
    counts = {}
    for i in range(sequence_length):
        pattern = evolve_pattern(tree)
        if pattern not in counts:
            counts[pattern] = float(1)
        else:
            counts[pattern] += 1
    probs = {}
    for k in sorted(
        counts.keys(), key=lambda p: [tree.state_space.index(c) for c in p]
    ):
        probs[k] = counts[k] / float(sequence_length)
    return probs


def draw_from_multinomial(pattern_probabilities, n):
    """Use a given table of probabilities from get_pattern_probabilities() and draw from its distribution"""
    patterns, probs = zip(*pattern_probabilities.items())
    data = np.random.multinomial(n, probs)
    data = [d / n for d in data]
    results = {}
    for p, pattern in enumerate(patterns):
        if data[p] != 0:
            results[pattern] = data[p]
    return


def get_pattern_probabilities(tree, model=None):
    """Returns a full table of site-pattern probabilities (binary character set)"""
    # Creating a table with binary labels and calling likelihood_start() to fill it in with probabilities
    combinations = list(
        itertools.product(
            "".join(s for s in tree.state_space), repeat=tree.get_num_taxa()
        )
    )
    combinations = ["".join(c) for c in combinations]
    emptyArray = {
        combination: __likelihood_start(combination) for combination in combinations
    }
    return emptyArray


def __likelihood(tree, n, likelihood_table):
    """Recursive part of the likelihood algorithm"""
    for b in range(tree.num_bases):
        L = 1
        for d in tree.get_descendants(n, return_iter=True):
            d_index = tree.node_index(d)
            M = (tree.nx_graph.nodes[d])["transition_matrix"]
            s = 0
            for b2 in range(tree.num_bases):
                if likelihood_table[d_index, b2] == None:
                    __likelihood(d, likelihood_table)
                s += M[b2, b] * likelihood_table[d_index, b2]
            L *= s
        likelihood_table[tree.node_index(n), b] = L
    if not tree.is_root(n):
        __likelihood(tree.get_parent(n), likelihood_table)


def __likelihood_start(tree, pattern):
    """Starts the likelihood algorithm.

    Starts the likelihood algorithm to determine the probability of seeing a given site pattern.

    Args:
        pattern: The site pattern to obtain the probability for.

    Returns:
        A probability (range 0-1) of observing the given site pattern given the tree.
    """

    def _to_int(p):
        return tree.state_space.index(p)

    pattern = [
        _to_int(p) for p in pattern
    ]  # A list of indices which correspond to taxa.
    taxa = tree.taxa  # The list of taxa.
    # Likelihood table for dynamic prog alg ~ lTable[node_index, character]
    likelihood_table = np.array(
        [[None for _ in range(tree.num_bases)] for _ in range(tree.num_nodes())]
    )

    # Encoding pattern into likelihood table
    for i, p in enumerate(pattern):
        likelihood_table[tree.node_index(taxa[i]), :] = [
            0 for _ in range(tree.num_bases)
        ]
        likelihood_table[tree.node_index(taxa[i]), p] = 1

    # Starting with node which has all descendants as leaves
    initNode = None
    for n in tree.nx_graph.nodes:
        desc = tree.get_descendants(n)
        if desc and all(
            tree.is_leaf(d) for d in desc
        ):  # If the node has descendants and they are all leaves
            initNode = n

    __likelihood(initNode, likelihood_table)  # Should fill in the entire table

    root_index = tree.node_index(tree.get_root(return_index=False))
    L = 0
    for i in range(tree.num_bases):
        L += tree.initDist[i] * likelihood_table[root_index, i]

    return L

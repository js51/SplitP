import itertools
import numpy as np
from warnings import warn
from random import choices
from networkx import dfs_tree
from splitp import constants


def evolve_pattern(tree, model=None):
    def __evolve_on_subtree(subtree, state):
        root_node = [n for n, d in subtree.in_degree() if d == 0][0]
        children = list(subtree.successors(root_node))
        if model == None: # Get the transition matrix from the phylogeny
                M = (tree.networkx_graph.nodes[root_node])["transition_matrix"]
        else:
            branch_length = (tree.networkx_graph.nodes[root_node])["branch_length"]
            M = model.transition_matrix(branch_length)
        probs = list(M[:, constants.DNA_state_space.index(state)])
        new_state = choices(constants.DNA_state_space, probs)[0]
        if len(children) == 0:
            return f"{str(root_node)}:{new_state}"
        else:
            subtrees = [dfs_tree(tree.networkx_graph, child) for child in children]
            return ",".join(
                __evolve_on_subtree(subtree, new_state) for subtree in subtrees
            )

    root_state = choices(constants.DNA_state_space, (1/4, 1/4, 1/4, 1/4))[0]
    root_node = [n for n, d in tree.networkx_graph.in_degree() if d == 0][0]
    children = list(tree.networkx_graph.successors(root_node))
    subtrees = [dfs_tree(tree.networkx_graph, child) for child in children]
    result_string = ",".join(
        __evolve_on_subtree(subtree, root_state) for subtree in subtrees
    )
    result = {
        pair.split(":")[0]: pair.split(":")[1] for pair in result_string.split(",")
    }
    taxa = tree.get_taxa()
    return "".join(result[k] for k in sorted(result.keys(), key=taxa.index))


def generate_alignment(tree, model, sequence_length):
    counts = {}
    for i in range(sequence_length):
        pattern = evolve_pattern(tree, model)
        if pattern not in counts:
            counts[pattern] = float(1)
        else:
            counts[pattern] += 1
    probs = {}
    for k in sorted(
        counts.keys(), key=lambda p: [constants.DNA_state_space.index(c) for c in p]
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
    return results


def get_pattern_probabilities(tree, model=None):
    """Returns a full table of site-pattern probabilities (binary character set)"""
    # Creating a table with binary labels and calling likelihood_start() to fill it in with probabilities
    combinations = list(
        itertools.product(
            "".join(s for s in constants.DNA_state_space), repeat=tree.get_num_taxa()
        )
    )
    combinations = ["".join(c) for c in combinations]
    emptyArray = {
        combination: __likelihood_start(tree, combination, model) for combination in combinations
    }
    return emptyArray

pattern_probabilities = get_pattern_probabilities


def __likelihood(tree, n, likelihood_table, model):
    """Recursive part of the likelihood algorithm"""
    for b in range(4):
        L = 1
        for d in tree.get_descendants(n, return_iter=True):
            d_index = tree.node_index(d)
            if model == None: # Get the transition matrix from the phylogeny
                M = (tree.networkx_graph.nodes[d])["transition_matrix"]
            else:
                branch_length = (tree.networkx_graph.nodes[d])["branch_length"]
                M = model.transition_matrix(branch_length)
            s = 0
            for b2 in range(4):
                if likelihood_table[d_index, b2] == None:
                    __likelihood(tree, d, likelihood_table, model)
                s += M[b2, b] * likelihood_table[d_index, b2]
            L *= s
        likelihood_table[tree.node_index(n), b] = L
    if not tree.is_root(n):
        __likelihood(tree, tree.get_parent(n), likelihood_table, model)


def __likelihood_start(tree, pattern, model):
    """Starts the likelihood algorithm.

    Starts the likelihood algorithm to determine the probability of seeing a given site pattern.

    Args:
        pattern: The site pattern to obtain the probability for.

    Returns:
        A probability (range 0-1) of observing the given site pattern given the tree.
    """

    def _to_int(p):
        return constants.DNA_state_space.index(p)

    pattern = [
        _to_int(p) for p in pattern
    ]  # A list of indices which correspond to taxa.
    taxa = tree.get_taxa()  # The list of taxa.
    # Likelihood table for dynamic prog alg ~ lTable[node_index, character]
    likelihood_table = np.array(
        [[None for _ in range(4)] for _ in range(tree.get_num_nodes())]
    )

    # Encoding pattern into likelihood table
    for i, p in enumerate(pattern):
        likelihood_table[tree.node_index(taxa[i]), :] = [
            0 for _ in range(4)
        ]
        likelihood_table[tree.node_index(taxa[i]), p] = 1

    # Starting with node which has all descendants as leaves
    initNode = None
    for n in tree.networkx_graph.nodes:
        desc = tree.get_descendants(n)
        if desc and all(
            tree.is_leaf(d) for d in desc
        ):  # If the node has descendants and they are all leaves
            initNode = n

    __likelihood(tree, initNode, likelihood_table, model)  # Should fill in the entire table

    root_index = tree.node_index(tree.root(return_index=False))
    L = 0
    for i in range(4):
        L += (1/4) * likelihood_table[root_index, i]

    return L

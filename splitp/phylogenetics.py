import splitp as sp
import numpy as np
from itertools import combinations
import networkx as nx
from networkx import dfs_postorder_nodes, bfs_successors
from splitp import matrix, splits
import scipy


def parsimony_score(self, pattern):
    """Calculate a parsimony score for a site pattern or split

    Args:
        pattern: A string representing a site pattern or split

    Returns:
        A parsimony score for the given split or site pattern
    """
    graph = self.nx_graph.copy()
    nodes = graph.nodes
    if "|" in pattern:
        pattern2 = [0 for x in range(len(pattern) - 1)]
        i = 0
        while pattern[i] != "|":
            pattern2[int(pattern[i])] = 1
            i += 1
        pattern = "".join(str(i) for i in pattern2)

    taxa = [t for t in self.taxa]
    for i, t in enumerate(taxa):
        nodes[t]["pars"] = set(pattern[i])
    score = 0
    for n in nodes:
        if not self.is_leaf(n) and "pars" not in nodes[n]:
            score += self.__parsimony(n, nodes=nodes)
    return score


def __parsimony(self, n, nodes=None):
    """Recursive step in Finch's Algorithm"""
    score = 0
    children = self.get_descendants(n)
    for c in children:
        if "pars" not in nodes[c]:
            score += self.__parsimony(c, nodes)

    nodes[n]["pars"] = nodes[children[0]]["pars"]  # children[0].pars
    for c in children:
        nodes[n]["pars"] = nodes[n]["pars"].intersection(nodes[c]["pars"])
    if nodes[n]["pars"] == set():
        nodes[n]["pars"] = nodes[children[0]]["pars"]
        for c in children:
            nodes[n]["pars"] = (nodes[n]["pars"]).union(nodes[c]["pars"])
        return score + 1
    return score


def hartigan_algorithm(self, pattern):
    score = 0
    graph = self.nx_graph.copy()
    nodes = graph.nodes
    taxa = [t for t in self.taxa]
    for i, t in enumerate(taxa):
        nodes[t]["S1"] = set(pattern[i])
    postorder_nodes = list(
        dfs_postorder_nodes(graph, source=self.get_root(return_index=False))
    )
    for n in postorder_nodes:
        if not self.is_leaf(n):
            children = self.get_descendants(n)
            k = {}
            for state in self.state_space:
                k[state] = len(
                    set(child for child in children if state in nodes[child]["S1"])
                )
            K = max(k.values())
            nodes[n]["S1"] = {state for state in self.state_space if k[state] == K}
            nodes[n]["S2"] = {state for state in self.state_space if k[state] == K - 1}
    # Now compute the score
    top_to_bottom_nodes = [
        x[0] for x in bfs_successors(graph, source=self.get_root(return_index=False))
    ] + taxa
    for n in top_to_bottom_nodes:
        if n == self.get_root(return_index=False):
            nodes[n]["hart_state"] = list(nodes[n]["S1"])[0]
        else:
            parent = nodes[list(graph.predecessors(n))[0]]
            if parent["hart_state"] not in nodes[n]["S1"]:
                nodes[n]["hart_state"] = list(nodes[n]["S1"])[0]
                score += 1
            else:
                nodes[n]["hart_state"] = parent["hart_state"]
    return score


def erickson_SVD(alignment, taxa=None, method=sp.Method.flattenings):
    all_scores = {}
    bundle = {}
    labels = {}

    def _erickstep(all_taxa, alignment):
        scores = {}
        for pair in combinations(taxa, 2):
            flat_pair = tuple(
                sorted(
                    element
                    for tup in pair
                    for element in (tup if not isinstance(tup, str) else (tup,))
                )
            )
            other = tuple(
                sorted(
                    element
                    for tup in all_taxa
                    for element in (tup if not isinstance(tup, str) else (tup,))
                    if element not in flat_pair
                )
            )
            split = (flat_pair, other)
            try:
                score = all_scores[split]
            except KeyError:
                if method == sp.Method.flattenings:
                    xflat = tree.reduced_sparse_flattening(split, alignment)
                elif method == sp.Method.subflattenings:
                    xflat = tree.fast_signed_sum_subflattening(
                        split, alignment, bundle=bundle, labels=labels
                    )
                score = tree.split_score(xflat)
                all_scores[split] = score
            scores[pair] = (pair, split, score)
        best_pair, best_split, best_score = min(scores.values(), key=lambda x: x[2])
        return best_pair, best_split, best_score

    num_taxa = len(list(alignment.keys())[0])  # Length of first pattern
    if taxa is None:
        taxa = [
            str(np.base_repr(i, base=max(i + 1, 2))) if num_taxa <= 36 else f"t{str(i)}"
            for i in range(num_taxa)
        ]
    tree = sp.NXTree.dummy_tree(taxa=taxa)
    true_splits = []
    while len(true_splits) < num_taxa - 2:
        best_pair, best_split, best_score = _erickstep(taxa, alignment)
        true_splits.append(tuple(sorted(best_split)))
        taxa = tuple(
            [
                element
                for element in taxa
                if (
                    element not in best_split[0]
                    and not set(element).issubset(best_split[0])
                )
            ]
            + [best_split[0]]
        )
    return true_splits


def newick_string_from_splits(splits):
    def _consolidate(tup, smaller_halves):
        if len(tup) == 2:
            return tup
        if len(tup) == 1:
            return tup[0]
        if isinstance(tup, str):
            return tup
        for smaller_half in smaller_halves:
            if set(smaller_half).issubset(tup) and len(smaller_half) < len(tup):
                # then the consolidation is made up of the smaller half and what is left over
                left_over = tuple(set(tup).difference(set(smaller_half)))
                try:
                    smaller_halves.remove(smaller_half)
                    smaller_halves.remove(left_over)
                except ValueError:
                    pass
                return tuple(
                    (
                        _consolidate(left_over, smaller_halves),
                        _consolidate(smaller_half, smaller_halves),
                    )
                )
    splits = sorted(splits, key=lambda x: min(len(x[0]), len(x[1])), reverse=True)
    if len(splits) == 1:
        return str(splits[0]).replace("'", "").replace(" ", "") + ";"
    if len(splits) == 0:
        return ";"
    splits = iter(splits)
    first_split = next(splits)
    smaller_halves = [min(split, key=len) for split in splits]
    consolidated_split = tuple(
        _consolidate(half, smaller_halves) for half in first_split
    )
    return str(consolidated_split).replace("'", "").replace(" ", "") + ";"


def tree_from_splits(splits):
    return sp.NXTree(newick_string_from_splits(splits))


def split_tree_parsimony(alignment, splits=None):
    if type(alignment) is dict:
        alignment_dict = alignment
    else:
        alignment_dict = {}
        for table_pattern, value in alignment.itertuples(index=False, name=None):
            alignment_dict[table_pattern] = value
    num_taxa = len(list(alignment_dict.keys())[0])  # Length of first pattern
    all_splits = list(splits.all_splits(num_taxa)) if splits is None else splits
    scores = {split: 0 for split in all_splits}
    for split in all_splits:
        newick_string = []
        for part in split.split("|"):
            newick_string.append(f'({",".join(c for c in part)})')
        newick_string = f"({newick_string[0]},{newick_string[1]});"
        split_tree = sp.NXTree(newick_string, taxa_ordering="sorted")
        for pattern, value in alignment_dict.items():
            scores[split] += value * (
                split_tree.hartigan_algorithm(pattern) / (num_taxa - 1)
            )
    return scores


def JC_corrected_distance_matrix(alignment):
    num_taxa = len(list(alignment.keys())[0])  # Length of first pattern
    taxa = [str(i) for i in range(num_taxa)]
    distance_matrix = [[0 for _ in range(num_taxa)] for _ in range(num_taxa)]
    for i in range(num_taxa):
        for j in range(i + 1, num_taxa):
            for pattern, value in alignment.items():
                if pattern[i] != pattern[j]:
                    distance_matrix[i][j] += value
                    distance_matrix[j][i] += value
    for i in range(num_taxa):
        for j in range(i + 1, num_taxa):
            distance_matrix[i][j] = (
                -3.0 / 4.0 * np.log(1 - 4.0 / 3.0 * distance_matrix[i][j])
            )
            distance_matrix[j][i] = distance_matrix[i][j]
    return distance_matrix


def euclidean_split_distance(alignment, splits):
    print("assuming sorted taxa")
    states = ("A", "G", "C", "T")
    alignment_dict = {}
    for table_pattern, value in alignment.itertuples(index=False, name=None):
        alignment_dict[table_pattern] = value
    num_taxa = len(list(alignment_dict.keys())[0])  # Length of first pattern
    all_splits = list(splits.all_splits(num_taxa)) if splits == None else splits
    scores = {split: 0 for split in all_splits}
    for split in all_splits:
        split_list = split.split("|")
        for pattern, value in alignment_dict.items():
            part_a = "".join(pattern[int(s, base=num_taxa + 1)] for s in split_list[0])
            part_b = "".join(pattern[int(s, base=num_taxa + 1)] for s in split_list[1])
            vec_a = np.array([part_a.count(state) for state in states])
            vec_b = np.array([part_b.count(state) for state in states])
            vec_a = vec_a / np.linalg.norm(vec_a)
            vec_b = vec_b / np.linalg.norm(vec_b)
            scores[split] += value * (2 - np.linalg.norm(vec_a - vec_b)) / 2
    return scores


def __dense_split_score(matrix, k=None, singularValues=False, force_frob_norm=False):
    singular_values = list(
        scipy.linalg.svd(
            np.array(matrix), full_matrices=False, check_finite=False, compute_uv=False
        )
    )
    if force_frob_norm:
        return (
            1
            - (sum(val**2 for val in singular_values[0:4]))
            / (matrix.frobrenius_norm(matrix) ** 2)
        ) ** (1 / 2)
    else:
        min_shape = min(matrix.shape)
        return (
            1
            - (
                sum(val**2 for val in singular_values[0:4])
                / sum(val**2 for val in singular_values[0:min_shape])
            )
        ) ** (1 / 2)


def __sparse_split_score(
    matrix, return_singular_values=False, data_table_for_frob_norm=None
):
    largest_four_singular_values = scipy.sparse.linalg.svds(
        matrix, 4, return_singular_vectors=False
    )
    squared_singular_values = [sigma**2 for sigma in largest_four_singular_values]
    norm = matrix.frobenius_norm(matrix, data_table=data_table_for_frob_norm)
    return (1 - (sum(squared_singular_values) / (norm**2))) ** (1 / 2)


def split_score(
    matrix,
    return_singular_values=False,
    force_frob_norm_on_dense=False,
    data_table_for_frob_norm=None,
):
    if matrix.is_sparse(matrix):
        return __sparse_split_score(
            matrix, return_singular_values, data_table_for_frob_norm
        )
    else:
        return __dense_split_score(
            matrix, return_singular_values, force_frob_norm_on_dense
        )


def flattening_rank_1_approximation(flattening, return_vectors=False):
    r = np.array([sum(flattening)])
    c = np.array([sum(flattening.T)])
    if return_vectors:
        return r.T @ c, r.tolist()[0], c.tolist()[0]
    else:
        return r.T @ c


def flattening_rank_1_approximation_divergence(flattening):
    _, r, c = flattening_rank_1_approximation(flattening, return_vectors=True)
    total = 0
    for x in range(len(c)):
        for y in range(len(r)):
            if flattening[x, y] != 0:
                total += flattening[x, y] * np.log(flattening[x, y] / (r[y] * c[x]))
    return total


def star_tree(num_leaves):
    root_index = -1
    G = nx.DiGraph()
    G.add_node(root_index)
    # add the leaves as nodes and edges to the central node
    for i in range(0, num_leaves):
        G.add_node(i)
        G.add_edge(root_index, i)
    return G


def join_nodes(T, i, j, new_node, root_index):
    # Add new node
    T.add_node(new_node)
    T.add_edge(root_index, new_node)
    # Join nodes
    T.add_edge(new_node, i)
    T.add_edge(new_node, j)
    # Remove old edges
    T.remove_edge(root_index, i)
    T.remove_edge(root_index, j)


def neighbour_joining(distance_matrix, labels=None, return_newick=False):

    # Initialise
    D = distance_matrix.copy().astype(float)
    n = D.shape[0]
    new_node = n
    root_index = -1
    ignore = {-1}
    T = star_tree(n)
    num_leaves = n
    if labels is not None:
        # Add a label to each node
        for i in range(n):
            T.nodes[i]['label'] = labels[i]

    # NJ Algorithm
    while num_leaves > 2:
        # Instantiate Q matrix
        Q = np.full((n, n), np.inf)
        for i in range(n):
            for j in range(n):
                if i != j and i not in ignore and j not in ignore:
                    Q[i, j] = (num_leaves - 2) * D[i, j] - sum(D[i, :]) - sum(D[:, j])

        # Smallest Q pair
        i, j = np.unravel_index(Q.argmin(), Q.shape)
        ignore.add(i)
        ignore.add(j)

        # Join new nodes
        join_nodes(T, i, j, new_node, root_index)

        # Estimate branch lengths
        dist_new_to_i = (1 / 2) * D[i, j] + (1 / (2 * (num_leaves - 2))) * (
            sum(D[i, :]) - sum(D[j, :])
        )
        T.edges[new_node, i]["weight"] = dist_new_to_i
        T.edges[new_node, j]["weight"] = D[i, j] - dist_new_to_i

        # Append new row and column to distance matrix
        D = np.append(D, np.zeros((1, D.shape[0])), axis=0)
        D = np.append(D, np.zeros((D.shape[0], 1)), axis=1)
        n = D.shape[0]

        # Compute distance from other leaves to new node
        for k in range(n - 1):
            if k not in ignore:
                D[k, new_node] = (1 / 2) * (D[i, k] + D[j, k] - D[i, j])
                D[new_node, k] = D[k, new_node]

        # 'Delete' i and j by setting row and column to an array of zero for i and j
        D[i, :] = np.zeros(n)
        D[:, i] = np.zeros(n)
        D[j, :] = np.zeros(n)
        D[:, j] = np.zeros(n)

        new_node += 1
        num_leaves -= 1

    # Join the last two nodes
    i = new_node - 2
    j = new_node - 1
    # Remove the root node and connected edges
    T.remove_node(root_index)
    # Join nodes
    T.add_edge(i, j)
    # Add branch length
    T.edges[i, j]["weight"] = D[i, j]

    return T

def distance_matrix(networkx_tree):
    """Distance matrix of a tree.

    Args:
        networkx_tree (networkx.DiGraph): A tree.
    
    Returns:
        numpy.ndarray: A distance matrix.
    """
    # Get all the leaves
    leaf_nodes = [node for node in networkx_tree.nodes if networkx_tree.out_degree(node) == 0]
    # Create the distance matrix
    distance_matrix = np.zeros((len(leaf_nodes), len(leaf_nodes)))
    for i in range(len(leaf_nodes)):
        for j in range(i + 1, len(leaf_nodes)):
            distance = nx.shortest_path_length(networkx_tree.to_undirected(), leaf_nodes[i], leaf_nodes[j], weight="weight")
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    return distance_matrix

def midpoint_rooting(networkx_tree, weight_label="weight"):
    """Midpoint rooting of a tree.

    Args:
        networkx_tree (networkx.DiGraph): A tree.
    
    Returns:
        networkx.DiGraph: A rooted tree.
    """
    # Get all the leaves
    leaf_nodes = [node for node in networkx_tree.nodes if networkx_tree.out_degree(node) == 0]
    # Get the distance matrix
    D = distance_matrix(networkx_tree)
    # Get the index of the largest distance
    max_dist = np.max(D)
    i, j = np.unravel_index(np.argmax(D, axis=None), D.shape)
    # Get the undirected version of the tree
    tree_undirected = networkx_tree.to_undirected()
    # Get the path between the two leaves
    path = nx.shortest_path(tree_undirected, leaf_nodes[i], leaf_nodes[j])
    midpoint_dist = max_dist / 2 
    # Travel along the path until the midpoint is reached. Then go back and add a new node
    current_dist = 0
    prev_dist = 0
    print(path)
    for k in range(len(path) - 1):
        prev_dist = current_dist
        current_dist += tree_undirected[path[k]][path[k + 1]][weight_label]
        if current_dist >= midpoint_dist:
            # Add a new node
            new_node = -1
            networkx_tree.add_node(new_node)
            # Add the edges
            networkx_tree.add_edge(new_node, path[k])
            networkx_tree.add_edge(new_node, path[k + 1])
            # Remove the old edges
            if (path[k], path[k + 1]) in networkx_tree.edges:
                networkx_tree.remove_edge(path[k], path[k + 1])
            elif (path[k + 1], path[k]) in networkx_tree.edges:
                networkx_tree.remove_edge(path[k + 1], path[k])
            else:
                raise ValueError("Edge not found. Is tree already rooted?")
            # Add the branch lengths
            networkx_tree.edges[new_node, path[k]][weight_label] = current_dist - midpoint_dist
            networkx_tree.edges[new_node, path[k + 1]][weight_label] = midpoint_dist - prev_dist
            break





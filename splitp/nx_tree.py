import copy as cp
import itertools
from warnings import warn
from math import exp
from networkx import exception
# Numpy
import numpy as np
# Pandas
import pandas as pd
# NetworkX
import networkx as nx
from networkx.readwrite import json_graph
from networkx.algorithms.traversal.depth_first_search import dfs_tree, dfs_postorder_nodes
from networkx.algorithms.traversal.breadth_first_search import bfs_successors
# Scipy
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse import dok_matrix
from scipy.sparse.linalg import svds
# Random
import random
from random import choices
# SplitP
from splitp import tree_helper_functions as hf
from splitp import parsers
from splitp.enums import Model

class NXTree:
    """A rooted phylogenetic tree.

    A rooted phylogenetic tree consisting of a collection of node objects and an adjacency matrix

    Attributes:
        num_nodes: The number of nodes in the tree
        num_bases: The size of the character set for the phylogeny
        nodes: A list of node objects representing the nodes in the tree
        adjM: An adjacency matrix describing the connections within the tree
        initDist: The initial distribution of states/characters at the root of the tree.
    """

    def __init__(self, newickString, name=None, numStates=4, taxa_ordering=None):
        """Initialises a tree with the number of nodes and size of the character set"""
        self.num_bases = numStates
        self.initDist = [1/4,1/4,1/4,1/4]
        if self.num_bases == 4:
            self.state_space = ('A', 'C', 'G', 'T')
            self.state_space_dict = {char : self.state_space.index(char) for char in self.state_space}
            self.special_state = 'T'
        json_tree = parsers.newick_to_json(newickString, generate_names=True)
        self.nx_graph = json_graph.tree_graph(json_tree)
        self.subflatLRMats = {}
        # Check if branch lengths have been assigned for every edge:
        if all(('branch_length' in self.nx_graph.nodes[n]) or self.is_root(n) for n in self.nx_graph.nodes):
            print("Branch lengths assigned")
            for n in self.nx_graph.nodes:
                if not self.is_root(n):
                    b = self.nx_graph.nodes[n]['branch_length']
                    self.nx_graph.nodes[n]['transition_matrix'] = self.build_JC_matrix(b)
        else:
            print("Branch lengths missing: none were assigned")

        leaves = [n for n in self.nx_graph.nodes if self.is_leaf(n)]
        if not taxa_ordering:
            taxa_ordering = leaves
        elif taxa_ordering == 'sorted':
            taxa_ordering = sorted(leaves)
        self.taxa = taxa_ordering
        for i, taxon_name in enumerate(self.taxa):
            self.nx_graph.nodes[taxon_name]['t_index'] = i

        for i, n in enumerate(self.nx_graph.nodes):
            self.nx_graph.nodes[n]['index'] = i
        self.name = name

    def __str__(self):
        return parsers.json_to_newick(json_graph.tree_data(self.nx_graph, self.get_root(return_index=False)))

    def all_splits(self, trivial=False, only_balance=None, randomise=False):
        k = only_balance
        n = self.get_num_taxa()
        taxa_string = "".join(self.taxa)
        r = 0 if trivial else 1
        loop_over = range(r, 2**(n-1) - r)
        if randomise:
            loop_over = [i for i in loop_over]
            random.shuffle(loop_over)
        for i in loop_over:
            template = format(i, f'0{n}b')
            if not only_balance or sum(int(b) for b in template) in [only_balance, n-only_balance]:
                if r < sum(int(b) for b in template) < n-r:
                    left  = ""
                    right = ""
                    for t, b in zip(taxa_string, template):
                        if b == '0':
                            left += t
                        else:
                            right += t
                    if self.taxa[0] in right:
                        left, right = right, left
                    yield f'{left}|{right}'

    def true_splits(self, include_trivial=False):
        """Returns set of all true splits in the tree."""
        all_taxa_string = ''.join(str(i) for i in self.taxa)
        splits = set()
        for node in list(self.nodes()):
            split = sorted((node, ''.join(i for i in all_taxa_string if i not in node)))
            if include_trivial or (len(split[0])>1 and len(split[1])>1):
                splits.add(f'{split[0]}|{split[1]}')
        for split in splits:
            yield split

    ## BROKEN because this function doesn't realise that all_splits is a generator (I think!)
    # TODO: add a warning or someting until it's fixed
    def false_splits(self, only_balance=None, randomise=False):
        warn("This function is completely broken! Sorry!")
        """Returns set of all false splits in the tree."""
        true_splits = self.true_splits(include_trivial=False)
        for split in hf.all_splits(self.get_num_taxa(), trivial=False, only_balance=only_balance, randomise=randomise):
            if split not in true_splits:
                yield split

    def reassign_all_transition_matrices(self, matrix):
        for n in self.nx_graph.nodes:
            self.nx_graph.nodes[n]['transition_matrix'] = matrix
            warn("branch lengths have not been recalculated.")
            if 'branch_length' in self.nx_graph.nodes[n]: 
                self.nx_graph.nodes[n].pop('branch_length') # TODO: recompute branch lengths instead

    def build_JC_matrix(self, l):
        matrix = [[0 for i in range(self.num_bases)] for n in range(self.num_bases)]
        for r, row in enumerate(matrix):
            for c, _ in enumerate(row):
                matrix[r][c] = (1 / 4) + (3 / 4) * exp((-4 * l) / 3) if r == c else (1 / 4) - (1 / 4) * exp(
                    (-4 * l) / 3)
        return np.array(matrix).T
    
    def build_K2ST_matrix(self, transition, transversion):
        if self.num_bases != 4:
            warn(f"K2ST matrices are 4x4 but your model has {self.num_bases} states!" )
        if transversion > transition: 
            warn(f"transitions are known to be more likely than transversions!")
        purines = ('A', 'G')
        pyrimidines = ('C', 'T')
        matrix = [[0 for i in range(self.num_bases)] for n in range(self.num_bases)]
        for r, row in enumerate(matrix):
            from_state = self.state_space[r]
            for c, _ in enumerate(row):
                to_state = self.state_space[c]
                if from_state == to_state:
                    # No change
                    matrix[r][c] = 1-(transition+2*transversion)
                elif from_state in purines and to_state in purines:
                    matrix[r][c] = transition 
                elif from_state in pyrimidines and to_state in pyrimidines:
                    matrix[r][c] = transition
                else:
                    matrix[r][c] = transversion
        return np.array(matrix).T

    def rate_matrix(self, model):
        def _JC_rate_matrix(mutation_rate=None):
            if mutation_rate:
                a = mutation_rate
                return [[-3*a, a, a, a],
                        [a, -3*a, a, a],
                        [a, a, -3*a, a],
                        [a, a, a, -3*a]]
        def _K2ST_rate_matrix(rate_transition=None, rate_transversion=None, ratio=None):
            if self.num_bases != 4:
                warn(f"K2ST matrices are 4x4 but your model has {self.num_bases} states!" )
            purines = ('A', 'G')
            pyrimidines = ('C', 'T')
            matrix = [[0 for i in range(self.num_bases)] for n in range(self.num_bases)]
            if rate_transition and rate_transversion:
                transition = rate_transition
                transversion = rate_transversion
                if transversion > transition: 
                    warn(f"transitions are known to be more likely than transversions!")
            elif ratio:
                transition = ratio
                transversion = 1
            for r, row in enumerate(matrix):
                from_state = self.state_space[r]
                for c, _ in enumerate(row):
                    to_state = self.state_space[c]
                    if from_state == to_state:
                        # No change
                        matrix[r][c] = -(transition+2*transversion)
                    elif from_state in purines and to_state in purines:
                        matrix[r][c] = transition 
                    elif from_state in pyrimidines and to_state in pyrimidines:
                        matrix[r][c] = transition
                    else:
                        matrix[r][c] = transversion
            return np.array(matrix).T
        if   model is Model.JC:   return _JC_rate_matrix
        elif model is Model.K2ST: return _K2ST_rate_matrix

    def scale_TR_rate_matrix(self, Q, return_scale_factor=False):
        scale_factor = 1/(-sum(self.initDist[i]*Q[i][i] for i in range(4)))
        if return_scale_factor:
            return scale_factor
        else:
            return scale_factor*np.array(Q)

    def adjacency_matrix(self):
        return np.array(nx.adjacency_matrix(self.nx_graph).todense())
    
    def add_node(self, n, branch_length=0, in_node=None):
        """Adds a new node to the tree

        Args:
            n: The node object to add into the tree.
            in_node: The name of the parent node, default is None and is used for the root.
        """
        warn("addNode is deprecated, trees should be instantiated from newick strings and never changed.", DeprecationWarning)
        self.nx_graph.add_node(n,
                               branch_length=branch_length,
                               transition_matrix=self.build_JC_matrix(branch_length)
                               )
        if in_node:
            self.nx_graph.add_edge(in_node, n)

    def get_num_taxa(self):
        """Returns the number of taxa/leaf-nodes in the tree"""
        num_leaves = 0
        for i, n in enumerate(self.nodes()):
            if self.is_leaf(i):
                num_leaves += 1
        return num_leaves

    def num_nodes(self):
        return len(list(self.nx_graph.nodes))

    def node_index(self, n):
        return self.nx_graph.nodes[n]['index']

    def nodes_list(self):
        return list(self.nx_graph.nodes)

    def nodes(self):
        return self.nx_graph.nodes

    def is_leaf(self, n_index_or_name):
        """Determines whether a node is a leaf node from it's index."""
        if type(n_index_or_name) == type(str()):
            return self.nx_graph.out_degree(n_index_or_name) == 0
        else:
            return self.nx_graph.out_degree(list(self.nx_graph.nodes)[n_index_or_name]) == 0

    def is_root(self, n_index_or_name):
        """Determines whether a node is a root node from it's index."""
        if type(n_index_or_name) == type(str()):
            return self.nx_graph.in_degree(n_index_or_name) == 0
        else:
            return self.nx_graph.in_degree(list(self.nx_graph.nodes)[n_index_or_name]) == 0

    def get_root(self, return_index=True):
        """Returns the root node"""
        for i, n in enumerate(self.nx_graph.nodes):
            if self.is_root(n):
                return i if return_index else n

    def index_of_node(self, node_name):
        return self.nodes_list().index(node_name)

    def node_name(self, index):
        return self.nodes_list()[index]

    def get_parent(self, n):
        """Returns the parent node for a given node"""
        return list(self.nx_graph.predecessors(n))[0]

    def get_descendants(self, n, return_iter=False):
        """Returns a list of children/descendents of a given node"""
        return list(self.nx_graph.successors(n)) if not return_iter else self.nx_graph.successors(n)

    def __likelihood(self, n, lTable):
        """Recursive part of the likelihood algorithm"""
        for b in range(self.num_bases):
            L = 1
            for d in self.get_descendants(n, return_iter=True):
                d_index = self.node_index(d)
                M = (self.nx_graph.nodes[d])['transition_matrix']
                s = 0
                for b2 in range(self.num_bases):
                    if lTable[d_index, b2] == None:
                        self.__likelihood(d, lTable)
                    s += M[b2, b] * lTable[d_index, b2]
                L *= s
            lTable[self.node_index(n), b] = L
        if not self.is_root(n):
            self.__likelihood(self.get_parent(n), lTable)

    def __likelihood_start(self, pattern):
        """Starts the likelihood algorithm.

        Starts the likelihood algorithm to determine the probability of seeing a given site pattern.

        Args:
            pattern: The site pattern to obtain the probability for.

        Returns:
            A probability (range 0-1) of observing the given site pattern given the tree.
        """

        def _to_int(p):
            return self.state_space.index(p)

        pattern = [_to_int(p) for p in pattern]  # A list of indices which correspond to taxa.
        taxa = self.taxa  # The list of taxa.
        # Likelihood table for dynamic prog alg ~ lTable[node_index, character]
        lTable = np.array([[None for _ in range(self.num_bases)] for _ in range(self.num_nodes())])

        # Encoding pattern into likelihood table
        for i, p in enumerate(pattern):
            lTable[self.node_index(taxa[i]), :] = [0 for _ in range(self.num_bases)]
            lTable[self.node_index(taxa[i]), p] = 1

        # Starting with node which has all descendants as leaves
        initNode = None
        for n in self.nx_graph.nodes:
            desc = self.get_descendants(n)
            if desc and all(self.is_leaf(d) for d in desc):  # If the node has descendants and they are all leaves
                initNode = n

        self.__likelihood(initNode, lTable)  # Should fill in the entire table

        root_index = self.node_index(self.get_root(return_index=False))
        L = 0
        for i in range(self.num_bases):
            L += self.initDist[i] * lTable[root_index, i]

        return L

    def svd_error(self, M):
        """"Returns the SVD for a given matrix (All but two/four largest SVs)"""
        warn("svd_error is deprecated, use split_score() for a much faster result", DeprecationWarning)
        sv = list(np.linalg.svd(np.array(M).astype(np.float64), compute_uv=False))
        sv[0] = 0
        sv[1] = 0
        if self.num_bases == 4:
            sv[2] = 0
            sv[3] = 0
        return sum(sv)

    def all_singular_values(self, M):
        warn("all_singular_values is deprecated, use split_score() for a much faster result", DeprecationWarning)
        return (list(np.linalg.svd(np.array(M).astype(np.float64), compute_uv=False)))
    
    def __dense_split_score(self, matrix, k=None, singularValues=False, force_frob_norm=False):
        singular_values = list(scipy.linalg.svd(np.array(matrix), full_matrices=False, check_finite=False, compute_uv=False))
        if force_frob_norm:
            return (1-(sum(val**2 for val in singular_values[0:4]))/(hf.frob_norm(matrix)**2))**(1/2)
        else:
            min_shape = min(matrix.shape)
            return (1-(sum(val**2 for val in singular_values[0:4])/sum(val**2 for val in singular_values[0:min_shape])))**(1/2)

    def __sparse_split_score(self, matrix, return_singular_values=False, data_table_for_frob_norm=None):
        largest_four_singular_values = svds(matrix, 4, return_singular_vectors=False)
        squared_singular_values = [sigma**2 for sigma in largest_four_singular_values]
        norm = hf.frob_norm(matrix, data_table=data_table_for_frob_norm)
        return (1-(sum(squared_singular_values)/(norm**2)))**(1/2)

    def split_score(self, matrix, return_singular_values=False, force_frob_norm_on_dense=False, data_table_for_frob_norm=None):
        if hf.is_sparse(matrix):
            return self.__sparse_split_score(matrix, return_singular_values, data_table_for_frob_norm)
        else:
            return self.__dense_split_score(matrix, return_singular_values, force_frob_norm_on_dense)

    def get_pattern_probabilities(self):
        """Returns a full table of site-pattern probabilities (binary character set)"""
        # Creating a table with binary labels and calling likelihood_start() to fill it in with probabilities
        if self.num_bases == 2:
            emptyArray = np.array(
                [
                    ['{:0{}b}'.format(i, self.get_num_taxa()),
                     self.__likelihood_start('{:0{}b}'.format(i, self.get_num_taxa()))]
                    for i in range(self.num_bases ** self.get_num_taxa())
                ]
            )
        elif self.num_bases == 4:
            combinations = list(itertools.product(''.join(s for s in self.state_space), repeat=self.get_num_taxa()))
            combinations = [''.join(c) for c in combinations]
            emptyArray = pd.DataFrame(
                [
                    [combinations[i], self.__likelihood_start(combinations[i])]
                    for i in range(len(combinations))
                ]
            )
        return emptyArray

    def evolve_pattern(self):
        def __evolve_on_subtree(subtree, state):
            root_node = [n for n,d in subtree.in_degree() if d==0][0]
            children = list(subtree.successors(root_node))
            probs = list(self.nx_graph.nodes[root_node]['transition_matrix'][:,self.state_space.index(state)])
            new_state = choices(self.state_space, probs)[0]
            if len(children) == 0:
                return f'{str(root_node)}:{new_state}'
            else:
                subtrees = [dfs_tree(self.nx_graph, child) for child in children]
                return ','.join(__evolve_on_subtree(subtree, new_state) for subtree in subtrees)
        root_state = choices(self.state_space, self.initDist)[0]
        root_node = [n for n,d in self.nx_graph.in_degree() if d==0][0]
        children = list(self.nx_graph.successors(root_node))
        subtrees = [dfs_tree(self.nx_graph, child) for child in children]
        result_string = ','.join(__evolve_on_subtree(subtree, root_state) for subtree in subtrees)
        result = { pair[0] : pair[2] for pair in result_string.split(',') }
        return ''.join(result[k] for k in sorted(result.keys(), key=self.taxa.index))

    def generate_alignment(self, sequence_length):
        counts = {}
        for i in range(sequence_length):
            pattern = self.evolve_pattern()
            if pattern not in counts:
                counts[pattern] = float(1)
            else:
                counts[pattern] += 1
        probs = {}
        for k in sorted(counts.keys(), key=lambda p: [self.state_space.index(c) for c in p]):
            probs[k] = counts[k]/float(sequence_length)
        return pd.DataFrame(probs.items())

    def draw_from_multinomial(self, LT, n):
        """Use a given table of probabilities from getLikelihoods() and draw from its distribution"""
        probs = list(LT.iloc[:, 1])
        probs = [float(p) for p in probs]
        data = np.random.multinomial(n, probs)
        data = [d / n for d in data]
        if self.num_bases == 2:
            return np.array(
                [['{:0{}b}'.format(i, self.get_num_taxa()), data[i]] for i in range(self.num_bases ** self.get_num_taxa())])
        elif self.num_bases == 4:
            patterns = list(LT.iloc[:, 0])
            results = []
            for i in range(len(patterns)):
                if data[i] != 0:
                    results.append([patterns[i], data[i]])
            return pd.DataFrame(results)

    def __index_of(self, string):
        string = reversed(string)
        index = 0
        for o, s in enumerate(string):
            index += (4**o)*self.state_space_dict[s]
        return index
        
    def __reconstruct_pattern(self, split, row_label, col_label):
        n = len(split[0]) + len(split[1])
        pattern = {}
        for splindex, loc in enumerate(split[0]):
            pattern[int(loc, n)] = row_label[splindex]
        for splindex, loc in enumerate(split[1]):
            pattern[int(loc, n)] = col_label[splindex]
        return "".join(pattern[i] for i in range(n))
    
    def __subflattening_labels(self, length):
        #TODO: make into a generator and don't iterate twice for templates list
        n = length
        templates = [
            "".join(self.special_state for _ in range(i)) + '*' + "".join(self.special_state for _ in range(n-i - 1))
            for i in range(n)
        ]
        patterns = []
        for template in templates:
            for c in self.state_space:
                if c != self.special_state:
                    patterns.append(template.replace('*', c))
        patterns.append("".join(self.special_state for _ in range(n)))
        patterns.sort(key=self.__index_of)
        return patterns

    def __subflattening_labels_generator(self, length):
        n = length
        other_states = self.state_space[0:-1]
        templates = (
            ("".join(self.special_state for _ in range(i)), "".join(self.special_state for _ in range(n-i - 1)))
            for i in range(n)
        )
        patterns = []
        for template in templates:
            for c in other_states:
                yield f'{template[0]}{c}{template[1]}'
        yield "".join(self.special_state for _ in range(n))

    def sparse_flattening(self, split, table, format='dok'):
        split = split.split('|')
        num_taxa = sum(len(part) for part in split)
        if format == 'coo':
            rows = []
            cols = []
            data = []
            for r in table.itertuples(index=False, name=None):
                if r[1] != 0:
                    pattern = r[0]
                    row = self.__index_of(''.join([str(pattern[int(s, num_taxa)]) for s in split[0]]))
                    col = self.__index_of(''.join([str(pattern[int(s, num_taxa)]) for s in split[1]]))
                    rows.append(row)
                    cols.append(col)
                    data.append(r[1])
            return coo_matrix((data, (rows, cols)), shape=(4**len(split[0]),4**len(split[1])))
        elif format == 'dok':
            flattening = dok_matrix((4**len(split[0]),4**len(split[1])))
            for r in table.itertuples(index=False, name=None):
                pattern = r[0]
                row = self.__index_of(''.join([str(pattern[int(s, num_taxa)]) for s in split[0]]))
                col = self.__index_of(''.join([str(pattern[int(s, num_taxa)]) for s in split[1]]))
                flattening[row, col] = r[1]
            return flattening

    def fast_sparse_flattening(self, split, table, format='dok'):
        """A faster version of signed sum subflattening. Requires a data dictionary and can be supplied with a bundle of 
        re-usable information to reduce the number of calls to the multiplications function.
        """
        split = split.split('|')
        num_taxa = sum(len(part) for part in split)
        if format == 'coo':
            rows = []
            cols = []
            data = []
            for r in table.itertuples(index=False, name=None):
                if r[1] != 0:
                    pattern = r[0]
                    row = self.__index_of(''.join([str(pattern[int(s, num_taxa)]) for s in split[0]]))
                    col = self.__index_of(''.join([str(pattern[int(s, num_taxa)]) for s in split[1]]))
                    rows.append(row)
                    cols.append(col)
                    data.append(r[1])
            return coo_matrix((data, (rows, cols)), shape=(4**len(split[0]),4**len(split[1])))
        elif format == 'dok':
            flattening = dok_matrix((4**len(split[0]),4**len(split[1])))
            for r in table.itertuples(index=False, name=None):
                pattern = r[0]
                row = self.__index_of(''.join([str(pattern[int(s, num_taxa)]) for s in split[0]]))
                col = self.__index_of(''.join([str(pattern[int(s, num_taxa)]) for s in split[1]]))
                flattening[row, col] = r[1]
            return flattening

    #def subflattening(self, split, data, build_from='flattening', return_sparse=False):
     
    def signed_sum_subflattening(self, split, data_table):
        data_dict = {}
        for table_pattern, value in data_table.itertuples(index=False, name=None):
            data_dict[table_pattern] = value
        split = split.split('|')
        num_taxa = sum(len(part) for part in split)
        subflattening = [[0 for i in range(3*len(split[1])+1)] for j in range(3*len(split[0])+1)]
        H = np.array([[1, -1], [1, 1]])
        S = np.kron(H, H)
        S = { (c1, c2) : S[self.state_space.index(c1)][self.state_space.index(c2)] for c1 in self.state_space for c2 in self.state_space}
        row_labels = self.__subflattening_labels(len(split[0]))
        col_labels = self.__subflattening_labels(len(split[1]))
        rec_pat = self.__reconstruct_pattern
        items = data_dict.items
        # SOME OF THIS CAN BE RE-USED!!!!
        for row in range(len(row_labels)):
            for col in range(len(col_labels)):
                pattern = rec_pat(split, row_labels[row], col_labels[col])
                signed_sum = 0
                for table_pattern, value in items():
                    product = value
                    for t in zip(pattern, table_pattern):
                        product *= S[t]
                    signed_sum += product
                subflattening[row][col] = signed_sum
        return np.array(subflattening)

    def fast_signed_sum_subflattening(self, split, data_dict, bundle=None, labels=None):
        """
        A faster version of signed sum subflattening. Requires a data dictionary and can be supplied with a bundle of 
        re-usable information to reduce the number of calls to the multiplications function.
        """
        split = split.split('|')
        sp1, sp2 = len(split[0]), len(split[1])
        subflattening = [[0 for _ in range(3*sp2+1)] for _ in range(3*sp1+1)]
        try:
            row_labels = labels[sp1]
        except KeyError:
            row_labels = list(self.__subflattening_labels_generator(sp1))
            labels[sp1] = row_labels
        try:
            col_labels = labels[sp2]
        except KeyError:
            col_labels = list(self.__subflattening_labels_generator(sp2))
            labels[sp2] = col_labels
        banned = {('C','C'), ('G','G'), ('A','T')} | {(x, 'A') for x in self.state_space} | {('T', x) for x in self.state_space}
        for r, row in enumerate(row_labels):
            for c, col in enumerate(col_labels):
                pattern = self.__reconstruct_pattern(split, row, col)
                signed_sum = 0
                for table_pattern, value in data_dict.items():
                    try:
                        product = bundle[(pattern, table_pattern)]
                    except KeyError:
                        product = 1
                        for t in zip(pattern, table_pattern):
                            if t not in banned:
                                product *= -1
                        bundle[(pattern, table_pattern)] = product
                    signed_sum += product*value
                subflattening[r][c] = signed_sum
        return np.array(subflattening)

    def sparse_subflattening(self, split, data_table, as_dense_array=False):
        split = split.split('|')
        num_taxa = len(split[0]) + len(split[1])
        H = np.array([[1, -1], [1, 1]])
        S = np.kron(H, H)
        rows = []
        cols = []
        data = []
        row_labels = self.__subflattening_labels(len(split[0])) #then these can go
        col_labels = self.__subflattening_labels(len(split[1]))
        for row_label in row_labels: # loop over generators instead of lists :)
            for col_label in col_labels: 
                pattern = self.__reconstruct_pattern(split, row_label, col_label) 
                rows.append(row_labels.index(row_label))
                cols.append(col_labels.index(col_label))
                signed_sum = 0
                for r in data_table.itertuples(index=True, name=None):
                    table_pattern = r[1]
                    value = r[2]
                    product = 1
                    for p, tp in zip(pattern, table_pattern):
                        product *= S[self.state_space.index(p)][self.state_space.index(tp)]
                    product *= value
                    signed_sum += product
                data.append(signed_sum)
        subflattening = coo_matrix((data, (rows, cols)), shape=(3*len(split[0])+1, 3*len(split[1])+1))
        if as_dense_array:
            return subflattening.toarray()
        else:
            return subflattening.todok()
        
    def flattening(self, split, table):
        """Build a flattening matrix from a split

        Args:
            split: A string representing the split to build the flattening from
            table: A pandas data-frame table with site pattern `probabilities'

        Returns:
            A flattening matrix as a data-frame (to provide labels)
        """
        split = split.split('|')

        if self.num_bases == 2:
            rowLables = ['{:0{}b}'.format(i, len(split[0])) for i in range(2 ** len(split[0]))]
            colLables = ['{:0{}b}'.format(i, len(split[1])) for i in range(2 ** len(split[1]))]
        elif self.num_bases == 4:
            rowLables = list(map(''.join, list(itertools.product(''.join(self.state_space), repeat=len(split[0])))))
            colLables = list(map(''.join, list(itertools.product(''.join(self.state_space), repeat=len(split[1])))))

        F = pd.DataFrame(0, index=rowLables, columns=colLables, dtype=np.float64)
        # Searching through data table and moving likelihoods to the F matrix
        for r in table.itertuples(index=True, name='Pandas'):
            if r[2] != '0.0':
                pattern = r[1]
                row = ''.join([str(pattern[int(s)]) for s in split[0]])
                col = ''.join([str(pattern[int(s)]) for s in split[1]])
                F.loc[row, col] = r[2]

        return F  # .astype(pd.SparseDtype("float64", 0))

    def __make_mat(self, F, S_hat, left=True):
        colLabels = list(F)
        rowLabels = F.index
        if (str(S_hat), left, F.shape) in self.subflatLRMats:
            return self.subflatLRMats[(str(S_hat), left, F.shape)]
        if left:
            M = np.empty((0, F.shape[0]))
            rows = len(rowLabels[0])
        else:
            M = np.empty((0, F.shape[1]))
            rows = len(colLabels[0])
        O = np.ones(len(S_hat) + 1)
        for r in range(rows):
            A = np.empty((0, 0))
            if r == 0:  # First row of M
                A = S_hat
                for i in range(rows - 1):
                    A = np.kron(A, O)
            else:
                A = O
                for i in range(0, r - 1):
                    A = np.kron(A, O)
                A = np.kron(A, S_hat)
                for i in range(rows - r - 1):
                    A = np.kron(A, O)

            M = np.append(M, A, axis=0)

        rowOfOnes = O
        for i in range(rows - 1):
            rowOfOnes = np.kron(rowOfOnes, O)
        rowOfOnes = rowOfOnes.reshape((1, -1))
        M = np.append(M, rowOfOnes, axis=0)
        self.subflatLRMats[(str(S_hat), left, F.shape)] = M
        return M

    def subflattening_alt(self, F, S=None, returnLRMats=False):
        if np.all(S) == None:
            S = np.array([[1, -1], [1, 1]])  # Swapped rows compared to other way
            S = np.kron(S, S)

        # Remove the constant row to obtain S
        S_hat = np.empty((0, self.num_bases))
        for row in S:
            if list(np.asarray(row)) != [1 for i in range(self.num_bases)]:
                S_hat = np.append(S_hat, [row], axis=0)

        L = self.__make_mat(F, S_hat, True)
        R = self.__make_mat(F, S_hat, False)
        F = np.asarray(F)  # , dtype=np.float64)
        if not returnLRMats:
            return L @ F @ R.T
        else:
            return (L @ F @ R.T, L, R)

    def transformed_flattening(self, F, S=None):
        """ Creates a transformed flattening from a flattening data frame """
        if np.all(S) == None:
            H = np.array([[1, -1], [1, 1]])  # Swapped rows compared to other way
            S = np.kron(H, H)

        colLabels = list(F)
        rowLabels = F.index
        F = np.array(F).astype(np.float64)
        L, R = S, S
        while L.shape[0] != F.shape[0]:
            L = np.kron(L, S)
        while R.shape[1] != F.shape[1]:
            R = np.kron(R, S)
        Ft = np.matmul(np.matmul(L, F), np.transpose(R))
        return pd.DataFrame(Ft, index=rowLabels, columns=colLabels, dtype=np.float64)

    def subflattening(self, Ft, specialState='T', type=(1, 1)):
        """ Creates a subflattening from a transformed flattening data frame """
        matrix = []
        self.special_state = specialState
        if self.num_bases == 2:
            for r in list(Ft.index):
                row = []
                for c in list(Ft):
                    if r.count('1') <= 1 and c.count('1') <= 1:
                        row.append(Ft.loc[r, c])
                if row != []: matrix.append(row)
        elif self.num_bases == 4:
            for r in list(Ft.index):
                row = []
                for c in list(Ft):
                    # Could be optimised
                    if len("".join([x for x in r if x != specialState])) <= type[0] and len(
                            "".join([x for x in c if x != specialState])) <= type[1]:
                        row.append(Ft.loc[r, c])
                if row: matrix.append(row)
        return np.asarray(matrix)

    def hartigan_algorithm(self, pattern):
        score = 0
        graph = self.nx_graph.copy()
        nodes = graph.nodes
        taxa = [t for t in self.taxa]
        for i, t in enumerate(taxa):
            nodes[t]['S1'] = set(pattern[i])
        postorder_nodes = list(dfs_postorder_nodes(graph, source=self.get_root(return_index=False)))
        for n in postorder_nodes:
            if not self.is_leaf(n):
                children = self.get_descendants(n)
                k={}
                for state in self.state_space:
                    k[state] = len(set(child for child in children if state in nodes[child]['S1']))
                K = max(k.values())
                nodes[n]['S1'] = {state for state in self.state_space if k[state] == K}
                nodes[n]['S2'] = {state for state in self.state_space if k[state] == K-1}
        # Now compute the score
        top_to_bottom_nodes = [x[0] for x in bfs_successors(graph, source=self.get_root(return_index=False))] + taxa
        for n in top_to_bottom_nodes:
            if n == self.get_root(return_index=False):
                nodes[n]['hart_state'] = list(nodes[n]['S1'])[0]
            else:
                parent = nodes[list(graph.predecessors(n))[0]]
                if parent['hart_state'] not in nodes[n]['S1']:
                    nodes[n]['hart_state'] = list(nodes[n]['S1'])[0]
                    score += 1
                else:
                    nodes[n]['hart_state'] = parent['hart_state']
        return score


    def parsimony_score(self, pattern):
        """Calculate a parsimony score for a site pattern or split

        Args:
            pattern: A string representing a site pattern or split

        Returns:
            A parsimony score for the given split or site pattern
        """
        graph = self.nx_graph.copy()
        nodes = graph.nodes
        if '|' in pattern:
            pattern2 = [0 for x in range(len(pattern) - 1)]
            i = 0
            while pattern[i] != '|':
                pattern2[int(pattern[i])] = 1
                i += 1
            pattern = "".join(str(i) for i in pattern2)

        taxa = [t for t in self.taxa]
        for i, t in enumerate(taxa):
            nodes[t]['pars'] = set(pattern[i])
        score = 0
        for n in nodes:
            if not self.is_leaf(n) and 'pars' not in nodes[n]:
                score += self.__parsimony(n, nodes=nodes)
        return score

    def __parsimony(self, n, nodes=None):
        """Recursive step in Finch's Algorithm"""
        score = 0
        children = self.get_descendants(n)
        for c in children:
            if 'pars' not in nodes[c]:
                score += self.__parsimony(c, nodes)

        nodes[n]['pars'] = nodes[children[0]]['pars'] # children[0].pars
        for c in children:
            nodes[n]['pars'] = nodes[n]['pars'].intersection(nodes[c]['pars'])
        if nodes[n]['pars'] == set():
            nodes[n]['pars'] = nodes[children[0]]['pars']
            for c in children:
                nodes[n]['pars'] = (nodes[n]['pars']).union(nodes[c]['pars'])
            return score + 1
        return score
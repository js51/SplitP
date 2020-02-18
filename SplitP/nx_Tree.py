import numpy as np
import copy as cp
import networkx as nx
from networkx.readwrite import json_graph
import pandas as pd
import itertools
from SplitP import TreeHelperFunctions as hf
from SplitP import Parsers


class nxtree:
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
        json_tree = Parsers.newick_to_json(newickString, generate_names=True)
        self.nx_graph = json_graph.tree_graph(json_tree)
        # Check if branch lengths have been assigned for every edge:

        if all(('branch_length' in self.nx_graph.nodes[n]) or self.isRoot(n) for n in self.nx_graph.nodes):
            print("Branch lengths assigned")
            for n in self.nx_graph.nodes:
                if not self.isRoot(n):
                    b = self.nx_graph.nodes[n]['branch_length']
                    self.nx_graph.nodes[n]['transition_matrix'] = self.build_JC_matrix(b)
        else:
            print("Branch lengths missing: none were assigned")

        if not taxa_ordering:
            taxa_ordering = [n for n in self.nx_graph.nodes if self.isLeaf(n)]
        self.taxa = taxa_ordering
        for i, taxon_name in enumerate(self.taxa):
            self.nx_graph.nodes[taxon_name]['t_index'] = i

        for i, n in enumerate(self.nx_graph.nodes):
            self.nx_graph.nodes[n]['index'] = i
        self.name = name

    def __str__(self):
        return Parsers.json_to_newick(json_graph.tree_data(self.nx_graph, self.getRoot(return_index=False)))

    def reassign_all_transition_matrices(self, matrix):
        for n in self.nx_graph.nodes:
            self.nx_graph.nodes[n]['transition_matrix'] = matrix
        # Recompute branch lengths ??

    def build_JC_matrix(self, l):
        from math import exp
        matrix = [[0 for i in range(self.num_bases)] for n in range(self.num_bases)]
        for r, row in enumerate(matrix):
            for c, _ in enumerate(row):
                matrix[r][c] = (1 / 4) + (3 / 4) * exp((-4 * l) / 3) if r == c else (1 / 4) - (1 / 4) * exp(
                    (-4 * l) / 3)
        return np.array(matrix).T

    def addNode(self, n, branch_length=0, in_node=None):
        """Adds a new node to the tree

        Args:
            n: The node object to add into the tree.
            in_node: The name of the parent node, default is None and is used for the root.
        """
        self.nx_graph.add_node(n,
                               branch_length=branch_length,
                               transition_matrix=self.build_JC_matrix(branch_length)
                               )
        if in_node:
            self.nx_graph.add_edge(in_node, n)

    def getNumTaxa(self):
        """Returns the number of taxa/leaf-nodes in the tree"""
        num_leaves = 0
        for i, n in enumerate(self.nodes_collection()):
            if self.isLeaf(i):
                num_leaves += 1
        return num_leaves

    def num_nodes(self):
        return len(list(self.nx_graph.nodes))

    def node_index(self, n):
        return self.nx_graph.nodes[n]['index']

    def nodes_list(self):
        return list(self.nx_graph.nodes)

    def nodes_collection(self):
        return self.nx_graph.nodes

    def isLeaf(self, n_index_or_name):
        """Determines whether a node is a leaf node from it's index."""
        if type(n_index_or_name) == type(str()):
            return self.nx_graph.out_degree(n_index_or_name) == 0
        else:
            return self.nx_graph.out_degree(list(self.nx_graph.nodes)[n_index_or_name]) == 0

    def isRoot(self, n_index_or_name):
        """Determines whether a node is a root node from it's index."""
        if type(n_index_or_name) == type(str()):
            return self.nx_graph.in_degree(n_index_or_name) == 0
        else:
            return self.nx_graph.in_degree(list(self.nx_graph.nodes)[n_index_or_name]) == 0

    def getRoot(self, return_index=True):
        """Returns the root node"""
        for i, n in enumerate(self.nodes_collection()):
            if self.isRoot(n):
                return i if return_index else n

    def index_of_node(self, node_name):
        return self.nodes_list().index(node_name)

    def node_name(self, index):
        return self.nodes_list()[index]

    def getParent(self, n):
        """Returns the parent node for a given node"""
        return list(self.nx_graph.predecessors(n))[0]

    def getDescendants(self, n, return_iter=False):
        """Returns a list of children/descendents of a given node"""
        return list(self.nx_graph.successors(n)) if not return_iter else self.nx_graph.successors(n)

    def __likelihood(self, n, lTable):
        """Recursive part of the likelihood algorithm"""
        for b in range(self.num_bases):
            L = 1
            for d in self.getDescendants(n, return_iter=True):
                d_index = self.node_index(d)
                M = (self.nx_graph.nodes[d])['transition_matrix']
                s = 0
                for b2 in range(self.num_bases):
                    if lTable[d_index, b2] == None:
                        self.__likelihood(d, lTable)
                    s += M[b2, b] * lTable[d_index, b2]
                L *= s
            lTable[self.node_index(n), b] = L
        if not self.isRoot(n):
            self.__likelihood(self.getParent(n), lTable)

    def __likelihood_start(self, pattern):

        """Starts the likelihood algorithm.

        Starts the likelihood algorithm to determine the probability of seeing a given site pattern.

        Args:
            pattern: The site pattern to obtain the probability for.

        Returns:
            A probability (range 0-1) of observing the given site pattern given the tree.
        """

        def toInt(p):
            return self.state_space.index(p)

        pattern = [toInt(p) for p in pattern]  # A list of indices which correspond to taxa.
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
            desc = self.getDescendants(n)
            if desc and all(self.isLeaf(d) for d in desc):  # If the node has descendants and they are all leaves
                initNode = n

        self.__likelihood(initNode, lTable)  # Should fill in the entire table

        root_index = self.node_index(self.getRoot(return_index=False))
        L = 0
        for i in range(self.num_bases):
            L += self.initDist[i] * lTable[root_index, i]

        return L

    def SVDError(self, M):
        """"Returns the SVD for a given matrix (All but two/four largest SVs)"""
        sv = list(np.linalg.svd(np.array(M).astype(np.float64), compute_uv=False))
        sv[0] = 0
        sv[1] = 0
        if self.num_bases == 4:
            sv[2] = 0
            sv[3] = 0
        return sum(sv)

    def allSVs(self, M):
        return (list(np.linalg.svd(np.array(M).astype(np.float64), compute_uv=False)))

    def splitScore(self, M, k=None, singularValues=False):
        """Returns the split score for a given flattening matrix"""
        sVals = list(np.linalg.svd(np.array(M).astype(np.float64), compute_uv=False))
        if singularValues: sing_vals = sVals.copy()
        sumSq = 0
        for i in range((self.num_bases if not k else k), min(M.shape)):
            sumSq += (sVals[i]) ** 2
        sumSq2 = sumSq
        for i in range(1, (self.num_bases if not k else k)):
            sumSq2 += (sVals[i]) ** 2

        score = np.sqrt(sumSq) / hf.fNorm(np.array(M).astype(np.float64))

        if singularValues:
            return score, sing_vals
        return score

    def getLikelihoods(self):
        """Returns a full table of site-pattern probabilities (binary character set)"""
        # Creating a table with binary labels and calling likelihood_start() to fill it in with probabilities
        if self.num_bases == 2:
            emptyArray = np.array(
                [
                    ['{:0{}b}'.format(i, self.getNumTaxa()),
                     self.__likelihood_start('{:0{}b}'.format(i, self.getNumTaxa()))]
                    for i in range(self.num_bases ** self.getNumTaxa())
                ]
            )
        elif self.num_bases == 4:
            combinations = list(itertools.product(''.join(s for s in self.state_space), repeat=self.getNumTaxa()))
            combinations = [''.join(c) for c in combinations]
            emptyArray = pd.DataFrame(
                [
                    [combinations[i], self.__likelihood_start(combinations[i])]
                    for i in range(len(combinations))
                ]
            )

        return emptyArray

    def drawFromMultinomial(self, LT, n):
        """Use a given table of probabilities from getLikelihoods() and draw from its distribution"""
        probs = list(LT.iloc[:, 1])
        probs = [float(p) for p in probs]
        data = np.random.multinomial(n, probs)
        data = [d / n for d in data]
        if self.num_bases == 2:
            return np.array(
                [['{:0{}b}'.format(i, self.getNumTaxa()), data[i]] for i in range(self.num_bases ** self.getNumTaxa())])
        elif self.num_bases == 4:
            patterns = list(LT.iloc[:, 0])
            results = []
            for i in range(len(patterns)):
                if data[i] != 0:
                    results.append([patterns[i], data[i]])
            return pd.DataFrame(results)

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
            rowLables = list(map(''.join, list(itertools.product('ACGT', repeat=len(split[0])))))
            colLables = list(map(''.join, list(itertools.product('ACGT', repeat=len(split[1])))))

        F = pd.DataFrame(0, index=rowLables, columns=colLables, dtype=np.float64)
        # Searching through data table and moving likelihoods to the F matrix
        for r in table.itertuples(index=True, name='Pandas'):
            if r[2] != '0.0':
                pattern = r[1]
                row = ''.join([str(pattern[int(s)]) for s in split[0]])
                col = ''.join([str(pattern[int(s)]) for s in split[1]])
                F.loc[row, col] = r[2]

        return F  # .astype(pd.SparseDtype("float64", 0))

    def toXML(self):
        """ Returns a string-representation of the tree """
        h = 0
        tree_text = '<?xml version="1.0" encoding="UTF-8"?>\n<tree>\n'
        r = self.getRoot()
        tree_text = tree_text + '<' + str(r) + '>\n'
        for c in self.getDescendants(r):
            tree_text = tree_text + self.__nodeXML(self.nodes[c], h + 1)
        tree_text = tree_text + '</' + str(r) + '>\n</tree>'
        return tree_text

    def __nodeXML(self, n, h):
        text = ''
        tabs = ''.join('\t' for x in range(h))
        children = self.getDescendants(n)
        if not children:
            text = text + tabs + '<' + str(n) + '/>\n'
        else:
            text = text + tabs + '<' + str(n) + '>\n'
            for c in self.getDescendants(n):
                text = text + self.__nodeXML(self.nodes[c], h + 1)
            text = text + tabs + '</' + str(n) + '>\n'
        return text

    def __makeMat(self, F, S_hat, left=True):
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

    def subFlatteningAlt(self, F, S=None, returnLRMats=False):
        if np.all(S) == None:
            S = np.matrix([[1, -1], [1, 1]])  # Swapped rows compared to other way
            S = np.kron(S, S)

        # Remove the constant row to obtain S
        S_hat = np.empty((0, self.num_bases))
        for row in S:
            if list(np.asarray(row)[0]) != [1 for i in range(self.num_bases)]:
                S_hat = np.append(S_hat, row, axis=0)

        L = self.__makeMat(F, S_hat, True)
        R = self.__makeMat(F, S_hat, False)
        F = np.asarray(F)  # , dtype=np.float64)
        if not returnLRMats:
            return L @ F @ R.T
        else:
            return (L @ F @ R.T, L, R)

    def transformedFlattening(self, F, S=None):
        """ Creates a transformed flattening from a flattening data frame """
        if np.all(S) == None:
            H = np.matrix([[1, -1], [1, 1]])  # Swapped rows compared to other way
            S = np.kron(H, H)

        colLabels = list(F)
        rowLabels = F.index
        F = np.matrix(F).astype(np.float64)
        L, R = S, S
        while L.shape[0] != F.shape[0]:
            L = np.kron(L, S)
        while R.shape[1] != F.shape[1]:
            R = np.kron(R, S)
        Ft = np.matmul(np.matmul(L, F), np.transpose(R))
        return pd.DataFrame(Ft, index=rowLabels, columns=colLabels, dtype=np.float64)

    def subFlattening(self, Ft, specialState='T', type=(1, 1)):
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
                if row != []: matrix.append(row)
        return np.asarray(matrix)

    def getParsimony(self, pattern):
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

        taxa = [n for n in nodes if self.isLeaf(n)]
        for i, t in enumerate(taxa):
            nodes[t]['pars'] = set(pattern[i])
        score = 0
        for n in nodes:
            if not self.isLeaf(n) and 'pars' not in nodes[n]:
                score += self.__parsimony(n, nodes=nodes)
        return score

    def __parsimony(self, n, nodes=None):
        """Recursive step in Finch's Algorithm"""
        score = 0
        children = self.getDescendants(n)
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

    def wronginess(self, split, verbose=False):

        split = split.split("|")
        S1 = set([int(i) for i in split[0]])
        S2 = set([int(i) for i in split[1]])

        partition = [S1, S2]
        matrices = []

        identityMatToIgnore = []

        for S in partition:
            T = cp.deepcopy(self)

            taxa = [n for n in T.nodes if T.isLeaf(n.index)]

            for s in S:
                for i in range(len(T.adjM[:, taxa[s].index])):
                    T.adjM[:, taxa[s].index] = 0
                if taxa[s] in T.nodes:
                    T.nodes.remove(taxa[s])
                    T.num_nodes -= 1

            new_taxa = [n for n in T.nodes if T.isLeaf(n.index)]
            new_taxa = [t for t in new_taxa if t not in taxa]

            while new_taxa != []:
                for t in new_taxa:
                    for i in range(len(T.adjM[:, t.index])):
                        T.adjM[:, t.index] = 0
                    if t in T.nodes:
                        T.nodes.remove(t)
                        T.num_nodes -= 1

                new_taxa = [n for n in T.nodes if T.isLeaf(n.index)]
                new_taxa = [t for t in new_taxa if t not in taxa]

            # After all new taxa are removed, we remove roots with only one descendants,
            # And in the case where the root has multiple descendants, we take note of the ones
            # which have the identity as their in-edge

            rootMightNeedRemoving = True  # placeholder
            while rootMightNeedRemoving:
                root = T.getRoot()
                rootDes = T.getDescendants(root)
                if len(rootDes) == 1:
                    for i in range(len(T.adjM[:, root.index])):
                        T.adjM[:, root.index] = 0
                        T.adjM[root.index, :] = 0
                    T.nodes = [n for n in T.nodes if n.index != root.index]
                    T.num_nodes -= 1
                else:
                    rootMightNeedRemoving = False

                root = T.getRoot()
                if len(T.getDescendants(root)) == 1:
                    rootMightNeedRemoving = True

            rootDes = T.getDescendants(T.getRoot())
            rootDesNodes = [n for n in T.nodes if (n.index in rootDes)]
            for d in rootDesNodes:
                if np.array_equal(d.inEdgeMat, np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])):
                    identityMatToIgnore.append(True)

            A = cp.copy(T.adjM)
            matrices.append(A)

        diffMat = matrices[0] * matrices[1]

        if len(identityMatToIgnore) == 2:
            sharedEdges = np.sum(diffMat) - 1
        else:
            sharedEdges = np.sum(diffMat)
        return sharedEdges
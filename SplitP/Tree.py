import numpy as np
import copy as cp
import pandas as pd
import itertools
import collections
from SplitP import TreeHelperFunctions as hf
from itertools import permutations

class tree:
    """A rooted phylogenetic tree.
    
    A rooted phylogenetic tree consisting of a collection of node objects and an adjacency matrix
    
    Attributes:
        num_nodes: The number of nodes in the tree
        num_bases: The size of the character set for the phylogeny
        nodes: A list of node objects representing the nodes in the tree
        adjM: An adjacency matrix describing the connections within the tree
        initDist: The initial distribution of states/characters at the root of the tree.
    """
    
    def __init__(self, n, b, name=""):
        """Initialises a tree with the number of nodes and size of the character set
        
        Args:
            n: The number of nodes which will be entered into the tree.
            b: The size of the character set for the tree.
        """
        self.num_nodes = n
        self.state_space = ('A', 'C', 'G', 'T')
        self.special_state = None
        self.num_bases = b
        self.nodes = [None for x in range(self.num_nodes)]
        self.adjM = np.array([[0 for x in range(self.num_nodes)] for x in range(self.num_nodes)])
        self.initDist = [None for x in range(self.num_bases)]
        self.name = name
        self.subflatLRMats = {}
        
    def __str__(self):
        return str(self.adjM)
        
    def addNode(self, n, in_node = None):
        """Adds a new node to the tree
        
        Args:
            n: The node object to add into the tree.
            in_node: The index of the parent node, default is None and is used for the root.
        """
        self.nodes[n.index] = n
        if in_node != None:
            self.adjM[in_node][n.index] = 1
    
    def getNumTaxa(self):
        """Returns the number of taxa/leaf-nodes in the tree"""
        o = 0
        for n in self.nodes:
            if self.isLeaf(n.index):
                o += 1
        return o
    
    def isLeaf(self, n_index):
        """Determines whether a node is a leaf node from it's index."""
        return np.all(self.adjM[n_index][:] == 0)
        
    def isRoot(self, n_index):
        """Determines whether a node is a root node from it's index."""
        return np.all(self.adjM[:,n_index] == 0) and not np.all(self.adjM[n_index,:] == 0) # THE AND PART IS UNTESTED
        
    def getRoot(self):
        """Returns the root node"""
        for n in self.nodes:
            if self.isRoot(n.index):
                return n

    def getParent(self, n):
        """Returns the parent node for a given node"""
        return self.nodes[np.nonzero(self.adjM[:,n.index])[0][0]]
    
    def getDescendants(self, n):
        """Returns a list of children/descendents of a given node"""
        return list(np.nonzero(self.adjM[n.index,:])[0])    
        
    def __likelihood(self, n, lTable):
        """Recursive part of the likelihood algorithm"""
        for b in range(self.num_bases):
            L = 1
            for d in self.getDescendants(n):
                M = (self.nodes[d]).inEdgeMat
                s = 0
                for b2 in range(self.num_bases):
                    if lTable[d, b2] == None:
                        self.__likelihood(self.nodes[d], lTable)
                    s += M[b2, b] * lTable[d, b2]
                L *= s
            lTable[n.index, b] = L
        if not self.isRoot(n.index):
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
            if p == 'A': return 0
            if p == 'C': return 1
            if p == 'G': return 2
            if p == 'T': return 3
            return int(p)
        pattern = [toInt(p) for p in pattern] # A list of indices which correspond to taxa.
        taxa = [n for n in self.nodes if self.isLeaf(n.index)] # The list of taxa.
        # Likelihood table for dynamic prog alg ~ lTable[node_index, character]
        lTable = np.array([[None for x in range(self.num_bases)] for x in range(self.num_nodes)])
        
        # Encoding pattern into likelihood table
        for i in range(len(pattern)):
            lTable[taxa[i].index, :] = [0 for x in range(self.num_bases)]
            lTable[taxa[i].index, pattern[i]] = 1
            
        # Starting with node which has all descendants as leaves
        initNode = None
        for n in self.nodes:
            desc = self.getDescendants(n)
            if desc and all(self.isLeaf(d) for d in desc):     # If the node has descendants and they are all leaves
                initNode = n

        self.__likelihood(initNode, lTable) # Should fill in the entire table
        
        root_index = self.getRoot().index
        L = 0
        for i in range(self.num_bases):
            L += self.initDist[i]*lTable[root_index, i]

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
        return(list(np.linalg.svd(np.array(M).astype(np.float64), compute_uv=False)))
    
    def splitScore(self, M, k = None, singularValues = False):
        """Returns the split score for a given flattening matrix"""
        sVals = list(np.linalg.svd(np.array(M).astype(np.float64), compute_uv=False))
        if singularValues: sing_vals = sVals.copy()
        sumSq = 0
        for i in range((self.num_bases if not k else k), min(M.shape)):
            sumSq += (sVals[i])**2
        sumSq2 = sumSq
        for i in range(1, (self.num_bases if not k else k)):
            sumSq2 += (sVals[i])**2
            
        score = np.sqrt(sumSq)/hf.fNorm(np.array(M).astype(np.float64))
        
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
                    for i in range(self.num_bases**self.getNumTaxa())
                ]
            )
    
        elif self.num_bases == 4:
            combinations = list(itertools.product('ACGT', repeat=self.getNumTaxa()))
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
        probs = list(LT.iloc[:,1])
        probs = [float(p) for p in probs]
        data = np.random.multinomial(n, probs)
        data = [d/n for d in data]
        if self.num_bases == 2:
            return np.array([['{:0{}b}'.format(i, self.getNumTaxa()),  data[i]] for i in range(self.num_bases**self.getNumTaxa())])
        elif self.num_bases == 4:
            patterns = list(itertools.product('ACGT', repeat=self.getNumTaxa()))
            patterns = [''.join(c) for c in patterns]
            return np.array([[patterns[i], data[i]] for i in range(self.num_bases**self.getNumTaxa())])

    def drawFromMultinomialFast(self, LT, n):
        """Use a given table of probabilities from getLikelihoods() and draw from its distribution"""
        probs = list(LT.iloc[:, 1])
        probs = [float(p) for p in probs]
        data = np.random.multinomial(n, probs)
        data = [d / n for d in data]
        if self.num_bases == 2:
            return np.array(
                [['{:0{}b}'.format(i, self.getNumTaxa()), data[i]] for i in range(self.num_bases ** self.getNumTaxa())])
        elif self.num_bases == 4:
            patterns = list(LT.iloc[:,0])
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
            rowLables = ['{:0{}b}'.format(i, len(split[0])) for i in range(2**len(split[0]))]
            colLables = ['{:0{}b}'.format(i, len(split[1])) for i in range(2**len(split[1]))]
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
            
        return F #.astype(pd.SparseDtype("float64", 0))
    
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
    
    
    def subFlatteningAlt(self, F, S = None, returnLRMats = False):
        def makeMat(S_hat, left = True):
            if left:
                M = np.empty((0, F.shape[0]))
                rows = len(rowLabels[0])
            else: 
                M = np.empty((0, F.shape[1]))
                rows = len(colLabels[0]) 
            O = np.ones(len(S_hat) + 1)
            for r in range(rows):
                A = np.empty((0,0))
                if r == 0: # First row of M
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
            return M
        # end makeMat()
        
        colLabels = list(F)
        rowLabels = F.index
        
        if np.all(S) == None:
            S = np.matrix([[1,-1],[1,1]]) # Swapped rows compared to other way
            S = np.kron(S,S)
        
        # Remove the constant row to obtain S
        S_hat = np.empty((0, 4))
        for row in S:
            if list(np.asarray(row)[0]) != [1,1,1,1]:
                S_hat = np.append(S_hat, row, axis=0)
        
        L = makeMat(S_hat, True)
        R = makeMat(S_hat, False)
        F = np.asarray(F, dtype=np.float64)
        if not returnLRMats: 
            return L @ F @ R.T
        else: 
            return (L @ F @ R.T, L, R)

    def subFlatteningAltFast(self, F, S=None, returnLRMats=False):

        def makeMat(S_hat, left=True):
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

        # end makeMat()

        colLabels = list(F)
        rowLabels = F.index

        if np.all(S) == None:
            S = np.matrix([[1, -1], [1, 1]])  # Swapped rows compared to other way
            S = np.kron(S, S)

        # Remove the constant row to obtain S
        S_hat = np.empty((0, self.num_bases))
        for row in S:
            if list(np.asarray(row)[0]) != [1 for i in range(self.num_bases)]:
                S_hat = np.append(S_hat, row, axis=0)

        L = makeMat(S_hat, True)
        R = makeMat(S_hat, False)
        F = np.asarray(F) #, dtype=np.float64)
        if not returnLRMats:
            return L @ F @ R.T
        else:
            return (L @ F @ R.T, L, R)
    
    def transformedFlattening(self, F, S = None):
        """ Creates a transformed flattening from a flattening data frame """
        if np.all(S) == None:
            H = np.matrix([[1,-1],[1,1]]) # Swapped rows compared to other way
            S = np.kron(H,H)
        
        colLabels = list(F)
        rowLabels = F.index
        F = np.matrix(F).astype(np.float64)
        L, R = S, S
        while L.shape[0] != F.shape[0]:
            L = np.kron(L,S)
        while R.shape[1] != F.shape[1]:
            R = np.kron(R,S)
        Ft = np.matmul(np.matmul(L,F), np.transpose(R))
        return pd.DataFrame(Ft, index=rowLabels, columns=colLabels, dtype=np.float64)

    def subFlattening(self, Ft, specialState='T', type=(1,1)):
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
                    if len("".join([x for x in r if x != specialState])) <= type[0] and len("".join([x for x in c if x != specialState])) <= type[1]:
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
        if '|' in pattern:
            pattern2 = [0 for x in range(len(pattern)-1)]
            i = 0
            while pattern[i] != '|':
                pattern2[int(pattern[i])] = 1
                i += 1
            pattern = "".join(str(i) for i in pattern2)
            
        taxa = [n for n in self.nodes if self.isLeaf(n.index)]
        for i in range(len(taxa)):
            taxa[i].pars = set(pattern[i])
        score = 0
        for n in self.nodes:
            if not self.isLeaf(n.index) and n.pars == None:
                score += self.__parsimony(n)
        for n in self.nodes:
            n.pars = None
        return score

    def __parsimony(self, n):
        """Recursive step in Finch's Algorithm"""
        score = 0
        children = self.getDescendants(n)
        children = [self.nodes[c] for c in children]
        for c in children:
            if c.pars == None:
                score += self.__parsimony(c)
                
        n.pars = children[0].pars
        for c in range(1, len(children)):
            n.pars = (n.pars).intersection(children[c].pars)
        if not n.pars:
            n.pars = children[0].pars
            for c in range(1, len(children)):
                n.pars = (n.pars).union(children[c].pars)
            return score + 1
        return score 
    
    
    def wronginess(self, split, verbose = False):
        
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
            
            rootMightNeedRemoving = True # placeholder
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
                if np.array_equal(d.inEdgeMat, np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])):
                    identityMatToIgnore.append(True)
            
            A = cp.copy(T.adjM)
            matrices.append(A)
    
        diffMat = matrices[0] * matrices[1]
    
        if len(identityMatToIgnore) == 2:
            sharedEdges = np.sum(diffMat) - 1
        else:
            sharedEdges = np.sum(diffMat)
        return sharedEdges
        
class node:
    """Represents a node, used by the tree class.
        
    Attributes:
        index: The number of nodes in the tree
        pars: The size of the character set for the phylogeny
        inEdgeMat: A list of node objects representing the nodes in the tree
    """
    
    def __init__(self, mat = None, index = None):
        self.index = index
        self.pars = None
        self.inEdgeMat = np.array(mat)
        
    def __str__(self):
         return "n_" + str(self.index)
     
    def checkMatrix(self, M):
        # Check that the transition matrix makes sense.
        pass
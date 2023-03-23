from splitp.parsers.newick import newick_to_json, json_to_newick
from networkx import json_graph
from dataclasses import dataclass
from networkx import draw
from splitp.enums import DrawFormat
import splitp.constants as constants
import numpy as np

class Phylogeny:
    __slots__ = (
        "name",
        "networkx_graph",
        "newick_string",
    )

    def __init__(
        self,
        newick_string,
        name=None,
    ):
        """A rooted phylogenetic tree.

        A rooted phylogenetic tree, stored as a wrapped networkx graph.

        Attributes:
            name: A name for the tree
            networkx_graph: the underlying networkx graph
        """

        # Set chosen tree properties
        self.name = name if name else "Unnamed Tree"
        self.newick_string = newick_string

        # Build networkx graph from newick string
        self.networkx_graph = json_graph.tree_graph(
            newick_to_json(self.newick_string, generate_names=True)
        )

    def __str__(self):
        """Return the tree in JSON format"""
        return json_to_newick(
            json_graph.tree_data(self.networkx_graph, self.root(return_index=False))
        )
    
    def reassign_transition_matrices(self, transition_matrix):
        """DEPRECATED: Reassign transition matrices to all nodes in the tree"""
        for node in self.networkx_graph.nodes:
            self.networkx_graph.nodes[node]["transition_matrix"] = np.array(transition_matrix)
    
    def get_num_nodes(self):
        return len(self.networkx_graph.nodes)
    
    def node_index(self, node):
        return list(self.networkx_graph.nodes).index(node)

    def __eq__(self, other):
        """Check equality of JSON objects"""
        return self.__str__() == other.__str__()

    def is_root(self, node):
        """Determines whether a node is a root node from its index."""
        if type(node) == type(str()):
            # A node name was supplied
            return self.networkx_graph.in_degree(node) == 0
        else:
            # An index was supplied
            return (
                self.networkx_graph.in_degree(list(self.networkx_graph.nodes)[node])
                == 0
            )

    def root(self, return_index=True):
        """Returns the root node"""
        for i, n in enumerate(self.networkx_graph.nodes):
            if self.is_root(n):
                return i if return_index else n

    def get_parent(self, n):
        """Returns the parent node for a given node"""
        return list(self.nx_graph.predecessors(n))[0]
    
    def nodes(self):
        """Returns a list of all nodes in the tree."""
        return list(self.networkx_graph.nodes)

    def get_taxa(self):
        return [n for n in self.networkx_graph.nodes if self.is_leaf(n)]
    
    def get_num_taxa(self):
        return len(self.get_taxa())

    def is_leaf(self, n_index_or_name):
        """Determines whether a node is a leaf node from it's index."""
        if type(n_index_or_name) == type(str()):
            return self.networkx_graph.out_degree(n_index_or_name) == 0
        else:
            return (
                self.networkx_graph.out_degree(
                    list(self.networkx_graph.nodes)[n_index_or_name]
                )
                == 0
            )

    def get_parent(self, n):
        """Returns the parent node for a given node"""
        return list(self.networkx_graph.predecessors(n))[0]

    def get_descendants(self, n, return_iter=False):
        """Returns a list of children/descendents of a given node"""
        return (
            list(self.networkx_graph.successors(n))
            if not return_iter
            else self.networkx_graph.successors(n)
        )

    def splits(self, include_trivial=False, as_strings=False):
        """Returns set of all splits displayed by the tree."""
        from networkx.algorithms.traversal.depth_first_search import dfs_tree
        all_taxa = [x for x in self.get_taxa()]
        splits = set()
        for node in list(self.nodes()):
            subtree = dfs_tree(self.networkx_graph, node)
            left = tuple(
                    sorted(
                        [
                            leaf
                            for leaf in subtree.nodes
                            if leaf in all_taxa
                        ],
                        key=all_taxa.index,
                    )
                )
            right = tuple(sorted((i for i in all_taxa if i not in left), key=all_taxa.index))
            split = (left, right)
            if all_taxa[0] not in split[0]:
                split = (split[1], split[0])
            if include_trivial or (len(split[0]) > 1 and len(split[1]) > 1):
                splits.add(split)
        for split in splits:
            yield split if not as_strings else self.format_split(split)

    def format_split(self, split):
        if isinstance(split, str):
            return split  # If already a string, just send it back
        if len(split[0]) + len(split[1]) > 35:
            raise ValueError(
                "Cannot produce string format for split with more than 35 taxa."
            )
        if all(len(taxon) == 1 for taxon in self.get_taxa()):
            return f'{"".join(split[0])}|{"".join(split[1])}'

    def draw(self, draw_format=DrawFormat.ASCII):
        try:
            from Bio.Phylo import draw_ascii, draw, read
            from io import StringIO
        except ModuleNotFoundError:
            print(
                "Additional dependencies required to draw tree. Please install BioPython and StringIO."
            )
            return
        tree = read(StringIO(str(self)), "newick")
        if draw_format == DrawFormat.ASCII:
            draw_ascii(tree)
        elif draw_format == DrawFormat.matplotlib:
            draw(tree, branch_labels=lambda c: c.branch_length)
        else:
            raise ValueError("Invalid format.")

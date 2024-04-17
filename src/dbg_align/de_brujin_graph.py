from functools import singledispatchmethod
from typing import List, Union

import cogent3
from cogent3.core.moltype import MolType

from enum import Enum

class AlignmentMethod(Enum):
    RAW = "RAW"
    DBG1D = "DBG1D"
    DBG2D = "DBG2D"

class DeBrujinGraph:
    """ A class to represent a de Bruijn graph for a set of sequences.

    Note: Indexes for sequences are 1 based (ie: start from 1).
    """
    def __init__(self, kmer_length: int, moltype: MolType = cogent3.DNA):
        self.kmer_length = kmer_length
        self.root = DeBrujinGraphNode(self, None)  # Root node of the graph
        self.graph = {}
        self.moltype = moltype
        self.sequence_names = []  # Index table to store names of sequences
        self.is_compressed = False

    @classmethod
    def generate_kmers(cls, sequence: str, k: int):
        for i in range(len(sequence) - k + 1):
            yield sequence[i:i + k]

    def compress_graph(self):
        """Compresses the graph by extending nodes until the next branch point or the end of a run."""
        visited = set()
        for node in list(self.graph.values()):
            if node not in visited:
                run_nodes = node.run()
                if len(run_nodes) > 1:
                    node.compress(run_nodes)  # Compress the entire run starting from 'node'
                    visited.update(run_nodes)  # Mark all nodes in the run as visited
        self.graph = {} 
        self.is_compressed = True

    def order_complexity(self, alignment_type: AlignmentMethod):
        """Returns the order complexity of aligninging the sequences."""
        if alignment_type == AlignmentMethod.RAW:
            product = 1
            for sequence in self:
                product *= len(sequence)
            return product
        elif alignment_type == AlignmentMethod.DBG1D or alignment_type == AlignmentMethod.DBG2D:
            return 0  #TODO: Implement the order complexity calculation for DBG1D and DBG2D
        else:
            raise ValueError("Unsupported alignment type")
        
    @singledispatchmethod
    def add_sequence(self, sequence, name=None):
        # Placeholder for unrecognized types
        raise TypeError("Unsupported sequence type")

    @add_sequence.register(str)
    def _(self, sequence: str, name=None):
        # check string for characcters in alphabet
        self.moltype.verify_sequence(sequence)
        if len(sequence) < self.kmer_length:
            raise ValueError("Sequence is shorter than kmer length")
        sequence_index = len(self.sequence_names)+1
        if not name:
            name = f"Sequence_{sequence_index}"
        self.sequence_names.append(name)  
        # Convert the sequence into kmers and add them to the graph
        current_node = self.root
        passage_index = 0
        for kmer in self.generate_kmers(sequence, self.kmer_length):
            next_node = self.graph.get(kmer)
            if next_node is None:
                next_node = DeBrujinGraphNode(self, kmer)
                self.graph[kmer] = next_node
            current_node.add_edge(target_node= next_node, 
                                  sequence_index=sequence_index,
                                  passage_index=passage_index)
            current_node = next_node
            passage_index += 1

    @add_sequence.register(cogent3.Sequence)
    def _(self, sequence: cogent3.Sequence, name=None):
        if sequence.moltype != self.moltype:
            raise ValueError("Sequence moltype does not match dBg moltype")
        name = name or sequence.name
        if not name:
            name = f"Sequence_{len(self.sequence_names)+1}"
            self.sequence_names.append(name)  # Store the name in the index table
        # Logic to add sequence goes here

    @add_sequence.register(cogent3.SequenceCollection)
    def _(self, sequences: cogent3.SequenceCollection):
        for sequence in sequences:
            self.add_sequence(sequence, sequence.name or None)

    @add_sequence.register(list)
    def _(self, sequences, names = None):
        if names and len(names) != len(sequences):
            raise ValueError("Names and sequences must have the same length")
        if names is None:
            names = [f"Sequence_{i+1}" for i in range(len(sequences))]
        for sequence, name in zip(sequences, names):
            self.add_sequence(sequence, name)

    @add_sequence.register
    def _(self, sequences: dict):
        for name, sequence in sequences.items():
            self.add_sequence(sequence, name)

    def names(self):
        """Returns an iterable collection of sequence names."""
        return self.sequence_names

    def has_cycles(self):
        """Returns True if the graph contains cycles."""
        for node in self.graph.values():
            if node.has_cycle():
                return True
        return False
    
    def has_bubbles(self)-> bool:
        """Returns True if the graph contains bubbles."""
        for node in self.graph.values():
            if node.has_bubble():
                return True
        return False
    
    def __len__(self):
        return len(self.sequence_names)
    
    @singledispatchmethod
    def __getitem__(self, index: Union[int, str]):
        raise TypeError("Index must be a string or an integer")

    @__getitem__.register
    def _(self, index: int):
        if index < 1 or index > len(self.sequence_names):
            raise IndexError("Sequence index out of range")
        # Start the sequence reconstruction from the root node
        sequence = self.root.get_sequence(index)
        return sequence

    @__getitem__.register
    def _(self, name: str):
        if name not in self.sequence_names:
            raise KeyError(f"Sequence name '{name}' not found")
        sequence_index = self.sequence_names.index(name) + 1
        # Start the sequence reconstruction from the root node
        sequence = self.root.get_sequence(sequence_index)
        return sequence

    def __repr__(self):
        return f"dbg k:{self.kmer_length}, mol:{self.moltype}, seq's:{len(self)})"
    
    def to_mermaid(self):
        """Generates a Mermaid graph description of the de Bruijn graph."""
        if not self.root:
            return "graph LR;"

        mermaid_str = "graph LR;\n"
        mermaid_str = self._mermaid_from_node(self.root, mermaid_str)
        return mermaid_str

    def _mermaid_from_node(self, root, mermaid_str):
        stack = [(root, 0)]  # Stack to keep nodes and their edge index
        visited = set()

        while stack:
            node, edge_index = stack[-1]  # Get the top item from the stack
            if node in visited:
                stack.pop()
                continue

            if edge_index < len(node.edges):
                edge = node.edges[edge_index]
                target = edge.target_node
                if node.kmer is None:
                    mermaid_str += f"{'Root'} -->|\"{len(edge.traversals)}\"| {target.kmer};\n"
                else:
                    mermaid_str += f"{node.kmer} -->|\"{len(edge.traversals)}\"| {target.kmer};\n"
                # Increment the current node's edge index for the next iteration
                stack[-1] = (node, edge_index + 1)
                # Push the target node onto the stack if it hasn't been visited
                if target not in visited:
                    stack.append((target, 0))
            else:
                # All edges for the current node have been processed
                visited.add(node)
                stack.pop()

        return mermaid_str
    
class DeBrujinGraphNode:
    def __init__(self, de_brujin_graph : DeBrujinGraph, kmer: str) -> None:
        self.kmer = kmer
        self.edges = []  # List of DeBrujinGraph_Edge objects
        self.parent_nodes = set()  # set of nodes that lead to this node
        
    @property 
    def parents(self):
        return self.parent_nodes
    
    @property
    def children(self):
        return {edge.target_node for edge in self.edges}

    @property
    def has_single_parent(self) -> bool:
        return len(self.parent_nodes) == 1
    
    @property
    def has_one_child(self) -> bool:
        return len(self.edges) == 1
    
    @property
    def first_child(self) -> "DeBrujinGraphNode":
        return self.edges[0].target_node
    
    @property
    def sole_child(self) -> "DeBrujinGraphNode":
        if not self.has_one_child:
            raise ValueError("Node does not have a single child")
        return self.edges[0].target_node

    def add_edge(self, target_node, sequence_index=None, passage_index=None):
        """Create a new edge to target_node if not already existing, or return existing one."""
        for edge in self.edges:
            if edge.target_node == target_node:
                edge.add_traversal(sequence_index, passage_index)
                return edge
        new_edge = DeBrujinGraphEdge(target_node)
        new_edge.add_traversal(sequence_index, passage_index)
        self.edges.append(new_edge)
        target_node.parent_nodes.add(self) # Add this node as a parent to the target node
        return new_edge
 
    def run(self):
        """Generates a list of nodes that form a linear run from this node."""
        run_nodes = [self]
        candidate = self
        while candidate.has_one_child \
                and candidate.sole_child.has_single_parent:
            child = candidate.first_child
            run_nodes.append(child)
            candidate = child
        return run_nodes
        
    def compress(self, run_nodes = None):
        if not run_nodes:
            run_nodes = self.run()
        if len(run_nodes) > 1:
            # Extend the base node's kmer and adopt the last node's edges
            last_node = self
            for node in run_nodes[1:]:
                self.kmer += node.kmer[-1]  # Assuming k-1 overlap
                last_node = node
            for child in last_node.children:
                if not last_node in child.parent_nodes:
                    pass
                child.parent_nodes.remove(last_node)
                child.parent_nodes.add(self)
            self.edges = last_node.edges

    def has_cycle(self, visited=None, ancestors=None):
        if visited is None:
            visited = set()
        if ancestors is None:
            ancestors = set()

        if self in ancestors:
            return True  # A cycle is detected.
        if self in visited:
            return False

        visited.add(self)
        ancestors.add(self)
        for edge in self.edges:
            if edge.target_node.has_cycle(visited, ancestors):
                return True
        ancestors.remove(self)
        return False
   
    def sequence_edges(self, sequence_index: int) -> List["DeBrujinGraphEdge"]:
        """Retrieve all edges that correspond to a specific sequence index."""
        if sequence_index is None:
            return self.edges
        filtered_edges = [edge for edge in self.edges 
                          if any(traversal.sequence_index == sequence_index 
                                 for traversal in edge.traversals)]
        return filtered_edges

    def sequence_passage_edges(self, sequence_index: int, passage_index: int) -> List["DeBrujinGraphEdge"]:
        """Retrieve all edges that correspond to a specific sequence and passage index."""
        filtered_edges = [edge for edge in self.edges 
                        if any(traversal.sequence_index == sequence_index 
                               and traversal.passage_index == passage_index 
                               for traversal in edge.traversals)]
        return filtered_edges

    def find(self, sequence_index: int, passage_index: int, visited=None):
        """return the node that contains an edge with the specified sequence and passage index."""
        if visited is None:
            visited = set()

        visited.add(self)  # mark the current node as visited

        for edge in self.sequence_edges(sequence_index):
            for traversal in edge.traversals:
                if traversal.passage_index == passage_index and traversal.sequence_index == sequence_index:
                    return self # found the target
            if edge.target_node not in visited:
                found_node = edge.target_node.find(sequence_index, passage_index, visited)
                if found_node:
                    return found_node  # return only if a node was found
        return None  # return None if no matching node was found

    def get_sequence(self, sequence_index: int, start_passage_index: int = 1, length: int = None) -> str:
        current_node = self.find(sequence_index, start_passage_index)
        if not current_node:
            raise ValueError("No starting node found for the specified sequence and passage index")
        current_passage_index = start_passage_index 
        sequence = current_node.kmer

        while current_node:
            # Retrieve the edge that corresponds to the current sequence and passage index
            edges = current_node.sequence_passage_edges(sequence_index, current_passage_index)
            if not edges:
                break  # No more edges match the criteria, stop the sequence construction

            # Assuming there's exactly one edge for each passage index (handle errors or warnings otherwise)
            if len(edges) > 1:
                # This is a simple error handling case where multiple edges for the same passage index are unexpected
                raise ValueError("Multiple edges found for the same passage index, expected only one.")
            
            next_edge = edges[0]
            
            # add the last character of the kmer
            sequence += next_edge.target_node.kmer[-1]
            
            # Update the character count and check if we've reached the desired length
            if length is not None and len(sequence) >= length:
                sequence = sequence[:length]  # Trim the sequence to the specified length if necessary
                break
            
            # Move to the next node and increment the passage index
            current_node = next_edge.target_node
            current_passage_index += 1

        return sequence

    def __getitem__(self, index):
        return self.edges[index].target_node
    
    def __len__(self):
        return len(self.edges)
    
    def __repr__(self):
        return f"Node:({self.kmer}) [{','.join([edge.target_node.kmer for edge in self.edges])}]"

class DeBrujinGraphTraversal:
    def __init__(self, sequence_index: int, passage_index: int) -> None:
        self.sequence_index = sequence_index
        self.passage_index = passage_index

    def __repr__(self):
        return f"({self.sequence_index}:{self.passage_index})"

class DeBrujinGraphEdge:
    def __init__(self, target_node: DeBrujinGraphNode) -> None:
        self.target_node = target_node
        self.traversals : List[DeBrujinGraphTraversal] = []  

    def add_traversal(self, sequence_index:int, passage_index:int):
        self.traversals.append(DeBrujinGraphTraversal(sequence_index, passage_index))
        return self

    def __repr__(self):
        return f"Edge to {self.target_node.kmer} ({len(self.traversals)} traversals)"
    
    def label(self):
        return f"{[(t.sequence_index, t.passage_index) for t in self.traversals]}"
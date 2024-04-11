from functools import singledispatchmethod
from typing import List, Union

import cogent3
from cogent3.core.moltype import MolType


class DeBrujinGraph:
    """ A class to represent a de Bruijn graph for a set of sequences.

    Note: Indexes for sequences are 1 based (ie: start from 1).
    """
    def __init__(self, kmer_length: int, moltype: MolType = cogent3.DNA):
        self.kmer_length = kmer_length
        self.root = DeBrujinGraph_Node(self, None)  # Root node of the graph
        self.graph = {}
        self.moltype = moltype
        self.sequence_names = []  # Index table to store names of sequences

    @classmethod
    def generate_kmers(cls, sequence: str, k: int):
        for i in range(len(sequence) - k + 1):
            yield sequence[i:i + k]

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
            passage_index += 1
            next_node = self.graph.get(kmer)
            if next_node is None:
                next_node = DeBrujinGraph_Node(self, kmer)
                self.graph[kmer] = next_node
            current_node.add_edge(target_node= next_node, 
                                  sequence_index=sequence_index,
                                  passage_index=passage_index)
            current_node = next_node

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

    def compress(self):
        from .partial_order_graph import PartialOrderGraph

        """Compresses the graph to a partial order graph"""
        return PartialOrderGraph(self)

    def has_cycles(self):
        """Returns True if the graph contains cycles."""
        for node in self.graph.values():
            if node.has_cycle():
                return True
        return False
    
    def has_bubbles(self):
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
    

class DeBrujinGraph_Node:
    def __init__(self, de_brujin_graph : DeBrujinGraph, kmer: str) -> None:
        self.kmer = kmer
        self.edges = []  # List of DeBrujinGraph_Edge objects
        self.de_brujin_graph = de_brujin_graph

    def add_edge(self, target_node, sequence_index, passage_index):
        new_edge = DeBrujinGraph_Edge(target_node, sequence_index, passage_index)
        self.edges.append(new_edge)
        return new_edge
 
    def has_cycle(self)-> bool:
        """Returns True if there are multiple edges for any sequence_index."""
        sequence_index_count = {}
        for edge in self.edges:
            if edge.sequence_index in sequence_index_count:
                return True  # Found a cycle
            sequence_index_count[edge.sequence_index] = 1
        return False

    def has_bubble(self) -> bool:
        """Returns True if there are at least two edges leading to different target nodes."""
        target_nodes = set()
        for edge in self.edges:
            target_nodes.add(edge.target_node)
            if len(target_nodes) > 1:
                return True
        return False

    def sequence_edges(self, sequence_index: int) -> List["DeBrujinGraph_Edge"]:
        """Retrieve all edges that correspond to a specific sequence index."""
        filtered_edges = [edge for edge in self.edges if edge.sequence_index == sequence_index]
        return filtered_edges

    def sequence_passage_edges(self, sequence_index: int, passage_index: int) -> List["DeBrujinGraph_Edge"]:
        """Retrieve all edges that correspond to a specific sequence and passage index."""
        filtered_edges = [edge for edge in self.edges 
                          if edge.sequence_index == sequence_index 
                          and edge.passage_count == passage_index]
        return filtered_edges
    
    def find(self, sequence_index: int, passage_index: int, visited=None):
        if visited is None:
            visited = set()

        visited.add(self)  # mark the current node as visited

        for edge in self.sequence_edges(sequence_index):
            if edge.passage_count == passage_index:
                return edge.target_node  # found the target

            if edge.target_node not in visited:
                found_node = edge.target_node.find(sequence_index, passage_index, visited)
                if found_node:
                    return found_node  # return only if a node was found
        return None  # return None if no matching node was found

    def get_sequence(self, sequence_index: int, start_passage_index: int = 1, length: int = None) -> str:
        sequence = ''
        current_node = self
        current_passage_index = start_passage_index
        char_count = 0  # Track the number of characters added to the sequence

        current_node = self.find(sequence_index, start_passage_index)
        if not current_node:
            raise ValueError("No starting node found for the specified sequence and passage index")

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
            
            # Decide whether to add the entire kmer or just the last character
            if not sequence:
                # If sequence is empty, start with the full kmer
                sequence += next_edge.target_node.kmer
            else:
                # Otherwise, just add the last character of the kmer
                sequence += next_edge.target_node.kmer[-1]
            
            # Update the character count and check if we've reached the desired length
            char_count += len(sequence[-len(next_edge.target_node.kmer):])
            if length is not None and char_count >= length:
                sequence = sequence[:length]  # Trim the sequence to the specified length if necessary
                break
            
            # Move to the next node and increment the passage index
            current_node = next_edge.target_node
            current_passage_index += 1

        return sequence

    def __repr__(self):
        if self.edges:
            return f"Node {self.kmer} edges:{len(self.edges)} sequence_ids:{[e.sequence_index for e in self.edges]}"
        else:
            return f"Node {self.kmer} edges:0"


class DeBrujinGraph_Edge:
    def __init__(self, target_node: DeBrujinGraph_Node, sequence_index: int, passage_count:int) -> None:
        self.sequence_index = sequence_index
        self.target_node = target_node
        self.passage_count = passage_count

    def __repr__(self):
        return f"Edge to {self.target_node.kmer} ({self.sequence_index}, {self.passage_count})"
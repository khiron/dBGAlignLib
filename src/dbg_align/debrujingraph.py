from collections import deque
from functools import singledispatchmethod
from typing import List, Union

import cogent3
from cogent3.core.moltype import MolType
from .dbg_edge import DBGNode
from .dbg_traversal import DBGTraversal
from enum import Enum
from graphviz import Digraph

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
        self.root = DBGNode(kmer = None, kmer_length = kmer_length)  # Root node of the graph
        self.graph = {}
        self.moltype = moltype
        self.sequence_names = {}  # dict keyed on sequence names, returns tuple containing index and lengths of the sequence
        self.is_compressed = False

    @classmethod
    def generate_kmers(cls, sequence: str, k: int):
        for i in range(len(sequence) - k + 1):
            yield sequence[i:i + k]

    def to_pog(self)->"DeBrujinGraph":
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
        return self

    def order_complexity(self, alignment_type: AlignmentMethod):
        """Returns the order complexity of aligninging the sequences."""
        if alignment_type == AlignmentMethod.RAW:
            product = 1
            for _, (_, length) in self.sequence_names.items():
                product *= length
            return product
        elif alignment_type == AlignmentMethod.DBG1D:
            if not self.is_compressed:
                raise ValueError("Graph must be compressed before calculating de Brujin graph order complexity")
                # for each level 1 bubble 
                # sum the product of lengths of each sequence in the bubble
            return 0 #TODO: Implement the order complexity calculation for DBG1D
        elif alignment_type == AlignmentMethod.DBG2D:
            if not self.is_compressed:
                raise ValueError("Graph must be compressed before calculating de Brujin graph order complexity")
            # for each bubble
            # sum the product of lengths of each sequence in the bubble
            return 0  #TODO: Implement the order complexity calculation for DBG2D
        else:
            raise ValueError("Unsupported alignment type")

    def compression_ratio(self) -> int:
        """Returns the mean kmer size of the POG over the original kmer length of the DBG that created the POG."""
        if not self.is_compressed:
            raise ValueError("Graph must be compressed before calculating compression ratio")
        total_kmer_length = 0
        node_count = 0

        for node in self.root.traverse_all():
            if node.kmer:
                total_kmer_length += len(node.kmer)
                node_count += 1

        if node_count > 0:
            return total_kmer_length / (node_count * self.kmer_length)
        else:
            return 0

    def kmers(self, kmer: str)-> List["DBGNode"]:
        """Returns the node for a given kmer."""
        result = []
        for node in self.root.traverse_all():
            if node.kmer == kmer:
                result.append(node)
        return result

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
        sequence_index = len(self)+1
        if not name:
            name = f"Sequence_{sequence_index}"
        self.sequence_names[name] = (sequence_index, len(sequence))
        # Convert the sequence into kmers and add them to the graph
        current_node = self.root
        passage_index = 0
        for kmer in self.generate_kmers(sequence, self.kmer_length):
            next_node = self.graph.get(kmer)
            if next_node is None:
                next_node = DBGNode(kmer, self.kmer_length)
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
            name = f"Sequence_{len(self)+1}"
        self.add_sequence(str(sequence), name)    

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
        return list(self.sequence_names.keys())

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
    
    def sequence_length(self, sequence_index: int) -> int:
        """Returns the length of a sequence in the graph."""
        for _, (index, length) in self.sequence_names.items():
            if index == sequence_index:
                return length
        raise KeyError(f"Sequence index '{sequence_index}' not found")
    
    def __len__(self):
        return len(self.sequence_names.keys())

    def __iter__(self):
        for index in range(1, len(self) + 1):
            yield self[index]
    
    @singledispatchmethod
    def __getitem__(self, index: Union[int, str]):
        raise TypeError("Index must be a string or an integer")

    @__getitem__.register
    def _(self, index: int):
        if index < 1 or index > len(self):
            raise IndexError("Sequence index out of range")
        # Start the sequence reconstruction from the root node
        sequence = self.root.get_sequence(index)
        return sequence

    def index_for_name(self, name: str)->int:
        """Returns the index for a sequence name."""
        seq = self.sequence_names[name]
        if not seq:
            raise KeyError(f"Sequence name '{name}' not found")
        return seq[0]
    
    def len_for_name(self, name: str)->int:
        """Returns the length for a sequence name."""
        seq = self.sequence_names[name]
        if not seq:
            raise KeyError(f"Sequence name '{name}' not found")
        return seq[1]

    @__getitem__.register
    def _(self, name: str):
        if name not in self.sequence_names:
            raise KeyError(f"Sequence name '{name}' not found")
        sequence_index = self.index_for_name(name)
        # Start the sequence reconstruction from the root node
        sequence = self.root.get_sequence(sequence_index)
        return sequence

    def __repr__(self):
        return f"dbg k:{self.kmer_length}, mol:{self.moltype}, seq's:{len(self)})"
    
    def to_mermaid(self, show_kmers: bool = True):
        """Generates a Mermaid graph description of the de Bruijn graph."""
        if not self.root:
            return "graph LR;"

        mermaid_str = "graph LR;\n"
        queue = deque([self.root])
        visited = set()
        
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            for edge in node.edges:
                target = edge.target_node
                if target not in visited:
                    queue.append(target)
                if not show_kmers:
                    mermaid_str += f'{node.kmer}(" ") --> {target.kmer}(" ");\n'
                else:
                    if node.kmer:
                        mermaid_str += f"{node.kmer} --> {target.kmer};\n"
                    else:
                        mermaid_str += f"{node.kmer} --> {target.kmer};\n"
        
        return mermaid_str

    def to_graphviz(self, show_kmers: bool = True):
        def sanitize_identifier(identifier):
            # Replace spaces and special characters with underscores
            return "".join(char if char.isalnum() else "_" for char in identifier)

        """Generates a Graphviz graph description of the de Bruijn graph."""
        if not self.root:
            return None

        dot = Digraph(comment='De Bruijn Graph')
        dot.attr('graph', rankdir='LR')  # Lay out the graph from left to right

        queue = deque([self.root])
        visited = set()
        
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            # Node label handling
            # Handle special case for the root node or nodes with None kmer
            if node.kmer is None:
                node_id = "Root"
                node_label = "Root"  # Always label the root node as "Root"
            else:
                node_id = sanitize_identifier(node.kmer)
                node_label = node.kmer if show_kmers else " "

            dot.node(node_id, label=node_label)

            for edge in node.edges:
                target = edge.target_node
                if target not in visited:
                    queue.append(target)
                    target_id = "Root" if target.kmer is None else sanitize_identifier(target.kmer)
                # Edge label handling, could be more sophisticated depending on your needs
                edge_label = str(len(edge.traversals)) if show_kmers else " "
                dot.edge(node_id, target_id)

        return dot
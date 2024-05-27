from collections import deque
from functools import singledispatchmethod
from typing import Union

import cogent3
from cogent3.core.moltype import MolType
from .constants import AlignmentMethod
from graphviz import Digraph

class DeBruijnGraph:
    """ A class to represent a de Bruijn graph for a set of sequences.

    Note: Indexes for sequences are 1 based (ie: start from 1).
    """
    def __init__(self, kmer_length: int, moltype: MolType = cogent3.DNA):
        from .dbg_node import DBGNode
        self.kmer_length = kmer_length
        self.root = DBGNode(kmer = None)  # Root node of the graph
        self.graph = {}
        self.moltype = moltype
        self.sequence_names = {}  # dict keyed on sequence names, returns tuple containing index and lengths of the sequence
        self.is_compressed = False

    @classmethod
    def generate_kmers(cls, sequence: str, k: int):
        for i in range(len(sequence) - k + 1):
            yield sequence[i:i + k]

    def to_pog(self)->"DeBruijnGraph":
        from .partialordergraph import PartialOrderGraph
        pog = PartialOrderGraph(self)
        return pog

    def expected_work(self, alignment_type: AlignmentMethod):
        """Returns the order complexity of aligninging the sequences."""
        if alignment_type == AlignmentMethod.EXACT:
            product = 1
            for _, (_, length) in self.sequence_names.items():
                product *= length
            return product
        elif alignment_type == AlignmentMethod.PROGRESSIVE:
            # Extract sequence lengths and sort them
            sequence_lengths = sorted(length for _, (_, length) in self.sequence_names.items())
            # Sum the product of each pair of consecutive lengths
            return sum(sequence_lengths[i] * sequence_lengths[i+1] for i in range(len(sequence_lengths) - 1))
        elif alignment_type in (AlignmentMethod.DBG_LENGTH, AlignmentMethod.DBG_LENGTH_NUMBER):
            raise ValueError("Graph must be transformed to PartialOrderGraph before calculating expected_work")
        else:
            raise ValueError("Unsupported alignment type")

    @singledispatchmethod
    def add_sequence(self, sequence, name=None):
        # Placeholder for unrecognized types
        raise TypeError("Unsupported sequence type")

    @add_sequence.register(str)
    def _(self, sequence: str, name=None):
        from .dbg_edge import DBGEdge, DBGNode
        
        # check string for characters in alphabet
        self.moltype.verify_sequence(sequence)
        if len(sequence) < self.kmer_length:
            raise ValueError("Sequence is shorter than kmer length")
        sequence_index = len(self)+1
        if not name:
            name = f"Sequence_{sequence_index}"
        self.sequence_names[name] = (sequence_index, len(sequence))
        # Convert the sequence into kmers and add them to the graph
        current_node = self.root
        for kmer in self.generate_kmers(sequence, self.kmer_length):
            next_node = self.graph.get(kmer)
            if not next_node: # Node doens't exist it, add it and make an edge for this sequence to it or connect a cycle edge to it
                next_node = DBGNode(kmer)
                self.graph[kmer] = next_node
                cycle_edge = current_node.get_cycle_edge(sequence_index)
                if cycle_edge: # it's a cycle we can close
                    cycle_edge.target_node = next_node
                else:    
                    current_node.edges.append(DBGEdge(target_node=next_node, sequence_index=sequence_index)) 
                current_node = next_node
            else: # Node already exists, check if we have an edge for this sequence
                if next_node.get_cycle_edge(sequence_index): # it's a cycle if the next node has an edge for this sequence
                    cycle_edge = current_node.get_cycle_edge(sequence_index)
                    cycle_edge.cycle += kmer[-1]
                    # keep current node the same
                else: # create a cycle_edge
                    if next_node.get_edge(sequence_index):# This sequence already passes through this node
                        current_node.edges.append(DBGEdge(target_node=None, sequence_index=sequence_index, cycle=kmer[-1]))
                        # keep current node the same
                    else:
                        cycle_edge = current_node.get_cycle_edge(sequence_index)
                        if cycle_edge: # it's a cycle we can close
                            cycle_edge.target_node = next_node
                        else:    
                            current_node.edges.append(DBGEdge(target_node=next_node, sequence_index=sequence_index)) 
                        current_node = next_node

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
            # if any node.edge has a non-empty cycle list then there is a cycle
            if any(edge.cycle for edge in node.edges):
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

        mermaid_str = "graph LR;\n "
        mermaid_str += "s(start);\n e(end);\n "
        queue = deque([self.root])
        visited = set()
        termini = []
        
        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            if not node:
                continue
            if not node.edges:
                termini.append(node)
            for edge in node.edges:
                target = edge.target_node
                if target not in visited:
                    queue.append(target)

                hide_kmers = '(" ")' if not show_kmers else ''    
                if node == self.root:
                    mermaid_str += f's --> {target.kmer};\n'
                else:
                    if edge.cycle:
                        if not target:
                            mermaid_str += f"{node.kmer}{hide_kmers} --{','.join(edge.cycle)}--> e;\n"
                        else:
                            mermaid_str += f"{node.kmer}{hide_kmers} --{','.join(edge.cycle)}--> {target.kmer}{hide_kmers};\n"
                    else:
                        if not target:
                            mermaid_str += f"{node.kmer}{hide_kmers} --> e;\n"
                        else:
                            mermaid_str += f"{node.kmer}{hide_kmers} --> {target.kmer}{hide_kmers};\n"
        for node in termini:
            mermaid_str += f"{node.kmer} --> e;\n"
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
from functools import singledispatchmethod
from typing import List, Set, Union

from .allignment_buffer import AlignmentBuffer
from .debruijngraph import DeBruijnGraph
from .dbg_node import DBGNode
from .pog_node import POG_Node
from .pog_bubble import POG_Bubble
from .constants import AlignmentMethod

class PartialOrderGraph:
    def __init__(self, debruijn_graph : DeBruijnGraph = None):
        if debruijn_graph is None:
            self.root = None
            self.sequence_names = {}  # dict keyed on sequence names, returns tuple containing index and lengths of the sequence
        else:
            self.sequence_names = debruijn_graph.sequence_names # dict keyed on sequence names, returns tuple containing index and lengths of the sequence
            self.transform_dbg_to_pog(debruijn_graph.root)

    def transform_dbg_to_pog(self, node : DBGNode):
        sequence_set = {value[0] for value in self.sequence_names.values()}
        self.root = POG_Node.from_dbg_node(node, sequence_set, read_full_kmer=True)

        # add a synthetic end node with no fragment to all final nodes
        end_node = POG_Node("", sequence_set)
        for sequence in sequence_set:
            node = self.root
            last_node = node
            while node:
                last_node = node
                node = node.get_next(sequence)
            end_node.add_node(last_node)

    def work(self, alignment_type: AlignmentMethod):
        """Returns the order complexity of aligninging the sequences."""
        if alignment_type == AlignmentMethod.EXACT:
            product = 1
            for length in [length for _, length in self.sequence_names.values()]:
                product *= length
            return product
        elif alignment_type == AlignmentMethod.PROGRESSIVE:
            # Extract sequence lengths and sort them
            sequence_lengths = sorted(length for _, (_, length) in self.sequence_names.items())
            # Sum the product of each length with the next one
            return sum(sequence_lengths[i] * sequence_lengths[i+1] for i in range(len(sequence_lengths) - 1))        
        elif alignment_type == AlignmentMethod.DEBRUIJNGRAPH:
            return sum(1 for _ in self.root)
        elif alignment_type == AlignmentMethod.BRAIDEDDEBRUIJGRAPH:
            return sum(1 for _ in self.root)
        else:
            raise ValueError("Unsupported alignment type")

    def __len__(self):
        return len(self.sequence_names)
    
    @singledispatchmethod
    def __getitem__(self, index: Union[int, str]):
        raise TypeError("Index must be a string or an integer")

    @__getitem__.register
    def _(self, index: int):
        if index < 1 or index > len(self):
            raise IndexError("Sequence index out of range")
        # Start the sequence reconstruction from the root node
        sequence = self.root.sequence(index)
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
        sequence = self.root.sequence(sequence_index)
        return sequence
    
    def bubbles(self)->List[POG_Bubble]:
        # if root.edges is empty then there are no bubbles - return an empty list
        if not self.root:
            return []
        else:
            bubbles = self.root.bubbles()
            # remove all leaf bubbles where the edge lengths are equal
            return bubbles
        
    def align(self, buffer : AlignmentBuffer):
        self.root.align(buffer)

    def names(self):
        """Returns an iterable collection of sequence names."""
        return list(self.sequence_names.keys())
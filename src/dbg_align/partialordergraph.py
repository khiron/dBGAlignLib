from functools import singledispatchmethod
from typing import List, Set

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
        elif alignment_type == AlignmentMethod.DBG_LENGTH:
            return sum(1 for _ in self.root)
        elif alignment_type == AlignmentMethod.DBG_LENGTH_NUMBER:
            return sum(1 for _ in self.root)
        elif alignment_type in (AlignmentMethod.DBG_LENGTH, AlignmentMethod.DBG_LENGTH_NUMBER):
            raise ValueError("Graph must be transformed to DirectedAcyclicGraph before calculating expected_work")
        else:
            raise ValueError("Unsupported alignment type")

    def __len__(self):
        return len(self.sequence_names)
    
    def __getitem__(self, index):
        return self.root.sequence(index)
    
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

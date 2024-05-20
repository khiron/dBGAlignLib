from functools import singledispatchmethod
from typing import List, Set
from .debruijngraph import DeBruijnGraph
from .dbg_node import DBGNode
from .dag_node import DAG_Node
from .dag_bubble import DAG_Bubble
from .constants import AlignmentMethod

class DirectedAcyclicGraph:
    def __init__(self, debruijn_graph : DeBruijnGraph = None):
        if debruijn_graph is None:
            self.root = None
            self.sequence_names = {}  # dict keyed on sequence names, returns tuple containing index and lengths of the sequence
        else:
            self.sequence_names = debruijn_graph.sequence_names
            self.root = self.transform_dbg_to_dag(debruijn_graph.root)

    def transform_dbg_to_dag(self, dbg_node : DBGNode):
        dag_node = DirectedAcyclicGraph.Node(dbg_node.sequence)
        for edge in dbg_node.edges:
            if edge.to_node is not None:
                dag_node.add_edge(self.transform_dbg_to_dag(edge.to_node))
        return dag_node
    
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
    
    def bubbles(self)->List[DAG_Bubble]:
        # if root.edges is empty then there are no bubbles - return an empty list
        if not self.root:
            return []
        else:
            bubbles = self.root.bubbles()
            # remove all leaf bubbles where the edge lengths are equal
            return bubbles
        

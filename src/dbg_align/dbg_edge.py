from typing import List
from .dbg_traversal import DBGTraversal
from .dbg_node import DBGNode


class DBGEdge:
    def __init__(self, target_node: DBGNode, cycle_count: int = 0) -> None:
        self.target_node = target_node
        self.cycle_count = cycle_count # a count of how many cycles preceded this edge
        self.sequences = set()
        self.traversals : List[DBGTraversal] = []  

    def add_sequence(self, sequence_index:int) -> "DBGEdge":
        self.sequences.add(sequence_index)
        return self
    
    def add_traversal(self, sequence_index:int, passage_index:int):
        self.traversals.append(DBGTraversal(sequence_index, passage_index))
        return self

    def __repr__(self):
        return f"Edge to {self.target_node.kmer} ({len(self.traversals)} traversals)"
    
    def label(self):
        return f"{[(t.sequence_index, t.passage_index) for t in self.traversals]}"
    
    def traversals_for_sequence(self, sequence_index:int):
        return [traversal for traversal in self.traversals if traversal.sequence_index == sequence_index]
    
    def next_node(self, sequence_index:int, passage_index:int):
        traversals = self.traversals_for_sequence(sequence_index)
        for traversal in traversals:
            if traversal.passage_index == passage_index:
                return self.target_node
        return None
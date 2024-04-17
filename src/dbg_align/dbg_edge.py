from typing import List
from .dbg_traversal import DBGTraversal
from .dbg_node import DBGNode


class DBGEdge:
    def __init__(self, target_node: DBGNode) -> None:
        self.target_node = target_node
        self.traversals : List[DBGTraversal] = []  

    def add_traversal(self, sequence_index:int, passage_index:int):
        self.traversals.append(DBGTraversal(sequence_index, passage_index))
        return self

    def __repr__(self):
        return f"Edge to {self.target_node.kmer} ({len(self.traversals)} traversals)"
    
    def label(self):
        return f"{[(t.sequence_index, t.passage_index) for t in self.traversals]}"
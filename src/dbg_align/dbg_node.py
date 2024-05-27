from collections import deque
from typing import List, Optional, Set, Tuple

class DBGNode:
    def __init__(self, kmer: str) -> None:
        self.kmer = kmer
        self.edges = []  # List of DeBrujinGraph_Edge objects

    def to_pog(self)->"POG_Node":
        current_node = self
        first_node = True
        while current_node.edges:
            for edge in current_node.edges:
                pass
    def get_sequence(self, sequence_index: int) -> Optional[str]:
        current_node = self
        sequence = ""
        if self.kmer:  # ie: this is not a root we are calling this from
            sequence += self.kmer
        while current_node:
            if current_node.kmer: # if we're not at the root
                if not sequence: # take the whole kmer
                    sequence += current_node.kmer
                else: # take the last character
                    sequence += current_node.kmer[-1]
            cycle = current_node.get_cycle(sequence_index)
            if cycle:
                sequence += "".join(current_node.get_cycle(sequence_index))
            next_node = current_node.get_next(sequence_index)
            if next_node is None:
                return sequence
            current_node = next_node
        return sequence
   
    def get_edge(self, sequence_id : int) -> Optional["DBGEdge"]:
        for edge in self.edges:
            if sequence_id == edge.sequence:
                return edge
        return None

    def get_cycle(self, sequence_id : int)-> List[str]:
        edge = self.get_edge(sequence_id)
        if edge and edge.cycle:
            return edge.cycle
        return None

    def get_cycle_edge(self, sequence_id : int)-> "DBGEdge":
        for edge in self.edges:
            if sequence_id == edge.sequence:
                if not edge.target_node:
                    return edge
        return None
    
    def get_next(self, sequence_id : int)->"DBGNode":
        edge = self.get_edge(sequence_id)
        if edge:
            return edge.target_node
        return None

    def __getitem__(self, index):
        return self.edges[index].target_node
    
    def __repr__(self):
        return f"Node:({self.kmer}) [{','.join([edge.target_node.kmer for edge in self.edges if edge.target_node is not None])}]"
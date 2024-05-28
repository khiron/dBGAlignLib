from __future__ import annotations # this is needed for forward references in type hints
from functools import singledispatchmethod
from typing import List, Set, Union

class POG_Node:
    def __init__(self, fragment : str = None, sequence_set: Set[int] = None):
        self.fragment = fragment
        self.sequence_set = sequence_set
        self.next = []

    def __add__(self, to_node: Union['POG_Node', List['POG_Node']]):
        # can't use @singledispatch with forward referenced types
        if isinstance(to_node, POG_Node):
            self.next.append(to_node)
        elif isinstance(to_node, list):
            self.next.extend(to_node)
        else:
            raise TypeError(f"Unsupported type {type(to_node)} for addition")
        return self

    @classmethod
    def from_dbg_node(cls, dbg_node : 'DBGNode', sequence_set: Set[int], read_full_kmer : bool = True):
        from .dbg_node import DBGNode
        if dbg_node.kmer == None: # is special case of root node
            instance = cls(None, sequence_set)

            unique_edges = {edge.target_node for edge in dbg_node.edges}
            for unique_edge in unique_edges:
                sequences = set()
                for edge in dbg_node.edges:
                    if edge.target_node == unique_edge:
                        sequences.add(edge.sequence) 
                instance += cls.from_dbg_node(unique_edge, sequences, True)

            for target, sequences in dbg_node.get_braids().items():
                child = cls.from_dbg_node(target, sequences, True)
                instance += child
        else:
            node = dbg_node
            kmer = node.kmer if read_full_kmer else node.kmer[-1]
            instance = cls(kmer, sequence_set)  

            while node and node.edges and node.edges_form_single_braid():# extend POG node until we reach a branch
                node = node.edges[0].target_node # in a single braid any edge will get you to the next node
                instance.fragment += node.kmer[-1] # add the last base of the kmer to the fragment
                if node.edges: # except if we are the terminal node
                    instance.fragment += node.edges[0].cycle #add any cycles to the fragment
            if not node:
                return instance
            
            for target, sequences in node.get_braids().items():
                child = cls.from_dbg_node(target, sequences, False)
                instance += child
        return instance

    def __getitem__(self, index: int):
        return self.next[index]
    
    def __len__(self):
        return len(self.next)
    
    def get_next(self, index: int) -> "POG_Node":
        # get the node in next that contains index in the sequence_numbers
        for node in self.next:
            if index in node.sequence_set:
                return node
        return None
    
    def sequence(self, index: int) -> str:
        if index not in self.sequence_set:
            return ''
        sequence = ''
        node = self
        while node is not None:
            fragment = node.fragment if node.fragment is not None else ''
            sequence += fragment
            node = node.get_next(index)
        return sequence

    def __repr__(self):
        return f"{self.sequence_set}:{self.fragment}:{len(self.next)}"

    def bubbles(self, depth : int = 0)->List['DAG_Bubble']:
        from .pog_bubble import POG_Bubble
        # if there are no edges then there are no bubbles
        if not self.next:
            return []
        else:
            # if there is only one edge then that's a weird case of an incomplete bubble
            if len(self.next) == 1:
                start = self # this is the start of the bubble
                # find the end node by following the edge
                end = self.next[0]
                inner_bubbles = []
                return [POG_Bubble(start, end, inner_bubbles, depth)]
            else:
                start = self
                # if there are multiple edges then we have at least 1 bubble splitting from this node
                # find the join node for this branch
                # follow the first edge until you get to a node with the same sequence_set
                end_candidate = start.next[0]
                first_child = end_candidate
                while not start.sequence_set.issubset(end_candidate.sequence_set): # if the sequence set starting a bubble is in a future node then that future node is the end of the bubble 
                    if not end_candidate.edges:
                        raise ValueError("Bubble never closes")
                    end_candidate = end_candidate.edges[0]
                inner_bubbles = first_child.bubbles(depth + 1)
                return [POG_Bubble(start, end_candidate, inner_bubbles, depth)]
        

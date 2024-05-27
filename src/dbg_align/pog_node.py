from __future__ import annotations # this is needed for forward references in type hints
from functools import singledispatchmethod
from typing import List, Set, Union

class POG_Node:
    def __init__(self, fragment : str = None, sequence_set: Set[int] = None):
        self.fragment = fragment
        self.sequence_set = sequence_set
        self.edges = []

    def __add__(self, to_node: Union['POG_Node', List['POG_Node']]):
        # can't use @singledispatch with forward referenced types
        if isinstance(to_node, POG_Node):
            self.edges.append(to_node)
        elif isinstance(to_node, list):
            self.edges.extend(to_node)
        else:
            raise TypeError(f"Unsupported type {type(to_node)} for addition")
        return self

    def __getitem__(self, index):
        return self.edges[index]
    
    def __len__(self):
        return len(self.edges)
    
    def sequence (self, index :int) ->str:
        # starting from root add all the fragments for each node along edges that include index in the sequence_numbers
        if index in self.sequence_set:
            return self.fragment + ''.join(edge.sequence(index) for edge in self.edges)
        else:
            return ''
        
    def __repr__(self):
        return f"{self.sequence_set}:{self.fragment}:{len(self.edges)}"

    def bubbles(self, depth : int = 0)->List['DAG_Bubble']:
        from .pog_bubble import POG_Bubble
        # if there are no edges then there are no bubbles
        if not self.edges:
            return []
        else:
            # if there is only one edge then that's a weird case of an incomplete bubble
            if len(self.edges) == 1:
                start = self # this is the start of the bubble
                # find the end node by following the edge
                end = self.edges[0]
                inner_bubbles = []
                return [POG_Bubble(start, end, inner_bubbles, depth)]
            else:
                start = self
                # if there are multiple edges then we have at least 1 bubble splitting from this node
                # find the join node for this branch
                # follow the first edge until you get to a node with the same sequence_set
                end_candidate = start.edges[0]
                first_child = end_candidate
                while not start.sequence_set.issubset(end_candidate.sequence_set): # if the sequence set starting a bubble is in a future node then that future node is the end of the bubble 
                    if not end_candidate.edges:
                        raise ValueError("Bubble never closes")
                    end_candidate = end_candidate.edges[0]
                inner_bubbles = first_child.bubbles(depth + 1)
                return [POG_Bubble(start, end_candidate, inner_bubbles, depth)]
        

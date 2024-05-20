from __future__ import annotations
from functools import singledispatchmethod
from typing import List, Set, Union

class DAG_Node:
    def __init__(self, fragment : str = None, sequence_set: Set[int] = None):
        self.fragment = fragment
        self.sequence_numbers = sequence_set
        self.edges = []

    def __add__(self, to_node: Union['DAG_Node', List['DAG_Node']]):
        # can't use @singledispatch with forward referenced types
        if isinstance(to_node, DAG_Node):
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
        if index in self.sequence_numbers:
            return self.fragment + ''.join(edge.sequence(index) for edge in self.edges)
        else:
            return ''
from __future__ import annotations # this is needed for forward references in type hints
from typing import List
from dbg_align.dag_node import DAG_Node


class DAG_Bubble:
    def __init__(self, split : DAG_Node, join : DAG_Node, inner_bubbles: List[DAG_Bubble] = None, depth : int = 0):
        self.split = split
        self.join = join
        self.depth = depth
        self.inner_bubbles = inner_bubbles
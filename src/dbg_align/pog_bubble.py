from __future__ import annotations # this is needed for forward references in type hints
from typing import List
from dbg_align.pog_node import POG_Node


class POG_Bubble:
    def __init__(self, start : POG_Node, end : POG_Node, inner_bubbles: List[POG_Bubble] = None, depth : int = 0):
        self.start = start
        self.end = end
        self.depth = depth
        self.inner_bubbles = inner_bubbles
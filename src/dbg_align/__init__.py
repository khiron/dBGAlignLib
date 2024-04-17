"""dbg_align: A package for aligning sequences using de Bruijn graphs."""

from .de_brujin_graph import DeBrujinGraph, AlignmentMethod
from .dbg_edge import DBGEdge
from .dbg_traversal import DBGTraversal
from .dbg_node import DBGNode


from .utils import display_mermaid_in_jupyter

__version__ = "0.1.0"
__author__ = "Richard Morris"

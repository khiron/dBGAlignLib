"""dBGAlign: A package for aligning sequences using de Bruijn graphs."""

from .de_brujin_graph import DeBrujinGraph
from .partial_order_graph import PartialOrderGraph

__version__ = "0.1.0"
__author__ = "Richard Morris"


def build_dBg(sequences, k) -> DeBrujinGraph:
    dBg = DeBrujinGraph(k)
    dBg.add_sequences(sequences)
    return dBg

def Build_POG(dBg: DeBrujinGraph) -> PartialOrderGraph:
    return PartialOrderGraph(dBg)
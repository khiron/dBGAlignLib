"""dBGAlign: A package for aligning sequences using de Bruijn graphs."""
from cogent3 import SequenceCollection
from .deBrujinGraph import DeBrujinGraph
__version__ = "0.1.0"
__author__ = "Richard Morris"
def build_dBg(sequences, k) -> DeBrujinGraph:
    dBg = DeBrujinGraph(k)
    dBg.add_sequences(sequences)
    return dBg
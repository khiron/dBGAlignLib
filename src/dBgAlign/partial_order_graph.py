from typing import List, Set

from .de_brujin_graph import DeBrujinGraph


class PartialOrderGraph:
    def __init__(self, debruijn_graph: DeBrujinGraph = None):
        if debruijn_graph:
            self.graph = self._compress(debruijn_graph)
        else:
            self.graph = {}
        self.original_kmer_length = (
            debruijn_graph.kmer_length if debruijn_graph else None
        )

    def _compress(self, debruijn_graph):
        self.original_kmer_length = debruijn_graph.kmer_length
        # Logic to compress the graph
        # from the starting dBg node
        # copy the kmer into the starting POG node's payload
        # if all edges in the dBg node are unique
        # then extend the POG node's payload to include the kmer for the next dBg node
        # if there are multiple unique edges in the dBg node
        # add a POG edge for each adding the dBg node's sequence index to the POG node sequences list
        pass

    def _uncompress(self, kmer_length: int = None) -> DeBrujinGraph:
        # Logic to uncompress the graph into a DeBruijnGraph
        kmer_length = kmer_length or self.original_kmer_length
        pass

    def uncompress(self) -> DeBrujinGraph:
        """Returns the uncompressed graph"""
        return self._uncompress(self)


class PartialOrderGraph_Node:
    def __init__(self, partial_order_graph: PartialOrderGraph, payload: str):
        self.payload = payload  # the kmer or sequence that this node represents
        self.edges = []  # 2 or more PartialOrderGraph_Edge objects
        self.partial_order_graph = partial_order_graph

    def Edges(self) -> List["PartialOrderGraph_Edge"]:
        return self.edges


class PartialOrderGraph_Edge:
    def __init__(self, target_node: PartialOrderGraph_Node, sequences: Set[int]):
        self.target_node = target_node  # the node that this edge points to
        self.sequences = (
            sequences  # the list of sequence indices that this edge represents
        )

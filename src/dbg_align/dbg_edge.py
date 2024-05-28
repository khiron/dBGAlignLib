from .dbg_node import DBGNode

class DBGEdge:
    def __init__(self, target_node: DBGNode, sequence_index : int = None, cycle : str = "") -> None:
        self.target_node = target_node
        self.sequence = sequence_index
        self.cycle = cycle

    def __repr__(self):
        return f"Edge ->{self.target_node.kmer} seq: ({self.sequences})"
    
    def label(self):
        return f"{','.join(self.sequences)}"
    

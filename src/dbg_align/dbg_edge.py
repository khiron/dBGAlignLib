from .dbg_node import DBGNode

class DBGEdge:
    def __init__(self, target_node: DBGNode, sequence_index : int = None, cycle : str = "") -> None:
        self.target_node = target_node
        self.sequences = set()
        if sequence_index is not None:
            self.sequences.add(sequence_index)
        self.cycle = cycle

    def __repr__(self):
        return f"Edge to {self.target_node.kmer} ({len(self.sequences)} sequences)"
    
    def label(self):
        return f"{','.join(self.sequences)}"
    

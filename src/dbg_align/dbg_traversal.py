class DBGTraversal:
    def __init__(self, sequence_index: int, passage_index: int) -> None:
        self.sequence_index = sequence_index
        self.passage_index = passage_index

    def __repr__(self):
        return f"({self.sequence_index}:{self.passage_index})"


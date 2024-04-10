from cogent3 import SequenceCollection, Sequence
class DeBrujinGraph:
    def __init__(self, kmer_length):
        self.kmer_length = kmer_length
        self.graph = {}
        self.start_nodes = {}  # Track the start node of each sequence

    def add_sequences(self, sequences : SequenceCollection):
        for seq_index, seq in enumerate(sequences.seqs):
            # Ensure each sequence has a unique start identifier
            start_node = f"start_{seq_index}"
            self.start_nodes[seq_index] = start_node
            previous_node = start_node

            for i in range(len(seq) - self.kmer_length + 1):
                kmer = seq[i:i+self.kmer_length]
                prefix = kmer[:-1]
                suffix = kmer[1:]

                if previous_node not in self.graph:
                    self.graph[previous_node] = {}
                if suffix not in self.graph[previous_node]:
                    self.graph[previous_node][suffix] = set()
                self.graph[previous_node][suffix].add(seq_index)

                previous_node = suffix

                # Ensure every suffix is represented in the graph, even if it has no outgoing edges
                if suffix not in self.graph:
                    self.graph[suffix] = {}

    def __getitem__(self, seq_index):
        if seq_index not in self.start_nodes:
            raise IndexError("Sequence index out of range")

        sequence = ''
        current_node = self.start_nodes[seq_index]
        while True:
            found = False
            for next_node, originating_sequences in self.graph[current_node].items():
                if seq_index in originating_sequences:
                    if 'start' not in current_node:  # Skip the initial start node
                        sequence += next_node[-1]  # Append the last character of the next node
                    current_node = next_node
                    found = True
                    break
            if not found:
                break

        return sequence

    def as_mermaid(self)->str:
        mermaid_str = "graph LR;\n"
        for node, edges in self.graph.items():
            for target, seq_indices in edges.items():
                label = ",".join(str(idx) for idx in seq_indices)
                node_label = node if not node.startswith('start') else node[:5] + node.split('_')[1]
                target_label = target
                mermaid_str += f'    {node_label} -->|"{label}"| {target_label};\n'
        return mermaid_str


from collections import deque
from typing import List, Optional, Set, Tuple

class DBGNode:
    def __init__(self, kmer: str, kmer_length: int) -> None:
        self.kmer = kmer
        self.kmer_length = kmer_length
        self.edges = []  # List of DeBrujinGraph_Edge objects
        self.parent_nodes = set()  # set of nodes that lead to this node

    def connect_node(self, node: "DBGNode", sequence_index: int, cycle_count: int):
        """Connect this node to another node and return the updated cycle_count."""
        from .dbg_edge import DBGEdge
        # if existing edge that does not include this sequence and it has a cycle count => cycle_count then add the sequence to that edge and return the edge's cycle_count
        for edge in self.edges:
            if sequence_index not in edge.sequences and edge.cycle_count >= cycle_count:
                edge.add_sequence(sequence_index)
                return edge.cycle_count
        # else add a new edge with the sequence and this cycle_count
        self.edges.append(DBGEdge(node, cycle_count).add_sequence(sequence_index))
        return cycle_count
            
    def next(self, sequences:Set[int], cycle_count:int) -> Tuple[Optional["DBGNode"], Optional[int]]:
        """Return the next node in the sequence for the specified set of sequences providing its cycle count does not decrease."""
        if self.edges:
            seq_edges = [edge for edge in self.edges if sequences in edge.sequences]
            if seq_edges:
                # Find the edge with the smallest edge.cyclecount >= cycle_count
                next_edge = min((edge for edge in seq_edges if edge.cyclecount >= cycle_count), 
                                key=lambda edge: edge.cyclecount, default=None)
                if next_edge:
                    return next_edge.target_node, next_edge.cyclecount
        return None, None

    def can_transform_to_pog(self) -> bool:
        """Check if the node can be transformed to the next node in the sequence for the specified sequence index and cycle count."""
        if not self.edges:
            return False # this is the last node
        if len(self.edges) == 1:
            return True
        # Check that all edges point to the same target_node
        if not all(edge.target_node == self.edges[0].target_node for edge in self.edges[1:]):
            return False
        # if no edge has a sequence in another edge it can be transformed
        return not set.intersection(*[edge.sequences for edge in self.edges])

    def transform_to_pog(self, cycle_count:int):
        current_node = self
        # find edge with the lowest cycle_count not below the current cycle_count
        next_edge = min((edge for edge in self.edges if edge.cycle_count >= cycle_count), 
                        key=lambda edge: edge.cycle_count, default=None)
        while next_edge:
            current_node = current_node.next(next_edge.sequences, next_edge.cycle_count)


            next_node = next_edge.target_node
            self.kmer += next_node.kmer[-1]
            cycle_count = next_edge.cycle_count
            # find edge on next_node with the lowest cycle_count not below the current cycle_count

            next_edge = next_edge.target_node.edges[0] if next_edge.target_node.edges else None

            self.edges = [next_edge]
                    



    @property 
    def parents(self):
        return self.parent_nodes
    
    @property
    def children(self):
        return {edge.target_node for edge in self.edges}

    @property
    def has_single_parent(self) -> bool:
        return len(self.parent_nodes) == 1
    
    @property
    def has_one_child(self) -> bool:
        return len(self.edges) == 1
    
    @property
    def first_child(self) -> "DBGNode":
        if not self.edges:
            raise ValueError("Node does not have any children")
        return self.edges[0].target_node
    
    @property
    def sole_child(self) -> "DBGNode":
        if not self.has_one_child:
            raise ValueError("Node does not have a single child")
        return self.edges[0].target_node
    
    @property
    def first_sequence_index(self) -> int:
        return self.first_child.traversals[0].sequence_index

    @property
    def starts_split(self) -> bool: 
        return len(self.edges) > 1
    
    def find_merge_node(self) -> "DBGNode":
        if not self.starts_split:
            raise ValueError(f"Node {self.kmer} does not start a split")
        # find the first node that has all the same sequences entering it as this node has exiting it
        current_node = self
        sequence_to_follow = self.first_sequence_index
        while current_node:
            if not all(edge.target_node.find(sequence_to_follow, 1) == current_node for edge in current_node.edges):
                break
            current_node = current_node.first_child
        return current_node

    def add_edge(self, target_node, sequence_index=None, passage_index=None):
        """Create a new edge to target_node if not already existing, or return existing one."""
        from dbg_align.dbg_edge import DBGEdge

        for edge in self.edges:
            if edge.target_node == target_node:
                if sequence_index is not None and any(traversal.sequence_index == sequence_index for traversal in edge.traversals):
                    continue # if there is already a traversal for this sequence then we are looping so create a new edge
                edge.add_traversal(sequence_index, passage_index)
                return edge
        new_edge = DBGEdge(target_node)
        new_edge.add_traversal(sequence_index, passage_index)
        self.edges.append(new_edge)
        target_node.parent_nodes.add(self) # Add this node as a parent to the target node
        return new_edge
 
    def can_compress(self):
        """Check if the node can be compressed."""
        return self.has_one_child and self.sole_child.has_single_parent

    def run(self):
        """Generates a list of nodes that form a linear run from this node."""
        run_nodes = [self]
        candidate = self
        while candidate.can_compress():
            child = candidate.first_child
            run_nodes.append(child)
            candidate = child
        return run_nodes
        
    def compress(self, run_nodes = None):
        if not run_nodes:
            run_nodes = self.run()
        if len(run_nodes) > 1:
            # Extend the base node's kmer and adopt the last node's edges
            last_node = self
            for node in run_nodes[1:]:
                self.kmer += node.kmer[-1]  # Assuming k-1 overlap
                last_node = node
            for child in last_node.children:
                if not last_node in child.parent_nodes:
                    pass
                child.parent_nodes.remove(last_node)
                child.parent_nodes.add(self)
            self.edges = last_node.edges

    def has_cycle(self, visited=None, ancestors=None):
        if visited is None:
            visited = set()
        if ancestors is None:
            ancestors = set()

        if self in ancestors:
            return True  # A cycle is detected.
        if self in visited:
            return False

        visited.add(self)
        ancestors.add(self)
        for edge in self.edges:
            if edge.target_node.has_cycle(visited, ancestors):
                return True
        ancestors.remove(self)
        return False
   
    def sequence_edges(self, sequence_index: int) -> List["DBGEdge"]:
        """Retrieve all edges that correspond to a specific sequence index."""
        if sequence_index is None:
            return self.edges
        filtered_edges = [edge for edge in self.edges 
                          if any(traversal.sequence_index == sequence_index 
                                 for traversal in edge.traversals)]
        return filtered_edges

    def sequence_passage_edges(self, sequence_index: int, passage_index: int) -> List["DBGEdge"]:
        """Retrieve all edges that correspond to a specific sequence and passage index."""
        filtered_edges = [edge for edge in self.edges 
                        if any(traversal.sequence_index == sequence_index 
                               and traversal.passage_index == passage_index 
                               for traversal in edge.traversals)]
        return filtered_edges

    def find(self, sequence_index: int, passage_index: int, visited=None):
        """return the node that contains an edge with the specified sequence and passage index."""
        if visited is None:
            visited = set()

        visited.add(self)  # mark the current node as visited

        for edge in self.sequence_edges(sequence_index):
            for traversal in edge.traversals:
                if traversal.passage_index == passage_index and traversal.sequence_index == sequence_index:
                    return self # found the target
            if edge.target_node not in visited:
                found_node = edge.target_node.find(sequence_index, passage_index, visited)
                if found_node:
                    return found_node  # return only if a node was found
        return None  # return None if no matching node was found

    def next_in_sequence(self, sequence_index: int, current_passage_index: int) -> "DBGNode":
        """Return the next node in the sequence for the specified sequence index 
           where the passage index is lowest but no lower than the current_passage_index."""
        valid_traversals = []
        for edge in self.edges:
            for traversal in edge.traversals:
                if traversal.sequence_index == sequence_index and traversal.passage_index >= current_passage_index:
                    valid_traversals.append((traversal.passage_index, edge.target_node))

        if valid_traversals:
            passage_index, target_node = min(valid_traversals, key=lambda x: x[0])
            return target_node, passage_index
        else:
            return None, current_passage_index

    def count_characters(self, sequence_index: int, current_passage_index: int = 0) -> int:
        """
        Count the number of characters remaining in the sequence from this node. 
        starting with the smallest passage index above the current passage index.
        """
        current_node = self
        character_count = len(current_node.kmer)
        passage_index = current_passage_index
        while current_node:
            current_node, passage_index = current_node.next_in_sequence(sequence_index, passage_index)
        return character_count
            
    def get_sequence(self, sequence_index: int) -> str:
        current_node = self
        current_passage_index = 0
        sequence = ""
        first_node = True
        overlap_length = self.kmer_length - 1

        while current_node:
            current_node, current_passage_index = current_node.next_in_sequence(sequence_index, current_passage_index)
            if current_node:
                if current_node.kmer: # nb: root node has no kmer
                    if first_node: # start with the full kmer
                        sequence += current_node.kmer
                        first_node = False
                    else:
                        sequence += current_node.kmer[overlap_length:] # in subsequent appends skip the overlap (eg: first 2 characters in 3mers)
        return sequence

    def traverse_all(self):
        """Traverse all nodes in the graph starting from this node."""
        visited = set()
        queue = deque([self])

        while queue:
            node = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            yield(node)
            for edge in node.edges:
                if edge.target_node not in visited:
                    queue.append(edge.target_node)

    def __getitem__(self, index):
        return self.edges[index].target_node
    
    def __repr__(self):
        return f"Node:({self.kmer}) [{','.join([edge.target_node.kmer for edge in self.edges])}]"


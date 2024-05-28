from typing import Any, List, Union, Tuple

from .alignment import AlignmentPlugin
from .alignment_operation import AlignmentOperation

class ProfileCostAlignment:
    def __init__(self, cost: Tuple[int, int], length_range: Tuple[int, int]):
        self.cost = cost
        self.length_range = length_range

    def __eq__(self, other):
        if not isinstance(other, ProfileCostAlignment):
            return NotImplemented
        return self.cost == other.cost and self.length_range == other.length_range

    def __repr__(self):
        return f"ProfileCostAlignment(cost={self.cost}, length_range={self.length_range})"        

class PluginCostAlignment(AlignmentPlugin):
    def __init__(self):
        self.records: list = []

    def align_sequences(self, seq1: str, seq2: str) -> ProfileCostAlignment:
        cost = len(seq1) * len(seq2)
        length_range = (max(len(seq1), len(seq2)), len(seq1) + len(seq2))
        self.records.append((AlignmentOperation.SEQUENCE_SEQUENCE, seq1, seq2, cost, length_range))
        return ProfileCostAlignment((cost, cost), length_range)

    def align_sequence_to_profile(self, seq: str, profile: ProfileCostAlignment) -> ProfileCostAlignment:
        cost_min = len(seq) * profile.cost[0]
        cost_max = len(seq) * profile.cost[1]
        length_range = (max(len(seq), profile.length_range[0]), len(seq) + profile.length_range[1])
        self.records.append((AlignmentOperation.SEQUENCE_PROFILE, seq, profile, (cost_min, cost_max), length_range))
        return ProfileCostAlignment((cost_min, cost_max), length_range)

    def align_profiles(self, profile1: ProfileCostAlignment, profile2: ProfileCostAlignment) -> ProfileCostAlignment:
        cost_min = profile1.length_range[0] * profile2.length_range[0]
        cost_max = profile1.length_range[1] * profile2.length_range[1]
        length_range = (max(profile1.length_range[0], profile2.length_range[0]), profile1.length_range[1] + profile2.length_range[1])
        self.records.append((AlignmentOperation.PROFILE_PROFILE, profile1, profile2, (cost_min, cost_max), length_range))
        return ProfileCostAlignment((cost_min, cost_max), length_range)

    def concatenate(self, elements: List[Union[str, ProfileCostAlignment]]) -> ProfileCostAlignment:
        total_min_len = 0
        total_max_len = 0
        total_min_cost = 0
        total_max_cost = 0

        for element in elements:
            if isinstance(element, str):
                total_min_len += len(element)
                total_max_len += len(element)
            elif isinstance(element, ProfileCostAlignment):
                total_min_len += element.length_range[0]
                total_max_len += element.length_range[1]

        self.records.append(('concatenate', elements, (total_min_cost, total_max_cost), (total_min_len, total_max_len)))
        return ProfileCostAlignment((total_min_cost, total_max_cost), (total_min_len, total_max_len))

    def get_records(self) -> list:
        return self.records

class TotalCost():
    from .allignment_buffer import AlignmentBuffer
    def __init__(self, buffer: AlignmentBuffer):
        if not isinstance(buffer.alignment_plugin, PluginCostAlignment):
            raise ValueError("Buffer must use PluginCostAlignment")
        self.buffer = buffer
        self.total_cost = (0,0)

    @property 
    def cost(self) -> Tuple[int,int]:
        if self.total_cost == (0,0):
            self.calculate_total_cost()
        return self.total_cost

    def calculate_total_cost(self) -> Tuple[int,int]:
        # add the cost tuple of each operation in the buffer
        for result in self.buffer.results.values():
            self.total_cost = (self.total_cost[0] + result.cost[0], self.total_cost[1] + result.cost[1])



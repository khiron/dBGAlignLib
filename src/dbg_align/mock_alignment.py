
from typing import List, Union
from .alignment import AlignmentPlugin
from .composite_alignment import CompositeAlignment


class MockAlignmentPlugin(AlignmentPlugin):
    def __init__(self):
        self.records: list = []

    def align_sequences(self, seq1: str, seq2: str) -> CompositeAlignment:
        profile = CompositeAlignment([seq1, seq2])
        self.records.append(('sequence', seq1, seq2))
        return profile

    def align_sequence_to_profile(self, seq: str, profile: CompositeAlignment) -> CompositeAlignment:
        new_profile = CompositeAlignment([seq, profile])
        self.records.append(('sequence_profile', seq, profile))
        return new_profile

    def align_profiles(self, profile1: CompositeAlignment, profile2: CompositeAlignment) -> CompositeAlignment:
        combined_profile = CompositeAlignment([profile1, profile2])
        self.records.append(('profile', profile1, profile2))
        return combined_profile

    def concatenate(self, elements: List[Union[str, CompositeAlignment]]) -> CompositeAlignment:
        concatenated_elements = []
        for element in elements:
            if isinstance(element, str):
                concatenated_elements.append(element)
            elif isinstance(element, CompositeAlignment):
                concatenated_elements.append(element.concatenate())
        concatenated_profile = CompositeAlignment(concatenated_elements)
        self.records.append(('concatenate', elements))
        return concatenated_profile

    def get_records(self) -> list:
        return self.records

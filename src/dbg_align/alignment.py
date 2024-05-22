from abc import ABC, abstractmethod
from typing import Any
from .alignment_operation import AlignmentOperation
from .composite_alignment import CompositeAlignment


from abc import ABC, abstractmethod

class AlignmentPlugin(ABC):
    @abstractmethod
    def align_sequences(self, seq1: str, seq2: str) -> Any:
        """
        Align two sequences and return a profile.

        Parameters
        ----------
        seq1 : str
            The first sequence to align.
        seq2 : str
            The second sequence to align.

        Returns
        -------
        Any
            The resulting profile.
        """
        pass

    @abstractmethod
    def align_sequence_to_profile(self, seq: str, profile: Any) -> Any:
        """
        Align a sequence to a previously performed profile.

        Parameters
        ----------
        seq : str
            The sequence to align.
        profile : Any
            The previously performed profile.

        Returns
        -------
        Any
            The resulting profile.
        """
        pass

    @abstractmethod
    def align_profiles(self, profile1: Any, profile2: Any) -> Any:
        """
        Align two previously performed profiles.

        Parameters
        ----------
        profile1 : Any
            The first profile to align.
        profile2 : Any
            The second profile to align.

        Returns
        -------
        Any
            The resulting profile.
        """
        pass

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

    def get_records(self) -> list:
        return self.records

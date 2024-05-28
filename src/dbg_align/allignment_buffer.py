from typing import Any, Dict, List, Tuple, Union
from .composite_alignment import AlignmentOperation, CompositeAlignment
from .alignment import AlignmentPlugin
from .alignment_operation import AlignmentOperation
from functools import singledispatch

class AlignmentBuffer:
    def __init__(self, alignment_plugin: callable):
        self.alignment_plugin: AlignmentPlugin = alignment_plugin
        self.operations: list = []
        self.results: Dict[int, Any] = {}
        self.next_index: int = 0

    def add_alignment(self, *args: Any) -> int:
        methods = {
            (str, str): self.add_sequences,
            (str, int): self.add_sequence_to_profile,
            (int, str): self.add_profile_to_sequence,
            (int, int): self.add_profiles,
        }
        arg_types = tuple(map(type, args))
        
        if len(args) == 2 and arg_types in methods:
            return methods[arg_types](*args)
        else:
            raise NotImplementedError("Unsupported argument types.")        
    def add_sequences(self, seq1: str, seq2: str) -> int:
        self.operations.append((AlignmentOperation.SEQUENCE_SEQUENCE, seq1, seq2))
        result = self.alignment_plugin.align_sequences(seq1, seq2)
        self.results[self.next_index] = result
        current_index = self.next_index
        self.next_index += 1
        return current_index

    def add_sequence_to_profile(self, seq: str, profile_index: int) -> int:
        self.operations.append((AlignmentOperation.SEQUENCE_PROFILE, seq, profile_index))
        result = self.alignment_plugin.align_sequence_to_profile(seq, self.results[profile_index])
        self.results[self.next_index] = result
        current_index = self.next_index
        self.next_index += 1
        return current_index

    def add_profile_to_sequence(self, profile_index: int, seq: str) -> int:
        self.operations.append((AlignmentOperation.SEQUENCE_PROFILE, profile_index, seq))
        result = self.alignment_plugin.align_sequence_to_profile(seq, self.results[profile_index])
        self.results[self.next_index] = result
        current_index = self.next_index
        self.next_index += 1
        return current_index

    def add_profiles(self, profile_index1: int, profile_index2: int) -> int:
        self.operations.append((AlignmentOperation.PROFILE_PROFILE, profile_index1, profile_index2))
        result = self.alignment_plugin.align_profiles(self.results[profile_index1], self.results[profile_index2])
        self.results[self.next_index] = result
        current_index = self.next_index
        self.next_index += 1
        return current_index

    def concatenate(self, elements: list) -> int:
        result = self.alignment_plugin.concatenate([self.results[e] if isinstance(e, int) else e for e in elements])
        self.results[self.next_index] = result
        current_index = self.next_index
        self.next_index += 1
        return current_index

    def build_structure(self) -> None:
        for index, op in enumerate(self.operations):
            op_type, *args = op
            if op_type == AlignmentOperation.SEQUENCE_SEQUENCE:
                seq1, seq2 = args
                profile = self.alignment_plugin.align_sequences(seq1, seq2)
            elif op_type == AlignmentOperation.SEQUENCE_PROFILE:
                seq, profile_index = args
                profile = self.alignment_plugin.align_sequence_to_profile(seq, self.results[profile_index])
            elif op_type == AlignmentOperation.PROFILE_PROFILE:
                profile_index1, profile_index2 = args
                profile = self.alignment_plugin.align_profiles(self.results[profile_index1], self.results[profile_index2])
            self.results[index] = profile

    def process_structure(self) -> Dict[int, Any]:
        final_results = {}
        for key, profile in self.results.items():
            elements = []
            if isinstance(key, int):
                elements.append(self.results[key])
            final_profile = self.alignment_plugin.concatenate(elements)
            final_results[key] = final_profile
        return final_results

    def clear(self) -> None:
        self.operations = []
        self.results = {}
        self.next_index = 0


from enum import Enum

class AlignmentOperation(Enum):
    SEQUENCE_SEQUENCE = 1
    SEQUENCE_PROFILE = 2
    PROFILE_PROFILE = 3

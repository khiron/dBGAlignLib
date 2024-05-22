from typing import List, Union, Any
from .alignment_operation import AlignmentOperation

class CompositeAlignment:
    def __init__(self, components: List[Union[str, 'CompositeAlignment']]):
        """
        Initialize a CompositeAlignment object.

        Parameters
        ----------
        components : list
            A list of components which can be either string fragments or other CompositeAlignments.
        """
        self.components = components

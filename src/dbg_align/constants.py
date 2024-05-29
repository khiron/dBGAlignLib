from enum import Enum


class AlignmentMethod(Enum):
    EXACT = "Exact"
    PROGRESSIVE = "Progressive"
    DEBRUIJNGRAPH = "de Bruijn Graph"
    BRAIDEDDEBRUIJGRAPH = "braided de Bruijn Graph"

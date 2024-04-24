import pytest
from dbg_align import DeBrujinGraph, AlignmentMethod
import cogent3

def test_analysis():
    dbg = DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCAT"
        })
    assert dbg.order_complexity(AlignmentMethod.RAW)==1404
    with pytest.raises(ValueError):
        assert dbg.order_complexity(AlignmentMethod.DBG1D)==360
        assert dbg.order_complexity(AlignmentMethod.DBG2D)==65
    dbg.to_pog()
    assert dbg.order_complexity(AlignmentMethod.RAW)==1404
    assert dbg.order_complexity(AlignmentMethod.DBG1D)==0
    assert dbg.order_complexity(AlignmentMethod.DBG2D)==0


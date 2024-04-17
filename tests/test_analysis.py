import pytest
from dbg_align import DeBrujinGraph, AlignmentMethod
import cogent3

@pytest.mark.xfail(reason="order_complexity is not yet implemented")
def test_analysis():
    dbg = DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCAT"
        })
    assert dbg.order_complexity(AlignmentMethod.RAW)==1404
    assert dbg.order_complexity(AlignmentMethod.DBG1D)==360
    assert dbg.order_complexity(AlignmentMethod.DBG2D)==65

import pytest
from dbg_align import DeBruijnGraph, AlignmentMethod
import cogent3

def test_analysis():
    dbg = DeBruijnGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCAT"
        })
    assert dbg.expected_work(AlignmentMethod.RAW)==1404
    with pytest.raises(ValueError):
        assert dbg.expected_work(AlignmentMethod.DBG1D)==360
        assert dbg.expected_work(AlignmentMethod.DBG2D)==65
    dbg.to_pog()
    assert dbg.expected_work(AlignmentMethod.RAW)==1404
    assert dbg.expected_work(AlignmentMethod.DBG1D)==0
    assert dbg.expected_work(AlignmentMethod.DBG2D)==0


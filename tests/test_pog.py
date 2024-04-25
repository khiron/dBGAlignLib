import cogent3
import dbg_align

def test_pog_reconstitutes_faithfully():
    dbg = dbg_align.DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCAT"
        })
    pog = dbg.to_pog()
    assert len(pog) == 3
    assert pog.names() == ["seq1", "seq2", "seq3"]
    assert pog["seq1"] == "ACAGTACGGCAT"
    assert pog["seq2"] == "ACAGTACTGGCAT"
    assert pog["seq3"] == "ACAGCGCAT"

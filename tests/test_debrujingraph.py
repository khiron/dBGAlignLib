import cogent3
import dBgAlign

def test_create_empty_dbg():
    dbg = dBgAlign.DeBrujinGraph(3)
    assert len(dbg) == 0
    assert dbg.names() == []
    assert dbg.moltype == cogent3.DNA

def test_create_dbg_from_string():
    dbg = dBgAlign.DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence("ACGT")
    assert len(dbg) == 1
    assert dbg.names() == ["Sequence_1"]

def test_has_cycles():
    dbg = dBgAlign.DeBrujinGraph(3)
    dbg.add_sequence("ACGTCATGCA")
    assert not dbg.has_cycles()
    dbg.add_sequence("ACATCATGCA")
    assert dbg.has_cycles()

def test_has_bubbles(output_dir):
    dbg = dBgAlign.DeBrujinGraph(3)
    dbg.add_sequence(["ACGTCATGCA", "ACGTCATGCA"])
    assert not dbg.has_bubbles()
    dbg.add_sequence(["ACGTCATGCA", "ACGTCATCATGCA"])
    assert dbg.has_bubbles()
    # write out the mermaid file
    with open(output_dir / "test_has_bubbles.mmd", "w") as f:
        f.write(dbg.to_mermaid())

def test_create_dbg_from_list():
    dbg = dBgAlign.DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence(["ACGT", "CGTA"])
    assert len(dbg) == 2
    assert dbg.names() == ["Sequence_1", "Sequence_2"]

def test_create_dbg_from_dict():
    dbg = dBgAlign.DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCAT"
        })
    assert len(dbg) == 3
    assert dbg.names() == ["seq1", "seq2", "seq3"]
    assert dbg[1] == "ACAGTACGGCAT"
    assert dbg[2] == "ACAGTACTGGCAT"
    assert dbg[3] == "ACAGCGCAT"
    assert dbg['seq1'] == "ACAGTACGGCAT"
    assert dbg['seq2'] == "ACAGTACTGGCAT"
    assert dbg['seq3'] == "ACAGCGCAT"

def test_sequence_reconstruction():
    dbg = dBgAlign.DeBrujinGraph(3, cogent3.DNA)
    sequences = {
        "seq1": "ACGTGAC",
        "seq2": "TACGTGA",
    }
    dbg.add_sequence(sequences)
    assert dbg["seq1"] == "ACGTGAC"
    assert dbg["seq2"] == "TACGTGA"
    # Test with specific start and length
    assert dbg.root.get_sequence(1, start_passage_index=2, length=4) == "CGTG"

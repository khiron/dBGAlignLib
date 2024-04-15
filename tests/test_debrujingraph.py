import cogent3
import dbg_align

def test_create_empty_dbg():
    dbg = dbg_align.DeBrujinGraph(3)
    assert len(dbg) == 0
    assert dbg.names() == []
    assert dbg.moltype == cogent3.DNA

def test_create_dbg_from_string():
    dbg = dbg_align.DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence("ACGT")
    assert len(dbg) == 1
    assert dbg.names() == ["Sequence_1"]
    assert dbg.root.__repr__() == "Node:(None) [ACG]"
    assert dbg.root[0].__repr__() == "Node:(ACG) [CGT]"
    assert dbg.root[0][0].kmer == "CGT"


def test_has_cycles():
    dbg = dbg_align.DeBrujinGraph(3)
    dbg.add_sequence("ACGTCATGCA")
    assert not dbg.has_cycles()
    dbg.add_sequence("ACATCATGCA")
    assert dbg.has_cycles()

def test_create_dbg_from_list():
    dbg = dbg_align.DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence(["ACGT", "CGTA"])
    assert len(dbg) == 2
    assert dbg.names() == ["Sequence_1", "Sequence_2"]
    assert dbg.root.__repr__() == "Node:(None) [ACG,CGT]"
    assert len(dbg.root) == 2
    assert dbg.root[0].kmer == "ACG"
    assert dbg.root[0][0].kmer == "CGT"
    assert dbg.root[1].kmer == "CGT"
    assert dbg.root[1][0].kmer == "GTA"

def test_create_dbg_from_dict():
    dbg = dbg_align.DeBrujinGraph(3,cogent3.DNA)
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
    dbg = dbg_align.DeBrujinGraph(3, cogent3.DNA)
    sequences = {
        "seq1": "ACGTGAC",
        "seq2": "TACGTGA",
    }
    dbg.add_sequence(sequences)
    assert dbg["seq1"] == "ACGTGAC"
    assert dbg["seq2"] == "TACGTGA"
    # Test with specific start and length
    assert dbg.root.get_sequence(1, start_passage_index=2, length=4) == "CGTG"

def test_compress(output_dir):
    dbg = dbg_align.DeBrujinGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCAT"
        })
    with open(output_dir / "nested_bubble.md", "w") as f:
        f.write("```mermaid\n")
        f.write(dbg.to_mermaid())
        f.write("```")    
    dbg.compress_graph()
    # write mermaid out to testout folder
    with open(output_dir / "nested_bubble_compressed.md", "w") as f:
        f.write("```mermaid\n")
        f.write(dbg.to_mermaid())
        f.write("```")
    
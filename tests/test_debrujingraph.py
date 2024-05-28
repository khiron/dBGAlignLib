import cogent3
import dbg_align


def test_create_kmers():
    seq = 'AGCT'
    kmers = dbg_align.DeBruijnGraph.generate_kmers(seq, 3)
    assert list(kmers) == ['AGC', 'GCT']

def test_DBGNode():
    node = dbg_align.DBGNode("ACG")
    assert node.kmer == "ACG"
    assert node.__repr__() == "Node:(ACG) []"
    assert node.__str__() == "ACG"

def test_braids():
    node1 = dbg_align.DBGNode("ACG")
    node2 = dbg_align.DBGNode("CGT")
    node3 = dbg_align.DBGNode("GTA")

    node1.add_edge(dbg_align.DBGEdge(node2, 1))
    node1.add_edge(dbg_align.DBGEdge(node2, 2))
    node1.add_edge(dbg_align.DBGEdge(node3, 1))
    assert len(node1.edges) == 3

    braids = node1.get_braids()
    assert len(braids) == 2
    assert braids == {node2: {1, 2}, node3: {1}}


def test_create_empty_dbg():
    dbg = dbg_align.DeBruijnGraph(3)
    assert len(dbg) == 0
    assert dbg.names() == []
    assert dbg.moltype == cogent3.DNA

def test_create_dbg_from_string():
    dbg = dbg_align.DeBruijnGraph(3,cogent3.DNA)
    dbg.add_sequence("ACGT")
    assert len(dbg) == 1
    assert dbg.names() == ["Sequence_1"]
    assert dbg.root.__repr__() == "Node:(None) [ACG]"
    assert dbg.root[0].__repr__() == "Node:(ACG) [CGT]"
    assert dbg.root[0][0].kmer == "CGT"

def test_has_cycles():
    dbg = dbg_align.DeBruijnGraph(3)
    dbg.add_sequence("ACGTCATGCA")
    assert not dbg.has_cycles()
    dbg.add_sequence("ACATCATGCA")
    assert dbg.has_cycles()
    assert dbg[2] == "ACATCATGCA" 

def test_create_dbg_from_list():
    dbg = dbg_align.DeBruijnGraph(3,cogent3.DNA)
    dbg.add_sequence(["ACGT", "CGTA"])
    assert len(dbg) == 2
    assert dbg.names() == ["Sequence_1", "Sequence_2"]
    assert dbg.root.__repr__() == "Node:(None) [ACG,CGT]"
    assert len(dbg.root.edges) == 2
    assert dbg.root[0].kmer == "ACG"
    assert dbg.root[0][0].kmer == "CGT"
    assert dbg.root[1].kmer == "CGT"
    assert dbg.root[1][0].kmer == "GTA"

def test_create_dbg_from_dict():
    dbg = dbg_align.DeBruijnGraph(3,cogent3.DNA)
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

def test_create_dbg_from_cogent3_sequences():
    dbg = dbg_align.DeBruijnGraph(3,cogent3.DNA)
    seq1 = cogent3.DNA.make_seq("ACAGTACGGCAT")
    seq2 = cogent3.DNA.make_seq("ACAGTACTGGCAT")
    seq3 = cogent3.DNA.make_seq("ACAGCGCAT")
    dbg.add_sequence({
        "seq1": seq1, 
        "seq2": seq2, 
        "seq3": seq3
        })
    assert len(dbg) == 3
    assert dbg.names() == ["seq1", "seq2", "seq3"]
    assert dbg[1] == "ACAGTACGGCAT"
    assert dbg[2] == "ACAGTACTGGCAT"
    assert dbg[3] == "ACAGCGCAT"

def test_cycle():
    dbg = dbg_align.DeBruijnGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCAT"
        })
    assert not dbg.has_cycles()
    dbg.add_sequence("ACATCATGCA")
    assert dbg.has_cycles()

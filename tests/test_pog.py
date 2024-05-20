from pathlib import Path
import cogent3
import dbg_align

def test_pog_reconstitutes_faithfully():
    dbg = dbg_align.DeBruijnGraph(3,cogent3.DNA)
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

def test_compress_to_pog(output_dir: Path):
    dbg = dbg_align.DeBruijnGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCAT"
        })
    with open(output_dir / "nested_bubble.md", "w") as f:
        f.write("```mermaid\n")
        f.write(dbg.to_mermaid())
        f.write("```")    
    dbg.to_pog()
    # write mermaid out to testout folder
    with open(output_dir / "nested_bubble_compressed.md", "w") as f:
        f.write("```mermaid\n")
        f.write(dbg.to_mermaid())
        f.write("```")

def test_pog_cycle(output_dir: Path):
    dbg = dbg_align.DeBruijnGraph(3,cogent3.DNA)
    dbg.add_sequence({
        "seq1": "ACAGTACGGCAT", 
        "seq2": "ACAGTACTGGCAT", 
        "seq3":"ACAGCGCGCAT" # contains cycle
        })
    with open(output_dir / "cycle.md", "w") as f:
        f.write("```mermaid\n")
        f.write(dbg.to_mermaid())
        f.write("```")   
    assert dbg.has_cycles()
    assert len(dbg) == 3
    assert dbg.names() == ["seq1", "seq2", "seq3"]
    assert dbg["seq1"] == "ACAGTACGGCAT"
    assert dbg["seq2"] == "ACAGTACTGGCAT"
    assert dbg["seq3"] == "ACAGCGCGCAT"
     
    dbg.to_pog()
    assert dbg["seq1"] == "ACAGTACGGCAT"
    assert dbg["seq2"] == "ACAGTACTGGCAT"
    assert dbg["seq3"] == "ACAGCGCGCAT"
    assert dbg.kmers("ACGG")
    assert dbg.kmers("ACTGG")
    assert dbg.kmers("AGCGCGC")
    # write mermaid out to testout folder
    with open(output_dir / "cycle_compressed.md", "w") as f:
        f.write("```mermaid\n")
        f.write(dbg.to_mermaid())
        f.write("```")


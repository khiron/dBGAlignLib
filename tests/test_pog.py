from dbg_align import PartialOrderGraph, POG_Node
from dbg_align.alignment_cost_plugin import PluginCostAlignment
from dbg_align.allignment_buffer import AlignmentBuffer
from dbg_align.constants import AlignmentMethod
import cogent3
import dbg_align
from pathlib import Path


def cost_alignment_plugin_factory() -> PluginCostAlignment:
    return PluginCostAlignment()

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
     
    pog = dbg.to_pog()
    assert pog["seq1"] == "ACAGTACGGCAT"
    assert pog["seq2"] == "ACAGTACTGGCAT"
    assert pog["seq3"] == "ACAGCGCGCAT"


def test_construct_pog():
    pog = PartialOrderGraph()
    pog.sequence_names = {'Sequence 1':(0,10),'Sequence 2':(1,10)}  # dict keyed on sequence names, returns tuple containing index and lengths of the sequence
    end = POG_Node("GCAT", {0,1}) 
    pog.root = POG_Node("AGT", {0, 1}) + [POG_Node("GCG", {0})+end, POG_Node("GTG",{1})+end]
    assert pog.root.fragment == "AGT"
    assert len(pog) == 2
    assert pog[0] == "AGTGCGGCAT"
    assert pog[1] == "AGTGTGGCAT"
    assert pog.expected_work(AlignmentMethod.EXACT) == 100
    assert pog.expected_work(AlignmentMethod.PROGRESSIVE) == 100
    # assert dag.expected_work(AlignmentMethod.DBG_LENGTH) == 9
    # assert dag.expected_work(AlignmentMethod.DBG_LENGTH_NUMBER) == 9
    bubbles = pog.bubbles()
    assert len(bubbles) == 1
    assert bubbles[0].start.fragment == "AGT"
    assert bubbles[0].end.fragment == "GCAT"
    assert bubbles[0].depth == 0
    assert bubbles[0].start.sequence_set == {0,1}
    assert bubbles[0].end.sequence_set == {0,1}
    assert len(bubbles[0].inner_bubbles) == 1

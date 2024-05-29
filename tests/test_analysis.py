from pathlib import Path
import sys
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

def test_primate_seven(data_dir: Path):
    filename = data_dir / "BRCA1" / "primates" / "brca1.fasta"

    graph = DeBruijnGraph(kmer_length=12, moltype=cogent3.DNA)
    sequences_collection = cogent3.load_unaligned_seqs(str(filename), format='FASTA', moltype=cogent3.DNA)
    for sequence in sequences_collection.iter_seqs():
        graph.add_sequence(sequence)
    pog = graph.to_pog()
    assert pog.work(AlignmentMethod.RAW)==0

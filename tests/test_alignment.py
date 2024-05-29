from dbg_align.alignment_cost_plugin import PluginCostAlignment, ProfileCostAlignment, TotalCost
from dbg_align.allignment_buffer import AlignmentBuffer
from dbg_align.partialordergraph import PartialOrderGraph
from dbg_align import partialordergraph

def test_cost_calculation():
    cost_alignment = PluginCostAlignment()
    
    # Test case 1: Align sequences "CGT" and "TACT"
    expected_output = ProfileCostAlignment((3*4, 3*4), (max(3, 4), 3 + 4))
    result = cost_alignment.align_sequences("CGT", "TACT")
    assert result == expected_output, f"Expected {expected_output}, but got {result}"
    
    # Test case 2: Align sequence "AGT" to profile ((12, 12), (4, 7))
    expected_output = ProfileCostAlignment((3*12, 3*12), (max(3, 4), 3 + 7))
    result = cost_alignment.align_sequence_to_profile("AGT", ProfileCostAlignment((12, 12), (4, 7)))
    assert result == expected_output, f"Expected {expected_output}, but got {result}"
    
    # Test case 3: Align profile ((12, 12), (4, 7)) to profile ((12, 12), (4, 7))
    expected_output = ProfileCostAlignment((12*12, 12*12), (max(4, 4), 7 + 7))
    result = cost_alignment.align_profiles(ProfileCostAlignment((12, 12), (4, 7)), ProfileCostAlignment((12, 12), (4, 7)))
    assert result == expected_output, f"Expected {expected_output}, but got {result}"

    # Test case 4: Concatenate ["AGT", ((12, 12), (4, 7)), "GCAT"]
    elements = ["AGT", ProfileCostAlignment((12, 12), (4, 7)), "GCAT"]
    expected_output = ProfileCostAlignment((12, 12), (3 + 4 + 4, 3 + 7 + 4))
    result = cost_alignment.concatenate(elements)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"

def test_alignment_buffer():
    cost_alignment = PluginCostAlignment()
    buffer = AlignmentBuffer(cost_alignment)

    # Test aligning two sequences
    index1 = buffer.add_alignment("CGT", "TACT")
    expected_result1 = ProfileCostAlignment((12, 12), (4, 7))
    assert buffer.results[index1] == expected_result1, f"Expected {expected_result1}, but got {buffer.results[index1]}"

    # Test aligning a sequence to a profile
    index2 = buffer.add_alignment("AGT", index1)
    expected_result2 = ProfileCostAlignment((36, 36), (4, 10))
    assert buffer.results[index2] == expected_result2, f"Expected {expected_result2}, but got {buffer.results[index2]}"

    # Test aligning two profiles
    index3 = buffer.add_alignment(index1, index2)
    expected_result3 = ProfileCostAlignment((16, 70), (4, 17))
    assert buffer.results[index3] == expected_result3, f"Expected {expected_result3}, but got {buffer.results[index3]}"

    # Test concatenation
    index4 = buffer.concatenate(["AGT", index3, "GCAT"])
    expected_result4 = ProfileCostAlignment((0, 0), (11, 24))
    assert buffer.results[index4] == expected_result4, f"Expected {expected_result4}, but got {buffer.results[index4]}"

    assert TotalCost(buffer).cost == (64, 118), f"Expected (64, 118), but got {TotalCost(buffer).cost}"

def test_dbg_to_cost_alignment():
    from dbg_align import DeBruijnGraph, AlignmentMethod
    from dbg_align.alignment_cost_plugin import PluginCostAlignment
    from dbg_align.allignment_buffer import AlignmentBuffer

    graph = DeBruijnGraph(3)
    graph.add_sequence({
        "seq1": "ACAGTACGGCAT", #12
        "seq2": "ACAGTACTGGCAT", #13
        "seq3": "ACAGCGCAT" #9
        })

    dag = PartialOrderGraph(graph)    
    assert dag.work(AlignmentMethod.PROGRESSIVE) == 9 * 12 + 12 * 13 
    assert dag.work(AlignmentMethod.BRAIDEDDEBRUIJGRAPH) == 2 * 3 + 7 * 3 




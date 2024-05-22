from dbg_align import DirectedAcyclicGraph, DAG_Node
from dbg_align.alignment_cost_plugin import PluginCostAlignment
from dbg_align.allignment_buffer import AlignmentBuffer
from dbg_align.constants import AlignmentMethod


def cost_alignment_plugin_factory() -> PluginCostAlignment:
    return PluginCostAlignment()

def test_construct_dag():
    dag = DirectedAcyclicGraph()
    dag.sequence_names = {'Sequence 1':(0,10),'Sequence 2':(1,10)}  # dict keyed on sequence names, returns tuple containing index and lengths of the sequence
    end = DAG_Node("GCAT", {0,1}) 
    dag.root = DAG_Node("AGT", {0, 1}) + [DAG_Node("GCG", {0})+end, DAG_Node("GTG",{1})+end]
    assert dag.root.fragment == "AGT"
    assert len(dag) == 2
    assert dag[0] == "AGTGCGGCAT"
    assert dag[1] == "AGTGTGGCAT"
    assert dag.expected_work(AlignmentMethod.EXACT) == 100
    assert dag.expected_work(AlignmentMethod.PROGRESSIVE) == 100
    # assert dag.expected_work(AlignmentMethod.DBG_LENGTH) == 9
    # assert dag.expected_work(AlignmentMethod.DBG_LENGTH_NUMBER) == 9
    bubbles = dag.bubbles()
    assert len(bubbles) == 1
    assert bubbles[0].start.fragment == "AGT"
    assert bubbles[0].end.fragment == "GCAT"
    assert bubbles[0].depth == 0
    assert bubbles[0].start.sequence_set == {0,1}
    assert bubbles[0].end.sequence_set == {0,1}
    assert len(bubbles[0].inner_bubbles) == 1
    cost_alignment = PluginCostAlignment()
    buffer = AlignmentBuffer(cost_alignment) 
    dag.align(buffer)
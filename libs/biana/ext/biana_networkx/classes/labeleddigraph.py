from digraph import DiGraph
from labeledgraph import LabeledGraph
from biana_networkx.exception import NetworkXException, NetworkXError
import biana_networkx.convert as convert

class LabeledDiGraph(LabeledGraph,DiGraph):
    pass  # just use the inherited classes

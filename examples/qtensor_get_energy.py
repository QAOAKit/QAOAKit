"""
QTensor is a fast tensor network simulator
See more: https://github.com/danlkv/QTensor
"""

import networkx as nx
from QAOAKit import opt_angles_for_graph, angles_to_qtensor_format
from qtensor import QAOA_energy
import numpy as np

# build graph
G = nx.star_graph(5)
# grab optimal angles
p = 3
angles = angles_to_qtensor_format(opt_angles_for_graph(G, p))
# use QTensor to find energy
E = QAOA_energy(G, angles["gamma"], angles["beta"])
print(E)

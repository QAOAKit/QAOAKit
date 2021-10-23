"""
Computing solutions from classical algorithms and comparing them to QAOA
"""

import networkx as nx
import numpy as np
from functools import partial
from qtensor import QAOA_energy

from QAOAKit import (
    opt_angles_for_graph,
    angles_to_qtensor_format,
)

from QAOAKit.utils import (
    maxcut_obj,
    get_adjacency_matrix,
)

from QAOAKit.classical import thompson_parekh_marwaha

from QAOAKit.qiskit_interface import goemans_williamson


import time

# get random graph
G = nx.random_regular_graph(3, 100, seed=42)
obj = partial(maxcut_obj, w=get_adjacency_matrix(G))

# solution 1: Goemans-Williamson
t0 = time.time()
soln = goemans_williamson(G, nsamples=100)
t1 = time.time()
print(f"GW finished in {t1-t0:.2f} sec")
print(f"GW cut {obj(soln[0])}")
scores = [obj(x) for x in soln]
print(f"GW mean cut {np.mean(scores):.2f}")
print(f"GW best cut {np.max(scores)}")

# solution 2: An explicit vector algorithm for high-girth MaxCut
# https://scirate.com/arxiv/2108.12477

t0 = time.time()
soln, exp_round = thompson_parekh_marwaha(G, nsamples=100, girth=10)
t1 = time.time()
print(f"TMP finished in {t1-t0:.2f} sec")
scores = [obj(x) for x in soln]

print(f"TPM best cut: {np.max(scores)}")
print(f"TPM mean cut {np.mean(scores):.2f}")
print(f"TPM expectation: {exp_round*G.number_of_edges():.2f}")

# solution 3: QAOA with p=3

# grab angles from https://scirate.com/arxiv/2107.00677
angles = angles_to_qtensor_format(opt_angles_for_graph(G, 3))
# use QTensor to find expected value of QAOA solution
E = QAOA_energy(G, angles["gamma"], angles["beta"])
print(f"Expected QAOA cut: {E[0]}")

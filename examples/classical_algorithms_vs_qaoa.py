"""
Computing solutions from classical algorithms and comparing them to QAOA
"""

import networkx as nx
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

from QAOAKit.classical import (
    goemans_williamson,
    thompson_parekh_marwaha,
)

# get random graph
G = nx.random_regular_graph(3, 100, seed=42)
obj = partial(maxcut_obj, w=get_adjacency_matrix(G))

# solution 1: Goemans-Williamson
sol_gw = goemans_williamson(G, nsamples=10)
print(f"GW cut {obj(sol_gw)}")

# solution 2: An explicit vector algorithm for high-girth MaxCut
# https://scirate.com/arxiv/2108.12477

soln, exp_round = thompson_parekh_marwaha(G, nsamples=100, girth=10)

# pick best of sample
best_score = float("-inf")
for x in soln:
    score = obj(x)
    if score > best_score:
        best_score = score
print(f"TPM cut: {best_score}")

# solution 3: QAOA with p=3

# grab angles from https://scirate.com/arxiv/2107.00677
angles = angles_to_qtensor_format(opt_angles_for_graph(G, 3))
# use QTensor to find expected value of QAOA solution
E = QAOA_energy(G, angles["gamma"], angles["beta"])
print(f"Expected QAOA cut: {E[0]}")

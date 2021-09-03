# Copyright 2021, Argonne National Laboratory

"""
Goemans-Williamson classical algorithm for MaxCut
"""

import networkx as nx
import numpy as np
import cvxgraphalgs as cvxgr
from copy import deepcopy


def goemans_williamson(G: nx.Graph, nsamples: int = 100) -> np.ndarray:
    best_score = float("-inf")
    best_cut = None
    for _ in range(nsamples):
        sdp_cut = cvxgr.algorithms.goemans_williamson_weighted(G)
        score = sdp_cut.evaluate_cut_size(G)
        if score > best_score:
            best_score = score
            best_cut = deepcopy(sdp_cut)
    x = np.zeros(G.number_of_nodes())
    for idx in best_cut.left:
        x[idx] = 1
    return x

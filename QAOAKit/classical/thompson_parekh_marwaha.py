import networkx as nx
import numpy as np
import scipy.sparse


def thompson_parekh_marwaha(G, nsamples=1, girth=0):
    """An explicit vector algorithm for high-girth MaxCut [1]
    Implementation courtesy of authors, refactored by Jonathan Wurtz

    input:
        G (nx.Graph):   graph on which to solve MaxCut
                          nodes must be labelled 0,..,|V|-1
        nsamples (int): number of samples to draw
        girth (int):    Assumed girth of the graph. >=2. Unknown behavior if
                         this value is set to be larger than the graph girth.
    returns:
        x (np.array):   binary string representing an approximate solution as +/-1
        exp (float):    Expected cut fraction over many samples
    [1] https://scirate.com/arxiv/2108.12477
    """

    # Determine some properties of the graph...
    connectivity = len(G[0])
    n = nx.number_of_nodes(G)
    m = nx.number_of_edges(G)

    # Translate variables
    if girth == 0:
        k = min(
            [len(q) for q in nx.algorithms.minimum_cycle_basis(G)]
        )  # This is expensive
    else:
        k = girth

    d = connectivity

    """
    Copied from email from J. Thompson 8/31/2021@11:15am "An explicit vector algorithm for high-girth MaxCut"
    Modifications are indistinguishable vectorizations of commented code
    """

    ## STEP 1: Create matrices representing optimization equation ##############

    a = 1.0 / np.sqrt(d)
    b = np.sqrt(d - 1) / d

    dim = k

    """
    M = np.zeros((dim, dim))
    M[0,1] = a
    M[1,0] = a

    for i in range(1, dim-1):
      M[i, i+1] = b
      M[i+1, i] = b


    ############################################################################


    ## STEP 2: Calculate eigenvalues + vectors #################################

    results = np.linalg.eigh(M)
    """
    # Replacement:
    mm = np.ones(dim - 1) * b
    mm[0] = a
    M = scipy.sparse.diags([mm] * 2, offsets=[1, -1]).tocsc()
    results = scipy.sparse.linalg.eigsh(M, k=1, which="SA")
    ###
    betas = results[1][:, 0]  # Or use min e-vector of B_k
    sigma1 = results[0][0]

    ## STEP 3 (optional): Testing on an actual graph ##########################

    ## Generate Vectors and Exp Matrix ###
    alphas = np.zeros(len(betas))
    alphas[0] = betas[0]
    """
    for i in range(1,len(betas)):
      alphas[i] = betas[i] / np.sqrt(d*(d-1.)**(i-1))
    """
    # Replacement:
    alphas[1::] = betas[1::] / np.sqrt(d * (d - 1.0) ** np.arange(len(alphas) - 1))
    ###

    V = np.zeros((n, n))
    """
    for i in range(n):
      for j in range(n):
        dist = nx.shortest_path_length(G, source=i, target=j)
        if dist < k:
          V[i,j] = alphas[dist]
    """
    # Replacement:
    for i, dist_ in nx.shortest_path_length(G):
        for j, dist in dist_.items():
            if dist < k:
                V[i, j] = alphas[dist]
    ###

    # Normalized covariance matrix
    V = V / np.linalg.norm(V, axis=1)
    W = V.T.dot(V)

    exp = 0.0
    exp_round = 0.0
    for e in G.edges:
        exp += 0.5 - 0.5 * W[e[0], e[1]]
        exp_round += np.arccos(W[e[0], e[1]])
    exp_round /= np.pi

    # Sample from a Gaussian with covariance V
    sample0 = np.random.normal(size=[n, nsamples])

    soln = (1 - np.sign(V.dot(sample0)).T) / 2
    return soln, exp_round / m


if __name__ == "__main__":
    from functools import partial
    from QAOAKit.utils import maxcut_obj, get_adjacency_matrix

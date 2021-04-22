import pickle
import pynauty
import networkx as nx
import numpy as np
from pathlib import Path


class LookupTableHandler:
    def __init__(self):
        self.graph2angles = None
        self.graph2pynauty = None
        self.folder = Path(__file__).parent

    def get_graph2angles(self):
        if self.graph2angles is None:
            self.graph2angles = pickle.load(open(Path(self.folder, f"../data/lookup_tables/graph2angles.p"),"rb"))
        return self.graph2angles

    def get_graph2pynauty(self):
        if self.graph2pynauty is None:
            self.graph2pynauty = pickle.load(open(Path(self.folder, f"../data/lookup_tables/graph2pynauty.p"),"rb"))
        return self.graph2pynauty

lookup_table_handler = LookupTableHandler()

def get_adjacency_dict(G):
    """Returns adjacency dictionary for G
    G must be a networkx graph
    Return format: { n : [n1,n2,...], ... }
    where [n1,n2,...] is a list of neighbors of n
    """
    adjacency_dict = {}
    for n, neigh_dict in G.adjacency():
        neigh_list = []
        for neigh, attr_dict in neigh_dict.items():
            assert(len(attr_dict) == 0)
            neigh_list.append(neigh)
        adjacency_dict[n] = neigh_list
    return adjacency_dict


def get_graph_id(G):
    graph2pynauty = lookup_table_handler.get_graph2pynauty()

    g = pynauty.Graph(number_of_vertices=G.number_of_nodes(), directed=nx.is_directed(G),
                adjacency_dict = get_adjacency_dict(G))
    cert = pynauty.certificate(g)

    return graph2pynauty[cert]


def opt_angles_for_graph(G, p):
    graph2angles = lookup_table_handler.get_graph2angles()
    graph_id = get_graph_id(G)
    return graph2angles[G.number_of_nodes()][p][graph_id]


def angles_to_qaoa_format(angles):
    """ Converts from format in graph2angles
    into the format used by qaoa.py
    get_maxcut_qaoa_circuit(G, angles['beta'], angles['gamma'])
    """
    angles['beta'] = np.pi*np.array(angles['beta'])
    angles['gamma'] = -np.pi*np.array(angles['gamma']) / 2
    return angles

def state_num2str(basis_state_as_num, nqubits):
    return '{0:b}'.format(basis_state_as_num).zfill(nqubits)

def state_str2num(basis_state_as_str):
    return int(basis_state_as_str, 2)

def state_reverse(basis_state_as_num, nqubits):
    basis_state_as_str = state_num2str(basis_state_as_num, nqubits)
    new_str = basis_state_as_str[::-1]
    return state_str2num(new_str)

def get_adjusted_state(state):
    nqubits = np.log2(state.shape[0])
    if nqubits % 1:
        raise ValueError("Input vector is not a valid statevector for qubits.")
    nqubits = int(nqubits)

    adjusted_state = np.zeros(2**nqubits, dtype=complex)
    for basis_state in range(2**nqubits):
         adjusted_state[state_reverse(basis_state, nqubits)] = state[basis_state]
    return adjusted_state

def state_to_ampl_counts(vec, eps=1e-15):
    """Converts a statevector to a dictionary
    of bitstrings and corresponding amplitudes
    """
    qubit_dims = np.log2(vec.shape[0])
    if qubit_dims % 1:
        raise ValueError("Input vector is not a valid statevector for qubits.")
    qubit_dims = int(qubit_dims)
    counts = {}
    str_format = '0{}b'.format(qubit_dims)
    for kk in range(vec.shape[0]):
        val = vec[kk]
        if val.real**2+val.imag**2 > eps:
            counts[format(kk, str_format)] = val
    return counts


def obj_from_statevector(sv, obj_f, precomputed=None):
    """Compute objective from Qiskit statevector
    For large number of qubits, this is slow. 
    To speed up for larger qubits, pass a vector of precomputed energies
    for QAOA, precomputed should be the same as the diagonal of the cost Hamiltonian
    """
    if precomputed is None:
        adj_sv = get_adjusted_state(sv)
        counts = state_to_ampl_counts(adj_sv)
        assert(np.isclose(sum(np.abs(v)**2 for v in counts.values()), 1))
        return sum(obj_f(np.array([int(x) for x in k])) * (np.abs(v)**2) for k, v in counts.items())
    else:
        return np.dot(precomputed, np.abs(sv)**2)

def maxcut_obj(x,G):
    """Compute the value of a cut.
    Args:
        x (numpy.ndarray): binary string as numpy array.
        w (numpy.ndarray): adjacency matrix.
    Returns:
        float: value of the cut.
    """
    cut = 0
    for e in G.edges():
        if x[e[0]] != x[e[1]]:
            cut += 1
    return cut

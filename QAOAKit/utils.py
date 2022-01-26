import pickle
import copy
import pynauty
import networkx as nx
import numpy as np
import pandas as pd
from pathlib import Path
from functools import partial
from qiskit.providers.aer import AerSimulator
import json
import re
import warnings

from QAOAKit.qaoa import get_maxcut_qaoa_circuit

utils_folder = Path(__file__).parent


class LookupTableHandler:
    """Singleton handling all the tables
    constructed from the qaoa-dataset-version1

    This object may consume a lot of memory!

    Attributes:
        graph2angles (dict): maps from n_qubits, p and graph_id to optimal parameters
                graph2angles[n_qubits][p][graph_id] = {'beta':optimal_beta, 'gamma':optimal_gamma}
        graph2pynauty (dict): maps from pynauty certificate to graph_id
        full_qaoa_dataset_table (pandas.DataFrame) : TBD
    """

    def __init__(self):
        self.graph2angles = None
        self.graph2pynauty = None
        # dictionary with mapping from nqubits to dictionary containing
        # 'graph_id2graph', 'graph_id2pynautycert', 'pynautycert2graph_id', 'pynautycert2graph' tables
        # example: large_graph_table[5]['graph_id2graph']
        self.large_graph_table = None
        self.full_qaoa_dataset_table = None
        self.three_reg_dataset_table = None
        self.fixed_angle_dataset_table = None
        self.full_weighted_qaoa_dataset_table = None

    def get_graph2angles(self):
        if self.graph2angles is None:
            self.graph2angles = pickle.load(
                open(Path(utils_folder, f"../data/lookup_tables/graph2angles.p"), "rb")
            )
        return self.graph2angles

    def get_graph2pynauty(self):
        if self.graph2pynauty is None:
            self.graph2pynauty = pickle.load(
                open(Path(utils_folder, f"../data/lookup_tables/graph2pynauty.p"), "rb")
            )
        return self.graph2pynauty

    def get_large_graph_table(self, nqubits):
        if self.large_graph_table is None:
            self.large_graph_table = {}
        if nqubits not in self.large_graph_table:
            self.large_graph_table[nqubits] = pickle.load(
                open(
                    Path(
                        utils_folder,
                        f"../data/lookup_tables/graph2pynauty_large_{nqubits}.p",
                    ),
                    "rb",
                )
            )
        return self.large_graph_table[nqubits]

    def get_full_qaoa_dataset_table(self):
        if self.full_qaoa_dataset_table is None:
            self.full_qaoa_dataset_table = pd.read_pickle(
                Path(utils_folder, "../data/lookup_tables/full_qaoa_dataset_table.p")
            ).set_index(["pynauty_cert", "p_max"])
        return self.full_qaoa_dataset_table

    def get_full_weighted_qaoa_dataset_table(self):
        if self.full_weighted_qaoa_dataset_table is None:
            self.full_weighted_qaoa_dataset_table = pd.read_json(
                Path(utils_folder, "../data/transfer_qaoa_weighted/all_transfer.zip")
            )
            self.full_weighted_qaoa_dataset_table[
                "G"
            ] = self.full_weighted_qaoa_dataset_table.apply(
                lambda row: nx.node_link_graph(row["G_json"]),
                axis=1,
            )
        return self.full_weighted_qaoa_dataset_table

    def get_3_reg_dataset_table(self):
        if self.three_reg_dataset_table is None:
            self.three_reg_dataset_table = pd.read_pickle(
                Path(utils_folder, "../data/lookup_tables/3_reg_dataset_table.p")
            ).set_index(["pynauty_cert", "p_max"])
        return self.three_reg_dataset_table

    def get_fixed_angle_dataset_table(self):
        # fixed angle dataset is small, no need to store pickle on disk
        if self.fixed_angle_dataset_table is None:
            with open(
                Path(
                    utils_folder,
                    "../data/fixed-angle-2021-08/angles_regular_graphs.json",
                )
            ) as json_file:
                data = json.load(json_file)

            lines = []
            for d in data.keys():
                for p in data[d]:
                    if int(p) < 12:  # remove bad value at p=12
                        gamma = [float(x) for x in data[d][p]["gamma"]]
                        beta = [float(x) for x in data[d][p]["beta"]]
                        line_d = {
                            "d": int(d),
                            "p": int(p),
                            "gamma": np.array(gamma) / np.pi,
                            "beta": np.array(beta) / np.pi,
                            "AR": float(data[d][p]["AR"]),
                        }
                        lines.append(line_d)
            self.fixed_angle_dataset_table = pd.DataFrame(
                lines, columns=lines[0].keys()
            )
        return self.fixed_angle_dataset_table


lookup_table_handler = LookupTableHandler()


def get_adjacency_dict(G):
    """Returns adjacency dictionary for G
    G must be a networkx graph
    Return format: { n : [n1,n2,...], ... }
    where [n1,n2,...] is a list of neighbors of n
    ignores all attributes
    """
    adjacency_dict = {}
    for n, neigh_dict in G.adjacency():
        neigh_list = []
        for neigh, attr_dict in neigh_dict.items():
            neigh_list.append(neigh)
        adjacency_dict[n] = neigh_list
    return adjacency_dict


def get_pynauty_certificate(G):
    """Get pynauty certificate for G

    Parameters
    ----------
    G : networkx.Graph
        Unweighted graph to compute certificate for

    Returns
    -------
    cert : binary
        isomorphism certificate for G
    """
    g = pynauty.Graph(
        number_of_vertices=G.number_of_nodes(),
        directed=nx.is_directed(G),
        adjacency_dict=get_adjacency_dict(G),
    )
    return pynauty.certificate(g)


def isomorphic(G1, G2):
    """Tests if two unweighted graphs are isomorphic using pynauty
    Ignores all attributes
    """
    g1 = pynauty.Graph(
        number_of_vertices=G1.number_of_nodes(),
        directed=nx.is_directed(G1),
        adjacency_dict=get_adjacency_dict(G1),
    )
    g2 = pynauty.Graph(
        number_of_vertices=G2.number_of_nodes(),
        directed=nx.is_directed(G2),
        adjacency_dict=get_adjacency_dict(G2),
    )
    return pynauty.isomorphic(g1, g2)


def get_graph_id(G):
    graph2pynauty = lookup_table_handler.get_graph2pynauty()

    g = pynauty.Graph(
        number_of_vertices=G.number_of_nodes(),
        directed=nx.is_directed(G),
        adjacency_dict=get_adjacency_dict(G),
    )
    cert = pynauty.certificate(g)

    return graph2pynauty[cert]


def get_graph_from_id(graph_id, nqubits):
    graph_id2graph = lookup_table_handler.get_large_graph_table(nqubits)[
        "graph_id2graph"
    ]
    return copy.deepcopy(graph_id2graph[graph_id])


def opt_angles_for_graph(G, p):
    if G.number_of_nodes() <= 9 and p <= 3:
        graph2angles = lookup_table_handler.get_graph2angles()
        graph_id = get_graph_id(G)
        return copy.deepcopy(graph2angles[G.number_of_nodes()][p][graph_id])
    elif nx.is_regular(G) and G.number_of_nodes() <= 16 and G.degree[0] == 3 and p <= 2:
        row = get_3_reg_dataset_table_row(G, p)
        return {"beta": row["beta"], "gamma": row["gamma"]}
    elif nx.is_regular(G) and p <= 11:
        d = G.degree[0]
        warnings.warn(
            f"Optimal angles not available, returning fixed angles for {d}-regular graphs"
        )
        return get_fixed_angles(d, p)
    elif p <= 11:
        warnings.warn("Optimal angles not available, returning closest fixed angles")
        d_ave = int(round(2 * G.number_of_edges() / G.number_of_nodes()))
        return get_fixed_angles(d_ave, p)
    else:
        raise NotImplementedError(
            f"Parameters for p > 11 are not available; requested p = {p}"
        )


def get_fixed_angles(d, p):
    row = get_fixed_angle_dataset_table_row(d, p)
    angles = {"beta": row.beta, "gamma": row.gamma}
    return angles


def get_full_qaoa_dataset_table():
    return lookup_table_handler.get_full_qaoa_dataset_table()


def get_full_weighted_qaoa_dataset_table():
    return lookup_table_handler.get_full_weighted_qaoa_dataset_table()


def get_3_reg_dataset_table():
    return lookup_table_handler.get_3_reg_dataset_table()


def get_fixed_angle_dataset_table():
    return lookup_table_handler.get_fixed_angle_dataset_table()


def get_fixed_angle_dataset_table_row(d, p):
    df = lookup_table_handler.get_fixed_angle_dataset_table()
    row = df[(df["d"] == d) & (df["p"] == p)]
    assert len(row) == 1
    row = row.squeeze()
    assert isinstance(row, pd.Series)
    return row


def get_full_qaoa_dataset_table_row(G, p):
    """Returns full table row for a given NetworkX graph"""
    full_qaoa_dataset_table = lookup_table_handler.get_full_qaoa_dataset_table()

    g = pynauty.Graph(
        number_of_vertices=G.number_of_nodes(),
        directed=nx.is_directed(G),
        adjacency_dict=get_adjacency_dict(G),
    )
    cert = pynauty.certificate(g)

    return full_qaoa_dataset_table.loc[(cert, p)]


def get_3_reg_dataset_table_row(G, p):
    """Returns full table row for a given NetworkX graph"""
    df = lookup_table_handler.get_3_reg_dataset_table()

    g = pynauty.Graph(
        number_of_vertices=G.number_of_nodes(),
        directed=nx.is_directed(G),
        adjacency_dict=get_adjacency_dict(G),
    )
    cert = pynauty.certificate(g)

    return df.loc[(cert, p)]


def angles_to_qaoa_format(angles):
    """Converts from format in graph2angles
    into the format used by qaoa.py
    get_maxcut_qaoa_circuit(G, angles['beta'], angles['gamma'])
    """
    res = copy.deepcopy(angles)
    res["beta"] = beta_to_qaoa_format(res["beta"])
    res["gamma"] = gamma_to_qaoa_format(res["gamma"])
    return res


def beta_to_qaoa_format(beta):
    """Converts from format in graph2angles
    into the format used by qaoa.py
    """
    return np.pi * np.array(beta)


def gamma_to_qaoa_format(gamma):
    """Converts from format in graph2angles
    into the format used by qaoa.py
    get_maxcut_qaoa_circuit(G, angles['beta'], angles['gamma'])
    """
    return -np.pi * np.array(gamma) / 2


def angles_to_qiskit_format(angles):
    """Converts from format in graph2angles
    into the format used by QAOAAnsatz
    """
    return np.concatenate(
        [[-np.pi * g, np.pi * b] for g, b in zip(angles["gamma"], angles["beta"])]
    )


def angles_from_qiskit_format(angles):
    """Converts from the format used by QAOAAnsatz
    into the format in graph2angles
    """
    res = {}
    assert len(angles) % 2 == 0
    res["gamma"] = list(x / (-np.pi) for x in angles[::2])
    res["beta"] = list(x / np.pi for x in angles[1::2])
    return res


def angles_to_qtensor_format(angles):
    """Converts from format in graph2angles
    into the format used by QTensor
    """
    return {"gamma": [-g / 2 for g in angles["gamma"]], "beta": angles["beta"]}


def read_graph_from_file(f, expected_nnodes=None):
    """Read a graph in a format used by qaoa-dataset-version1

    Parameters
    ----------
    f : file-like object
        Handler for the file to read from
    expected_nnodes : int, default None
        Number of nodes to expect
        If passed, a check will be performed
        to confirm that the actual number of nodes
        matches the expectation

    Returns
    -------
    G : networkx.Graph
    graph_id : int
        ID of the graph
    """
    f.readline(-1)  # first line is blank
    line_with_id = f.readline(-1)  # second line has graph number and order
    graph_id, graph_order = [
        int(x) for x in re.split(" |, |. |.\n", line_with_id) if x.isdigit()
    ]
    if expected_nnodes is not None:
        assert graph_order == expected_nnodes
    G = nx.Graph()
    for n in range(graph_order):
        G.add_nodes_from([n])
    edge_id = 0
    # third line is first row of upper triangle of adjacency matrix (without the diagonal element)
    for n in range(graph_order - 1):
        adj_str = f.readline(-1)
        for m in range(graph_order - 1 - n):
            q_num = n + m + 1
            if adj_str[m] == "1":
                G.add_edge(n, q_num, edge_id=edge_id)
                edge_id += 1
    return G, graph_id


def load_results_file_into_dataframe(n_qubits, p):
    """Loads one file from ../data/qaoa-dataset-version1/Results/ into a pandas.DataFrame
    Column names are from ../data/qaoa-dataset-version1/Results/How_to_read_data_columns.txt
    Columns added:
    p_max : maximal p allowed; this is to differentiate from p in the original dataset, which can be lower due to achieving optimal solution
    beta : concatenated beta_i
    gamma : concatenated gamma_i
    """
    colnames = ["graph_id", "C_{true opt}", "C_init", "C_opt", "pr(max)", "p"]
    for i in range(p):
        colnames.append(f"beta_{i}/pi")
    for i in range(p):
        colnames.append(f"gamma_{i}/pi")
    df = pd.read_csv(
        Path(
            utils_folder,
            f"../data/qaoa-dataset-version1/Results/p={p}/n={n_qubits}_p={p}.txt",
        ),
        delim_whitespace=True,
        names=colnames,
        header=None,
    )
    df[
        "p_max"
    ] = p  # maximal p allowed; this is to differentiate from p in the original dataset, which can be lower due to achieving optimal solution
    df["beta"] = df.apply(
        lambda row: np.array([row[f"beta_{i}/pi"] for i in range(p)]), axis=1
    )
    df["gamma"] = df.apply(
        lambda row: np.array([row[f"gamma_{i}/pi"] for i in range(p)]), axis=1
    )
    assert (df["p_max"] >= df["p"]).all()
    return df[
        [
            "graph_id",
            "C_{true opt}",
            "C_init",
            "C_opt",
            "pr(max)",
            "p",
            "p_max",
            "beta",
            "gamma",
        ]
    ].set_index("graph_id")


def get_graph_and_assign_weights(graph_id, weight_id, nqubits, df_weights, graphs_dict):
    """Retrieves a graph from graph_id and nqubits and assigns weights
    from df_weights
    """
    weights = list(
        df_weights[
            (df_weights["graph_id"] == graph_id)
            & (df_weights["weight_id"] == weight_id)
        ]["weights"]
    )
    if len(weights) != 1:
        raise ValueError(
            f"For graph_id={graph_id}, weight_id={weight_id} found unexpected number ({len(weights)}) weights"
        )
    weights = weights[0]
    if graphs_dict is None:
        assert (
            nqubits <= 9
        ), "If the number of nodes is greater than 9, must pass graphs file path"
        G = get_graph_from_id(graph_id, nqubits)
    else:
        G = graphs_dict[graph_id]
    for u, v, attr_dict in G.edges(data=True):
        G[u][v]["weight"] = weights[attr_dict["edge_id"]]
    return copy.deepcopy(G)


def load_weighted_results_into_dataframe(
    folder_path, p, nqubits, df_weights, graphs_file_path=None
):
    """Loads all result files from ../data/weighted_angle_dat/{p}
    if no path to file describing graphs (graphs_file_path) is provided, will presume that
    the graphs should be loaded from the full_qaoa_dataset_table
    df_weights is a dataframe mapping graph_id and weight_id to list of weights
    The column names and conventions are described in ../data/weighted_angle_dat/Readme.txt
    One column is added:
    p_max : maximal p allowed; this is to differentiate from p in the original dataset, which can be lower due to achieving optimal solution
    """
    colnames = [
        "graph_id",
        "weight_id",
        "C_{true opt}",
        "C_init",
        "C_opt",
        "pr(max)",
        "p",
    ]
    for i in range(p):
        colnames.append(f"beta_{i}/pi")
    for i in range(p):
        colnames.append(f"gamma_{i}/pi")
    colnames.append("mean(weight)")
    colnames.append("std(weight)")
    dfs = []
    for fname in folder_path.glob("QAOA_dat_weighted_*"):
        dfs.append(
            pd.read_csv(
                fname,
                delim_whitespace=True,
                usecols=list(range(len(colnames))),
                names=colnames,
                header=None,
            )
        )
    df = pd.concat(dfs)
    df[
        "p_max"
    ] = p  # maximal p allowed; this is to differentiate from p in the original dataset, which can be lower due to achieving optimal solution
    df["beta"] = df.apply(lambda row: [row[f"beta_{i}/pi"] for i in range(p)], axis=1)
    df["gamma"] = df.apply(lambda row: [row[f"gamma_{i}/pi"] for i in range(p)], axis=1)
    # load graphs if needed
    if graphs_file_path is not None:
        graphs_dict = {}
        with open(graphs_file_path, "r") as f:
            graph_ids = []
            while True:
                try:
                    G, graph_id = read_graph_from_file(f)
                    graphs_dict[graph_id] = G
                except ValueError:
                    break
    else:
        graphs_dict = None
    df["G"] = df.apply(
        lambda row: get_graph_and_assign_weights(
            row["graph_id"], row["weight_id"], nqubits, df_weights, graphs_dict
        ),
        axis=1,
    )
    assert np.all(
        np.isclose(
            df.apply(
                lambda row: np.mean(
                    [x[2]["weight"] for x in row["G"].edges(data=True)]
                ),
                axis=1,
            ),
            df["mean(weight)"],
            atol=1e-07,
        )
        | np.isnan(df["std(weight)"])
    )

    assert np.all(
        np.isclose(
            df.apply(
                lambda row: np.std([x[2]["weight"] for x in row["G"].edges(data=True)]),
                axis=1,
            ),
            df["std(weight)"],
            atol=1e-07,
        )
        | np.isnan(df["std(weight)"])
    )

    if nqubits <= 10:
        number_of_rows_to_check = 20
    elif nqubits <= 15:
        number_of_rows_to_check = 10
    else:
        number_of_rows_to_check = 1

    assert np.all(
        np.isclose(
            df.head(number_of_rows_to_check).apply(
                lambda row: qaoa_maxcut_energy(
                    row["G"],
                    beta_to_qaoa_format(row["beta"]),
                    gamma_to_qaoa_format(row["gamma"]),
                ),
                axis=1,
            ),
            df.head(number_of_rows_to_check)["C_opt"],
        )
    )

    return df


def load_weights_into_dataframe(folder_path):
    """Loads all weight files from the folder passed
    the first column is the graph number
    The next E columns are the weights of the E edges.
    Other columns beyond E+1 are meaningless
    These weights were drawn uniformly at random from 1,2
    returns pandas.DataFrame with columns
       -'weight_id' (int)
       -'graph_id' (int)
       -'weights' (list of floating point weights)
    """
    lines = []

    for fname in folder_path.glob("weights_*"):
        with open(fname, "r") as f:
            for line_num, line in enumerate(f.readlines()):
                line = line.strip().split()
                line_d = {
                    "weight_id": line_num + 1,
                    "graph_id": int(line[0]),
                    "weights": [float(x) for x in line[1:]],
                }
                lines.append(line_d)
    return pd.DataFrame(lines, columns=lines[0].keys())


def brute_force(obj_f, num_variables, minimize=False):
    """Get the maximum of a function by complete enumeration
    Returns the maximum value and the extremizing bit string
    """
    if minimize:
        best_cost_brute = float("inf")
        compare = lambda x, y: x < y
    else:
        best_cost_brute = float("-inf")
        compare = lambda x, y: x > y
    bit_strings = (
        (
            (
                np.array(range(2 ** num_variables))[:, None]
                & (1 << np.arange(num_variables))
            )
        )
        > 0
    ).astype(int)
    for x in bit_strings:
        cost = obj_f(np.array(x))
        if compare(cost, best_cost_brute):
            best_cost_brute = cost
            xbest_brute = x
    return best_cost_brute, xbest_brute


#############################
# QAOA utils
############################


def state_num2str(basis_state_as_num, nqubits):
    return "{0:b}".format(basis_state_as_num).zfill(nqubits)


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

    adjusted_state = np.zeros(2 ** nqubits, dtype=complex)
    for basis_state in range(2 ** nqubits):
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
    str_format = "0{}b".format(qubit_dims)
    for kk in range(vec.shape[0]):
        val = vec[kk]
        if val.real ** 2 + val.imag ** 2 > eps:
            counts[format(kk, str_format)] = val
    return counts


def precompute_energies(obj_f, nbits):
    """
    Precomputed a vector of objective function values
    that accelerates the energy computation in obj_from_statevector
    """
    bit_strings = (
        ((np.array(range(2 ** nbits))[:, None] & (1 << np.arange(nbits)))) > 0
    ).astype(int)

    return np.array([obj_f(x) for x in bit_strings])


def obj_from_statevector(sv, obj_f, precomputed_energies=None):
    """Compute objective from Qiskit statevector
    For large number of qubits, this is slow.
    """
    if precomputed_energies is None:
        qubit_dims = np.log2(sv.shape[0])
        if qubit_dims % 1:
            raise ValueError("Input vector is not a valid statevector for qubits.")
        qubit_dims = int(qubit_dims)
        # get bit strings for each element of the state vector
        # https://stackoverflow.com/questions/22227595/convert-integer-to-binary-array-with-suitable-padding
        bit_strings = (
            ((np.array(range(sv.shape[0]))[:, None] & (1 << np.arange(qubit_dims)))) > 0
        ).astype(int)

        return sum(
            obj_f(bit_strings[kk]) * (np.abs(sv[kk]) ** 2) for kk in range(sv.shape[0])
        )
    else:
        amplitudes = np.array([np.abs(sv[kk]) ** 2 for kk in range(sv.shape[0])])
        return precomputed_energies.dot(amplitudes)


def maxcut_obj(x, w):
    """Compute the value of a cut.

    Args:
        x (numpy.ndarray): binary string as numpy array.
        w (numpy.ndarray): adjacency matrix.

    Returns:
        float: value of the cut.
    """
    X = np.outer(x, (1 - x))
    return np.sum(w * X)


def get_adjacency_matrix(G):
    n = G.number_of_nodes()
    w = np.zeros([n, n])

    for e in G.edges():
        if nx.is_weighted(G):
            w[e[0], e[1]] = G[e[0]][e[1]]["weight"]
            w[e[1], e[0]] = G[e[0]][e[1]]["weight"]
        else:
            w[e[0], e[1]] = 1
            w[e[1], e[0]] = 1
    return w


def qaoa_maxcut_energy(G, beta, gamma, precomputed_energies=None):
    """Computes MaxCut QAOA energy for graph G
    qaoa format (`angles_to_qaoa_format`) used for beta, gamma
    """
    if precomputed_energies is None:
        obj = partial(maxcut_obj, w=get_adjacency_matrix(G))
    else:
        obj = None
    qc = get_maxcut_qaoa_circuit(G, beta, gamma)
    backend = AerSimulator(method="statevector")
    sv = backend.run(qc).result().get_statevector()
    return obj_from_statevector(sv, obj, precomputed_energies=precomputed_energies)

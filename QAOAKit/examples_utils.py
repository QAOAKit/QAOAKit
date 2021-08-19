import networkx as nx
import copy
from pathlib import Path
import pandas as pd
import numpy as np

from .utils import (
    beta_to_qaoa_format,
    gamma_to_qaoa_format,
    qaoa_maxcut_energy,
)

example_utils_folder = Path(__file__).parent


def get_20_node_erdos_renyi_graphs():
    # First, load the graphs

    n_graphs = 10
    n_qubits = 20

    graphs = {}

    with open(
        Path(
            example_utils_folder,
            "../data/optimal_parameters_n_20/Erdos_Renyi_n_20_density_9pt5_ngraphs_10.txt",
        )
    ) as f:
        for _ in range(1, n_graphs + 1):
            f.readline(-1)  # first line is blank
            line_with_id = f.readline(-1)  # second line has graph number and order
            graph_id = int(line_with_id.strip().split()[1])
            edges = []
            G = nx.Graph()
            for n in range(n_qubits):
                G.add_nodes_from([n])
            # third line is first row of upper triangle of adjacency matrix (without the diagonal element)
            for n in range(n_qubits - 1):
                adj_str = f.readline(-1)
                for m in range(n_qubits - 1 - n):
                    q_num = n + m + 1
                    if adj_str[m] == "1":
                        edges.append([n, q_num])
                        G.add_edge(n, q_num)
            graphs[graph_id] = copy.deepcopy(G)

    graph_id2graph = pd.Series(graphs, name="G")

    # Second, load optimal parameters

    def load_results_file_into_dataframe(p):
        colnames = ["graph_id", "C_{true opt}", "C_init", "C_opt", "pr(max)", "p"]
        for i in range(p):
            colnames.append(f"beta_{i}/pi")
        for i in range(p):
            colnames.append(f"gamma_{i}/pi")
        colnames.append("Dummy 1")
        colnames.append("Dummy 2")
        colnames.append("Dummy 3")

        df = pd.read_csv(
            Path(
                example_utils_folder,
                f"../data/optimal_parameters_n_20/QAOA_dat_p={p}",
            ),
            delim_whitespace=True,
            names=colnames,
            header=None,
            index_col="graph_id",
        )
        df[
            "p_max"
        ] = p  # maximal p allowed; this is to differentiate from p in the original dataset, which can be lower due to achieving optimal solution
        assert (df["p_max"] >= df["p"]).all()
        return df

    df = pd.DataFrame()

    for p in range(1, 4):
        df_orig = load_results_file_into_dataframe(p).merge(
            graph_id2graph, how="outer", right_index=True, left_index=True
        )
        assert len(df_orig) == n_graphs
        df_orig["beta"] = df_orig.apply(
            lambda row: [row[f"beta_{i}/pi"] for i in range(p)], axis=1
        )
        df_orig["gamma"] = df_orig.apply(
            lambda row: [row[f"gamma_{i}/pi"] for i in range(p)], axis=1
        )
        df_orig["theta"] = df_orig.apply(lambda row: row["gamma"] + row["beta"], axis=1)
        df_orig["n"] = df_orig.apply(lambda row: row["G"].number_of_nodes(), axis=1)
        assert (df_orig["n"] == n_qubits).all()
        df = df.append(df_orig.reset_index())
    return df.rename(
        columns={
            "C_opt": "QAOA energy with optimized parameters",
            "C_{true opt}": "Optimal MaxCut solution from brute force",
        }
    )

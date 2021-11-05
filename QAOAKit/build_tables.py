import numpy as np
import pandas as pd
import networkx as nx
import json
import gc
import copy
import pynauty
import pickle
from pathlib import Path
from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile
import shutil
from tqdm import tqdm
from multiprocessing import Process

from .utils import (
    load_results_file_into_dataframe,
    get_adjacency_dict,
    read_graph_from_file,
)

build_tables_folder = Path(__file__).parent

n_graphs = {3: 2, 4: 6, 5: 21, 6: 112, 7: 853, 8: 11117, 9: 261080}


def load_data():
    if not Path(build_tables_folder, "../data/qaoa-dataset-version1/").is_dir():
        print("Loading data, this may take a while")
        data_folder = Path(build_tables_folder, "../data/")
        zipurl = "https://github.com/QAOAKit/data/zipball/master"
        with urlopen(zipurl) as zipresp:
            with ZipFile(BytesIO(zipresp.read())) as zfile:
                zfile.extractall(data_folder)
        # remove top folder
        folder = list(Path(build_tables_folder, "../data/").glob("QAOAKit-data-*"))
        assert len(folder) == 1
        folder = folder[0]
        for subfolder in folder.glob("*"):
            if subfolder.stem == ".gitignore":
                continue
            subfolder = subfolder.absolute()
            parent_dir = subfolder.parents[1]
            subfolder.rename(parent_dir / subfolder.name)
        shutil.rmtree(folder)


def build_graph2angles():
    print(f"Building graph2angles table...")

    tables = {}

    for n_qubits in tqdm(range(3, 10)):
        tables[n_qubits] = {}

        for p in tqdm(range(1, 4), leave=False):
            tables[n_qubits][p] = {}

            df = load_results_file_into_dataframe(n_qubits, p)

            for index, data in df.iterrows():
                if data["p"] != p:
                    assert np.isclose(data["C_{true opt}"], data["C_opt"])
                tables[n_qubits][p][int(index)] = {
                    "beta": data["beta"],
                    "gamma": data["gamma"],
                }

    pickle.dump(
        tables,
        open(Path(build_tables_folder, f"../data/lookup_tables/graph2angles.p"), "wb"),
    )
    del tables
    print(f"Done building graph2angles table.")


def build_graph2pynauty():
    print(f"Building graph2pynauty table...")

    tables = {}

    for n_qubits in tqdm(range(3, 10)):

        with open(
            Path(
                build_tables_folder,
                "../data/qaoa-dataset-version1/Graphs/graph" + str(n_qubits) + "c.txt",
            )
        ) as f:
            for _ in tqdm(range(n_graphs[n_qubits]), leave=False):
                G, graph_id = read_graph_from_file(f, expected_nnodes=n_qubits)
                g = pynauty.Graph(
                    number_of_vertices=G.number_of_nodes(),
                    directed=nx.is_directed(G),
                    adjacency_dict=get_adjacency_dict(G),
                )
                cert = pynauty.certificate(g)
                tables[cert] = graph_id

    pickle.dump(
        tables,
        open(Path(build_tables_folder, f"../data/lookup_tables/graph2pynauty.p"), "wb"),
    )
    del tables
    print(f"Done building graph2pynauty table.")


def build_graph2pynauty_large():
    print(f"Building graph2pynauty_large tables...")

    for n_qubits in tqdm(range(3, 10)):
        table = {}

        table["graph_id2pynautycert"] = {}
        table["graph_id2graph"] = {}
        table["pynautycert2graph_id"] = {}
        table["pynautycert2graph"] = {}

        with open(
            Path(
                build_tables_folder,
                "../data/qaoa-dataset-version1/Graphs/graph" + str(n_qubits) + "c.txt",
            )
        ) as f:
            for _ in tqdm(range(n_graphs[n_qubits]), leave=False):
                G, graph_id = read_graph_from_file(f, expected_nnodes=n_qubits)
                g = pynauty.Graph(
                    number_of_vertices=G.number_of_nodes(),
                    directed=nx.is_directed(G),
                    adjacency_dict=get_adjacency_dict(G),
                )
                cert = pynauty.certificate(g)
                table["graph_id2graph"][graph_id] = copy.deepcopy(G)
                table["graph_id2pynautycert"][graph_id] = cert
                table["pynautycert2graph_id"][cert] = graph_id
                table["pynautycert2graph"][cert] = copy.deepcopy(G)

        assert len(table["graph_id2pynautycert"]) == n_graphs[n_qubits]
        assert len(table["graph_id2graph"]) == n_graphs[n_qubits]
        assert len(table["pynautycert2graph_id"]) == n_graphs[n_qubits]
        assert len(table["pynautycert2graph"]) == n_graphs[n_qubits]
        pickle.dump(
            table,
            open(
                Path(
                    build_tables_folder,
                    f"../data/lookup_tables/graph2pynauty_large_{n_qubits}.p",
                ),
                "wb",
            ),
        )
        del table
        gc.collect()
    print(f"Done building graph2pynauty_large tables.")


def build_full_qaoa_dataset():
    print(f"Building full_qaoa_dataset table...")
    """
    Specifications:
        index: 'pynauty_cert'+p (??)
        columns: [
            'pynauty_cert','graph_id','#nodes','C_{true opt}','C_init','C_opt','pr(max)','p','beta','gamma',
            'p_max', # maximal p allowed; this is to differentiate from p in the original dataset, which can be lower due to achieving optimal solution
            ]
    """
    df = pd.DataFrame()

    for n_qubits in tqdm(range(3, 10)):
        graph_table = pickle.load(
            open(
                Path(
                    build_tables_folder,
                    f"../data/lookup_tables/graph2pynauty_large_{n_qubits}.p",
                ),
                "rb",
            )
        )
        graph_id2pynauty = pd.Series(
            graph_table["graph_id2pynautycert"], name="pynauty_cert"
        )
        graph_id2graph = pd.Series(graph_table["graph_id2graph"], name="G")
        for p in tqdm(range(1, 4), leave=False):
            df_orig = (
                load_results_file_into_dataframe(n_qubits, p)
                .merge(graph_id2pynauty, how="outer", right_index=True, left_index=True)
                .merge(graph_id2graph, how="outer", right_index=True, left_index=True)
            )
            assert len(df_orig) == n_graphs[n_qubits]
            df_orig["n"] = df_orig.apply(lambda row: row["G"].number_of_nodes(), axis=1)
            assert (df_orig["n"] == n_qubits).all()
            df = df.append(df_orig.reset_index())
    assert len(df) == sum(v * 3 for v in n_graphs.values())
    df.to_pickle(
        Path(build_tables_folder, "../data/lookup_tables/full_qaoa_dataset_table.p")
    )
    del df
    print(f"Done building full_qaoa_dataset table.")


def build_3_reg_dataset():
    print("Building  3-regular graph table")
    with open(
        Path(build_tables_folder, "../data/3_regular/3r_WURTZ_ensemble.json")
    ) as json_file:
        data = json.load(json_file)

    rows = []
    for row in tqdm(data):
        G = nx.Graph()
        G.add_edges_from(row["edges"])
        g = pynauty.Graph(
            number_of_vertices=G.number_of_nodes(),
            directed=nx.is_directed(G),
            adjacency_dict=get_adjacency_dict(G),
        )
        cert = pynauty.certificate(g)
        c_true_opt = row["0"]["MaxCut"]

        for p in range(1, 11):
            d = {}
            d["G"] = copy.deepcopy(G)
            d["n"] = G.number_of_nodes()
            d["p_max"] = p
            d["pynauty_cert"] = cert
            d["C_{true opt}"] = c_true_opt
            d["C_fixed"] = row[str(p)]["fixed_val"]
            if "optimized_val" in row[str(p)]:
                d["C_opt"] = row[str(p)]["optimized_val"]
                d["beta"] = np.array(row[str(p)]["angles"][0]["beta"]) / np.pi
                d["gamma"] = -2 * np.array(row[str(p)]["angles"][0]["gamma"]) / np.pi
                d["all beta (degenerate optima)"] = [
                    np.array(x["beta"]) / np.pi for x in row[str(p)]["angles"]
                ]
                d["all gamma (degenerate optima)"] = [
                    -2 * np.array(x["gamma"]) / np.pi for x in row[str(p)]["angles"]
                ]
                d["theta"] = np.hstack([d["gamma"], d["beta"]])
            rows.append(copy.deepcopy(d))

    df = pd.DataFrame(rows, columns=rows[0].keys())
    df.to_pickle(
        Path(build_tables_folder, "../data/lookup_tables/3_reg_dataset_table.p")
    )
    print("Done building  3-regular graph table")


if __name__ == "__main__":
    print("Building tables (this may take a few minutes)...")

    # each function is run in its own process to force system to free memory

    p = Process(target=load_data)
    p.start()
    p.join()
    del p
    gc.collect()

    p1 = Process(target=build_graph2angles)
    p1.start()

    p2 = Process(target=build_graph2pynauty)
    p2.start()

    p1.join()
    p2.join()
    del p1
    del p2
    gc.collect()

    p = Process(target=build_graph2pynauty_large)
    p.start()
    p.join()
    del p
    gc.collect()

    p = Process(target=build_full_qaoa_dataset)
    p.start()
    p.join()
    del p
    gc.collect()

    p = Process(target=build_3_reg_dataset)
    p.start()
    p.join()
    del p
    gc.collect()

    print("All done")

import numpy as np
import pandas as pd
import copy
import pickle
from pathlib import Path

from utils import load_results_file_into_dataframe

build_table_graph2angles_folder = Path(__file__).parent

tables = {}

for n_qubits in range(3, 10):
    tables[n_qubits] = {}

    for p in range(1, 4):
        tables[n_qubits][p] = {}

        df = load_results_file_into_dataframe(n_qubits, p)

        for index, data in df.iterrows():
            if data["p"] != p:
                assert np.isclose(data["C_{true opt}"], data["C_opt"])
            beta = [data[f"beta_{i}/pi"] for i in range(p)]
            gamma = [data[f"gamma_{i}/pi"] for i in range(p)]
            tables[n_qubits][p][int(index)] = {
                "beta": copy.deepcopy(beta),
                "gamma": copy.deepcopy(gamma),
            }
        print(f"Done with n={n_qubits}, p={p}")

pickle.dump(
    tables,
    open(
        Path(build_table_graph2angles_folder, f"../data/lookup_tables/graph2angles.p"),
        "wb",
    ),
)

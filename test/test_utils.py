import networkx as nx
import numpy as np
import pandas as pd
from qiskit import execute, Aer
from functools import partial
import pickle

from QAOA_parameters import opt_angles_for_graph, get_graph_id, angles_to_qaoa_format, get_full_qaoa_dataset_table_row, get_full_qaoa_dataset_table
from QAOA_parameters.utils import obj_from_statevector, maxcut_obj
from QAOA_parameters.qaoa import get_maxcut_qaoa_circuit

def test_retrieval():
    for n in range(3,10):
        G = nx.complete_graph(n)
        for p in range(1,4):
            angles = opt_angles_for_graph(G,p)
            assert(len(angles['beta']) == p)
            assert(len(angles['gamma']) == p)
            assert(isinstance(get_full_qaoa_dataset_table_row(G,p)[['C_opt','graph_id','C_{true opt}']], pd.Series))

def test_angle_conversion():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3,5]:
        for p in [2,3]:
            df = full_qaoa_dataset_table.reset_index()
            df = df[(df['n'] == n_qubits) & (df['p_max'] == p)]
            for _, row in df.iterrows():
                obj = partial(maxcut_obj, G=row['G'])
                opt_cut = row['C_opt']
                angles = angles_to_qaoa_format(opt_angles_for_graph(row['G'], row['p_max']))
                qc = get_maxcut_qaoa_circuit(row['G'], angles['beta'], angles['gamma'])
                backend = Aer.get_backend('statevector_simulator')
                res = backend.run(qc).result()
                sv = res.get_statevector()
                obj_val = obj_from_statevector(sv, obj)
                assert(np.isclose(opt_cut, obj_val))

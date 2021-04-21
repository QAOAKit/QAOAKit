import networkx as nx
import numpy as np
import pandas as pd
from qiskit import execute, Aer
from functools import partial
import pickle

from utils import opt_angles_for_graph, obj_from_statevector, maxcut_obj, get_graph_id, angles_to_qaoa_format
from qaoa import get_maxcut_qaoa_circuit

graph_table = pickle.load(open(f"../data/lookup_tables/graph2pynauty_large.p", "rb"))

for n_qubits in range(3,10):
    for p in range(1,4):

        colnames=['graph_id','C_{true opt}','C_init','C_opt','pr(max)','p']
        for i in range(p):
            colnames.append(f"beta_{i}/pi")
        for i in range(p):
            colnames.append(f"gamma_{i}/pi")
        df = pd.read_csv(f"../data/qaoa-dataset-version1/Results/p={p}/n={n_qubits}_p={p}.txt", delim_whitespace=True, names=colnames, header=None, index_col='graph_id')

        for graph_id, G in graph_table[n_qubits]['graph_id2graph'].items():
            obj = partial(maxcut_obj, G=G)
            
            opt_cut = df.loc[graph_id]['C_opt']
            
            angles = angles_to_qaoa_format(opt_angles_for_graph(G, p))
            
            qc = get_maxcut_qaoa_circuit(G, angles['beta'], angles['gamma'])
            
            backend = Aer.get_backend('statevector_simulator')
            
            res = backend.run(qc).result()
            sv = res.get_statevector()
            
            obj_val = obj_from_statevector(sv, obj)
            try:
                assert(np.isclose(opt_cut, obj_val))
            except AssertionError as e:
                print(graph_id, n_qubits, opt_cut, obj_val)
                raise e
        print(f"Done with p={p}, n_qubits={n_qubits}")

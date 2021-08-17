import pandas as pd
import pickle
from pathlib import Path
import json
import pynauty
import copy
import networkx as nx
import numpy as np

from .utils import load_results_file_into_dataframe, get_adjacency_dict

build_full_qaoa_dataset_table_folder = Path(__file__).parent

n_graphs={
        3: 2,
        4: 6,
        5: 21,
        6: 112,
        7: 853,
        8: 11117,
        9: 261080,
}

"""
Specifications:
    index: 'pynauty_cert'+p (??)
    columns: [
        'pynauty_cert','graph_id','#nodes','C_{true opt}','C_init','C_opt','pr(max)','p','beta','gamma',
        'theta', # concatenated gamma,beta: theta[idx] = gamma_idx, theta[p+idx] = beta_idx
        'p_max', # maximal p allowed; this is to differentiate from p in the original dataset, which can be lower due to achieving optimal solution
        ]
"""

df = pd.DataFrame()

for n_qubits in range(3,10):
    print('n_qubits:', n_qubits)
    graph_table = pickle.load(open(Path(build_full_qaoa_dataset_table_folder, f"../data/lookup_tables/graph2pynauty_large_{n_qubits}.p"), "rb"))
    graph_id2pynauty = pd.Series(graph_table['graph_id2pynautycert'], name='pynauty_cert')
    graph_id2graph = pd.Series(graph_table['graph_id2graph'], name='G')
    for p in range(1,4):
        df_orig = (
                load_results_file_into_dataframe(n_qubits,p).merge(graph_id2pynauty, how='outer', right_index=True, left_index=True)
                                                            .merge(graph_id2graph, how='outer', right_index=True, left_index=True)
                  )
        assert(len(df_orig) == n_graphs[n_qubits])
        df_orig['beta'] = df_orig.apply(lambda row: [row[f"beta_{i}/pi"] for i in range(p)], axis=1)
        df_orig['gamma'] = df_orig.apply(lambda row: [row[f"gamma_{i}/pi"] for i in range(p)], axis=1)
        df_orig['theta'] = df_orig.apply(lambda row: row['gamma'] + row['beta'], axis=1)
        df_orig['n'] = df_orig.apply(lambda row: row['G'].number_of_nodes(), axis=1)
        assert((df_orig['n'] == n_qubits).all())
        df = df.append(df_orig.reset_index())
assert(len(df) == sum(v * 3 for v in n_graphs.values()))
df.to_pickle(Path(build_full_qaoa_dataset_table_folder, "../data/lookup_tables/full_qaoa_dataset_table.p"))

print("Loading the results for 3-regular graphs...")
with open(Path(build_full_qaoa_dataset_table_folder, "../data/3_regular/3r_WURTZ_ensemble.json")) as json_file:
     data = json.load(json_file)

rows = []
for row in data:
    G = nx.Graph()
    G.add_edges_from(row['edges'])
    g = pynauty.Graph(number_of_vertices=G.number_of_nodes(), directed=nx.is_directed(G),
            adjacency_dict = get_adjacency_dict(G))
    cert = pynauty.certificate(g)
    c_true_opt = row['0']['MaxCut']
    
    for p in range(1,11):
        d = {}
        d['G'] = copy.deepcopy(G)
        d['n'] = G.number_of_nodes()
        d['p_max'] = p
        d['pynauty_cert'] = cert
        d['C_{true opt}'] = c_true_opt
        d['C_fixed'] = row[str(p)]['fixed_val']
        if 'optimized_val' in row[str(p)]:
            d['C_opt'] = row[str(p)]['optimized_val']
            d['beta'] = np.array(row[str(p)]['angles'][0]['beta']) / np.pi
            d['gamma'] = -2 * np.array(row[str(p)]['angles'][0]['gamma']) / np.pi
            d['theta'] = np.hstack([d['gamma'], d['beta']])
        rows.append(copy.deepcopy(d))

df = pd.DataFrame(rows, columns=rows[0].keys())
df.to_pickle(Path(build_full_qaoa_dataset_table_folder, "../data/lookup_tables/3_reg_dataset_table.p"))

import numpy as np
import pandas as pd
import copy
import pickle

n_qubits = 4
p = 3

tables = {}

for n_qubits in range(3,10):
    tables[n_qubits] = {}
    
    for p in range(1,4):
        tables[n_qubits][p] = {}
        
        colnames=['graph_id','C_{true opt}','C_init','C_opt','pr(max)','p']
        for i in range(p):
            colnames.append(f"beta_{i}/pi")
        for i in range(p):
            colnames.append(f"gamma_{i}/pi")
        df = pd.read_csv(f"../data/qaoa-dataset-version1/Results/p={p}/n={n_qubits}_p={p}.txt", delim_whitespace=True, names=colnames, header=None, index_col='graph_id')
        
        for index, data in df.iterrows():
            if data['p'] != p:
                assert(np.isclose(data['C_{true opt}'], data['C_opt']))
            beta = [data[f"beta_{i}/pi"] for i in range(p)]
            gamma = [data[f"gamma_{i}/pi"] for i in range(p)]
            tables[n_qubits][p][int(index)] = {'beta':copy.deepcopy(beta), 'gamma':copy.deepcopy(gamma)} 
        print(f"Done with n={n_qubits}, p={p}")

pickle.dump(tables, open(f"../data/lookup_tables/graph2angles.p", "wb"))

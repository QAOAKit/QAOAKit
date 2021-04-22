from QAOA_parameters import opt_angles_for_graph, get_full_qaoa_dataset_table_row

import networkx as nx

for n in range(3,10):
    G = nx.complete_graph(n)
    for p in range(1,4):
        print(opt_angles_for_graph(G,p))
        print(get_full_qaoa_dataset_table_row(G,p)[['C_opt','graph_id','C_{true opt}']])


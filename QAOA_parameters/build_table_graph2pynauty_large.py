import networkx as nx
import re
import copy
import pynauty
import pickle
from pathlib import Path

from utils import get_adjacency_dict

build_table_graph2pynauty_folder = Path(__file__).parent

n_graphs={}

n_graphs[3]=2
n_graphs[4]=6
n_graphs[5]=21
n_graphs[6]=112
n_graphs[7]=853
n_graphs[8]=11117
n_graphs[9]=261080

for n_qubits in range(3,10):
    table = {}

    table['graph_id2pynautycert'] = {}
    table['graph_id2graph'] = {}
    table['pynautycert2graph_id'] = {}
    table['pynautycert2graph'] = {}

    with open(Path(build_table_graph2pynauty_folder, "../data/qaoa-dataset-version1/Graphs/graph"+str(n_qubits)+"c.txt")) as f:
        for graph in range(1,n_graphs[n_qubits]+1):
            f.readline(-1)#first line is blank
            line_with_id = f.readline(-1) #second line has graph number and order
            graph_id, graph_order = [int(x) for x in re.split(' |, |. |.\n', line_with_id) if x.isdigit()]
            assert(graph_order == n_qubits)
            G=nx.Graph()
            edge_id = 0
            for n in range(n_qubits):
                G.add_nodes_from([n])
            #third line is first row of upper triangle of adjacency matrix (without the diagonal element)
            for n in range(n_qubits-1):
                adj_str = f.readline(-1)
                for m in range(n_qubits-1-n):
                    q_num=n+m+1
                    if adj_str[m]=='1':
                        G.add_edge(n,q_num,edge_id=edge_id)
                        edge_id += 1
            g = pynauty.Graph(number_of_vertices=G.number_of_nodes(), directed=nx.is_directed(G),
                        adjacency_dict = get_adjacency_dict(G))
            cert = pynauty.certificate(g)
            table['graph_id2graph'][graph_id] = copy.deepcopy(G)
            table['graph_id2pynautycert'][graph_id] = cert
            table['pynautycert2graph_id'][cert] = graph_id
            table['pynautycert2graph'][cert] = copy.deepcopy(G)

    assert(len(table['graph_id2pynautycert']) == n_graphs[n_qubits])
    assert(len(table['graph_id2graph']) == n_graphs[n_qubits])
    assert(len(table['pynautycert2graph_id']) == n_graphs[n_qubits])
    assert(len(table['pynautycert2graph']) == n_graphs[n_qubits])
    print(f"Done with n={n_qubits}")
    pickle.dump(table, open(Path(build_table_graph2pynauty_folder, f"../data/lookup_tables/graph2pynauty_large_{n_qubits}.p"), "wb"))

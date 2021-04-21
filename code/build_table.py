import networkx as nx
import re
import copy
import pynauty
import pickle

def get_adjacency_dict(G):
    """Returns adjacency dictionary for G
    G must be a networkx graph
    Return format: { n : [n1,n2,...], ... }
    where [n1,n2,...] is a list of neighbors of n
    """
    adjacency_dict = {}
    for n, neigh_dict in G.adjacency():
        neigh_list = []
        for neigh, attr_dict in neigh_dict.items():
            assert(len(attr_dict) == 0)
            neigh_list.append(neigh)
        adjacency_dict[n] = neigh_list
    return adjacency_dict

p=3

n_graphs={}

n_graphs[3]=2
n_graphs[4]=6
n_graphs[5]=21
n_graphs[6]=112
n_graphs[7]=853
n_graphs[8]=11117
n_graphs[9]=261080

tables = {}

for n_qubits in range(3,10):
    tables[n_qubits] = {}

    tables[n_qubits]['graph_id2pynautycert'] = {}
    tables[n_qubits]['graph_id2graph'] = {}
    tables[n_qubits]['pynautycert2graph_id'] = {}
    tables[n_qubits]['pynautycert2graph'] = {}
    
    with open("../data/qaoa-dataset-version1/Graphs/graph"+str(n_qubits)+"c.txt") as f:
        for graph in range(1,n_graphs[n_qubits]+1):
            f.readline(-1)#first line is blank
            line_with_id = f.readline(-1) #second line has graph number and order
            graph_id, graph_order = [int(x) for x in re.split(' |, |. |.\n', line_with_id) if x.isdigit()]
            assert(graph_order == n_qubits)
            edges=[]
            G=nx.Graph()
            for n in range(n_qubits):
                G.add_nodes_from([n])
            #third line is first row of upper triangle of adjacency matrix (without the diagonal element)
            for n in range(n_qubits-1):
                adj_str = f.readline(-1)
                for m in range(n_qubits-1-n):
                    q_num=n+m+1
                    if adj_str[m]=='1':
                        edges.append([n,q_num])
                        G.add_edge(n,q_num)
            g = pynauty.Graph(number_of_vertices=G.number_of_nodes(), directed=nx.is_directed(G),
                        adjacency_dict = get_adjacency_dict(G))
            cert = pynauty.certificate(g)
            tables[n_qubits]['graph_id2graph'][graph_id] = copy.deepcopy(G) 
            tables[n_qubits]['graph_id2pynautycert'][graph_id] = cert
            tables[n_qubits]['pynautycert2graph_id'][cert] = graph_id
            tables[n_qubits]['pynautycert2graph'][cert] = copy.deepcopy(G)
    
    assert(len(tables[n_qubits]['graph_id2pynautycert']) == n_graphs[n_qubits])
    assert(len(tables[n_qubits]['graph_id2graph']) == n_graphs[n_qubits])
    assert(len(tables[n_qubits]['pynautycert2graph_id']) == n_graphs[n_qubits])
    assert(len(tables[n_qubits]['pynautycert2graph']) == n_graphs[n_qubits])
    print(f"Done with n={n_qubits}")

pickle.dump(tables, open(f"../data/lookup_tables/graph2pynauty.p", "wb"))

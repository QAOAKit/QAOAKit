import networkx as nx
import re
import copy


n_qubits=7
p=3

n_graphs={}

n_graphs[3]=2
n_graphs[4]=6
n_graphs[5]=21
n_graphs[6]=112
n_graphs[7]=853
n_graphs[8]=11117
n_graphs[9]=261080

bipartite_count=0
planar_count=0
eulerian_count=0

graphs = {}

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
        graphs[graph_id] = copy.deepcopy(G) 

print(graphs)

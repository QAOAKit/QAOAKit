import networkx as nx
import numpy as np
import pandas as pd
from qiskit import execute, Aer
from functools import partial
import pickle

from qiskit.optimization.applications.ising.max_cut import get_operator
from qiskit.aqua.algorithms.minimum_eigen_solvers.qaoa.var_form import QAOAVarForm
from qiskit.quantum_info import Statevector

from QAOA_parameters import opt_angles_for_graph, get_graph_id, angles_to_qaoa_format, beta_to_qaoa_format, gamma_to_qaoa_format, angles_to_qiskit_format, get_full_qaoa_dataset_table_row, get_full_qaoa_dataset_table, qaoa_maxcut_energy
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
                assert(np.allclose(angles['beta'], beta_to_qaoa_format(row['beta'])))
                assert(np.allclose(angles['gamma'], gamma_to_qaoa_format(row['gamma'])))
                qc = get_maxcut_qaoa_circuit(row['G'], angles['beta'], angles['gamma'])
                backend = Aer.get_backend('statevector_simulator')
                res = backend.run(qc).result()
                sv = res.get_statevector()
                obj_val = obj_from_statevector(sv, obj)
                assert(np.isclose(opt_cut, obj_val))

def test_qaoa_maxcut_energy():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3,4]:
        p = 3
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df['n'] == n_qubits) & (df['p_max'] == p)]
        for _, row in df.iterrows():
            assert(np.isclose(row['C_opt'], qaoa_maxcut_energy(row['G'], beta_to_qaoa_format(row['beta']), gamma_to_qaoa_format(row['gamma']))))

def test_qiskit_angle_conversion():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3,4]:
        p = 3
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df['n'] == n_qubits) & (df['p_max'] == p)]
        for _, row in df.iterrows():
            G = row['G']
            C, offset = get_operator(nx.adjacency_matrix(G))
            vf = QAOAVarForm(C.to_opflow(), p)
            angles = angles_to_qiskit_format(opt_angles_for_graph(row['G'], row['p_max']))
            qc = vf.construct_circuit(angles)
            backend = Aer.get_backend('statevector_simulator')
            sv = Statevector(backend.run(qc).result().get_statevector())
            obj_val = -(sv.expectation_value(C.to_opflow()) + offset)
            opt_cut = row['C_opt']
            assert(np.isclose(opt_cut, obj_val))

def test_qiskit_qaoa_circuit():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3,4]:
        p = 1
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df['n'] == n_qubits) & (df['p_max'] == p)]
        for _, row in df.iterrows():
            G = row['G']
            backend = Aer.get_backend('statevector_simulator')
            C, _ = get_operator(nx.adjacency_matrix(G))
            vf = QAOAVarForm(C.to_opflow(), p)
            angles1 = angles_to_qiskit_format(opt_angles_for_graph(row['G'], row['p_max']))
            qc1 = vf.construct_circuit(angles1)
            sv1 = Statevector(backend.run(qc1).result().get_statevector())
            angles2 = angles_to_qaoa_format(opt_angles_for_graph(row['G'], row['p_max']))
            qc2 = get_maxcut_qaoa_circuit(row['G'], angles2['beta'], angles2['gamma'])
            sv2 = Statevector(backend.run(qc2).result().get_statevector())
            assert(sv1.equiv(sv2))

def test_weighted_qaoa():
    G = nx.petersen_graph()
    np.random.seed(42)
    for (u,v) in G.edges():
        G[u][v]["weight"] = np.random.uniform(low=1, high=5)
    assert(nx.is_weighted(G))
    backend = Aer.get_backend('statevector_simulator')
    p=3
    angles = {'beta' : np.random.uniform(low=0, high=np.pi/2, size=3),
              'gamma': np.random.uniform(low=0, high=np.pi, size=3)}
    angles1 = angles_to_qiskit_format(angles)
    C, _ = get_operator(nx.adjacency_matrix(G))
    vf = QAOAVarForm(C.to_opflow(), p)
    qc1 = vf.construct_circuit(angles1)
    sv1 = Statevector(backend.run(qc1).result().get_statevector())
    angles2 = angles_to_qaoa_format(angles)
    qc2 = get_maxcut_qaoa_circuit(G, angles2['beta'], angles2['gamma'])
    sv2 = Statevector(backend.run(qc2).result().get_statevector())
    assert(sv1.equiv(sv2))

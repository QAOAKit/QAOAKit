import networkx as nx
import numpy as np
import pandas as pd
from qiskit.providers.aer import AerSimulator
from functools import partial
from pathlib import Path
import qtensor

from qiskit.quantum_info import Statevector

from QAOAKit import opt_angles_for_graph, get_graph_id, get_graph_from_id, angles_to_qaoa_format, beta_to_qaoa_format, gamma_to_qaoa_format, angles_to_qiskit_format, angles_to_qtensor_format, get_full_qaoa_dataset_table_row, get_full_qaoa_dataset_table, qaoa_maxcut_energy
from QAOAKit.utils import obj_from_statevector, maxcut_obj, isomorphic, load_weights_into_dataframe, load_weighted_results_into_dataframe
from QAOAKit.qaoa import get_maxcut_qaoa_circuit, get_maxcut_qaoa_qiskit_circuit
from qiskit_optimization import QuadraticProgram
from qiskit.algorithms.minimum_eigen_solvers.qaoa import QAOAAnsatz

test_utils_folder = Path(__file__).parent

def test_retrieval():
    for n in range(3,10):
        G = nx.complete_graph(n)
        for p in range(1,4):
            angles = opt_angles_for_graph(G,p)
            assert(len(angles['beta']) == p)
            assert(len(angles['gamma']) == p)
            assert(isinstance(get_full_qaoa_dataset_table_row(G,p)[['C_opt','graph_id','C_{true opt}']], pd.Series))

def test_tables_consistency():
    p = 1
    G1 = nx.complete_graph(5)
    row = get_full_qaoa_dataset_table_row(G1,p)
    G2 = get_graph_from_id(row['graph_id'], 5)
    assert(isomorphic(G1,G2))

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
                backend = AerSimulator(method="statevector")
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
            angles = opt_angles_for_graph(row['G'], row['p_max'])
            G = row['G']
            qc, C, offset = get_maxcut_qaoa_qiskit_circuit(G, p, angles_to_qiskit_format(angles))
            backend = AerSimulator(method="statevector")
            sv = Statevector(backend.run(qc).result().get_statevector())
            obj_val = -(sv.expectation_value(C) + offset)
            opt_cut = row['C_opt']
            assert(np.isclose(opt_cut, obj_val))

def test_qiskit_qaoa_circuit():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3,4]:
        p = 1
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df['n'] == n_qubits) & (df['p_max'] == p)]
        for _, row in df.iterrows():
            angles1 = opt_angles_for_graph(row['G'], row['p_max'])
            G = row['G']
            qc1, C, offset = get_maxcut_qaoa_qiskit_circuit(G, p, angles_to_qiskit_format(angles1))
            backend = AerSimulator(method="statevector")
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
    backend = AerSimulator(method="statevector")
    p=3
    angles = {'beta' : np.random.uniform(low=0, high=np.pi/2, size=3),
              'gamma': np.random.uniform(low=0, high=np.pi, size=3)}
    qc1, C, offset = get_maxcut_qaoa_qiskit_circuit(G, p, angles_to_qiskit_format(angles))
    sv1 = Statevector(backend.run(qc1).result().get_statevector())
    angles2 = angles_to_qaoa_format(angles)
    qc2 = get_maxcut_qaoa_circuit(G, angles2['beta'], angles2['gamma'])
    sv2 = Statevector(backend.run(qc2).result().get_statevector())
    assert(sv1.equiv(sv2))

def test_load_weighted_results():
    p = 2
    n = 5

    folder_path = Path(test_utils_folder, f"../data/weighted_angle_dat/p={p}/")
    df_weights = load_weights_into_dataframe(folder_path)
    df = load_weighted_results_into_dataframe(folder_path, p, n, df_weights)
    assert(np.all(
        np.isclose(
            df.head(100).apply(
                lambda row: qaoa_maxcut_energy(row['G'], beta_to_qaoa_format(row['beta']), gamma_to_qaoa_format(row['gamma'])),
                axis=1
            ),
            df.head(100)['C_opt']
        )
    ))


def test_qtensor_angle_conversion():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3,4]:
        p = 3
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df['n'] == n_qubits) & (df['p_max'] == p)]
        for _, row in df.iterrows():
            angles = angles_to_qtensor_format(opt_angles_for_graph(row['G'], row['p_max']))
            G = row['G']
            obj_val = qtensor.QAOA_energy(G, angles['gamma'], angles['beta'])
            opt_cut = row['C_opt']
            assert(np.isclose(opt_cut, obj_val))


import networkx as nx
import numpy as np
import pandas as pd
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute
from qiskit.providers.aer import AerSimulator
from functools import partial
from pathlib import Path
import copy
import pytest
from itertools import groupby
import timeit

from qiskit.quantum_info import Statevector

from QAOAKit import (
    opt_angles_for_graph,
    get_fixed_angles,
    get_graph_id,
    get_graph_from_id,
    angles_to_qaoa_format,
    beta_to_qaoa_format,
    gamma_to_qaoa_format,
    angles_to_qiskit_format,
    angles_to_qtensor_format,
    get_3_reg_dataset_table,
    get_3_reg_dataset_table_row,
    get_full_qaoa_dataset_table_row,
    get_full_qaoa_dataset_table,
    get_fixed_angle_dataset_table,
    get_fixed_angle_dataset_table_row,
    qaoa_maxcut_energy,
    angles_from_qiskit_format,
)
from QAOAKit.utils import (
    obj_from_statevector,
    precompute_energies,
    maxcut_obj,
    isomorphic,
    load_weights_into_dataframe,
    load_weighted_results_into_dataframe,
    get_adjacency_matrix,
    brute_force,
    get_pynauty_certificate,
    get_full_weighted_qaoa_dataset_table,
)
from QAOAKit.classical import thompson_parekh_marwaha
from QAOAKit.qaoa import get_maxcut_qaoa_circuit
from QAOAKit.qiskit_interface import get_maxcut_qaoa_qiskit_circuit, goemans_williamson
from QAOAKit.examples_utils import get_20_node_erdos_renyi_graphs
from QAOAKit.parameter_optimization import get_median_pre_trained_kde

from qiskit_optimization import QuadraticProgram
from qiskit.algorithms.minimum_eigen_solvers.qaoa import QAOAAnsatz

test_utils_folder = Path(__file__).parent


def test_retrieval():
    for n in range(3, 10):
        G = nx.complete_graph(n)
        for p in range(1, 4):
            angles = opt_angles_for_graph(G, p)
            assert len(angles["beta"]) == p
            assert len(angles["gamma"]) == p
            assert isinstance(
                get_full_qaoa_dataset_table_row(G, p)[
                    ["C_opt", "graph_id", "C_{true opt}"]
                ],
                pd.Series,
            )


def test_fixed_angle_retrieval():
    for d in range(3, 5):
        for p in range(1, 4):
            angles = get_fixed_angles(d, p)
            assert type(angles) is dict
            assert len(angles["beta"]) == p
            assert len(angles["gamma"]) == p


def test_tables_consistency():
    p = 1
    G1 = nx.complete_graph(5)
    row = get_full_qaoa_dataset_table_row(G1, p)
    G2 = get_graph_from_id(row["graph_id"], 5)
    assert isomorphic(G1, G2)


def test_angle_conversion():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3, 5]:
        for p in [2, 3]:
            df = full_qaoa_dataset_table.reset_index()
            df = df[(df["n"] == n_qubits) & (df["p_max"] == p)]
            for _, row in df.iterrows():
                obj = partial(maxcut_obj, w=get_adjacency_matrix(row["G"]))
                opt_cut = row["C_opt"]
                angles = angles_to_qaoa_format(
                    opt_angles_for_graph(row["G"], row["p_max"])
                )
                assert np.allclose(angles["beta"], beta_to_qaoa_format(row["beta"]))
                assert np.allclose(angles["gamma"], gamma_to_qaoa_format(row["gamma"]))
                qc = get_maxcut_qaoa_circuit(row["G"], angles["beta"], angles["gamma"])
                backend = AerSimulator(method="statevector")
                res = backend.run(qc).result()
                sv = res.get_statevector()
                obj_val = obj_from_statevector(sv, obj)
                assert np.isclose(opt_cut, obj_val)


def test_qaoa_maxcut_energy():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3, 4]:
        p = 3
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df["n"] == n_qubits) & (df["p_max"] == p)]
        for _, row in df.iterrows():
            assert np.isclose(
                row["C_opt"],
                qaoa_maxcut_energy(
                    row["G"],
                    beta_to_qaoa_format(row["beta"]),
                    gamma_to_qaoa_format(row["gamma"]),
                ),
            )


def test_qiskit_angle_conversion():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3, 4]:
        p = 3
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df["n"] == n_qubits) & (df["p_max"] == p)]
        for _, row in df.iterrows():
            angles = opt_angles_for_graph(row["G"], row["p_max"])
            G = row["G"]
            qc, C, offset = get_maxcut_qaoa_qiskit_circuit(
                G, p, angles_to_qiskit_format(angles)
            )
            backend = AerSimulator(method="statevector")
            sv = Statevector(backend.run(qc).result().get_statevector())
            obj_val = -(sv.expectation_value(C) + offset)
            opt_cut = row["C_opt"]
            assert np.isclose(opt_cut, obj_val)


def test_qiskit_qaoa_circuit():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3, 4]:
        p = 1
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df["n"] == n_qubits) & (df["p_max"] == p)]
        for _, row in df.iterrows():
            angles1 = opt_angles_for_graph(row["G"], row["p_max"])
            G = row["G"]
            qc1, C, offset = get_maxcut_qaoa_qiskit_circuit(
                G, p, angles_to_qiskit_format(angles1)
            )
            backend = AerSimulator(method="statevector")
            sv1 = Statevector(backend.run(qc1).result().get_statevector())
            angles2 = angles_to_qaoa_format(
                opt_angles_for_graph(row["G"], row["p_max"])
            )
            qc2 = get_maxcut_qaoa_circuit(row["G"], angles2["beta"], angles2["gamma"])
            sv2 = Statevector(backend.run(qc2).result().get_statevector())
            assert sv1.equiv(sv2)


def test_weighted_qaoa():
    G = nx.petersen_graph()
    np.random.seed(42)
    for (u, v) in G.edges():
        G[u][v]["weight"] = np.random.uniform(low=1, high=5)
    assert nx.is_weighted(G)
    backend = AerSimulator(method="statevector")
    p = 3
    angles = {
        "beta": np.random.uniform(low=0, high=np.pi / 2, size=3),
        "gamma": np.random.uniform(low=0, high=np.pi, size=3),
    }
    qc1, C, offset = get_maxcut_qaoa_qiskit_circuit(
        G, p, angles_to_qiskit_format(angles)
    )
    sv1 = Statevector(backend.run(qc1).result().get_statevector())
    angles2 = angles_to_qaoa_format(angles)
    qc2 = get_maxcut_qaoa_circuit(G, angles2["beta"], angles2["gamma"])
    sv2 = Statevector(backend.run(qc2).result().get_statevector())
    assert sv1.equiv(sv2)


def test_load_weighted_results():
    p = 2
    n = 5

    folder_path = Path(test_utils_folder, f"../data/weighted_angle_dat/p={p}/")
    df_weights = load_weights_into_dataframe(folder_path)
    df = load_weighted_results_into_dataframe(folder_path, p, n, df_weights)
    assert np.all(
        np.isclose(
            df.head(100).apply(
                lambda row: qaoa_maxcut_energy(
                    row["G"],
                    beta_to_qaoa_format(row["beta"]),
                    gamma_to_qaoa_format(row["gamma"]),
                ),
                axis=1,
            ),
            df.head(100)["C_opt"],
        )
    )


def test_qtensor_angle_conversion():
    qtensor = pytest.importorskip("qtensor")
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3, 4]:
        p = 3
        df = full_qaoa_dataset_table.reset_index()
        df = df[(df["n"] == n_qubits) & (df["p_max"] == p)]
        for _, row in df.iterrows():
            angles = angles_to_qtensor_format(
                opt_angles_for_graph(row["G"], row["p_max"])
            )
            G = row["G"]
            obj_val = qtensor.QAOA_energy(G, angles["gamma"], angles["beta"])
            opt_cut = row["C_opt"]
            assert np.isclose(opt_cut, obj_val)


def test_from_qiskit_conversion():
    angles = {
        "beta": np.random.uniform(low=0, high=np.pi / 2, size=3),
        "gamma": np.random.uniform(low=0, high=np.pi, size=3),
    }
    angles2 = angles_from_qiskit_format(angles_to_qiskit_format(copy.deepcopy(angles)))
    assert np.all(np.isclose(angles["beta"], angles2["beta"]))
    assert np.all(np.isclose(angles["gamma"], angles2["gamma"]))


def test_3_reg_table():
    df = get_3_reg_dataset_table().reset_index()
    assert len(df) >= 1000

    # test that the angles are correct
    for _, row in df[(df["p_max"] == 1) | (df["p_max"] == 2)].head(50).iterrows():
        G = row["G"]
        angles = angles_to_qaoa_format({"beta": row["beta"], "gamma": row["gamma"]})
        assert np.isclose(
            qaoa_maxcut_energy(G, angles["beta"], angles["gamma"]), row["C_opt"]
        )

    # test that the optima match the full qaoa table
    for _, row in (
        df[((df["p_max"] == 1) | (df["p_max"] == 2)) & (df["n"] < 9)]
        .head(50)
        .iterrows()
    ):
        full_row = get_full_qaoa_dataset_table_row(row["G"], row["p_max"])
        assert np.isclose(full_row["C_opt"], row["C_opt"])


def test_3_reg_degenerate_optima():
    df = get_3_reg_dataset_table().reset_index()
    for p in [1, 2]:
        for _, row in df[df["p_max"] == p].head(10).iterrows():
            for beta, gamma in zip(
                row["all beta (degenerate optima)"],
                row["all gamma (degenerate optima)"],
            ):
                angles = angles_to_qaoa_format({"beta": beta, "gamma": gamma})
                assert np.isclose(
                    qaoa_maxcut_energy(row["G"], angles["beta"], angles["gamma"]),
                    row["C_opt"],
                )


def test_fixed_angles_3_reg():
    df = get_3_reg_dataset_table().sample(n=10).reset_index()

    for _, row in df.iterrows():
        angles = angles_to_qaoa_format(get_fixed_angles(3, row["p_max"]))
        assert np.isclose(
            qaoa_maxcut_energy(row["G"], angles["beta"], angles["gamma"]),
            row["C_fixed"],
            rtol=1e-4,
        )


def test_fixed_angles_all():
    df = get_fixed_angle_dataset_table()
    n = 12

    for _, row in df.iterrows():
        beta = beta_to_qaoa_format(row["beta"])
        gamma = gamma_to_qaoa_format(row["gamma"])
        for s in range(5):
            G = nx.random_regular_graph(row["d"], n, seed=s)
            obj = partial(maxcut_obj, w=get_adjacency_matrix(G))
            opt_en = brute_force(obj, n)[0]
            fixed_angle_en = qaoa_maxcut_energy(G, beta, gamma)
            assert fixed_angle_en / opt_en >= row["AR"]


def test_get_opt_angles_large_3_reg():
    for s in range(10):
        G = nx.random_regular_graph(3, 10, seed=s)
        for p in [1, 2]:
            angles = angles_to_qaoa_format(opt_angles_for_graph(G, p))
            row = get_3_reg_dataset_table_row(G, p)
            assert np.isclose(
                row["C_opt"], qaoa_maxcut_energy(G, angles["beta"], angles["gamma"])
            )


def test_brute_force():
    n = 8
    df = get_3_reg_dataset_table().reset_index()
    df = df[(df["n"] == n) & (df["p_max"] == 1)]
    for _, row in df.iterrows():
        obj = partial(maxcut_obj, w=get_adjacency_matrix(row["G"]))
        assert np.isclose(row["C_{true opt}"], brute_force(obj, n)[0])


def test_brute_force_minimize_unweighted():
    def all_equal(iterable):
        g = groupby(iterable)
        return next(g, True) and not next(g, False)

    n = 4
    G = nx.complete_graph(n)
    obj = partial(maxcut_obj, w=get_adjacency_matrix(G))
    fval, soln = brute_force(obj, n, minimize=True)
    assert np.isclose(fval, 0)
    assert all_equal(soln)


def test_brute_force_minimize_neg_weighted():
    n = 4
    G = nx.complete_graph(n)
    for u, v in G.edges():
        G[u][v]["weight"] = -1
    G[0][1]["weight"] = 5
    G[2][3]["weight"] = 5
    obj = partial(maxcut_obj, w=get_adjacency_matrix(G))
    fval, _ = brute_force(obj, n, minimize=True)
    assert np.isclose(fval, -4)


def test_get_opt_angles_k_reg():
    n = 10
    for d, max_p in [(3, 11), (4, 4), (5, 3)]:
        for p in range(3, max_p + 1):
            AR = get_fixed_angle_dataset_table_row(d, p).AR
            for s in range(5):
                G = nx.random_regular_graph(d, n, seed=s)
                obj = partial(maxcut_obj, w=get_adjacency_matrix(G))
                opt_en = brute_force(obj, n)[0]
                with pytest.warns(
                    Warning,
                    match="Optimal angles not available, returning fixed angles",
                ):
                    angles = angles_to_qaoa_format(opt_angles_for_graph(G, p))
                fixed_angle_en = qaoa_maxcut_energy(G, angles["beta"], angles["gamma"])
                assert fixed_angle_en / opt_en >= AR


def test_get_opt_angles_large_non_reg():
    G = nx.lollipop_graph(5, 10)
    assert G.number_of_nodes() > 10
    with pytest.warns(
        Warning,
        match="Optimal angles not available, returning closest fixed angles",
    ):
        angles = angles_to_qaoa_format(opt_angles_for_graph(G, 2))


def test_example_in_README():
    # build graph
    G = nx.star_graph(5)
    # grab optimal angles
    p = 3
    angles = angles_to_qaoa_format(opt_angles_for_graph(G, p))
    # build circuit
    qc = get_maxcut_qaoa_circuit(G, angles["beta"], angles["gamma"])
    qc.measure_all()
    # run circuit
    backend = AerSimulator()
    assert isinstance(backend.run(qc).result().get_counts(), dict)


def test_no_save_state():
    n = 6
    qc = get_maxcut_qaoa_circuit(nx.star_graph(n - 1), [0.1], [0.1], save_state=False)
    creg = ClassicalRegister(n)
    creg_ancilla = ClassicalRegister(1)
    qreg_ancilla = QuantumRegister(1)
    qc.add_register(creg)
    qc.add_register(creg_ancilla)
    qc.add_register(qreg_ancilla)
    qc.h(qreg_ancilla)
    for i in range(n):
        qc.cx(qreg_ancilla, i)
    qc.h(qreg_ancilla)

    qc.measure(range(6), creg)
    qc.measure(qreg_ancilla, creg_ancilla)

    counts = (
        execute(qc, AerSimulator(), basis_gates=["u1", "u2", "u3", "cx"])
        .result()
        .get_counts()
    )
    assert isinstance(counts, dict)


def test_examples_erdos_renyi():
    df = get_20_node_erdos_renyi_graphs()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 30


def test_gw():
    df = get_3_reg_dataset_table().reset_index()
    df = df[df["p_max"] == 1]

    def check_gw(G, QAOA_val, true_opt_val):
        nsamples = 100
        soln = goemans_williamson(G, nsamples=nsamples)
        obj = partial(maxcut_obj, w=get_adjacency_matrix(G))
        scores = [obj(x) for x in soln]
        assert len(soln) == nsamples
        assert np.max(scores) >= QAOA_val
        assert np.max(scores) <= true_opt_val

    df.head(5).apply(
        lambda row: check_gw(row["G"], row["C_opt"], row["C_{true opt}"]), axis=1
    )


def test_tpm():
    G = nx.random_regular_graph(3, 1000, seed=42)

    soln, exp_round = thompson_parekh_marwaha(G, nsamples=100, girth=100)
    obj = partial(maxcut_obj, w=get_adjacency_matrix(G))
    vals = [obj(x) / G.number_of_edges() for x in soln]
    assert np.isclose(np.average(vals), exp_round, atol=0.01)


def test_pass_quantum_register_to_qaoa_circuit_generator():
    # Use case: appending a project circuit for error mitigation to the generated QAOA circuit
    def get_bit_flip_projector(G, qr):
        """N is the number of qubits
        returns projector circuit (qiskit.QuantumCircuit) for bit-flip symmetry
        """
        N = G.number_of_nodes()
        qc = QuantumCircuit(qr)
        qc.h(N)
        for i in range(0, N):
            qc.cx(N, i)
        qc.h(N)
        return qc

    G = nx.star_graph(4)
    N = G.number_of_nodes()
    qr = QuantumRegister(N + 1)
    cr = ClassicalRegister(N + 1)

    p = 2
    angles = angles_to_qaoa_format(opt_angles_for_graph(G, p))
    qc = get_maxcut_qaoa_circuit(
        G, angles["beta"], angles["gamma"], transpile_to_basis=True, qr=qr, cr=cr
    )
    qc.compose(get_bit_flip_projector(G, qr), inplace=True)
    qc.measure(qr[N], cr[N])

    print(qc.draw())


def test_get_median_pre_trained_kde():
    n = 8
    n_kde_samples = 1
    for p in [1, 2]:
        median, kde = get_median_pre_trained_kde(p)
        new_data = kde.sample(n_kde_samples, random_state=0)
        angles = np.vstack([np.atleast_2d(median), new_data])
        assert angles.shape == (n_kde_samples + 1, p * 2)
        converted_angles = []
        for angle in angles:
            converted_angles.append(
                np.hstack(
                    [gamma_to_qaoa_format(angle[:p]), beta_to_qaoa_format(angle[p:])]
                )
            )
        angles = np.vstack(converted_angles)
        assert angles.shape == (n_kde_samples + 1, p * 2)

        df = get_3_reg_dataset_table().reset_index()
        df = df[(df["n"] == n) & (df["p_max"] == p)]
        for _, row in df.head(10).iterrows():
            ave_d = 2 * row["G"].number_of_edges() / row["G"].number_of_nodes()
            if p == 1:
                scaling_factor = np.arctan(1 / np.sqrt(ave_d - 1))
            else:
                scaling_factor = 1 / np.sqrt(ave_d)
            en_transf = max(
                qaoa_maxcut_energy(row["G"], angle[p:], angle[:p] * scaling_factor)
                for angle in angles
            )
            assert en_transf >= 0.9 * row["C_opt"]


def test_get_pynauty_certificate():
    elist = [[0, 1], [1, 2], [2, 3]]
    G1 = nx.Graph()
    G1.add_edges_from(elist)
    G2 = nx.Graph(elist[::-1])

    assert get_pynauty_certificate(G1) == get_pynauty_certificate(G2)


def test_obj_from_statevector():
    full_qaoa_dataset_table = get_full_qaoa_dataset_table()
    for n_qubits in [3, 5]:
        for p in [2, 3]:
            df = full_qaoa_dataset_table.reset_index()
            df = df[(df["n"] == n_qubits) & (df["p_max"] == p)]
            for _, row in df.iterrows():
                obj = partial(maxcut_obj, w=get_adjacency_matrix(row["G"]))
                opt_cut = row["C_opt"]
                angles = angles_to_qaoa_format(
                    opt_angles_for_graph(row["G"], row["p_max"])
                )
                qc = get_maxcut_qaoa_circuit(row["G"], angles["beta"], angles["gamma"])
                backend = AerSimulator(method="statevector")
                res = backend.run(qc).result()
                sv = res.get_statevector()
                obj_val = obj_from_statevector(sv, obj)
                assert np.isclose(opt_cut, obj_val)
                precomputed_energies = precompute_energies(obj, n_qubits)
                obj_val2 = obj_from_statevector(
                    sv, obj, precomputed_energies=precomputed_energies
                )
                assert np.isclose(opt_cut, obj_val2)


def test_precomputed_energies_fast():
    n_qubits = 20
    G = nx.random_regular_graph(3, n_qubits, seed=42)
    obj = partial(maxcut_obj, w=get_adjacency_matrix(G))
    precomputed_energies = precompute_energies(obj, n_qubits)
    beta = np.random.uniform(0, np.pi, 2)
    gamma = np.random.uniform(0, np.pi, 2)
    t1 = timeit.timeit(lambda: qaoa_maxcut_energy(G, beta, gamma), number=5)
    t2 = timeit.timeit(
        lambda: qaoa_maxcut_energy(
            G, beta, gamma, precomputed_energies=precomputed_energies
        ),
        number=5,
    )
    assert 2 * t2 < t1


def test_weighted_table():
    df = get_full_weighted_qaoa_dataset_table()

    # test that the angles are correct
    for _, row in df[(df["n"] == 8) | (df["p_max"] == 2)].head(50).iterrows():
        angles = angles_to_qaoa_format({"beta": row["beta"], "gamma": row["gamma"]})
        assert np.isclose(
            qaoa_maxcut_energy(row["G"], angles["beta"], angles["gamma"]), row["C_opt"]
        )

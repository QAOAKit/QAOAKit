import numpy as np
import networkx as nx
from qiskit_optimization import QuadraticProgram
from qiskit.algorithms.minimum_eigen_solvers.qaoa import QAOAAnsatz


def get_maxcut_quadratic_problem(G):
    """
    Construct Qiskit QuadraticProgram for MaxCut on graph G
    """
    n_qubits = G.number_of_nodes()
    problem = QuadraticProgram()
    _ = [problem.binary_var("x{}".format(i)) for i in range(n_qubits)]
    problem.maximize(
        linear=nx.adjacency_matrix(G).dot(np.ones(n_qubits)),
        quadratic=-nx.adjacency_matrix(G),
    )
    return problem


def get_maxcut_qaoa_qiskit_circuit(G, p, qiskit_angles):
    """
    Constructs a max cut QAOA circuit with Qiskit.

    Returns (quantum circuit, cost operator, cut size offset).
    On average, the circuit gives max cut size proportional to -(<C> + offset)
    """
    problem = get_maxcut_quadratic_problem(G)
    C, offset = problem.to_ising()
    ansatz = QAOAAnsatz(C, p).decompose()
    qc = ansatz.bind_parameters(qiskit_angles)
    qc.save_state()
    return qc, C, offset

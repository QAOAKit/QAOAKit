import numpy as np
import networkx as nx
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import GoemansWilliamsonOptimizer
from qiskit.algorithms.minimum_eigen_solvers.qaoa import QAOAAnsatz

from .utils import get_adjacency_matrix


def get_maxcut_quadratic_problem(G):
    """
    Construct Qiskit QuadraticProgram for MaxCut on graph G
    """
    n_qubits = G.number_of_nodes()
    problem = QuadraticProgram()
    _ = [problem.binary_var("x{}".format(i)) for i in range(n_qubits)]
    problem.maximize(
        linear=get_adjacency_matrix(G).dot(np.ones(n_qubits)),
        quadratic=-get_adjacency_matrix(G),
    )
    return problem


def goemans_williamson(G, nsamples):
    """
    Get GoemansWilliamson solutions using qiskit implementation
    """
    problem = get_maxcut_quadratic_problem(G)
    algorithm = GoemansWilliamsonOptimizer(num_cuts=nsamples, unique_cuts=False)
    return [z.x for z in algorithm.solve(problem).samples]


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

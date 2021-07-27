# QAOA circuit for MAXCUT

import networkx as nx
from qiskit import QuantumCircuit, Aer, execute
from qiskit.compiler import transpile
from qiskit_optimization import QuadraticProgram
from qiskit.algorithms.minimum_eigen_solvers.qaoa import QAOAAnsatz
import numpy as np

def append_zz_term(qc, q1, q2, gamma):
    qc.cx(q1,q2)
    qc.rz(2*gamma, q2)
    qc.cx(q1,q2)

def get_maxcut_cost_operator_circuit(G, gamma):
    N = G.number_of_nodes()
    qc = QuantumCircuit(N)
    for i, j in G.edges():
        if nx.is_weighted(G):
            append_zz_term(qc, i, j, gamma*G[i][j]["weight"])
        else:
            append_zz_term(qc, i, j, gamma)
    return qc

def append_x_term(qc, q1, beta):
    qc.h(q1)
    qc.rz(2*beta, q1)
    qc.h(q1)

def get_mixer_operator_circuit(G, beta):
    N = G.number_of_nodes()
    qc = QuantumCircuit(N)
    for n in G.nodes():
        append_x_term(qc, n, beta)
    return qc

def get_maxcut_qaoa_circuit(G, beta, gamma, transpile_to_basis=True):
    assert(len(beta) == len(gamma))
    p = len(beta) # infering number of QAOA steps from the parameters passed
    N = G.number_of_nodes()
    qc = QuantumCircuit(N)
    # first, apply a layer of Hadamards
    qc.h(range(N))
    # second, apply p alternating operators
    for i in range(p):
        qc = qc.compose(get_maxcut_cost_operator_circuit(G,gamma[i]))
        qc = qc.compose(get_mixer_operator_circuit(G,beta[i]))
    if transpile_to_basis:
        qc = transpile(qc, optimization_level=0,basis_gates=['u1', 'u2', 'u3', 'cx'])
    qc.save_state()
    return qc


def get_maxcut_qaoa_qiskit_circuit(G, p, qiskit_angles):
    """
    Constructs a max cut QAOA circuit with Qiskit.

    Returns (quantum circuit, cost operator, cut size offset).
    On average, the circuit gives max cut size proportional to -(<C> + offset)
    """
    n_qubits = len(G.nodes())
    problem = QuadraticProgram()
    _ = [problem.binary_var('x{}'.format(i)) for i in range(n_qubits)]
    problem.maximize(linear=nx.adjacency_matrix(G).dot(np.ones(n_qubits)), quadratic=-nx.adjacency_matrix(G))
    C, offset = problem.to_ising()
    ansatz = QAOAAnsatz(C, p).decompose()
    qc = ansatz.bind_parameters(qiskit_angles)
    qc.save_state()
    return qc, C, offset

# `QAOAKit`: A Toolkit for Reproducible Application and Verification of QAOA

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![Tests](https://github.com/QAOAKit/QAOAKit/actions/workflows/python-package-conda.yml/badge.svg)

### Installation

Recommended: create an Anaconda environment

```
conda create -n qaoa python=3
conda activate qaoa
```

Note that current implementation requires significant amounts of RAM (~5GB) as it loads the entire dataset into memory. Linux and macOS are currently supported.

```
pip install QAOAKit
python -m QAOAKit.build_tables
```

### Example

```python
import networkx as nx
from qiskit.providers.aer import AerSimulator
from QAOAKit import opt_angles_for_graph, angles_to_qaoa_format
from QAOAKit.qaoa import get_maxcut_qaoa_circuit

# build graph
G = nx.star_graph(5)
# grab optimal angles
p = 3
angles = angles_to_qaoa_format(opt_angles_for_graph(G,p))
# build circuit
qc = get_maxcut_qaoa_circuit(G, angles['beta'], angles['gamma'])
qc.measure_all()
# run circuit
backend = AerSimulator()
print(backend.run(qc).result().get_counts())
```

Almost all counts you get should correspond to one of the two optimal MaxCut solutions for star graph: `000001` or `111110`.

For graphs where no pre-optimized angles are available, the angles from "The fixed angle conjecture for QAOA on regular MaxCut graphs" ([arXiv:2107.00677](https://scirate.com/arxiv/2107.00677)) will be returned.

### Advanced usage

More advanced examples are available in `examples` folder:

- Using optimal parameters in state-of-the-art tensor network QAOA simulator [QTensor](https://github.com/danlkv/QTensor): [`examples/qtensor_get_energy.py`](https://github.com/QAOAKit/QAOAKit/blob/master/examples/qtensor_get_energy.py)
- Transferring parameters to large unseen instances: [`examples/Transferability to unseen instances.ipynb`](https://github.com/QAOAKit/QAOAKit/blob/master/examples/Transferability%20to%20unseen%20instances.ipynb)
- Exploring the role of problem structure in QAOA performance [`examples/Tackling open problems.ipynb`](https://github.com/QAOAKit/QAOAKit/blob/master/examples/QAOA%20symmetry%20and%20performance.ipynb)
- Exploring the performance of QAOA on small graphs as a function of average degree: [`examples/performance.ipynb`](https://github.com/QAOAKit/QAOAKit/blob/master/examples/performance.ipynb)
- Running classical algorithms for MaxCut: [`examples/classical_algorithms.py`](https://github.com/QAOAKit/QAOAKit/blob/master/examples/classical_algorithms_vs_qaoa.py)
- Comparing QAOA with classical algorithms for MaxCut: [`examples/classical_vs_quantum.ipynb`](https://github.com/QAOAKit/QAOAKit/blob/master/examples/classical_vs_quantum.ipynb)
- Clustering degenerate QAOA angles: [`examples/degenerate_optima_in_angle_space.py`](https://github.com/QAOAKit/QAOAKit/blob/master/examples/degenerate_optima_in_angle_space.py)

### Install from source

```
git clone https://github.com/QAOAKit/QAOAKit.git
cd QAOAKit
pip install -e .
python -m QAOAKit.build_tables
pytest
```

If you have an issue like "Illegal Instruction (core dumped)", you may have to force pip to recompile Nauty binaries (`pip install --no-binary pynauty pynauty`) or install Nauty separately: https://pallini.di.uniroma1.it/

#### Contributing

You should set up the linter to run before every commit.
```
pip install pre-commit
pre-commit install
```
Note that linter checks passing is a necessary condition for your contribution to be reviewed.

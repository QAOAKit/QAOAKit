# Set of tools for manipulating the ORNL qaoa dataset

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

### Conversion to QTensor

It is straightforward to use QTensor from these predefined angles.
```python
import networkx as nx
from QAOA_parameters import opt_angles_for_graph, angles_to_qaoa_format
from qtensor import QAOA_energy
import numpy as np

# build graph
G = nx.star_graph(5)
# grab optimal angles
p = 3
angles = angles_to_qaoa_format(opt_angles_for_graph(G,p))
# use QTensor to find energy
E = QAOA_energy(G, angles['gamma'], angles['beta'])
print(E)
```

### Installation

Optional: create an Anaconda environment

```
conda create -n qaoa python=3
conda activate qaoa
```

Note that current implementation requires significant amounts of RAM (~5GB) as it loads the entire dataset into memory.

```
git clone https://github.com/QAOAKit/QAOAKit.git
cd QAOAKit 
pip install -e .
python -m QAOAKit.build_table_graph2pynauty_large
python -m QAOAKit.build_full_qaoa_dataset_table
pytest
```

If you have an issue like "Illegal Instruction (core dumped)", you may have to force pip to recompile Nauty binaries (`pip install --no-binary pynauty pynauty`) or install Nauty separately: https://pallini.di.uniroma1.it/


### TODO

- [ ] Add angles from this recent paper: https://scirate.com/arxiv/2107.00677
- [ ] Add Kernel Density Estimation to generate initial points (https://doi.org/10.1609/aaai.v34i03.5616)
- [ ] Add optimal parameters for a few larger instances to showcase the power of transferability (add `examples` folder?)
- [ ] Libensemble example with multiprocessing that uses sampled points as initial guesses

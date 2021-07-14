# Set of tools for manipulating the ORNL qaoa dataset

### Example

```python
import networkx as nx
from qiskit import Aer
from QAOA_parameters import opt_angles_for_graph, angles_to_qaoa_format 
from QAOA_parameters.qaoa import get_maxcut_qaoa_circuit 

# build graph
G = nx.star_graph(5)
# grab optimal angles
p = 3
angles = angles_to_qaoa_format(opt_angles_for_graph(G,p))
# build circuit
qc = get_maxcut_qaoa_circuit(G, angles['beta'], angles['gamma'])
qc.measure_all()
# run circuit
backend = Aer.get_backend('qasm_simulator')
print(backend.run(qc).result().get_counts())
```

Almost all counts you get should correspond to one of the two optimal MaxCut solutions for star graph: `000001` or `111110`.

### Installation

Optional: create an Anaconda environment

```
conda create -n qaoa python=3
conda activate qaoa
```

Note that current implementation requires significant amounts of RAM (~5GB) as it loads the entire dataset into memory.

```
git clone https://github.com/rsln-s/QAOA_parameters.git
cd QAOA_parameters
pip install -e .
python QAOA_parameters/build_table_graph2pynauty_large.py
python QAOA_parameters/build_full_qaoa_dataset_table.py
pytest
```

If you have an issue like "Illegal Instruction (core dumped)", you may have to force pip to recompile Nauty binaries (`pip install --no-binary pynauty pynauty`) or install Nauty separately: https://pallini.di.uniroma1.it/


### TODO

- [ ] Qtensor parameter conversion
- [ ] Add angles from this recent paper: https://scirate.com/arxiv/2107.00677
- [ ] Update `setup.py` (what is `zip_safe`? what else should be there?)

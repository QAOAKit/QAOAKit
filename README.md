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

```
git clone https://github.com/rsln-s/QAOA_parameters.git
cd QAOA_parameters
pip install -e .
python QAOA_parameters/build_table_graph2pynauty_large.py
python QAOA_parameters/build_full_qaoa_dataset_table.py
pytest
```

### TODO

- [ ] Qtensor parameter conversion
- [ ] Update `setup.py` (what is `zip_safe`? what else should be there?)
- [ ] Continuous integration (eventually, if the goal is to make public)
- [ ] Add `pynauty-nice` to `setup.py` to ease installation (https://pypi.org/project/pynauty-nice/) or `pynauty` (https://github.com/michaelRadigan/pynauty-nice/issues/2, https://pypi.org/project/pynauty/). Maybe you need `pip install --no-binary pynauty pynauty`

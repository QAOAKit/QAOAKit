# Set of tools for manipulating the ORNL qaoa dataset


### Retrieving the dataset

The dataset is available here: https://code.ornl.gov/qci/qaoa-dataset-version1

```
git clone https://code.ornl.gov/qci/qaoa-dataset-version1.git
```

### Installation

```
pip install -e .
python QAOA_parameters/build_table_graph2pynauty_large.py
python QAOA_parameters/build_full_qaoa_dataset_table.py
pytest
```

### TODO

- [ ] Qtensor parameter conversion
- [ ] Update `setup.py` (what is `zip_safe`? what else should be there?)
- [ ] Examples (in tests?)

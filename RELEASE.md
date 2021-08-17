### Releasing

Get packages
```
python -m pip install --user --upgrade setuptools wheel
python -m pip install --user --upgrade twine
```

Build archive (don't forget to update the version in setup.py!)
```
rm dist/*
python setup.py sdist bdist_wheel
```
Upload to testpypi
```
python -m twine upload --repository testpypi dist/*
```
Test install
```
cd /tmp
conda deactivate
conda env remove -n test_qaoa
conda create -y -n test_qaoa python=3
conda activate test_qaoa
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple QAOAKit
python -m QAOAKit.build_table_graph2pynauty_large
python -m QAOAKit.build_full_qaoa_dataset_table
python -c 'from QAOAKit import opt_angles_for_graph; import networkx as nx; print(opt_angles_for_graph(nx.star_graph(5), 2))'
```

To publish to real PyPI (be careful!)
```
python -m twine upload dist/*
```

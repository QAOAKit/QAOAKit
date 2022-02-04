### Releasing

Get packages
```
python -m pip install --user --upgrade setuptools wheel
python -m pip install --user --upgrade twine
```

Build archive (don't forget to update the version in setup.py!)
```
rm -r build *.egg-info dist
python setup.py sdist bdist_wheel
```

Tag the version
```
git tag v999.999
git push origin --tags
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
conda create -y -n test_qaoa python=3.8
conda activate test_qaoa
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple QAOAKit
python -m QAOAKit.build_tables
python -c 'from QAOAKit import opt_angles_for_graph; import networkx as nx; print(opt_angles_for_graph(nx.star_graph(5), 2))'
conda deactivate
```

To publish to real PyPI (be careful!)
```
python -m twine upload dist/*
```

Test real PyPI install
```
cd /tmp
conda deactivate
conda env remove -n test_qaoa
conda create -y -n test_qaoa python=3.8
conda activate test_qaoa
pip install QAOAKit
python -m QAOAKit.build_tables
python -c 'from QAOAKit import opt_angles_for_graph; import networkx as nx; print(opt_angles_for_graph(nx.star_graph(5), 2))'
conda deactivate
```

Don't forget to get it as a release!

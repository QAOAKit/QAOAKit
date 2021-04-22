from setuptools import setup

setup(name='QAOA_parameters',
	description='Tools for using the pre-optimized parameters from ORNL dataset (and more)',
	author='Ruslan Shaydulin',
	author_email='rshaydu@g.clemson.edu',
	packages=['QAOA_parameters'],
    install_requires=['qiskit','pandas','networkx','numpy','pytest'],
	zip_safe=False)

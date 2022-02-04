from setuptools import setup

# read the contents of README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="QAOAKit",
    version="0.1.12",
    description="A Toolkit for Reproducible Study, Application and Verification of QAOA",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Ruslan Shaydulin",
    author_email="ruslan@shaydul.in",
    url="https://github.com/QAOAKit/QAOAKit",
    python_requires=">=3, <4",
    packages=["QAOAKit"],
    install_requires=[
        "qiskit==0.29.0",
        "pynauty==1.0.0",
        "qiskit-optimization",
        "pandas",
        "networkx",
        "numpy",
        "pytest",
        "tqdm",
        "cvxgraphalgs",
        "cvxopt",
        "scikit-learn==1.0",
        "notebook",
        "matplotlib",
        "seaborn",
    ],
    zip_safe=True,
)

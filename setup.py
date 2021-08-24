from setuptools import setup

setup(
    name="QAOAKit",
    version="0.1.7",
    description="A Toolkit for Reproducible Study, Application and Verification of QAOA",
    author="Ruslan Shaydulin",
    author_email="ruslan@shaydul.in",
    python_requires=">=3, <4",
    packages=["QAOAKit"],
    install_requires=[
        "qiskit>=0.28.0",
        "pynauty==1.0.0",
        "qiskit-optimization",
        "pandas",
        "networkx",
        "numpy",
        "pytest",
        "tqdm",
    ],
    zip_safe=True,
)

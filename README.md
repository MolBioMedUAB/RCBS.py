# RCBS.py
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/dynamicsUAB/RCBS.py)](https://github.com/dynamicsUAB/RCBS.py)
[![PyPI](https://img.shields.io/pypi/v/RCBS.py)](https://pypi.org/project/RCBS.py/)
![GitHub](https://img.shields.io/github/license/dynamicsUAB/RCBS.py)
[![MDAnalysis](https://img.shields.io/badge/Powered%20by-MDAnalysis-lightgray.svg)](https://www.mdanalysis.org)

## Description
RCBS.py (Reactivity of Chemical and Biochemical Systems) is a Python package that contains several scripts and functions which try to simplify the analysis of MD simulations and, from these, to prepare the system to carry out QM/MM simulations using ChemShell. MD simulations analysis has been built in top of the MDAnalysis python package.

The package is divided into two main modules: *md_analyser* and *qmmm_setup*.

### md_analyser
This module form RCBS.py is useful for performing analysis of MD trajectories in a simple way.

## Installation

### Using pip

```bash
pip install RCBS.py
```

### From source code

1. Clone the repository in you local machine

```bash
git clone https://github.com/dynamicsUAB/RCBS.py.git
```

2. Move to the folder

```bash
cd RCBS.py
```

3. Install the package using pip or pip3

```bash
pip install .
```


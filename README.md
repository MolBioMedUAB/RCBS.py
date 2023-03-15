# RCBS.py

[![DOI](https://zenodo.org/badge/266803510.svg)](https://zenodo.org/badge/latestdoi/266803510)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/dynamicsUAB/RCBS.py)](https://github.com/dynamicsUAB/RCBS.py)
[![PyPI](https://img.shields.io/pypi/v/RCBS.py)](https://pypi.org/project/RCBS.py/)
![GitHub](https://img.shields.io/github/license/dynamicsUAB/RCBS.py)
[![MDAnalysis](https://img.shields.io/badge/Powered%20by-MDAnalysis-lightgray.svg)](https://www.mdanalysis.org)

## Description

RCBS.py (Reactivity of Chemical and Biochemical Systems) is a Python package that contains several scripts, functions and classes that simplify the analysis of chemical and biochemical simulations. It is built on top of [MDAnalysis](https://www.mdanalysis.com) package. 

The package currently is divided in one main module, `md_analyser`. This module is mainly built in top of MDAnalysis' and is designed in order to symplify the analysis of MD simulations. In the future, capabilities for building QM/MM models from simulations will be added in a new module.

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
    pip install RCBS.py
    ```

## Usage

### md_analyser

`md_analyser` is made of two principal scripts: `measurements.py` and `calculators.py`.

#### `measurements.py`

This is the main file that has to be imported. In it a class called `Measurement` is created. This class has one only argument: a preloaded MDAnalysis' `Universe`. Once the object is created, measurements such as distances or angles can be created using the `add_*` functions from the class. They require previous creation of atom groups using the `select_atoms` function form `Universe` or the `selections.selection` function from RCBS.py. The syntax for the addition of a measurement is the following (shown using a distance as an example):

```python
from RCBS.md_analyser import measurements

m = measurements.Measurement(u)
m.add_distance('name', sel1, sel2)
```

More examples are available in the [example Notebook](https://github.com/dynamicsUAB/RCBS.py-examples/blob/eddcf8fe39f558717f4a86474f3938872b93adbf/md_analyser/trajectory_analysis_and_frame_extraction.ipynb) from [RCBS.py-examples](https://github.com/dynamicsUAB/RCBS.py-examples).

Once all measurements are added, `m.run_measure()` has to be used in order to run the measurements. The code works in such a way that each frame is loaded, each of the measurements is calculated and the next frame is loaded. A progress bar is shown during the process.

Once the measures are ran, `m.run_boolean()` can be used to do simple statistical studies. Examples can be found in the [example Notebook](https://github.com/dynamicsUAB/RCBS.py-examples/blob/eddcf8fe39f558717f4a86474f3938872b93adbf/md_analyser/trajectory_analysis_and_frame_extraction.ipynb).

#### `calculators.py`

In this file all the calculators are stored (work is being done to move all the calculators from `measurements.py`). Calculators can do very different things --in fact, anything you can think and code-- from the given selections. Calculators have to be implemented so they do the measure over only one structure (or frame) and they return a single variable.

`calculators.py` is loaded from `measurements.py`. Each calculator has to be incorporated in `measurements.py`, both as a `add_*` function and in the `run_measurement` function. 

## Future additions

- Create a script for automatically add new calculators to `measurements.py`.
- Add more calculators.
- Create the `qmmm` module.



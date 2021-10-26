# spin-statistics
Calculations of spin-statistical factors for triplet-triplet annihilation upconversion.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/davidbossanyi/spin-statistics/0b72942209d8d2a895c619287da9095e50f39c33)

## about
The spin-statistical factor &eta; for triplet-triplet annihilation upconversion gives the probability that an annihilation event between two triplet excitons will result in a singlet exciton. The value &eta; depends in a non-trivial way on a lot of different parameters. For the full story, please refer to [our paper recently published open access in JACS Au](https://pubs.acs.org/doi/10.1021/jacsau.1c00322). This repository contains code that allows you to perform calculations of &eta;.

If you use any of the code contained in this repository, please cite our paper [Spin statistics for triplet-triplet annihilation upconversion: exchange coupling, intermolecular orientation and reverse intersystem crossing](https://pubs.acs.org/doi/10.1021/jacsau.1c00322) by Bossanyi et al., recently published open access in JACS Au.

## spinstats package
#### installation
The code is structured as a python package, which you can install locally on your machine by navigating to the folder where you downloaded the repository and running
```sh
python setup.py install
```
You can then use the spin Hamiltonian and rate model classes in your own code:
```python
from spinstats import SpinHamiltonian
from spinstats.models import Model2 as Model
```
#### examples
The examples folder contains several scripts that illustrate how to use spinstats as a python package.
#### documentation
Documentation is provided in the docs folder. To view the docs, open `doc/build/html/index.html` in your browser. Alternatively, you can build the documentation from source by running
```sh
cd doc
make html
```
which requires you to have [sphinx](https://www.sphinx-doc.org/en/master/) installed.

## jupyter notebook
The code can also be accessed through the jupyter notebook contained in the notebooks folder. You can interact with the code directly by launching the notebook in binder using the badge at the top of this page.

## graphical user interface
If writing and editing code is not for you, there is also a graphical user interface (GUI) that you can use to explore the effects of various parameters on the value of &eta;. To use the GUI you will need a python installation on your machine, as well as a few packages. The best way to achieve this is by installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html), then using a virtual environment. In a command prompt or terminal that recognises `conda` and `python` as commands (e.g. the anaconda prompt), you would navigate to the GUI folder
```commandline
cd GUI
```
and install the virtual environment:
```commandline
conda env create -f environment.yml
```
next, activate the environment:
```commandline
conda activate spin-stats-gui
```
and launch the application:
```commandline
python run.py
```

# spin-statistics
Calculations of spin-statistical factors for triplet-triplet annihilation upconversion.

## about
The spin-statistical factor &eta; for triplet-triplet annihilation upconversion gives the probability that an annihilation event between two triplet excitons will result in a singlet exciton. The value &eta; depends in a non-trivial way on a lot of different parameters. For the full story, please refer to [our paper recently uploaded to ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/60dae5946bdaa41e623d4a33). This repository contains code that allows you to perform calculations of &eta;.

If you use any of the code contained in this repository, please cite our paper [Spin statistics for triplet-triplet annihilation upconversion: exchange coupling, intermolecular orientation and reverse intersystem crossing](https://chemrxiv.org/engage/chemrxiv/article-details/60dae5946bdaa41e623d4a33) by Bossanyi et al., recently uploaded to ChemRxiv.

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
Documentation is provided in the docs folder. To view the docs, open docs/build/html/index.html in your browser.

## jupyter notebook
The code can also be accessed through the jupyter notebook contained in the notebooks folder. You can interact with the code directly by launching the notebook in binder using the badge at the top of this page (still need to add this badge!).

## graphical user interface
If writing and editing code is not for you, there is also a graphical user interface (GUI) that you can use to explore the effects of various parameters on the value of &eta;. To use the GUI you will need a python installation on your machine, as well as a few packages. For example, with [Miniconda](https://docs.conda.io/en/latest/miniconda.html), you would need
```sh
conda install numpy pyqt pyqtgraph pandas scipy
```
Note that numpy, scipy and pandas are installed alongside spinstats if you decide to install it as a package on your local machine.

To launch the GUI, in a command prompt that recognises `python` as a command (for example the anaconda prompt), enter
```sh
cd GUI
python run.py
```

# Python Minima Hopping
Python Minima Hopping code for global strucutre optimization. This code can simulate both, cluster and bulk systems. Moreover can any ASE calculator be coupled to the code. For periodic systems variable cell shape MD, softening and optimization is implemented. 

## Installation
Python Minima Hopping requires:
* Python >= 3.7 

To install:
```bash 
git clone https://gitlab.com/goedeckergroup/ase_mh.git
cd ase_mh
pip install .
```

## Documentation & Tutorial
The documentation as well as tutorials for the usage of Minima Hopping can be found on https://python-minima-hopping.readthedocs.io/en/latest/.

## Usage
An example of how to use Python Minima Hopping can be found in the example folder on GitHub. It can be executed by
```bash
python example/mh_na13.py
```
An example of the pre-processing is given in the script:
``` bash
python example/fp_adjust.py
```
Further examples, tutorials as well as description of the input parameters can be found in the code documentation (https://python-minima-hopping.readthedocs.io/en/latest/). 


## References & Citing
1. https://arxiv.org/abs/2309.08418
2. https://doi.org/10.1016/j.jcpx.2023.100131
3. https://doi.org/10.1063/1.1724816


## Authors
* Marco Krummenacher
* Moritz Gubler
* Jonas Finkler
* Hannes Huber
* Martin Sommer-Jörgensen
* Stefan Goedecker






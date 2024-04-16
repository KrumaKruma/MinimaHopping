# Python Minima Hopping
Python Minima Hopping code for global strucutre optimization. This code can simulate both, cluster and bulk systems. Moreover can any ASE calculator be coupled to the code. For periodic systems variable cell shape MD, softening and optimization are implemented. 

## Installation
Python Minima Hopping requires:
* Python >= 3.7

For the usage of the MPI parallelization mpi4py (https://mpi4py.readthedocs.io/en/stable/) must be installed correctly and to generate graphs in the post-processing graphviz and PyGraphviz (https://pygraphviz.github.io/documentation/stable/) have to be installed.  

The code is installed by using Python pip:
```bash 
git clone https://gitlab.com/goedeckergroup/ase_mh.git
cd ase_mh
pip install .
```

Our code is also available on GitHub (https://github.com/KrumaKruma/MinimaHopping.git) and can be installed the same way:
```bash
git clone https://github.com/KrumaKruma/MinimaHopping.git
cd MinimaHopping
pip install .
```

## Documentation & Tutorial
The documentation including a description of the parameters as well as the output and several tutorials for the usage of Minima Hopping can be found on https://python-minima-hopping.readthedocs.io/en/latest/.

## Usage
An example of the pre-processing is given in the script in the example/clusters folder:
``` bash
python mh_na13_preprocess.py
```

An example of how to use Python Minima Hopping can be found in the example/clusters folder on GitHub. It can be executed by
```bash
python mh_na13.py
```

The MPI parallelization can be tested by executing the following script:
```bash
mpirun python mh_na13_mpi.py
```

Note: if you ran first the single thread simulation we recommend to delete the output of these simulation before starting with the MPI parallelized simulation.

Further examples, tutorials as well as description of the input parameters can be found in the code documentation (https://python-minima-hopping.readthedocs.io/en/latest/).
In the example folder are also examples including periodic boundary calculations or running Minima Hopping with two different calculators. 



## References & Citing
1. Reference for this code: https://doi.org/10.1016/j.softx.2024.101632
2. Reference to the original implementation: https://doi.org/10.1063/1.1724816
3. Reference for the local geometry optimizer: https://doi.org/10.1016/j.jcpx.2023.100131
4. Reference for the Overlap Matrix Fingerprint: https://doi.org/10.1063/1.4940026 and https://doi.org/10.1063/1.4828704

## Authors
* Marco Krummenacher
* Moritz Gubler
* Jonas Finkler
* Hannes Huber
* Martin Sommer-Jörgensen
* Stefan Goedecker






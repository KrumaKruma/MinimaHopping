# ASE_MH

Minima Hopping (MH) code in python coupled with the ASE library. This code can simulate both, cluster and bulk systems. The code can be used with any calcualtor which is aviable in ase where a stress tensor is present (if periodic systems are simulated). For periodic systems variable cell shape MD, softening and optimization is implemented. 

## Requirements
* Python 3.7 or later
* ase (atomic simulation envrionment)
* numpy
* scipy

## Usage
Currently the package cannot be installed via pip. However, once the repository is cloned it can be tested by executing the main.py. 
```bash
python main.py
```
This will start a minima hopping run with the LJ potential for the double funnel LJ38 cluster. After the simulation is finished several output files are written. The lowest n (parameter n_poslow in the Minimahopping class) minima are written in a directory minima and can be visualized by using v_sim (or any other software which can visualize .xyz and .ascii files). 
In the min.extxyz file all uniquely found minima are stored and the corresponding fingerprints are stored in fp.dat. The acc.extxyz file stores all the accepted minima and history.dat is a log file which shows all the useful information about the found minima. 

### Restart 
If at least the min.extxyz and the history.dat file is detected in the folder, a restart of the previous run is performed. If the acc.extxyz file contains any structure the last accepted structure is taken as input structure for the restart else the last found minimum in min.extxyz is taken as starting structure. If a fp.dat file is aviable the fingerprints are read else they are calculated in the begining of the restart. 


## ToDo
* Implementation of active learning with Gaussian process (FLARE++)
* setup so that it can be installed with pip

## Contribution
Several people have contributed to this code either with direct python implementations of their work or reference implementations in Fortran:
* Marco Krummenacher (marco.krummenacher@unibas.ch) --> main development
* Moritz Gubler --> variable cell shape [SQNM optimizer](https://github.com/moritzgubler/vc-sqnm)
* Jonas Finkler --> overlap matrix fingerprint
* Hannes Huber --> variable cell shape softening (Fortran reference implementation)
* Martin Sommer --> variable cell shape MD (Fortran reference implementation)






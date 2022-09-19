# ASE_MH

Minima Hopping (MH) code in python coupled with the ASE library. This code can simulate both, cluster and bulk systems. The code can be used with any calcualtor which is aviable in ase where a stress tensor is present (if periodic systems are simulated). For periodic systems variable cell shape MD, softening and optimization is implemented. 


## ToDo
* Implementation of active learning with Gaussian process (FLARE)
* Automatic restart mechanism
* Testing on several NN potential generate with NequIP


## Contribution
Several people have contributed to this code either with direct python implementations of their work or reference implementations in Fortran:
* Marco Krummenacher (marco.krummenacher@unibas.ch) --> main development
* Moritz Gubler --> variable cell shape SQNM optimizer
* Jonas Finkler --> overlap matrix fingerprint
* Hannes Huber --> variable cell shape softening
* Martin Sommer --> variable cell shape MD






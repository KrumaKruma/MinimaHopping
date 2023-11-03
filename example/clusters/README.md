# Example: Minima Hopping with a Na13 cluster

In this example the global optimization of an Na13 cluster is performed using the EAM ASE calculator (https://wiki.fysik.dtu.dk/ase/ase/calculators/eam.html#module-ase.calculators.eam). Please make sure that the MinimaHopping package is installed properly before doing this example. Further information of how to install the package can be found on the GitLab https://gitlab.com/goedeckergroup/ase_mh or in the documentation (https://python-minima-hopping.readthedocs.io/en/latest/).

## Distinguishing Minima
Before running Minima Hopping it is crucial to determine the fingerprint distance of structures which can be considered the same and structures which are different. Several MD trajectories are perfomed and relaxed back to the same local minimum. To do this you can run the following:

```bash
python mh_na13_preprocess.py
```

## Minima Hopping
In this example the Minima Hopping for Na13 is performed. First an initial structure is constructed where the atoms are construcuted as a sting. The global minimum which is to be found is an icosahedral strucutre. This example can be run:
```bash
python mh_na13.py
```

It is also possible to run multiple Minima Hopping processes sharing a single database with MPI. To do that the mpi4py library must be installed correctly on your system. Running the follwing command will perfom such a simulation:
```bash
mpirun python mh_na13_mpi.py
```




# Example for constraint Minima Hopping
It is also possible to include the constraints implemented in ASE in Minima Hopping. In this example two atoms of a Na13 are fixed so that they do not move neither in the MD nor in the geometry optimization. An EAM potential is used parametrized for sodium. The example can be started by
```bash
python mh_constrained.py
```
Because the optimizer checkes if all atoms are moved a warining is printed by the optimizer. In general we recommend to carefully use the constraints and check the output produced by Minima Hopping carefully. Not all constraints provided by ASE and all combinations of them have been tested. 




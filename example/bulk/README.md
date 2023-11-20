# Example: Minima Hopping with a Si4 crystal
In this example the global optimization of a Si4 crystal is performed using the openKim ASE calculator (https://wiki.fysik.dtu.dk/ase/ase/calculators/kim.html#module-ase.calculators.kim). To use this calculator the KIM API package as well as kimpy is required. The codes are hosted in GitHub (https://github.com/openkim) and is explained on https://openkim.org/kim-api/. Please follow the installation instructions on these website to install the calculator properly. For this example any silicon force field can be used, however, in the python script by default the Sillinger-Weber (SW_StillingerWeber_1985_Si__MO_405512056662_005) potential is used. After having installed the KIM API properly the potential can be installed via the following command:

```bash
kim-api-collections-management install user SW_StillingerWeber_1985_Si__MO_405512056662_005
```

Please be aware that we have only tested this example this the given Stillinger-Weber force field. In case you use another silicon potential the result can differ.

## Distinguishing Minima
Before running Minima Hopping it is crucial to determine the fingerprint distance of structures which can be considered the same and structures which are different. Several MD trajectories are perfomed and relaxed back to the same local minimum. To do this you can run the following:

```bash
python mh_si4_preprocess.py
```

## Minima Hopping
In this example the Minima Hopping for Na13 is performed. First an initial structure is constructed where the atoms are construcuted as a sting. The global minimum which is to be found is an icosahedral strucutre. This example can be run:
```bash
python mh_si4.py
```



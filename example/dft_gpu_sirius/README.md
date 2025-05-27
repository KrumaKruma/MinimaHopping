Install the [SIRIUS](https://github.com/electronic-structure/SIRIUS) electronic structure library and enable the Python bindings. At the time of writing, these SPACK install options worked just fine:
``` bash
spack install "sirius@7.5.2%gcc@12.3.0 build_type=Release +python +cuda cuda_arch=86 ^openblas threads=openmp ^openmpi ^libxc@6.2.2"
```

Install the [MPI-enabled SIRIUS ASE calculator](https://github.com/moritzgubler/sirius-python-interface)
``` bash
pip install 'git@https://github.com/moritzgubler/sirius-python-interface'
```

This repository contains and installs a `siriusMH` command line utility for starting Minima Hopping simulations with SIRIUS. Start a minima hopping simulation with two minima hopping workers, each using two processes for DFT calculations and one database process for communication (2 * 2 + 1 = 5 processes in total).
``` bash
mpirun -np 5 python3 siriusMH mhparams.json params.json initial.extxyz --slaveRankSize=2 --numberOfSlaves=2
```
`params.json` contains the DFT options for SIRIUS. All keywords can be found here: [https://github.com/electronic-structure/SIRIUS/blob/develop/src/context/input_schema.json](https://github.com/electronic-structure/SIRIUS/blob/develop/src/context/input_schema.json).
`mhparamns.json` contains all Minima Hopping parameters.

This example runs GPU accelerated Minima Hopping on the DFT level.
# VCS SLAB Example
In this example two copper atoms on a platinum surface are simulated. In x- and y-direction are periodic boundary conditions and not in z-direction. In this example also variable cell shape MD and optimization is allowed in two directions. For those calulations it is important that the calculator has no stress in z-direction. An example of such a calculator can be found in the vcsSlabCalculator which is then combined with ASE's EMT calculator. This example can be started with:
```bash
python vcs_slab_ptcu.py
```

If you look at the output you notice that the cell in z-direction is 1.0. This is due to the definitions used in the variable cell shape geometry optimization. If you are interested in a more detailed insight into the underlying theory we recommend reading the corresponding paper of the geometry optimizer https://doi.org/10.1016/j.jcpx.2023.100131 section *5.1 Transformations of variables for surface boundary conditions*.

## Note
### Calulator
- The calculator presented in this example is tailored to this example and we recommend to carfully create or select your own calculator for variable cell shape slab simulations. 
- Make sure that in z-direction no periodic images are included in the calculation
- The stress can depend on the cell volume. Be aware that the stress is only dependent on the surface and not on the cell volume in that case.

### Minima Hopping
- The variable cell shape slab only works for inputs where the z-direction is non-peridoic. Please ensure that before you start the simulation. Otherwise the simulation will be aborded. 



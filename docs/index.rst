.. Minima Hopping documentation master file, created by
   sphinx-quickstart on Fri Oct 28 16:04:54 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Minima Hopping: a python implementation compatible with the atomic simulation environment!
===========================================================================================
The minimahopping method efficiently finds the global minimum of complex condensed matter systems.
This version of the code is based on the `atomic simulation environment (ASE) <https://wiki.fysik.dtu.dk/ase/>`_ library.
Any calculator implemented in there can be used for the calculation of the energy and forces.

Code development
-----------------
The MH python code is developed through our GitLab repository:

https://gitlab.com/goedeckergroup/ase_mh.git

If you detect any bug please write an issue in this repository. If you can solve the bug yourself or you would like to
add any code please fork the code and make a pull-request to the "develop" branch. Every now and then the develop branch
is pushed to the "master" branch.


Current list of contributors:
------------------------------

* Marco Krummenacher (University of Basel)
* Moritz Gubler (University of Basel)
* Jonas A. Finkler (University of Basel)


.. toctree::
   :maxdepth: 2
   :caption: MINIMA HOPPING DOCUMENTATION:

   GETTING_STARTED/install
   GETTING_STARTED/theory
   GETTING_STARTED/parameters
   GETTING_STARTED/output
   EXAMPLES/clusters
   EXAMPLES/bulk
   EXAMPLES/foundation_models
   EXAMPLES/dft_gpu
   EXAMPLES/md_calculator
   EXAMPLES/slab
   EXAMPLES/constraints

Citing MinimaHopping
--------------------
This implementation

   M. Krummenacher, M. Gubler, J. A. Finkler, H. Huber, M. Sommer-Jörgensen, S.Goedecker, `Performing highly efficient Minima Hopping strucutre predictions using the Atomic Simulation Environment <https://doi.org/10.1016/j.softx.2024.101632>`_, SoftwareX Volume 25, February 2024, 101632


Geometry Optimizer:

   M. Gubler, M. Krummenacher, H. Huber, S. Goedecker, `Efficient variable cell shape geometry optimization <https://doi.org/10.1016/j.jcpx.2023.100131>`_, Journal of Computational Physics: X Volume 17, November 2023, 100131


Minima Hopping method:

   S. Goedecker, `Minima hopping: An efficient search method for the global minimum of the potential energy surface of complex molecular systems <https://aip.scitation.org/doi/10.1063/1.1724816>`_,The Journal of chemical physics 120 (21) (2004) 9911-9917

Fingerprint:

   A. Sadeghi, S. A. Ghasemi, B. Schaefer, S. Mohr, M. A. Lill, S. Goedecker, `Metrics for measuring distances in configuration spaces <https://aip.scitation.org/doi/full/10.1063/1.4828704>`_, The Journal of chemical physics 139 (18) (2013) 184118

   L. Zhu, M. Amsler, T. Fuhrer, B. Schaefer, S. Faraji, S. Rostami, S. A. Ghasemi, A. Sadeghi, M. Grauzinyte, C. Wolverton, et al., `A fingerprint based metric for measuring similarities of crystalline structures <https://aip.scitation.org/doi/full/10.1063/1.4940026>`_, The Journal of chemical physics 144 (3) (2016) 034203.


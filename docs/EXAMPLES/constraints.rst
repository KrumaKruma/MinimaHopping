Example: Constraints
+++++++++++++++++++++++++++++++++++++++
It is also possible to include the constraints implemented in ASE in Minima Hopping.
In this example two atoms of a 13 atoms sodium cluster are fixed so that they do not move neither in the MD nor in the geometry optimization. 
An EAM potential is used parametrized for `sodium <https://wiki.fysik.dtu.dk/ase/ase/calculators/eam.html#module-ase.calculators.eam>`_. 
To use this calculator a parameter file is needed. This file can be downloaded the following way:

.. code-block::

    wget https://www.ctcms.nist.gov/potentials/Download/2016--Nichol-A-Ackland-G-J--Na/3/Na_v2.eam.fs

First the imports needed are done.

.. code-block:: python

    from minimahopping.minhop import Minimahopping
    from ase.cluster.wulff import wulff_construction
    from ase.calculators.eam import EAM
    from ase.constraints import FixAtoms

Next the cluster is constructed. We start with a string of sodium atoms. Next also the EAM calculator is attachted to the atoms object.

.. code-block:: python

    atoms = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.15],
                            size=13, # maximum number of atoms
                            structure='bcc', rounding='above')


    calculator = EAM(potential="Na_v2.eam.fs")
    atoms.calc = calculator

In the next step we use the ASE's constraint method to fix the first two atoms.

.. code-block:: python

    constraints = [FixAtoms(indices=[0,1])]

In the next step the Minima Hopping is set up and run for 50 steps.

.. code-block:: python

    with Minimahopping(atoms, 
                       constraints=constraints, 
                       mdmin=2,
                       initial_step_size=1e-3,
                       verbose_output=True, 
                       T0=1000, 
                       dt0=0.1, 
                       use_MPI=False) as mh:
        mh(totalsteps=50)


.. warning::

    Please be aware that not all constraints and all combination of them implemented in ASE have been tested. Furthermore, the use of constraints in combination with variable cell shape 
    features can lead to problems with the cell in geometry optimization. 

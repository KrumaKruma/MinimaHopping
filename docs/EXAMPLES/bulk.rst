Example: Bulk
++++++++++++++++++++++++++++++++
In this example the global optimization of a Si4 crystal is performed using the `openKim ASE calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/kim.html#module-ase.calculators.kim>`_. 
To use this calculator the KIM API package as well as kimpy is required. The codes are hosted in `GitHub <https://github.com/openkim>`_ and is explained on official `openKim Webpage <https://openkim.org/kim-api/>`_. 
Please follow the installation instructions on these website to install the calculator properly. 
For this example any silicon force field can be used, however, in the python script by default the Sillinger-Weber (SW_StillingerWeber_1985_Si__MO_405512056662_005) potential is used. 
After having installed the KIM API properly the potential can be installed via the following command:
.. code-block:: 

    kim-api-collections-management install user SW_StillingerWeber_1985_Si__MO_405512056662_005


Exercice 1: Distinguishing minima
---------------------------------
The aim of this tutorial is to determine the fingerprint distance of structures which can be considered to be the same
and structures which are different. The idea is to perform various MD trajectories and relax them back to the same
local minimum with a certain force norm. To do this task the package includes a class `adjust_fp` which is doing this
task automatically. First we need all the imports which are required:

.. code-block:: python

    from minimahopping.adjust_fp import adjust_fp
    import logging
    from ase.lattice.cubic import FaceCenteredCubic
    from ase.calculators.kim.kim import KIM

Now we create a structure and set up a calculator. In this example we create a silicon fcc structure and
we set up the openKim calculator:

.. code-block:: python

    logging.INFO
    atoms = FaceCenteredCubic(symbol='Si', latticeconstant=4., size=(1,1,1))
    calculator = KIM("SW_StillingerWeber_1985_Si__MO_405512056662_005")
    atoms.calc = calculator

Afterwards we set up the adjustment class.

.. code-block:: python

    fnrm =  0.001
    adjust = adjust_fp(initial_configuration=atoms, 
                       iterations=10, 
                       T0=100, 
                       dt0=0.01, 
                       mdmin=1, 
                       n_S_orbitals=1, 
                       n_P_orbitals=1, 
                       width_cutoff=4, 
                       fmax=fnrm, 
                       write_graph_output=False)

Now we run 100 mds followed by a geometry optimization.

.. code-block:: python

    outdict = adjust.run()

A dictionairy containing all the information is returned and can be printed the following way:

.. code-block:: python

    msg = "\n=======================FINGERPRINT================================\n"
    msg += '\nMaximal fingerprint distance between the same local minima:\n' + str(fp_max)
    msg += '\n Mean fingerprint distance between the same local minima:\n' + str(fp_mean)
    msg += '\n Standard deviaton of the fingerprint distances:\n' + str(fp_std)
    msg += '\n Suggested minimum threshold (mean + 3 * std):\n' + str(fp_mean + 3 * fp_std)
    msg += "\n==================================================================\n"
    print(msg)

    e_max = outdict['energy']['max']
    e_mean = outdict['energy']['mean']
    e_std = outdict['energy']['std']
    msg = "\n=========================ENERGIES=================================\n"
    msg += '\nMaximal difference between the same local minima:\n' + str(e_max)
    msg += '\n Mean energy difference between the same local minima:\n' + str(e_mean)
    msg += '\n Standard deviaton of the energy differences:\n' + str(e_std)
    msg += '\n Suggested energy threshold (mean + 3 * std):\n' + str(e_mean + 3 * e_std)
    msg += "\n==================================================================\n"
    print(msg)

.. note::
    Be aware that it is very important to use the same parameters for the calculation of the energy and force, the OMFP and the local geometry optimization in the Minima Hopping method.


Exercise 2: Starting Minimahopping
----------------------------------

The aim of this tutorial is to start the minima hopping algorithm with the given default settings. If you want to use
different parameters you can find a detailed description of them :doc:`here <parameters>`. First all the required
classes are imported:

.. code-block:: python

    from ase.lattice.cubic import FaceCenteredCubic
    from ase.calculators.kim.kim import KIM
    from minimahopping.minhop import Minimahopping

Now we create a structure and set up a calculator. As in exercise 1 we create the structure of a 4 atoms Silicon crystal and
we set up the corresponding calculator:

.. code-block:: python

    initial_configuration = FaceCenteredCubic(symbol='Si', latticeconstant=5.25, size=(1,1,1))

In a next step we set up the openKim calculator

.. code-block:: python

    calculator = KIM("SW_StillingerWeber_1985_Si__MO_405512056662_005")
    initial_configuration.calc = calculator

Now we can set up the minima hopping class and run it:

.. code-block:: python

    with Minimahopping(initial_configuration, 
                       verbose_output=True, 
                       T0=2000, 
                       dt0=0.1, 
                       use_MPI=False) as mh:
        mh(totalsteps=50)

The minima hopping algorithm cycles now through 50 escape loops.

.. note::
    If a second calculator is desired this can easily be done by setting up a second md calculator and give it as an argument to the ```MinimaHopping``` class.
    
    .. code-block:: python

        calculator = SOME_ASE_CALCULATOR
        md_calculator = SOME_OTHER_ASE_CALCULATOR

        with Minimahopping(initial_configuration,
                           md_calculator = md_calculator
                           verbose_output=True,
                           T0=2000, 
                           dt0=0.1,
                           use_MPI=False) as mh:

        mh(totalsteps=50)


.. caution::
    Be aware that in case you want to examine periodic systems your calculator needs the stress property included so
    that variable cell shape md and geometry optimization is possible.

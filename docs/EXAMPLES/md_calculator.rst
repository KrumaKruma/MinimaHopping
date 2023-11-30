Example: MD Calculator
+++++++++++++++++++++++++++++++++++++++++++++++++
In this example two calculators are used in the Minima Hopping algorithm. In this example the global minimum of bulk 4 atoms silicon is searched.
The MD and a pre-optimization is performed using a force field using `openKim <https://openkim.org/kim-api/>`_ and a Stillinger-Weber potential is used (SW_StillingerWeber_1985_Si__MO_405512056662_005).
If the openKim api is installed correctly the potential can be installed with

.. code-block::
    
    kim-api-collections-management install user SW_StillingerWeber_1985_Si__MO_405512056662_005

The actual geometry optimization and, hence, the database is performed using `Quantum Espresso <https://www.quantum-espresso.org/>`_. Quantum espresso can easily be installed with sudo.

.. code-block::

    sudo apt install quantum-espresso

However, if you are using a high perfomance computing facility and you would like to perfom long Minima Hopping runs we highly recommend to compile quantum-espresso yourself to get the best performance.
First we import the modules nessecairy for running the example:

.. code-block:: python

    from ase.lattice.cubic import FaceCenteredCubic
    from minimahopping.minhop import Minimahopping
    from ase.calculators.espresso import Espresso
    from ase.calculators.kim.kim import KIM
    from ase.io import read
    from ase.build import bulk
    from ase.io import read, write

Next we read the structures and set up the two calculators.

.. code-block:: python

    # read the initial strucutre
    initial_configuration = read("input.extxyz")

    # give the pseudo potential for silicon
    pseudopotentials = {'Si': 'Si.pbe-nl-rrkjus_psl.1.0.0.UPF'}

    # set parameters for quantum espresso dft calculation
    input_data = {
        'system': {
            'ecutwfc': 60.,
            'ecutrho': 125.,
            'ibrav'  : 0,
            'nosym' : True},
        'disk_io': 'none',
        'electrons': {
            'electron_maxstep': 500,
            'mixing_mode': 'local-TF',
            'mixing_beta': 0.7},
        'rsim':{
            'smear3d': 2.}}

    # set up quantum espresso calculator
    calculator = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(3, 3, 3), koffset=(0, 0, 0),input_data=input_data)
    
    # set up the md calculator in this case openKim with a Stillinger-Weber potential
    md_calculator = KIM("SW_StillingerWeber_1985_Si__MO_405512056662_005")


Next we attach the quantum-espresso calculator to the input structure. The geometry optimization and, hence, the database will consist of quantum-espresso energies.

.. code-block:: python

    initial_configuration.calc = calculator

Next we set up the Minima Hopping and run it for 50 steps. In this example the extra MD calculator is given as an input to Minima Hopping.

.. code-block:: python

    with Minimahopping(initial_configuration, md_calculator=md_calculator, mdmin=2,verbose_output=True, T0=4000, dt0=0.01, use_MPI=False) as mh:
        mh(totalsteps=50)





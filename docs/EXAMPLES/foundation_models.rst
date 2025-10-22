Example: Foundation Models
++++++++++++++++++++++++++
Machine learned interatomic potentials can be used to predict energies and forces of atomic structures with almost DFT accuracy at a fraction of the computational cost. 
Recently, so called foundation models have become available which are trained on large datasets covering large parts of the periodic table. 
These models can serve as a great entry point to structure search with the Minima Hopping method. 

In this example we will use the `GRACE potential <https://doi.org/10.1103/PhysRevX.14.021036>`_. 
An overview of available foundation models can be found `here <https://gracemaker.readthedocs.io/en/latest/gracemaker/foundation/>`_, but you are of course free to use your favorite foundation model. 

After installing the required packages, we simply load the calculator: 

.. code-block:: python

   from tensorpotential.calculator import grace_fm
   calculator = grace_fm('GRACE-2L-OAM')

An example script to run a Minima Hopping simulation is provided under `examples/foundation_models`.
A small script to generate an input structure from a given composition is also provided. 
It simply creates places atoms in a cell while avoiding too close overlap. 
Running Minima Hopping for a new composition is thus as easy as entering the right composition and starting the script. 
The provided example usually finds the perovskite structure of :math:`\mathrm{CaTiO_3}` within a few hours using a single GPU.
The provided example also shows how to use MPI parallelization with multiple GPUs.

.. code-block:: python

   from minimahopping.minhop import Minimahopping
   from make_structure import get_random_packed
   from tensorpotential.calculator import grace_fm

   composition = f'(CaTiO3)4'
   init_structure = get_random_packed(composition, scale=1.5)
   init_structure.calc = grace_fm('GRACE-2L-OAM')
   with Minimahopping(init_structure, 
                      symprec=1e-5,
                      verbose_output=False, 
                      T0=500, 
                      dt0=0.1, 
                      use_MPI=False, 
                      mdmin=3,
                      fixed_cell_simulation=False,
                      fmax=0.002,
                      energy_threshold=0.002,
                      fingerprint_threshold=0.005,
                      write_graph_output=False,
                      alpha_reject=1.05,
                      alpha_accept=0.7,
                      beta_increase=1.1,
                      beta_decrease=0.8,
                      enhanced_feedback=True,
                      collect_md_data=False,
     ) as mh:
     mh(totalsteps=100000)



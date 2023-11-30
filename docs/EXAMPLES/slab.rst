
Example: Slab Boundary Conditions
+++++++++++++++++++++++++++++++++++++++++++++++++++
This example presents two ways how to deal with slab boundary conditions or surfaces. 
In the first case a fixed unit cell is introduced where the z-direction pbc is set False and the cell vector in z-direction is very large so that not periodic image is taken into account in the calculator.
In the second case also variable cell shape features are activated and the cell is moved in x-y-direction. 

Example 1: Fixed cell slab
--------------------------
In this example we simulate a platinum surface with two copper atoms attached to it. 
The underlying energy and force calculation is performed using ASE's EMT calculator.
First all the nessecairy modules are imported:

.. code-block:: python

    from ase import Atom, Atoms
    from ase.build import fcc110
    from ase.calculators.emt import EMT
    from minimahopping.minhop import Minimahopping

In a next step the surface with the two copper atoms on it is constructed.

.. code-block:: python

    # Make the Pt 110 slab.
    atoms = fcc110('Pt', (2, 2, 2), vacuum=7.)

    # Add the Cu2 adsorbate.
    adsorbate = Atoms([Atom('Cu', atoms[7].position + (0., 0., 2.5)),
                   Atom('Cu', atoms[7].position + (0., 0., 5.0))])
    atoms.extend(adsorbate)

.. note::
    Please make sure that the cell in z-direction is large enough so that the potential is not interacting with its periodic images in this direction. 
    Furthermore, also the periodic boundary condition is z-direction should be False.

In the last step the calculator is attached to the atoms object and the Minima Hopping is started.

.. code-block:: python

    # Set the calculator.
    atoms.calc = EMT()

    with Minimahopping(atoms, mdmin=3, fixed_cell_simulation = True,verbose_output=True, T0=500, dt0=0.01, use_MPI=False) as mh:
        mh(totalsteps=50)

Please be aware that the `fixed_cell_simulation` parameter is set to True. This parameter is fixing the cell and hence no variable cell shape features will be applied. 

.. note::
    In case the unit cell of the global minimum of a periodic system is known in advance or a the global minimum of a specific cell shape is searched for, it is possible to use this feature also for full periodic boundary conditions.


Example 2: Variable cell shape slab 
-----------------------------------
It is also possible to run variable cell shape surface simulations where the cell is only moved in x-y-direction. 
To do that the cell should be formed in the follwoing way:
    
    .. math::

        M = \begin{bmatrix}
                a & b & 0 \\
                c & d & 0 \\
                0 & 0 & 1
            \end{bmatrix}

In general there is no further adjustment needed regarding the Minima Hopping to perform variable cell shape surface simulations. 
In the follwoing example again two copper atoms on a platinum surface are simulated using an EMT calculator.
First all the nessecairy imports are done.

.. code-block:: python

    from ase import Atom, Atoms
    from ase.build import fcc110
    from ase.calculators.emt import EMT
    from vcsSlabCalculator import SlabCalculator
    from minimahopping.minhop import Minimahopping

In a next step the surface with the two copper atoms on it is constructed.

.. code-block:: python

    # Make the Pt 110 slab.
    atoms = fcc110('Pt', (2, 2, 2), vacuum=7.)

    # Add the Cu2 adsorbate.
    adsorbate = Atoms([Atom('Cu', atoms[7].position + (0., 0., 2.5)),
                   Atom('Cu', atoms[7].position + (0., 0., 5.0))])
    atoms.extend(adsorbate)

Then the cell is changed to the right form 

.. code-block:: python

    cell = atoms.get_cell()
    cell[3,3] = 1.0
    atoms.set_cell(cell)

The other vecotrs are in this case already in the right shape and ,therefore, only the z-direction has to be adjusted.
Next the calculator is set up and the Minima Hopping is started.

.. code-block:: python

    # Set the calculator.
    calc = SlabCalculator(EMT())
    atoms.calc = calc

    with Minimahopping(atoms, mdmin=3,verbose_output=True, T0=500, dt0=0.01, use_MPI=False) as mh:
        mh(totalsteps=50)

.. caution::
    Here we have implemented a slab calculator for this example.
    In case you aim to do similar simulations carefully choose your calculator. 
    It is very important that there is no stress in the z-direction given back from the calculator.
    Please make sure that your calculator is adjusted accordingly.
    Otherwise the cell might also be moved in z-direction during the simulation.





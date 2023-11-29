from ase import Atom, Atoms
from ase.build import fcc110
from ase.calculators.emt import EMT
from vcsSlabCalculator import SlabCalculator
from minimahopping.minhop import Minimahopping


def main():
    # Make the Pt 110 slab.
    initial_configuration = fcc110('Pt', (2, 2, 2), vacuum=7.)

    # Add the Cu2 adsorbate.
    adsorbate = Atoms([Atom('Cu', initial_configuration[7].position + (0., 0., 2.5)),
                   Atom('Cu', initial_configuration[7].position + (0., 0., 5.0))])
    initial_configuration.extend(adsorbate)

    # Set the calculator.
    calc = SlabCalculator(EMT())
    initial_configuration.calc = calc


    with Minimahopping(initial_configuration, mdmin=3,verbose_output=True, T0=2000, dt0=0.01, use_MPI=False) as mh:
        mh(totalsteps=50)


if __name__ == '__main__':
    main()
    quit()

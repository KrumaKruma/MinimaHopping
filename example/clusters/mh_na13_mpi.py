#!/usr/bin/env python3
from ase.calculators.eam import EAM
from minimahopping.minhop import Minimahopping
from ase.cluster.wulff import wulff_construction


def main():
    # construct a chain with 13 atoms:
    # initial_configuration = Icosahedron('Na', 2, latticeconstant=None)

    # initial_configuration = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.25],
    #                         size=55, # maximum number of atoms
    #                         structure='bcc', rounding='above')

    initial_configuration = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.15],
                            size=13, # maximum number of atoms
                            structure='bcc', rounding='above')


    calculator = EAM(potential="Na_v2.eam.fs")
    initial_configuration.calc = calculator
    with Minimahopping(initial_configuration, verbose_output=True, T0=2000, dt0=0.1, use_MPI=True) as mh:
        mh(totalsteps=50)


if __name__ == '__main__':
    main()
    quit()


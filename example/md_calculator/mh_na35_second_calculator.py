#!/usr/bin/env python3
from ase.calculators.eam import EAM
from ase.calculators.lj import LennardJones
from minimahopping.minhop import Minimahopping
from ase.cluster.wulff import wulff_construction

def main():

    initial_configuration = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 2, 0),(0, 0, 1)], energies=[0.501, 0.501, 0.5015],
                            size=13, # maximum number of atoms
                            structure='bcc', rounding='above')

    calculator = EAM(potential="Na_v2.eam.fs")
    calculator2 = LennardJones(sigma=3.5)
    initial_configuration.calc = calculator
    with Minimahopping(initial_configuration, calculator2,verbose_output=True, mdmin=10,T0=2000, dt0=0.1, use_MPI=False) as mh:

        mh(totalsteps=500)


if __name__ == '__main__':
    main()
    quit()


#!/usr/bin/env python3
from ase.lattice.cubic import FaceCenteredCubic
from ase.calculators.kim.kim import KIM
from minimahopping.minhop import Minimahopping

def main():
    
    initial_configuration = FaceCenteredCubic(symbol='Si', latticeconstant=5.25, size=(1,1,1))
    calculator = KIM("SW_StillingerWeber_1985_Si__MO_405512056662_005")
    initial_configuration.calc = calculator

    with Minimahopping(initial_configuration, verbose_output=True, T0=2000, dt0=0.1, use_MPI=False) as mh:
    # or using mpi in the minima hopping simulation:
    # with Minimahopping(initial_configuration, verbose_output=True, T0=2000, dt=0.1, use_MPI=True) as mh:
        mh(totalsteps=50)


if __name__ == '__main__':
    main()
    quit()


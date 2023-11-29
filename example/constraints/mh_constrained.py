from minimahopping.minhop import Minimahopping
from ase.cluster.wulff import wulff_construction
from ase.calculators.eam import EAM
from ase.constraints import FixAtoms


def main():
    # construct a chain with 13 atoms:
    initial_configuration = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.15],
                            size=13, # maximum number of atoms
                            structure='bcc', rounding='above')

    # set up eam calculator and attach it to atoms object
    calculator = EAM(potential="Na_v2.eam.fs")
    initial_configuration.calc = calculator

    # set up constraints and fixing atom 0 and 1
    constraints = [FixAtoms(indices=[0,1])]

    initial_configuration.calc = calculator

    with Minimahopping(initial_configuration, mdmin=2, constraints=constraints,initial_step_size=1e-3,verbose_output=True, T0=1000, dt0=0.1, use_MPI=False) as mh:
        mh(totalsteps=50)

if __name__ == '__main__':
    main()
    quit()

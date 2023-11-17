from minimahopping.minhop import Minimahopping
from ase.cluster.wulff import wulff_construction
from ase.calculators.eam import EAM
from ase.constraints import FixAtoms


def main():
    atoms = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.15],
                            size=13, # maximum number of atoms
                            structure='bcc', rounding='above')


    calculator = EAM(potential="Na_v2.eam.fs")
    atoms.calc = calculator
    constraints = [FixAtoms(indices=[0,1])]

    atoms.calc = calculator

    with Minimahopping(atoms, mdmin=2, constraints=constraints,initial_step_size=1e-3,verbose_output=True, T0=1000, dt0=0.1, use_MPI=False) as mh:
        mh(totalsteps=50)

if __name__ == '__main__':
    main()
    quit()


from ase.calculators.eam import EAM
from mh import Minimahopping
from ase.cluster.wulff import wulff_construction
from ase.cluster import Icosahedron


def main():
    # construct a chain with 13 atoms:
    # atoms = Icosahedron('Na', 2, latticeconstant=None)

    atoms = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.15],
                            size=13, # maximum number of atoms
                            structure='bcc', rounding='above')

    calculator = EAM(potential="Na_v2.eam.fs")
    atoms.calc = calculator

    mh = Minimahopping(atoms, verbose=True)
    mh(totalsteps=1000)






if __name__ == '__main__':
    main()


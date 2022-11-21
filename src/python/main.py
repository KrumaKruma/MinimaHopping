from ase.io import read, write
# from ase.calculators.lj import LennardJones
from mh import Minimahopping
from ase.calculators.eam import EAM
from eam_sym_calc import EAMSymmetryCalculator


def main():
    filename = "./posinp.xyz"
    atoms = read(filename)
    calculator = EAM(potential="/kernph/hubhan00/hannes-python/ASE_git/ase_mh/src/python/Na_v2.eam.fs")
    atoms.calc = calculator

    #calc2 = EAM(potential="Na_v2.eam.fs")# EAMSymmetryCalculator(width_cutoff=5.0)
    calc2 = EAMSymmetryCalculator(width_cutoff=4.0)

    mh = Minimahopping(atoms, calc2, verbose=True)
    mh(totalsteps=5000)



if __name__ == '__main__':
    main()

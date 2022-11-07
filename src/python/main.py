from ase.io import read, write
from ase.calculators.lj import LennardJones
from mh import Minimahopping
from bazant_calc import BazantCalculator
from bazant_sym_calc import BazantSymmetryCalculator


def main():
    filename = "../../data/Si_in_non_gm.extxyz"
    atoms = read(filename)
    calculator = BazantCalculator()
    atoms.calc = calculator

    calc2 = BazantSymmetryCalculator()

    mh = Minimahopping(atoms, calc2, verbose=True)
    mh(totalsteps=100)






if __name__ == '__main__':
    main()

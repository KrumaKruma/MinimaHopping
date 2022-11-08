from ase.io import read, write
from ase.calculators.eam import EAM
from mh import Minimahopping
from bazant_calc import BazantCalculator
def main():
    filename = "../../data/SiC_in.extxyz"
    atoms = read(filename)
    # calculator = LennardJones()
    # calculator.parameters.epsilon = 1.0
    # calculator.parameters.sigma = 1.0
    # calculator.parameters.rc = 6.0
    #calculator = EAM(potential="Na_v2.eam.fs")
    calculator = BazantCalculator()
    atoms.calc = calculator

    mh = Minimahopping(atoms, verbose=True)
    mh(totalsteps=100)






if __name__ == '__main__':
    main()


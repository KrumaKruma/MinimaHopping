from ase.io import read, write
from ase.calculators.lj import LennardJones
from mh import Minimahopping


def main():
    filename = "../../data/LJ38.xyz"
    atoms = read(filename)
    calculator = LennardJones()
    calculator.parameters.epsilon = 1.0
    calculator.parameters.sigma = 1.0
    calculator.parameters.rc = 6.0
    atoms.calc = calculator

    mh = Minimahopping(atoms, verbose=True)
    mh(totalsteps=100)






if __name__ == '__main__':
    main()


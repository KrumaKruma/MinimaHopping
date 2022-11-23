from ase.io import read, write
# from ase.calculators.lj import LennardJones
from mh import Minimahopping
from ase.calculators.eam import EAM
from eam_sym_calc import EAMSymmetryCalculator

import sys
sys.path.append('.') # append local directory to module directory paths
import params


def main():
    # print(params.bias_parameter_list)
    filename = "./posinp.xyz"
    atoms = read(filename)
    potential_file = '/kernph/hubhan00/hannes-python/ASE_git/ase_mh/src/python/Na_v2.eam.fs'
    calculator = EAM(potential = potential_file)
    atoms.calc = calculator

    if (params.use_bias):
        print("Biased MH is used with params:", params.bias_parameter_list)
        calc2 = EAMSymmetryCalculator(**params.bias_parameter_list)
    else:
        print("Non biased MH is used")
        calc2 = EAM(potential = potential_file)

    mh = Minimahopping(atoms, calc2, verbose=True)
    mh(totalsteps = params.totalsteps) # **params.mh_parameter_list ?Bug?



if __name__ == '__main__':
    main()

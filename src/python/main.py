
from ase.io import read, write
# from ase.calculators.lj import LennardJones
from mh import Minimahopping
from ase.calculators.eam import EAM
from eam_sym_calc import EAMSymmetryCalculator
# from ase.cluster.wulff import wulff_construction
# from ase.cluster import Icosahedron
import sys
sys.path.append('.') # append local directory to module directory paths
import params



def main():

    # print(params.use_bias)
    # print(params.totalsteps)
    # print(params.mh_parameter_list)
    # print(params.bias_parameter_list)
    # construct a chain with 13 atoms:
    # atoms = Icosahedron('Na', 2, latticeconstant=None)

    # atoms = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.15],
    #                         size=13, # maximum number of atoms
    #                         structure='bcc', rounding='above')

    # print(params.bias_parameter_list)
    filename = "./posinp.xyz"
    atoms = read(filename)
    potential_file = '/kernph/hubhan00/hannes-python/ASE_git/ase_mh/src/python/Na_v2.eam.fs'
    calculator = EAM(potential = potential_file)
    atoms.calc = calculator

    if (params.use_bias):
        print("Biased MH is used with params:", params.bias_parameter_list)
        calculator2 = EAMSymmetryCalculator(**params.bias_parameter_list)
    else:
        print("Non biased MH is used")
        calculator2 = EAM(potential = potential_file)


    mh = Minimahopping(atoms, calc2=calculator2, switch_calc=params.use_bias, verbose=True, **params.mh_parameter_list) # **params.mh_parameter_list
    mh(totalsteps = params.totalsteps)



if __name__ == '__main__':
    main()

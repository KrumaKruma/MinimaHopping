#!/usr/bin/env python3
from minimahopping.adjust_fp import adjust_fp
import logging
from ase.lattice.cubic import FaceCenteredCubic
from ase.calculators.kim.kim import KIM

def main():

    logging.INFO
    
    
    
    atoms = FaceCenteredCubic(symbol='Si', latticeconstant=4., size=(1,1,1))
    calculator = KIM("SW_StillingerWeber_1985_Si__MO_405512056662_005")
    atoms.calc = calculator


    fnrm =  0.001
    adjust = adjust_fp(initial_configuration=atoms, iterations=10, T0=100, dt0=0.01, mdmin=1, n_S_orbitals=1, n_P_orbitals=1, width_cutoff=4, fmax=fnrm, write_graph_output=False)

    outdict = adjust.run()

    fp_max = outdict['fp']['max']
    fp_mean = outdict['fp']['mean']
    fp_std = outdict['fp']['std']
    
    msg = "\n=======================FINGERPRINT================================\n"
    msg += '\nMaximal fingerprint distance between the same local minima:\n' + str(fp_max)
    msg += '\n Mean fingerprint distance between the same local minima:\n' + str(fp_mean)
    msg += '\n Standard deviaton of the fingerprint distances:\n' + str(fp_std)
    msg += '\n Suggested minimum threshold (mean + 3 * std):\n' + str(fp_mean + 3 * fp_std)
    msg += "\n==================================================================\n"
    print(msg)

    e_max = outdict['energy']['max']
    e_mean = outdict['energy']['mean']
    e_std = outdict['energy']['std']
    msg = "\n=========================ENERGIES=================================\n"
    msg += '\nMaximal difference between the same local minima:\n' + str(e_max)
    msg += '\n Mean energy difference between the same local minima:\n' + str(e_mean)
    msg += '\n Standard deviaton of the energy differences:\n' + str(e_std)
    msg += '\n Suggested energy threshold (mean + 3 * std):\n' + str(e_mean + 3 * e_std)
    msg += "\n==================================================================\n"
    print(msg)




if __name__  == '__main__':
    main()






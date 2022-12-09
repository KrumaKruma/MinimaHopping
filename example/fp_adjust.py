from ase.calculators.eam import EAM
from minimahopping.adjust_fp import adjust_fp
from ase.cluster.wulff import wulff_construction
from ase.cluster import Icosahedron
from ase.io import read, write

def main():

    atoms = Icosahedron('Na', 3, latticeconstant=None)
    calculator = EAM(potential="Na_v2.eam.fs")
    atoms.calc = calculator


    fnrm =  0.001
    adjust = adjust_fp(fmax=fnrm, 
                    iterations=100, 
                    temperature=500, 
                    dt=0.1, 
                    md_min=1, 
                    ns_orb=1, 
                    np_orb=1, 
                    width_cutoff=3.5, 
                    exclude=[])

    outdict = adjust.run(atoms=atoms)

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






import numpy as np
import sympenlib
from ase.io import read, write
from ase.calculators.emt import EMT

def symmetry_bias(atoms, width_cutoff=4.0, natx_sphere=50, nums=1, nump=1, lengthfp=200, num_cat=1, nex_cutoff=2):
    import numpy as np
    # Calculation of Symmetry bias via f2py
    # subroutine "sympen_f2py_interface" needs atomic units as input => transform input angstroem to bohr
    # subroutine "sympen_f2py_interface" gives atomic units as output => transform output to si units
    Ha_eV = 27.211399
    Bohr_Ang = 0.529177249
    Ang_Bohr = 1.0/0.529177249

    nat = len(atoms)
    positions = atoms.get_positions().T * Ang_Bohr
    atom_element_number_array = atoms.get_atomic_numbers()


    # get lattice:
    pbcs = atoms.pbc
    pbc = set(pbcs)
    pbc = {True}
    assert len(pbc) == 1, "mixed boundery conditions in symmetry penalty"
    if True in pbc:
        lattice = atoms.get_cell().T * Ang_Bohr #3x3 for pbc else array of 3
    else:
        # in case of non PBC => big cell => no ghost atoms in atom sphere
        lattice = np.array( [[100.0, 0.0, 0.0 ], [0.0, 100.0, 0.0 ], [0.0, 0.0, 100.0 ]] )

    alat = lattice
    rxyz = positions

    nat_sphere_current_max, dpenaldr, penalty, dpenaldalat = sympenlib.sympen_f2py_interface(alat,rxyz,atom_element_number_array,natx_sphere,nums,nump,width_cutoff,num_cat,nex_cutoff,lengthfp,nat)


    bias = penalty * Ha_eV
    dbiasdr = dpenaldr.T * (Ha_eV/Bohr_Ang)
    dbiasdalat = dpenaldalat * (Ha_eV/Bohr_Ang)

    # print(nat_sphere_current_max, bias, "nat_sphere, bias(pweight=1.0)")

    return bias, dbiasdr, dbiasdalat, nat_sphere_current_max



def add_bias_2_pes(atoms, e_pot, forces, deralat, pweight, bias, dbiasdr, dbiasdalat):

    positions = atoms.get_positions()#*Ang_Bohr
    nat = positions.shape[0]

    forces = forces.T
    dbiasdr = dbiasdr.T
    # deralat = deralat.T
    # dbiasdalat = dbiasdalat.T

    b_e_pot = e_pot + pweight * bias

    #add biased forces to forces:
    for iat in range(nat):
        forces[0,iat] = forces[0,iat] - pweight * dbiasdr[0,iat]
        forces[1,iat] = forces[1,iat] - pweight * dbiasdr[1,iat]
        forces[2,iat] = forces[2,iat] - pweight * dbiasdr[2,iat]
        # print(dbiasdr.item((0,iat)))
        # print(dbiasdr.item((1,iat)))
        # print(dbiasdr.item((2,iat)))


    b_deralat = deralat + pweight * dbiasdalat

    b_forces = forces.T
    # b_deralat = deralat.T

    return b_e_pot, b_forces, b_deralat




def main():
    # filename = "Si_in2.extxyz"
    # filename = "N_copy.extxyz"
    # filename = "goal_struc_1.ascii"
    filename = "poscur.ascii"

    # Read local minimum input file
    atoms = read(filename)

    # positions = atoms.get_positions()#*Ang_Bohr
    # cell = atoms.get_cell(complete=False).T#*Ang_Bohr
    # nat = positions.shape[0]

    # e_pot, forces, deralat, stress_tensor = bazant.energyandforces_bazant(cell,positions.T,nat)


    bias, dbiasdr, dbiasdalat = symmetry_bias(atoms)

    print(bias)

    pweight = 1.0

    # atoms.calc = EMT()
    # e_pot = atoms.get_potential_energy()
    # forces = atoms.get_forces()
    # deralat = np.zeros((3, 3))
    # b_e_pot, b_forces, b_deralat = add_bias_2_pes(atoms, e_pot, forces, deralat, pweight, bias, dbiasdr, dbiasdalat)


if __name__ == '__main__':
    main()

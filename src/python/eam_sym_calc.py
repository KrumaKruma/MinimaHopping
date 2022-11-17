import numpy as np
from ase.calculators.calculator import Calculator, all_changes
# import EMT
from ase.calculators.eam import EAM
import sym_bias_wrapper
from time import perf_counter

class EAMSymmetryCalculator(EAM):


    implemented_properties = ['energy', 'forces', 'free_energy']
    implemented_properties += ['stress']  # bulk properties
    nolabel = True

    # def __init__(self,):
    def __init__(self, pweight = 1.0, width_cutoff = 4.5, natx_sphere = 50, nums = 1, nump = 1, lengthfp = 200, num_cat = 1, nex_cutoff = 2):
        self.pweight = pweight
        self.width_cutoff = width_cutoff
        self.natx_sphere = natx_sphere
        self.nums = nums
        self.nump = nump
        self.lengthfp = lengthfp
        self.num_cat = num_cat
        self.nex_cutoff = nex_cutoff
        EAM.__init__(self, potential="Na_v2.eam.fs")


    def lattice_derivative(self, stress_tensor, cell):
        """
        Calculation of the lattice derivative from the stress tensor. This function cannot be used or has to be changed
        if the stress tensor is not included in the calculator used
        Input:
            stress_tensor: stress tensor from ase atoms object
                stress tensor from ase atoms object (atoms.get_stress(voigt=False,apply_constraint=False))
            cell: cell from ase atoms object (atoms.get_cell(complete=False))
        Return:
            deralat: np array
                numpy array containing the lattice derivatives
        """
        assert stress_tensor.shape == (3,3), 'Stress tensor is not a 3x3 array'
        assert cell.shape == (3,3), 'Cell is not a 3x3 array'

        inv_cell = np.linalg.inv(cell)
        prefact = np.linalg.det(cell)
        deralat = - prefact * np.matmul(stress_tensor, inv_cell)
        return deralat

    def calculate(
            self,
            atoms=None,
            properties=None,
            system_changes=all_changes,
    ):

        if properties is None:
            properties = self.implemented_properties




        #Set base_calculator:
        EAM.calculate(self, atoms, properties, system_changes)
        self.calc = EAM(potential="Na_v2.eam.fs")
        self.atoms.calc = self.calc

        #Get energy, forces, deralat from calculator
        time1 = perf_counter()
        energy = self.atoms.get_potential_energy()
        forces = self.atoms.get_forces()
        time2 = perf_counter()


        # no lattice, no stress, no deralat
        if self.atoms.cell.rank == 3:
            stress = self.atoms.get_stress(voigt=False,apply_constraint=False)
            cell = self.atoms.get_cell(complete=False)
            deralat = self.lattice_derivative(stress, cell)
        # large cell, zero deralat
        else:
            cell = np.array( [[100.0, 0.0, 0.0 ], [0.0, 100.0, 0.0 ], [0.0, 0.0, 100.0 ]] )
            deralat = np.array( [[0.0, 0.0, 0.0 ], [0.0, 0.0, 0.0 ], [0.0, 0.0, 0.0 ]] )


        time3 = perf_counter()
        bias, dbiasdr, dbiasdalat = sym_bias_wrapper.symmetry_bias(atoms, width_cutoff = self.width_cutoff, natx_sphere = self.natx_sphere, nums = self.nums, nump = self.nump, lengthfp = self.lengthfp, num_cat = self.num_cat, nex_cutoff = self.nex_cutoff)
        # b_e_pot = bias
        # b_forces = -dbiasdr
        # b_deralat = dbiasdalat
        b_e_pot, b_forces, b_deralat = sym_bias_wrapper.add_bias_2_pes(atoms, energy, forces, deralat, self.pweight, bias, dbiasdr, dbiasdalat)
        time4 = perf_counter()

        # print(bias, energy, time2 - time1, time4 - time3)

        # no lattice, no stress
        if self.atoms.cell.rank == 3:
            b_stress = - np.matmul(b_deralat, cell) / np.linalg.det(cell)
            stress_out = np.zeros((6,))
            stress_out[0] = b_stress[0,0]
            stress_out[1] = b_stress[1,1]
            stress_out[2] = b_stress[2,2]
            stress_out[3] = b_stress[1,2]
            stress_out[4] = b_stress[0,2]
            stress_out[5] = b_stress[0,1]
            self.results['stress'] = stress_out

        self.results['energy'] = b_e_pot

        self.results['free_energy'] = b_e_pot

        self.results['forces'] = b_forces #.T

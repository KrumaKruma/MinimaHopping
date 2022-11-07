import numpy as np
from ase.calculators.calculator import Calculator, all_changes
import bazant

class BazantCalculator(Calculator):


    implemented_properties = ['energy', 'forces', 'free_energy']
    implemented_properties += ['stress']  # bulk properties
    nolabel = True

    def __init__(self,):

        Calculator.__init__(self,)


    def calculate(
            self,
            atoms=None,
            properties=None,
            system_changes=all_changes,
    ):

        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)
        positions = self.atoms.positions
        cell = self.atoms.cell

        energy, forces, deralat, stress_intern = bazant.energyandforces_bazant(cell.T, positions.T, natoms)
        stress = - np.matmul(deralat,cell) / np.linalg.det(cell)

        # no lattice, no stress
        if self.atoms.cell.rank == 3:
            stress_out = np.zeros((6,))
            stress_out[0] = stress[0,0]
            stress_out[1] = stress[1,1]
            stress_out[2] = stress[2,2]
            stress_out[3] = stress[1,2]
            stress_out[4] = stress[0,2]
            stress_out[5] = stress[0,1]
            self.results['stress'] = stress_out

        self.results['energy'] = energy

        self.results['free_energy'] = energy

        self.results['forces'] = forces.T



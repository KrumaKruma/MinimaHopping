import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.emt import EMT
from ase import atoms


class SlabCalculator(Calculator):

    implemented_properties = ['energy', 'forces', 'stress', 'bandgap', 'fermienergy', 'chargedensity', 'chargedensityandgrid', 'charges']
    def __init__(self, original_calculator) -> None:
        self.original_calculator = original_calculator
        super().__init__()

    def calculate(
    self,
    atoms: atoms.Atoms = None,
    properties = None,
    system_changes = all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties
        calc = atoms.calc
        # print(system_changese properties)
        atoms.calc = self.original_calculator
        # super().calculate(atoms, properties, system_changes)
        orig_cell = atoms.cell.copy()
        atoms.cell[2,2] = 100
        volume = atoms.get_volume()
        surface = np.linalg.norm(np.cross(atoms.get_cell()[0,:], atoms.get_cell()[1,:]))
        #print(surface, volume)
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        stress = (atoms.get_stress(voigt=False)*volume)/surface

        super().calculate(atoms, properties, system_changes)
                # no lattice, no stress
        for i in np.where(atoms.pbc==False):
            ind = i[0]
            stress[ind,:] = 0.
            stress[:,ind] = 0.

        stress_out = np.zeros((6,))
        stress_out[0] = stress[0,0]
        stress_out[1] = stress[1,1]
        stress_out[2] = stress[2,2]
        stress_out[3] = stress[1,2]
        stress_out[4] = stress[0,2]
        stress_out[5] = stress[0,1]
        self.results['stress'] = stress_out
        #print("STRESS:   ", stress)
        self.results['energy'] = energy

        self.results['free_energy'] = energy

        self.results['forces'] = forces

        atoms.cell = orig_cell
        atoms.calc = calc
        
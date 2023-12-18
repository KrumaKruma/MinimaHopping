from ase.calculators.calculator import Calculator, all_changes, all_properties
from ase import atoms
import numpy as np
from ase import units

class EnthalpyCalculator(Calculator):
    """
    Calculates p*V and adds it to the potential energy of the base calculator. The stress is also adjusted. Forces are also implemented.
    >>> atoms.calc = EnthalpyCalculator(normalCalculator, pressure_giga_pascale = 5)
    The new calculator will now calculate H=E+pV and adjust the stress accordingly.
    """
    implemented_properties = ['energy', 'forces', 'stress']
    default_parameters = {}
    nolabel = True

    def __init__(self, base_calculator: Calculator, pressure_giga_pascale: float = 0.0):
        super().__init__()
        self.base_calculator = base_calculator
        # self.pressure is in ase units.
        self.pressure = pressure_giga_pascale * units.GPa
        

    def calculate(
        self,
        atoms: atoms.Atoms = None,
        properties = None,
        system_changes = all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties

        super().calculate(atoms, properties, system_changes)
        if False in atoms.pbc:
            print('External pressure can only be applied when atoms.pbc = [True, True, True.]')
            print('atoms.pbc currently holds: ', atoms.pbc)

            if 'stress' in properties:
                self.results['stress'] = np.zeros(6)
            if 'energy' in properties:
                self.results['energy'] = self.base_calculator.get_potential_energy(atoms)
            if 'forces' in properties:
                self.results['forces'] = self.base_calculator.get_forces(atoms)

        else:
            if 'stress' in properties:
                self.results['stress'] = self.base_calculator.get_stress(atoms) + np.array([self.pressure, self.pressure, self.pressure, .0, .0, .0])
            if 'energy' in properties:
                self.results['energy'] = self.base_calculator.get_potential_energy(atoms) +  self.pressure * atoms.get_volume()
            if 'forces' in properties:
                self.results['forces'] = self.base_calculator.get_forces(atoms)

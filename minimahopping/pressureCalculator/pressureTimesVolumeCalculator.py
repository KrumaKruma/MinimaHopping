from ase.calculators.calculator import Calculator, all_changes, all_properties
from ase import atoms
import numpy as np
from ase import units

class SIRIUS(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']
    default_parameters = {}
    nolabel = True

    def __init__(self, pressure_giga_pascale: float = 0.0):

        super().__init__()
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
                self.results['energy'] = 0.0
            if 'forces' in properties:
                self.results['forces'] = np.zeros((len(atoms, 3)))


        elif not system_changes == []:
            if 'stress' in properties:
                self.results['stress'] = np.array([self.pressure, self.pressure, self.pressure, .0, .0, .0])
            if 'energy' in properties:
                self.results['energy'] = self.pressure * atoms.get_volume()
            if 'forces' in properties:
                self.results['forces'] = np.zeros((len(atoms, 3)))

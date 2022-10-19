import numpy as np
from ase import units

class Cell_atom:
    def __init__(self, positions, mass=None, velocities=None):
        self.positions = positions
        self.masses = np.array([mass, mass, mass])
        self.velocities = velocities

    def set_velocities_boltzmann(self, temperature):
        """
        Set the velocity of the cell atoms for the MD part accorting to a temperature and the boltzmann distribution
        Input:
            temperature: float
                Temperature of the cell atoms
        Return:
            velocities of the cell atoms
        """
        xi = np.random.standard_normal((len(self.masses), 3))
        temp = units.kB * temperature
        self.velocities = xi * np.sqrt(self.masses * temp)[:, np.newaxis]
        self.velocities /= self.masses

        return None
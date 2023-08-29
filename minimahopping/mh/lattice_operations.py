import numpy as np
import spglib
import minimahopping.logging.logger as logging
from ase.atoms import Atoms

def lattice_derivative(stress_tensor: np.array, cell: np.array):
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
    deralat = (- prefact * np.matmul(stress_tensor, inv_cell)).T
    return deralat


def reshape_cell(atoms: Atoms, symprec: float):
    lattice, scaled_positions, numbers = spglib.standardize_cell(atoms, to_primitive=False, no_idealize=False, symprec=symprec, angle_tolerance=-1.0)
    if lattice is None:
        msg = "cell reshape did not work"
        logging.logger.warning(msg)
    else:
        atoms.set_cell(lattice)
        atoms.set_scaled_positions(scaled_positions)



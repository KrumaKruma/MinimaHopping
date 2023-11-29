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


def check_boundary_conditions(atoms):
    """
    Function that returns checks and returns the type periodic boundary conditions.
    0: non-periodic/cluster simulation
    1: periodicity in one direction (error not implemented)
    2: 2-d periodic boundaries check if initialized correctly
    3: fully periodic simulation
    """
    
    number_of_periodic_axes = sum(atoms.pbc)

    if number_of_periodic_axes == 0:
        boundary_type = 0
    elif number_of_periodic_axes == 1:
        # 1-dimensional pbc is not implemented
        logging.logger.info("1-D boundary conditions are not implemented")
        quit()
    elif number_of_periodic_axes == 2:
        # check if lattice is correct
        index = np.where(atoms.pbc==False)[0]
        sum_offset_zcell = np.sum(atoms.cell[2,:2] + atoms.cell[:2,2])
        if sum_offset_zcell > 1e-10:
            logging.logger.info("Lattice is not correctly adjusted so that variable cell shape slab can be used.")
            quit()
        elif index != 2:
            logging.logger.info("Abort simulation: For slab simulations z-dimension has to be non-periodic")
            quit()
        else:
            boundary_type = 2
    elif number_of_periodic_axes == 3:
        boundary_type = 3

    return boundary_type



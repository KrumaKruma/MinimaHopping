import numpy as np
import spglib
import minimahopping.logging.logger as logging
from ase.atoms import Atoms
from ase.io import read, write
import ase.atom

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
    nat_in = len(atoms)
    cell = atoms.get_cell(complete=True)
    numbers = atoms.get_atomic_numbers()
    positions = atoms.get_scaled_positions()
    
    spglib_cell = (cell, positions, numbers)
    lattice, scaled_positions, numbers = spglib.standardize_cell(spglib_cell, to_primitive=False, no_idealize=False, symprec=symprec, angle_tolerance=-1.0)
    if lattice is None:
        msg = "cell reshape did not work"
        logging.logger.warning(msg)
    else:
        nat_out = scaled_positions.shape[0]
        if nat_in == nat_out:
            atoms.set_cell(lattice)
            atoms.set_scaled_positions(scaled_positions)
        else:
            logging.logger.warn("standardize_cell operation did not work use lower symprec value")


def transform_deralat(atoms: ase.atom.Atom, forces: np.ndarray, deralat: np.ndarray):
    """
    Transformation of the lattice derivative (see Gubler et al., J.Chem.Phys X, 17, 2023)
    Input:
        atoms: ASE atoms object
        forces: forces of the atoms
        deralat: lattice/cell derivative
    Output:
        new_deralat: transformed lattice derivative
    """
    index = len(np.where(atoms.pbc==True)[0])
    new_deralat = deralat.copy()
    reduced_positions = atoms.get_scaled_positions(wrap=False)
    sumsum = np.zeros((index,index))
    sumsum = np.dot(forces.T[:index,:], reduced_positions[:,:index])
    new_deralat[:index,:index] = deralat[:index,:index] - sumsum.T
    return new_deralat



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



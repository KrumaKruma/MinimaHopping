import numpy as np



def lattice_derivative(stress_tensor, cell):
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



def cart2frac(positions, cell):
    """
    Conversiton of the cartesian coordinates to fractional coordinates
    Input:
        positons: np array (3,nat)
            numpy array containing the postions
        cell: np array (3,3)
            numpy array containing the cell vectors
    Return:
        reduced_positions: np array
            numpy array containing the reduced coordinates
    """

    assert positions.shape[1] == 3, 'Atom postions have not 3 coordinates x,y and z'
    assert cell.shape == (3, 3), 'Cell is not a 3x3 array'

    inv_cell = np.linalg.inv(cell)
    reduced_positions = np.zeros(positions.shape)
    for i, at in enumerate(positions):
        reduced_positions[i, :] = np.matmul(inv_cell, at)

    return reduced_positions



def frac2cart(reduced_positions, cell):
    """
    Conversion of the fractional coordinates to cartesian coordinates
    Input:
        reduced_positions: np array (3,nat)
            numpy array containing the reduced positions
        cell: np array (3,3)
            numpy array containing the lattice vectors
    Return:
        position: np array
            numpy array containing the cartesian coordinates
    """

    assert reduced_positions.shape[1] == 3, 'Atom postions have not 3 coordinates x,y and z'
    assert cell.shape == (3, 3), 'Cell is not a 3x3 array'

    positions = np.zeros(reduced_positions.shape)

    for i, at in enumerate(reduced_positions):
        positions[i, :] = np.matmul(cell, at)

    return positions





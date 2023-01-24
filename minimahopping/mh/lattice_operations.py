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
    deralat = (- prefact * np.matmul(stress_tensor, inv_cell)).T
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


def reshape_cell2(atoms, imax):
    """
    Function that reshapes the cell so that the cell is as cubic as possible
    Input:
        atoms: ASE atoms object
            atoms object containing the lattice parameters. Cell will be reshaped in place
        imax: int
            maximum of the lattice expansion to try
    """
    lattice_in = atoms.get_cell().T
    positions_in = atoms.get_positions().T

    positions, lattice = alat2ascii(positions_in, lattice_in)
    lattice_temp, success = reshape_cell_ascii(lattice, imax)
    if success:
        lattice[:, :] = lattice_temp[:, :]

    permutation_matrix = np.zeros((3, 3))
    permutation_matrix[0, 1] = 1.
    permutation_matrix[1, 2] = 1.
    permutation_matrix[2, 0] = 1.
    inverse_permutation_matrix = np.linalg.inv(permutation_matrix)
    lattice_temp = np.matmul(lattice, permutation_matrix)
    lattice_temp, success = reshape_cell_ascii(lattice_temp, imax)
    lattice_temp = np.matmul(lattice_temp, inverse_permutation_matrix)

    if not success:
        positions[:, :] = positions_in[:, :]
    else:
        atoms.set_cell(lattice_temp.T, scale_atoms=False, apply_constraint=False)
        atoms.set_positions(positions.T)
        positions = atoms.get_positions(wrap=True).T
        positions_in[:, :] = positions[:, :]
        lattice[:, :] = lattice_temp[:, :]

    # positions_temp = positions
    permutation_matrix = np.zeros((3, 3))
    permutation_matrix[0, 2] = 1.
    permutation_matrix[1, 0] = 1.
    permutation_matrix[2, 1] = 1.
    inverse_permutation_matrix = np.linalg.inv(permutation_matrix)
    lattice_temp = np.matmul(lattice, permutation_matrix)
    lattice_temp, success = reshape_cell_ascii(lattice_temp, imax)
    lattice_temp = np.matmul(lattice_temp, inverse_permutation_matrix)

    if not success:
        positions[:, :] = positions_in[:, :]
    else:
        atoms.set_cell(lattice_temp.T, scale_atoms=False, apply_constraint=False)
        atoms.set_positions(positions.T)
        positions = atoms.get_positions(wrap=True).T
        lattice[:, :] = lattice_temp[:, :]


    atoms.set_positions(positions.T)
    atoms.set_cell(lattice.T, scale_atoms=False, apply_constraint=False)



def reshape_cell_ascii(lattice_in, imax):
    """
    Reshapes the cell which is in ascii format
    Input:
        lattice_in: numpy array
            numpy array containing the lattice vectors
        imax: int
            maximum iterations for trying new cell combinations
    Return:
        lattice_min: np array
            numpy array containing the best lattice found
        success: bool
            True if a new cell is found, False otherwise
    """
    area_min = get_area2(lattice_in)
    lattice = np.zeros((3, 3))
    lattice_min = np.zeros((3, 3))
    lattice[:, 0] = lattice_in[:, 0].copy()
    success = False
    lattice_min[:,:] = lattice[:,:]
    for iba in range(-imax, imax, 1):
        lattice[:, 1] = lattice_in[:, 1] + iba * lattice_in[:, 0]
        for ica in range(-imax, imax, 1):
            for icb in range(-imax, imax, 1):
                lattice[:, 2] = lattice_in[:, 2] + ica * lattice_in[:, 0] + icb * lattice_in[:, 1]
                area = get_area2(lattice)
                if area < area_min:
                    area_min = area
                    lattice_min[:, :] = lattice[:, :].copy()
                    success = True
                    #print(area)
                    #print(lattice_min)

    # quit()
    return lattice_min, success


def alat2ascii(pos, alat_in):
    """
    function that rotates the cell so that the lattice is in the ascii format
    input:
        pos: np array
            numpy array which contains the atom positions
        alat_in: np array
            numpy array containing the lattice parameters
    Retrun:
        pos_out: np array
            numpy array containing the positions after the rotation
        alat_out: np array
            numpy array containing the lattice vecotr after the rotation
    """

    alat_out = np.zeros((3, 3))

    r1 = np.linalg.norm(alat_in[:, 0])
    r2 = np.linalg.norm(alat_in[:, 1])
    r3 = np.linalg.norm(alat_in[:, 2])

    alat_out[0, 0] = r1
    alat_out[0, 1] = np.dot(alat_in[:, 0], alat_in[:, 1]) / r1
    alat_out[0, 2] = np.dot(alat_in[:, 0], alat_in[:, 2]) / r1
    alat_out[1, 1] = np.sqrt(r2 ** 2 - alat_out[0, 1] ** 2)
    alat_out[1, 2] = (np.dot(alat_in[:, 1], alat_in[:, 2]) - alat_out[0, 1] * alat_out[0, 2]) / alat_out[1, 1]
    alat_out[2, 2] = np.sqrt(r3 ** 2 - alat_out[0, 2] ** 2 - alat_out[1, 2] ** 2)

    inv_alat_out = np.linalg.inv(alat_in)
    t = np.matmul(alat_out, inv_alat_out)
    pos_out = np.matmul(t, pos)
    return pos_out, alat_out

def get_area2(lattice):
    """
    Calculates the area of the lattice
    Input:
        lattice: np array
            numpy array containing the lattice vectors
    Return:
        area: float
            area of the lattice
    :return:
    """

    area = 0.
    area += np.linalg.norm(np.cross(lattice[0, :], lattice[1, :]))
    area += np.linalg.norm(np.cross(lattice[0, :], lattice[2, :]))
    area += np.linalg.norm(np.cross(lattice[1, :], lattice[2, :]))

    return area



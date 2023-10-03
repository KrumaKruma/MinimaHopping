import numpy as np
import scipy.linalg.lapack as lapack
from minimahopping.mh.periodictable import getRcov_s



def split_bond_forces(positions, atomnames,lattice, forces):

    nat = positions.shape[0]

    rcov = [getRcov_s(x) for x in atomnames]
    rcov = np.array(rcov)
    rcov *= 1.8897161646320724

    forces = forces.T.flatten()

    number_of_bonds, is_bonded = make_bonds(positions, lattice, rcov)
    print(number_of_bonds)

    bond_matrix = get_bond_matrix(number_of_bonds, is_bonded, positions)

    a_solve = np.matmul(np.transpose(bond_matrix), bond_matrix)
    b_solve = np.matmul(np.transpose(bond_matrix), forces)

    lu, piv, coefficients, info = lapack.dgesv(a_solve, b_solve)

    forces_covalent = np.sum(coefficients*bond_matrix, axis=1)

    forces_rest = forces - forces_covalent

    forces_covalent = forces_covalent.reshape(3,nat)
    forces_rest = forces_rest.reshape(3,nat)


    return forces_covalent.T, forces_rest.T


def get_bond_matrix(number_of_bonds, is_bonded, positions):
    nat = is_bonded.shape[0]
    bond_matrix = np.zeros((3*nat, number_of_bonds))
    i_bond = 0
    for iat in range(nat):
        for jat in range(iat+1, nat, 1):
            is_bond = is_bonded[iat,jat]
            if is_bond:
                bond_matrix[(3*iat):(3*iat+3), i_bond] = positions[jat,:] - positions[iat,:]
                bond_matrix[(3*jat):(3*jat+3), i_bond] = - bond_matrix[(3*iat):(3*iat+3), i_bond]
                i_bond += 1
    return bond_matrix


def make_bonds(positions, lattice, rcov):
    ixyzmax = 1
    stretch = 1.2

    distances = get_distances(positions, lattice, ixyzmax)
    number_of_bonds, is_bonded = get_bonds(distances, rcov, stretch) 

    return number_of_bonds, is_bonded


def get_distances(positions, lattice, ixyzmax):
    nat = positions.shape[0]
    distances = np.zeros((nat, nat))

    for iat in range(nat):
        for jat in range(iat+1, nat, 1):
            minimal_distance = 1e100
            for i3 in range(-ixyzmax,ixyzmax, 1):
                for i2 in range(-ixyzmax,ixyzmax, 1):
                    for i1 in range(-ixyzmax,ixyzmax, 1):
                        dx = positions[iat, 0] - positions[jat, 0] - i1*lattice[0,0] - i2*lattice[1,0] - i3*lattice[2,0]
                        dy = positions[iat, 1] - positions[jat, 1] - i1*lattice[0,1] - i2*lattice[1,1] - i3*lattice[2,1]
                        dz = positions[iat, 2] - positions[jat, 2] - i1*lattice[0,2] - i2*lattice[1,2] - i3*lattice[2,2]
                        distance = dx**2 + dy**2 + dz**2
                        minimal_distance = np.minimum(distance, minimal_distance)
            minimal_distance = np.sqrt(minimal_distance)
            distances[iat,jat] = minimal_distance
            distances[jat,iat] = minimal_distance
    return distances


def get_bonds(distances, rcov, stretch):
    nat = distances.shape[0]
    number_of_bonds = 0
    is_bonded = np.zeros((nat,nat), dtype=bool)
    for iat in range(nat):
        for jat in range(iat+1, nat, 1):
            distance = distances[iat,jat]
            bond_distance = stretch * (rcov[iat] * rcov[jat])
            if distance < bond_distance:
                is_bonded[iat, jat] = True
                is_bonded[jat, iat] = True
                number_of_bonds += 1
    return number_of_bonds, is_bonded




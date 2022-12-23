import numpy as np
from copy import deepcopy
import minimahopping.mh.lattice_operations as lat_opt
import warnings
import minimahopping.biomode.biomode as split_forces 



def soften(atoms, calculator, nsoft, alpha_pos = 1e-3, cell_atoms = None, alpha_lat = None):
    atoms = atoms.copy()
    atoms.calc = calculator
    eps_dd = 1e-2

    # initialization and normalization of velocities
    positions_in, cell_positions_in, normed_velocities, normed_cell_velocities, e_pot_in, norm_const = initialize(atoms, cell_atoms, eps_dd)

    # first softening step to get the initial residum
    res_initial, curve, new_normed_velocities, new_normed_cell_velocities = update_velocities(atoms, 
                                                                                        cell_atoms, 
                                                                                        positions_in, 
                                                                                        cell_positions_in, 
                                                                                        normed_velocities, 
                                                                                        normed_cell_velocities, 
                                                                                        e_pot_in, 
                                                                                        eps_dd, 
                                                                                        alpha_pos, 
                                                                                        alpha_lat)
    normed_velocities = new_normed_velocities
    normed_cell_velocities = new_normed_cell_velocities

    # Soften loop
    for i in range(nsoft):
        res, curve, new_normed_velocities, new_normed_cell_velocities = update_velocities(atoms, 
                                                                                cell_atoms, 
                                                                                positions_in, 
                                                                                cell_positions_in, 
                                                                                normed_velocities, 
                                                                                normed_cell_velocities, 
                                                                                e_pot_in, 
                                                                                eps_dd, 
                                                                                alpha_pos, 
                                                                                alpha_lat)
        # Criterion for early stopping of softening
        if res < (curve * eps_dd * 0.5):
            break
        
        # Update velocities
        normed_velocities = new_normed_velocities
        normed_cell_velocities = new_normed_cell_velocities

    # Warning if softening has diverged and not converged
    if res_initial < res:
        warning_msg = "Softening diverged"
        warnings.warn(warning_msg, UserWarning)

    # Renormalize velocities
    velocities = normed_velocities / norm_const

    # Return output
    if cell_atoms is not None:
        cell_velocities = normed_cell_velocities / norm_const
        return velocities, cell_velocities
    else:
        return velocities

    

def initialize(atoms, cell_atoms, eps_dd):
    '''
    Initialization before the iterative part of the softening
    '''

    # Get the initial positions and velocities
    positions_in = atoms.get_positions()
    e_pot_in = atoms.get_potential_energy()
    velocities = atoms.get_velocities()

    # If periodic then also get initial cell positions
    if cell_atoms is not None:
        cell_positions_in = cell_atoms.positions
        cell_velocities = cell_atoms.velocities
    else:
        cell_positions_in = None
        cell_velocities = None

    # Get the normalization constant
    norm_const = get_norm_constant(velocities, cell_velocities, eps_dd)
    if cell_atoms is not None:
        normed_velocities = velocities * norm_const
        normed_cell_velocities = cell_velocities * norm_const
    else:
        normed_velocities = velocities * norm_const
        normed_cell_velocities = None

    return positions_in, cell_positions_in, normed_velocities, normed_cell_velocities, e_pot_in, norm_const


def get_norm_constant(velocities, cell_velocities, eps_dd):
    '''
    get the normalization constant for the velocities
    '''
    if cell_velocities is None:
        norm_const = eps_dd / np.sqrt(np.sum(velocities ** 2))
    else:
        norm_const = eps_dd / np.sqrt(np.sum(velocities ** 2) + np.sum(cell_velocities ** 2))

    return norm_const


def update_velocities(atoms, cell_atoms, positions_in, cell_positions_in, normed_velocities, normed_cell_velocities, e_pot_in, eps_dd, alpha_pos, alpha_lat):
    '''
    Performing one softening steps of the velocities
    '''
    positions = positions_in + normed_velocities
    atoms.set_positions(positions)

    if cell_atoms is not None:
        cell_positions = cell_positions_in + normed_cell_velocities
        reduced_positions = lat_opt.cart2frac(atoms.get_positions(), atoms.get_cell())
        atoms.set_cell(cell_positions, scale_atoms=False, apply_constraint=False)
        positions = lat_opt.frac2cart(reduced_positions, atoms.get_cell())
        atoms.set_positions(positions)


    e_pot = atoms.get_potential_energy()
    forces = atoms.get_forces()

    # Only a parameter for check
    fd2 = 2* (e_pot - e_pot_in) / eps_dd ** 2




    sdf = np.sum(normed_velocities * forces)
    sdd = np.sum(normed_velocities * normed_velocities)

    if cell_atoms is not None:
        stress_tensor = atoms.get_stress(voigt=False,apply_constraint=False)
        cell = atoms.get_cell()
        deralat = lat_opt.lattice_derivative(stress_tensor, cell)
        sdf += np.sum(normed_cell_velocities * deralat)
        sdd += np.sum(normed_cell_velocities * normed_cell_velocities)


    curve = -sdf/sdd

    tt = np.sum(forces * forces) # only for debugging reasons
    forces += curve * normed_velocities
    res = np.sum(forces * forces)


    if cell_atoms is not None:
        tt += np.sum(deralat*deralat)
        deralat = deralat + curve * normed_cell_velocities
        res = res + np.sum(deralat*deralat)

    tt = np.sqrt(tt)
    res = np.sqrt(res)

    # print("SOFTEN:  ", tt, res, curve , fd2, e_pot - e_pot_in)

    positions = positions + alpha_pos * forces
    normed_velocities = positions - positions_in

    if cell_atoms is not None:
        #self._cell_positions = self._cell_positions + self._alpha_lat*_deralat
        cell_positions[0, 0] = cell_positions[0, 0] + alpha_lat * deralat[0, 0]
        cell_positions[0, 1] = cell_positions[0, 1] + alpha_lat * deralat[0, 1]
        cell_positions[0, 2] = cell_positions[0, 2] + alpha_lat * deralat[0, 2]
        cell_positions[1, 0] = cell_positions[1, 0] + alpha_lat * deralat[1, 0]
        cell_positions[1, 1] = cell_positions[1, 1] + alpha_lat * deralat[1, 1]
        cell_positions[2, 0] = cell_positions[2, 0] + alpha_lat * deralat[2, 0]

        normed_cell_velocities = cell_positions - cell_positions_in
        normed_velocities = elim_moment(normed_velocities)
        normed_cell_velocities = elim_torque(normed_cell_velocities, cell_positions, cell_atoms.masses)
        divisor = np.sqrt(np.sum(normed_velocities ** 2) + np.sum(normed_cell_velocities ** 2))
    else:
        normed_velocities = elim_moment(normed_velocities)
        normed_velocities = elim_torque(normed_velocities, positions, atoms.get_masses())
        divisor = np.sqrt(np.sum(normed_velocities ** 2))

    sdd = eps_dd / divisor

    normed_velocities = normed_velocities * sdd
    if cell_atoms is not None:
        normed_cell_velocities = normed_cell_velocities * sdd    

    return res, curve, normed_velocities, normed_cell_velocities


def elim_moment(velocities):
    """
    Elimination of the momentum in the velocities
    """
    # eliminiation of momentum
    _s = np.sum(velocities, axis=0) / velocities.shape[0]
    no_moment_velocities = velocities - _s
    return no_moment_velocities


def elim_torque(velocities, positions, masses):
    """
    Elimination of the torque in the velocites
    """
    # elimination of torque
    # calculate center of mass and subtracti it from positions
    total_mass = np.sum(masses)
    masses_3d = np.vstack([masses] * 3).T
    weighted_positions = positions * masses_3d
    cm = np.sum(weighted_positions, axis=0)
    cm /= total_mass
    weighted_positions -= cm

    evaleria, teneria = moment_of_inertia(positions, masses)

    vrot = np.zeros((positions.shape[0], 3, 3))
    for iat, at in enumerate(positions):
        vrot[iat, :, 0] = np.cross(teneria[:, 0], at)
        vrot[iat, :, 1] = np.cross(teneria[:, 1], at)
        vrot[iat, :, 2] = np.cross(teneria[:, 2], at)

    velocities = velocities.flatten()
    vrot = vrot.reshape((positions.shape[0] * 3, 3), order="C")

    for i, vec in enumerate(vrot.T):
        vrot[:, i] = normalize(vec)

    weighted_positions += cm

    for i, eval in enumerate(evaleria):
        if abs(eval) > 1e-10:
            alpha = np.dot(vrot[:, i], velocities)
            velocities -= alpha * vrot[:, i]

    velocities = velocities.reshape((positions.shape[0], 3))

    # For debugging reasons this can be switched on to controle if torque is eliminated
    # get_torque(weighted_positions, velocities, masses)
    return velocities


def moment_of_inertia(positions, masses):
    '''
    Calcualtion of the eigenvalues and eigenvectors of the inertia tensor
    '''
    inertia_tensor = np.zeros((3, 3))
    for at, mass in zip(positions, masses):
        inertia_tensor[0, 0] += mass * (at[1] ** 2 + at[2] ** 2)
        inertia_tensor[1, 1] += mass * (at[0] ** 2 + at[2] ** 2)
        inertia_tensor[2, 2] += mass * (at[0] ** 2 + at[1] ** 2)
        inertia_tensor[0, 1] -= mass * (at[0] * at[1])
        inertia_tensor[0, 2] -= mass * (at[0] * at[2])
        inertia_tensor[1, 2] -= mass * (at[1] * at[2])

    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
    return eigenvalues, eigenvectors


def normalize(v):
    """
    Function that normalized a vector of arbitrary length
    """
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

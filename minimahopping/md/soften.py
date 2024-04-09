import numpy as np
import minimahopping.mh.lattice_operations as lat_opt
import minimahopping.biomode.biomode as split_forces 
import minimahopping.logging.logger as logging
import ase.atom
from minimahopping.mh.cell_atom import Cell_atom


def soften(atoms: ase.atom.Atom, calculator: ase.calculators.calculator.Calculator, nsoft: int, alpha_pos: float = 1e-3, cell_atoms: Cell_atom = None, alpha_lat: float = None):
    atoms = atoms.copy()
    atoms.calc = calculator
    eps_dd = 1e-2


    # Check if softenig steps is larger than zero
    if nsoft > 0:
        # initialization and normalization of velocities
        positions_in, cell_positions_in, normed_velocities, normed_cell_velocities, e_pot_in, norm_const = initialize(atoms = atoms, cell_atoms = cell_atoms, eps_dd = eps_dd)

        # first softening step to get the initial residum
        res_initial, curve, new_normed_velocities, new_normed_cell_velocities = update_velocities(atoms = atoms, 
                                                                                                  cell_atoms = cell_atoms, 
                                                                                                  positions_in = positions_in, 
                                                                                                  cell_positions_in = cell_positions_in, 
                                                                                                  normed_velocities = normed_velocities, 
                                                                                                  normed_cell_velocities = normed_cell_velocities, 
                                                                                                  e_pot_in = e_pot_in, 
                                                                                                  eps_dd = eps_dd, 
                                                                                                  alpha_pos = alpha_pos, 
                                                                                                  alpha_lat = alpha_lat)
        normed_velocities = new_normed_velocities
        normed_cell_velocities = new_normed_cell_velocities

        # Soften loop
        for i in range(nsoft):
            res, curve, new_normed_velocities, new_normed_cell_velocities = update_velocities(atoms = atoms, 
                                                                                cell_atoms = cell_atoms, 
                                                                                positions_in = positions_in, 
                                                                                cell_positions_in = cell_positions_in, 
                                                                                normed_velocities = normed_velocities, 
                                                                                normed_cell_velocities = normed_cell_velocities, 
                                                                                e_pot_in = e_pot_in, 
                                                                                eps_dd = eps_dd, 
                                                                                alpha_pos = alpha_pos, 
                                                                                alpha_lat = alpha_lat)
            # Criterion for early stopping of softening
            if res < (curve * eps_dd * 0.5):
                break
        
            # Update velocities
            normed_velocities = new_normed_velocities
            normed_cell_velocities = new_normed_cell_velocities

        # Warning if softening has diverged and not converged
        if res_initial < res:
            warning_msg = "Softening diverged"
            logging.logger.warning(warning_msg)

        # Renormalize velocities
        velocities = normed_velocities / norm_const

        # Return output
        if cell_atoms is not None:
            cell_velocities = normed_cell_velocities / norm_const
        else:
            cell_velocities = None
    else:
        # if there is no softening remove momentum and torque
        velocities = elim_moment(velocities=atoms.get_velocities())
        velocities = elim_torque(velocities=velocities, positions=atoms.get_positions(), masses=atoms.get_masses())
        if cell_atoms is not None:
            cell_velocities = elim_moment(velocities = cell_atoms.velocities)
            cell_velocities = elim_torque(velocities = cell_velocities,positions = cell_atoms.positions,masses = cell_atoms.masses)
        else:
            cell_velocities = None 
    
    return velocities, cell_velocities


    

def initialize(atoms: ase.atom.Atom, cell_atoms: Cell_atom, eps_dd: float):
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


def update_velocities(atoms: ase.atom.Atom, 
                      cell_atoms: Cell_atom, 
                      positions_in: np.ndarray, 
                      cell_positions_in: np.ndarray, 
                      normed_velocities: np.ndarray, 
                      normed_cell_velocities: np.ndarray, 
                      e_pot_in: float, 
                      eps_dd: float, 
                      alpha_pos: float, 
                      alpha_lat: float):
    '''
    Performing one softening steps of the velocities
    '''
    positions = positions_in + normed_velocities
    atoms.set_positions(positions)

    if cell_atoms is not None:
        cell_positions = cell_positions_in + normed_cell_velocities
        atoms.set_cell(cell_positions, scale_atoms=True, apply_constraint=False)
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
        deralat = lat_opt.lattice_derivative(stress_tensor = stress_tensor, cell = cell)
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

    debug_msg = "SOFTEN:   {:1.5f}    {:1.5f}    {:1.5f}    {:1.5f}    {:1.5f}".format(tt, res, curve , fd2, e_pot - e_pot_in)
    logging.logger.debug(debug_msg)

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
        normed_velocities = elim_moment(velocities = normed_velocities)
        normed_cell_velocities = elim_torque(velocities = normed_cell_velocities, positions = cell_positions, masses = cell_atoms.masses)
        divisor = np.sqrt(np.sum(normed_velocities ** 2) + np.sum(normed_cell_velocities ** 2))
    else:
        normed_velocities = elim_moment(velocities = normed_velocities)
        normed_velocities = elim_torque(velocities = normed_velocities, positions = positions,masses = atoms.get_masses())
        divisor = np.sqrt(np.sum(normed_velocities ** 2))

    sdd = eps_dd / divisor

    normed_velocities = normed_velocities * sdd
    if cell_atoms is not None:
        normed_cell_velocities = normed_cell_velocities * sdd    

    return res, curve, normed_velocities, normed_cell_velocities


def elim_moment(velocities: np.ndarray):
    """
    Elimination of the momentum in the velocities
    """
    # eliminiation of momentum
    _s = np.sum(velocities, axis=0) / velocities.shape[0]
    no_moment_velocities = velocities - _s
    return no_moment_velocities


def elim_torque(velocities: np.ndarray, positions: np.ndarray, masses: np.ndarray):
    """
    Elimination of the torque in the velocites
    """

    # Calculate center of mass and subtract it from positions
    total_mass = np.sum(masses)
    masses_3d = np.vstack([masses] * 3).T
    weighted_positions = positions * masses_3d

    cm = np.sum(weighted_positions, axis=0) / total_mass
    weighted_positions -= cm

    # Calculate moments of inertia with both functions
    evaleria, teneria = moment_of_inertia(positions = positions,masses = masses)

    # New vrot calculation: Vectorized operation replacing the loop
    teneria_reshaped = teneria.T.reshape(3, 1, 3)
    vrot = np.cross(teneria_reshaped, positions[None, :, :]).transpose(1, 2, 0)

    # flatten velocities and reshape vrot to match dimensions
    velocities = velocities.flatten()
    vrot = vrot.reshape((positions.shape[0] * 3, 3), order="C")

    # normalize vrot using a loop
    for i, vec in enumerate(vrot.T):
        vrot[:, i] = normalize(v = vec)
        weighted_positions += cm

    # New Implementation: Vectorized operation replacing the above for loop
    # mask for elements of evaleria that are greater than 1e-10
    mask = np.abs(evaleria) > 1e-10

    # calculate alpha and update velocities using np.einsum for dot product and updating velocities
    alpha = np.einsum('ij,i->j', vrot[:, mask], velocities)
    velocities -= np.einsum('ij,j->i', vrot[:, mask], alpha)

    # reshape velocities back to the original shape
    velocities = velocities.reshape((positions.shape[0], 3))


    return velocities

def moment_of_inertia(positions: np.ndarray, masses: np.ndarray):
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


def normalize(v: np.ndarray):
    """
    Function that normalized a vector of arbitrary length
    """
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

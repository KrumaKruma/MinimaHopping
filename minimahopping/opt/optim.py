import numpy as np
import warnings
from ase.io import write
import minimahopping.opt.periodic_sqnm as periodic_sqnm
import minimahopping.opt.free_or_fixed_cell_sqnm as free_or_fixed_cell_sqnm
import minimahopping.mh.lattice_operations as lat_opt



def optimization(atoms, calculator, max_force_threshold, outpath, initial_step_size=None, nhist_max=10, lattice_weight=2, alpha_min=1e-3, eps_subsp=1e-3, verbose=True):
    # copy the atoms object and attach calculator to it
    atoms = atoms.copy()
    atoms.calc = calculator

    # If verbose then open the files for writing
    if verbose:
        write(outpath + "geometry_optimization_trajectory.extxyz", atoms, parallel=False)
        f = open(outpath + "geometry_optimization_log.dat", "w")
        msg = 'STEP      ETOT              MAX_FORCE       GAIN_RATIO       STEPSIZE           DIM_SUPSP         MAX_DISP\n'
        f.write(msg)
        f.close()

    # Run geometry optimization
    trajectory, optimizer, number_of_steps = geometry_optimization(atoms, max_force_threshold, outpath,initial_step_size, nhist_max, lattice_weight, alpha_min, eps_subsp, verbose)
    positions_out = atoms.get_positions()
    lattice_out = atoms.get_cell()
    noise = optimizer.lower_bound()

    return positions_out, lattice_out, noise, trajectory, number_of_steps



def geometry_optimization(atoms, max_force_threshold, outpath,initial_step_size, nhist_max, lattice_weight, alpha_min, eps_subsp, verbose):
    # check if periodic boundary condition and assert that either fully periodic or non-periodic
    _pbc = list(set(atoms.pbc))
    assert len(_pbc) == 1, "mixed boundary conditions"
    if True in _pbc:
        trajectory, optimizer, number_of_steps = vcs_geometry_optimization(atoms, max_force_threshold, outpath,initial_step_size, nhist_max, lattice_weight, alpha_min, eps_subsp, verbose)
    else:
        trajectory, optimizer, number_of_steps = free_geometry_optimization(atoms, max_force_threshold, outpath,initial_step_size, nhist_max, alpha_min, eps_subsp, verbose)

    return trajectory, optimizer, number_of_steps



def vcs_geometry_optimization(atoms, max_force_threshold, outpath,initial_step_size, nhist_max, lattice_weight, alpha_min, eps_subsp, verbose):
    '''
    variable cell shape geometry optimization
    '''
    trajectory = []
    nat = len(atoms)
    i_step = 0
    max_disp = 0
    positions_old = atoms.get_positions()

    if initial_step_size is None:
        initial_step_size = -0.001

    initial_lattice = atoms.get_cell().T
    optimizer = periodic_sqnm.periodic_sqnm(nat, initial_lattice, initial_step_size, nhist_max, lattice_weight, alpha_min, eps_subsp)
    max_force_comp = 100
    while max_force_comp > max_force_threshold:
        pos_in = atoms.get_positions()
        max_force_comp = vcs_optimizer_step(atoms, optimizer)
        i_step += 1

        if verbose:
            write_log(atoms, optimizer, outpath, i_step, max_force_comp, max_disp)
        
        is_append_trajectory, positions_current = check_coordinate_shift(atoms, positions_old)
        if is_append_trajectory:
            trajectory.append(atoms.copy())

        is_aboard = check(i_step)
        if is_aboard:
            break

        pos_out = atoms.get_positions()
        max_disp = get_max_disp(pos_in, pos_out)
        positions_old = positions_current

    return trajectory, optimizer, i_step


def vcs_optimizer_step(atoms, optimizer):
    '''
    variable cell shape geometry optimization step
    '''
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress_tensor = atoms.get_stress(voigt=False, apply_constraint=False)
    lattice = atoms.get_cell()
    deralat = lat_opt.lattice_derivative(stress_tensor, lattice)

    max_force_comp = np.max(np.abs(forces))
    max_deralat_comp = np.max(np.abs(deralat))
    max_force_comp = np.maximum(max_force_comp, max_deralat_comp)

    positions = atoms.get_positions().T
    lattice = atoms.get_cell().T

    positions_new, lattice_new = optimizer.optimizer_step(positions, lattice, energy, forces.T, deralat.T)

    atoms.set_positions(positions_new.T)
    atoms.set_cell(lattice_new.T)

    return max_force_comp


def free_geometry_optimization(atoms, max_force_threshold, outpath,initial_step_size, nhist_max, alpha_min, eps_subsp, verbose):
    '''
    Cluster or fixed cell geometry optimization
    '''
    trajectory = []
    nat = len(atoms)
    i_step = 0
    max_disp = 0
    positions_old = atoms.get_positions()

    if initial_step_size is None:
        initial_step_size = -0.001

    optimizer = free_or_fixed_cell_sqnm.free_sqnm(nat, initial_step_size, nhist_max, alpha_min, eps_subsp)
    max_force_comp = 100
    while max_force_comp > max_force_threshold:
        pos_in = atoms.get_positions()
        max_force_comp = free_optimizer_step(atoms, optimizer)
        i_step += 1

        if verbose:
            write_log(atoms, optimizer, outpath, i_step, max_force_comp, max_disp)

        is_append_trajectory, positions_current = check_coordinate_shift(atoms, positions_old)
        if is_append_trajectory:
            trajectory.append(atoms.copy())

        is_aboard = check(i_step)
        if is_aboard:
            break

        pos_out = atoms.get_positions()
        max_disp = get_max_disp(pos_in, pos_out)
        positions_old = positions_current

    return trajectory, optimizer, i_step


def free_optimizer_step(atoms, optimizer):
    '''
    Cluster or fixed cell geometry optimization step
    '''
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    max_force_comp = np.max(np.abs(forces))

    positions = atoms.get_positions()

    positions_new = optimizer.optimizer_step(positions.T, energy, forces.T)

    atoms.set_positions(positions_new.T)

    return max_force_comp


def write_log(atoms, optimizer, outpath, i_step, max_force_comp, max_disp):
    '''
    If verbose is True each optimization step is written to a file and energy and the max force component is
    printed
    '''
    
    energy = atoms.get_potential_energy()
    opt_msg = "{:4d}     {:1.8f}       {:1.5e}     {:1.5e}      {:1.5e}        {:1.5e}       {:1.5e}\n".format(i_step,
                                                                                                            energy,
                                                                                                            max_force_comp,
                                                                                                            optimizer.optimizer.gainratio,
                                                                                                            optimizer.optimizer.alpha,
                                                                                                            optimizer.optimizer.dim_subsp,
                                                                                                            max_disp)
    f = open(outpath+"geometry_optimization_log.dat", "a")
    f.write(opt_msg)
    f.close()
    write(outpath + "geometry_optimization_trajectory.extxyz", atoms, append=True, parallel=False)


def check_coordinate_shift(atoms, positions_old):
    """
    checks if maximal coordinate shift is larger than 0.1
    """
    # Get the maximal coordinate shift
    positions_cur = atoms.get_positions()
    pos_diff = np.abs(positions_cur-positions_old)
    max_diff = np.max(pos_diff)
    # if maximal shift is larger than 0.1 -> write to trajectory and update current position
    if max_diff > 0.1:
        append_traj = True
        positions_current = positions_cur
    else:
        positions_current = positions_old
        append_traj = False
    return append_traj, positions_current


def check(i_step):
    '''
    Check if the geometry optimization has reached the limit of 10000 optimization steps.
    '''
    if i_step > 10000:
        warning_msg = "Geometry did not converge in {:d} optimizations steps".format(i_step)
        warnings.warn(warning_msg, UserWarning)
        is_aboard = True
    else:
        is_aboard = False
    return is_aboard


def get_max_disp(pos_in, pos_out):
    displacements = pos_in-pos_out
    max_disp = np.max(displacements)
    return max_disp






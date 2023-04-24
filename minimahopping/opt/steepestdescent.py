import numpy as np
import warnings
from ase.io import write
# import minimahopping.mh.lattice_operations as lat_opt


def steepestdescent(atoms, max_force_threshold,initial_step_size, verbose, optimization_trajectory_file, optimization_log_file):
    # check if periodic boundary condition and assert that either fully periodic or non-periodic
    '''
    geometry optimization
    '''
    # Initializations
    trajectory = []
    i_step = 0
    max_disp = 0
    positions_old = atoms.get_positions()
    max_force_comp = 100
    epot_max = -1e10
    
    # Assert that no mixed boundary conditions
    _pbc = list(set(atoms.pbc))
    assert len(_pbc) == 1, "mixed boundary conditions"

    # Negative value to get automatic initialize step size
    if initial_step_size is None:
        initial_step_size = -0.001
        # TODO: Rise an error is initial step size is None
    
    # Set if free relax or vcs relax
    if True in _pbc:
        vc_relax = True
    else:
        vc_relax = False

    #TODO: assert that only clusters are optimized with this function

    # GD initializations 
    forces = atoms.get_forces()
    step_size = initial_step_size
    # Optimization loop
    while max_force_comp > max_force_threshold:
        pos_in = atoms.get_positions()
        previous_forces = forces    

        # make an steepest GD
        updated_positions = atoms.get_positions() + step_size * forces
        atoms.set_positions(updated_positions)

        # stepsize feedback
        forces = atoms.get_forces()
        cosangle = np.sum(forces * previous_forces) / (np.linalg.norm(forces) * np.linalg.norm(previous_forces))
        print("DEBUGGG:   ", cosangle,np.sum(forces * previous_forces), np.linalg.norm(forces), np.linalg.norm(previous_forces),  step_size, max_force_comp)
        if cosangle < 0.5:
            step_size /= 2.0
        else:
            step_size *= 1.05

        max_force_comp = np.max(forces)        

        i_step += 1

        energy = atoms.get_potential_energy()
        if energy > epot_max:
            epot_max = energy

        if verbose:
            write_log(atoms, i_step, step_size, max_force_comp, max_disp, optimization_trajectory_file, optimization_log_file)
        
        is_append_trajectory, positions_current = check_coordinate_shift(atoms, positions_old)
        if is_append_trajectory:
            trajectory.append(atoms.copy())

        is_aboard = check(i_step)
        if is_aboard:
            break

        pos_out = atoms.get_positions()
        max_disp = get_max_disp(pos_in, pos_out)
        positions_old = positions_current

    return trajectory, i_step, epot_max


def write_log(atoms, i_step, step_size, max_force_comp, max_disp, optimization_trajectory_file, optimization_log_file):
    '''
    If verbose is True each optimization step is written to a file and energy and the max force component is
    printed
    '''
    
    energy = atoms.get_potential_energy()
    opt_msg = "{:4d}     {:1.8f}       {:1.5e}     {:1.5e}      {:1.5e}     \n".format(i_step,
                                                                                                            energy,
                                                                                                            max_force_comp,
                                                                                                            step_size,
                                                                                                            max_disp)
    
    optimization_log_file.write(opt_msg)
    optimization_log_file.flush()
    write(optimization_trajectory_file, atoms, parallel=False)
    optimization_trajectory_file.flush()

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

import numpy as np
import warnings
from ase.io import write
import minimahopping.mh.lattice_operations as lattice_operations
from sqnm.vcsqnm_for_ase import aseOptimizer
import minimahopping.md.md
import ase.atom
import typing
from ase.calculators.singlepoint import SinglePointCalculator


def optimization(atoms: ase.atom.Atom, 
                 calculator: ase.calculators.calculator.Calculator, 
                 max_force_threshold: float, 
                 outpath: str, 
                 fixed_cell_simulation: bool = False, 
                 initial_step_size: float = None, 
                 nhist_max: int = 10, 
                 lattice_weight:float = 2, 
                 alpha_min: float = 1e-3, 
                 eps_subsp: float = 1e-3,
                 opt_max_steps: float = 10000, 
                 verbose: bool = True):
    # copy the atoms object and attach calculator to it
    atoms = atoms.copy()
    atoms.calc = calculator

    # If verbose then open the files for writing
    if verbose:
        optimization_trajectory_file = open(outpath + "geometry_optimization_trajectory.extxyz", "w")
        write(filename = optimization_trajectory_file, images = atoms, parallel=False)
        optimization_log_file = open(outpath + "geometry_optimization_log.dat", "w")
        msg = 'STEP      ETOT              MAX_FORCE       GAIN_RATIO       STEPSIZE           DIM_SUPSP         MAX_DISP\n'
        optimization_log_file.write(msg)
    else:
        optimization_trajectory_file = None
        optimization_log_file = None

    try:
        # Run geometry optimization
        trajectory, optimizer, number_of_steps, epot_max = geometry_optimization(atoms = atoms,
                                                                    fixed_cell_simulation = fixed_cell_simulation, 
                                                                    max_force_threshold = max_force_threshold, 
                                                                    initial_step_size = initial_step_size, 
                                                                    nhist_max = nhist_max, 
                                                                    lattice_weight = lattice_weight, 
                                                                    alpha_min = alpha_min, 
                                                                    eps_subsp = eps_subsp,
                                                                    opt_max_steps = opt_max_steps, 
                                                                    verbose = verbose, 
                                                                    optimization_trajectory_file = optimization_trajectory_file,
                                                                    optimization_log_file = optimization_log_file)
        positions_out = atoms.get_positions()
        lattice_out = atoms.get_cell()
        noise = optimizer.optimizer.lower_bound()
    finally:
        # Close files
        if verbose:
            optimization_trajectory_file.close()
            optimization_log_file.close()
    return positions_out, lattice_out, noise, trajectory, number_of_steps, epot_max



def geometry_optimization(atoms: ase.atom.Atom, 
                          fixed_cell_simulation: bool,
                          max_force_threshold: float,
                          initial_step_size: float, 
                          nhist_max: int, 
                          lattice_weight: float, 
                          alpha_min: float, 
                          eps_subsp: float, 
                          opt_max_steps: float,
                          verbose: bool, 
                          optimization_trajectory_file: typing.IO, 
                          optimization_log_file: typing.IO):
    # check if periodic boundary condition and assert that either fully periodic or non-periodic
    '''
    geometry optimization
    '''
    # Initializations
    trajectory = [atoms.copy()]
    trajectory[0].calc = SinglePointCalculator(atoms, energy=atoms.get_potential_energy())
    trajectory[0].info.pop('label', None)
    i_step = 0
    max_disp = 0
    positions_old = atoms.get_positions()
    lattice_old = atoms.get_cell()
    max_force_comp = 100
    epot_max = -1e10
    
    # Assert that no mixed boundary conditions
    
    #assert sum(atoms.pbc) == 0 or sum(atoms.pbc) == 3, "mixed boundary conditions"

    # Negative value to get automatic initialize step size
    if initial_step_size is None:
        initial_step_size = -0.001
    
    # Set if free relax or vcs relax
    periodic_type = lattice_operations.check_boundary_conditions(atoms = atoms)
    if periodic_type != 0 and not fixed_cell_simulation:
        vc_relax = True
    else:
        vc_relax = False

    # Initialize optimizer
    optimizer = aseOptimizer(initial_structure=atoms, 
                             vc_relax=vc_relax, 
                             initial_step_size=initial_step_size, 
                             nhist_max=nhist_max, 
                             lattice_weigth=lattice_weight, 
                             alpha_min=alpha_min, 
                             eps_subsp=eps_subsp)
    
    # Optimization loop
    while max_force_comp > max_force_threshold:
        pos_in = atoms.get_positions()
        
        optimizer.step(atoms = atoms)
        max_force_comp = optimizer._getDerivativeNorm()
        
        i_step += 1
        
        energy = atoms.get_potential_energy()
        if energy > epot_max:
            epot_max = energy

        if verbose:
            write_log(atoms = atoms, 
                      optimizer = optimizer, 
                      i_step = i_step, 
                      max_force_comp = max_force_comp, 
                      max_disp = max_disp, 
                      optimization_trajectory_file = optimization_trajectory_file, 
                      optimization_log_file = optimization_log_file)
        
        positions_old, lattice_old, is_append_trajectory = minimahopping.md.md.check_coordinate_shift(atoms = atoms, 
                                                                                                      positions_old = positions_old, 
                                                                                                      lattice_old = lattice_old)
        if is_append_trajectory:
            trajectory.append(atoms.copy())
            trajectory[-1].calc = SinglePointCalculator(atoms, energy=energy)
            trajectory[-1].info.pop('label', None)

        is_aboard = check(i_step=i_step, opt_max_steps=opt_max_steps)
        if is_aboard:
            break

        pos_out = atoms.get_positions()
        max_disp = get_max_disp(pos_in = pos_in, pos_out = pos_out)

    return trajectory, optimizer, i_step, epot_max


def write_log(atoms: ase.atom.Atom, 
              optimizer: aseOptimizer, 
              i_step: int, 
              max_force_comp: float, 
              max_disp: float, 
              optimization_trajectory_file: typing.IO, 
              optimization_log_file: typing.IO):
    '''
    If verbose is True each optimization step is written to a file and energy and the max force component is
    printed
    '''
    
    energy = atoms.get_potential_energy()
    opt_msg = "{:4d}     {:1.8f}       {:1.5e}     {:1.5e}      {:1.5e}        {:1.5e}       {:1.5e}\n".format(i_step,
                                                                                                            energy,
                                                                                                            max_force_comp,
                                                                                                            optimizer.optimizer.optimizer.gainratio,
                                                                                                            optimizer.optimizer.optimizer.alpha,
                                                                                                            optimizer.optimizer.optimizer.dim_subspace,
                                                                                                            max_disp)
     
    forces = atoms.get_forces()
    optimization_log_file.write(opt_msg)
    optimization_log_file.flush()
    write(filename = optimization_trajectory_file, images = atoms, parallel=False)
    optimization_trajectory_file.flush()

def check(i_step: int, opt_max_steps: int):
    '''
    Check if the geometry optimization has reached the limit of opt_max_steps optimization steps.
    '''
    if i_step > opt_max_steps:
        warning_msg = "Geometry did not converge in {:d} optimizations steps".format(i_step)
        warnings.warn(warning_msg, UserWarning)
        is_aboard = True
    else:
        is_aboard = False
    return is_aboard


def get_max_disp(pos_in: np.ndarray, pos_out: np.ndarray):
    displacements = pos_in-pos_out
    max_disp = np.max(displacements)
    return max_disp

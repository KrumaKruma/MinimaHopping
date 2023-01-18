import numpy as np
import minimahopping.mh.lattice_operations as lat_opt
import warnings
from ase.io import write
import minimahopping.md.dbscan as dbscan




def md(atoms, calculator, outpath, cell_atoms = None, dt = 0.001, n_max = 3, verbose = True):
    """ 
    Performs an MD which is visiting n_max minima
    Input:
        atoms: ASE atoms object
        calculator: ASE calculator class
        outpath: path where the output is written
        cell_atoms: object containing the postions and the velocity of the cell atoms
        dt: time step
        n_max: number of minima that are visited by the md
        verbose: if True and MD_log and the MD positions are written
    """
    # Make a copy of the atoms object and attach a calculator
    atoms = atoms.copy()
    atoms.calc = calculator
    # Get the initial positions to get later the deviation from the calculated positions
    positions_old = atoms.get_positions()

    # Initialization of the MD. Get the energy and forces for the first MD step
    e_pot, forces, lattice_force = initialize(atoms, cell_atoms, outpath, verbose)
    # Run the MD until n_max minima have been visited
    etot_max, etot_min, e_pot_max, e_pot_min, trajectory, i_steps = run(atoms, cell_atoms, dt, forces, lattice_force, positions_old, e_pot, n_max, verbose, outpath)

    # adjust the time step for the next MD
    new_dt = adjust_dt(etot_max, etot_min, e_pot_max, e_pot_min, dt)
    # attach the last structure to the MD trajectory
    temp = atoms.copy()
    trajectory.append(temp.copy())
    return atoms.get_positions(), atoms.get_cell(), new_dt, trajectory, e_pot_max, i_steps



def initialize(atoms, cell_atoms, outpath, verbose):
    '''
    Initialization of the MD before the iterative part starts
    '''
    # If verbosity is true make the output files
    if verbose:
        write(outpath + "MD_trajectory.extxyz", atoms, parallel=False)
        f = open(outpath + "MD_log.dat", "w")
        msg = 'STEP      EPOT          EKIN          ETOT               DT\n'
        f.write(msg)
        f.close()

    # Get the energy and forces
    forces = atoms.get_forces()
    e_pot = atoms.get_potential_energy()

    # Get the lattice derivatives if the structure is periodic
    lattice_force = 0.
    if cell_atoms is not None:
        stress_tensor = atoms.get_stress(voigt=False,apply_constraint=False)
        lattice = atoms.get_cell()
        lattice_force = lat_opt.lattice_derivative(stress_tensor, lattice)

    return e_pot, forces, lattice_force


def calc_etot_and_ekin(atoms, cell_atoms, e_pot, etot_max, etot_min):
    """
    Calculation of the total and kinetic energy
    """
    # get the masses of the system (the mass is always 1 here)
    masses = get_masses(atoms=atoms)
    # calculate the kinetic energy
    e_kin = 0.5 * np.sum(masses * atoms.get_velocities() * atoms.get_velocities())
    # calculate the kinetic energy contribution of the cell atoms if the system is periodic
    if cell_atoms is not None:
        cell_masses = get_cell_masses(cell_atoms=cell_atoms)
        e_kin = e_kin + 0.5 * np.sum(cell_masses * cell_atoms.velocities * cell_atoms.velocities)
    e_tot = e_kin + e_pot

    return e_pot, e_kin, e_tot


def update_etot_minmax(e_tot, etot_min, etot_max):
    """
    Function that updates the maximal and minimal total energy
    """
    if e_tot > etot_max:
        etot_max = e_tot

    if e_tot < etot_min:
        etot_min = e_tot
    return etot_min, etot_max


def get_masses(atoms):
    """
    get the masses of the system in the shape (3,nat)
    """
    masses = atoms.get_masses()[:, np.newaxis]/ atoms.get_masses()[:, np.newaxis]
    return masses


def get_cell_masses(cell_atoms):
    """
    get the masses of the cell atoms in the shape (3,3)
    """
    cell_masses = cell_atoms.masses[:, np.newaxis]
    return cell_masses


def run(atoms, cell_atoms, dt, forces, lattice_force, positions_old, e_pot, n_max, verbose, outpath):
    '''
    Running the MD over n_max maxima. If this is not reached after 10'000 steps the MD stops
    '''
    # Initialization of some variables
    i_steps = 0
    sign_old = -1
    n_change = 0
    i_max = 0
    trajectory = []
    e_pot_old = atoms.get_potential_energy()
    is_one_cluster = True
    epot_max = -1e10
    epot_min = 1e10
    etot_max = -1e10
    etot_min = 1e10

    for i in range(10000):
        # perform velocity verlet step
        forces_new, lattice_force_new = verlet_step(atoms, cell_atoms, dt, forces, lattice_force)
        i_steps += 1
        # update the trajectory
        positions_new, trajectory = update_trajectory(atoms, positions_old, trajectory)
        # update the maximal/minimal potential energy
        epot_min_new, epot_max_new = update_epot_minmax(atoms, epot_min, epot_max)

        # check if a new minimum was found
        e_pot_new, sign_new, n_change, i_max = check(atoms, e_pot_old, i_steps, i_max, n_max, n_change, sign_old)
        # calculate the kinetic and total energy
        e_pot, e_kin, e_tot = calc_etot_and_ekin(atoms, cell_atoms, e_pot_new, etot_max, etot_min)
        # update the maximal/minimal total energy
        e_tot_min_new, e_tot_max_new = update_etot_minmax(e_tot, etot_min, etot_max)

        # TODO: Write here a function after the if statement
        if cell_atoms is None:
            if i_steps%5 == 0:
                is_one_cluster = check_and_fix_fragmentation(atoms)



        # Write a log file if verbosity is True
        if verbose:
            write_log(atoms, e_pot, e_kin, e_tot, i_steps, dt, outpath, forces_new)

        # Update energy and forces
        forces = forces_new
        lattice_force = lattice_force_new
        e_pot_old = e_pot_new
        positions_old = positions_new
        epot_max = epot_max_new
        epot_min = epot_min_new
        etot_max = e_tot_max_new
        etot_min = e_tot_min_new
        sign_old = sign_new

        # Stop MD if n_max minim have been visited
        if i_max > n_max:
            if is_one_cluster:
                break

    return etot_max, etot_min, epot_max, epot_min, trajectory, i_steps


def check_and_fix_fragmentation(atoms):
    # check if the cluster has not fragmented
    positions = atoms.get_positions()
    elements = atoms.get_atomic_numbers()
    is_one_cluster = dbscan.one_cluster(positions, elements)
    if not is_one_cluster:
        velocities = atoms.get_velocities()
        masses = atoms.get_masses()
        velocities = dbscan.adjust_velocities(positions, velocities, elements, masses)
        atoms.set_velocities(velocities)
    return is_one_cluster



def update_epot_minmax(atoms, epot_min, epot_max):
    """
    Function that updates the maximal and minimal potential energy
    """
    e_pot = atoms.get_potential_energy()
    if e_pot > epot_max:
        epot_max = e_pot
    if e_pot < epot_min:
        epot_min = e_pot
    return epot_min, epot_max


def update_trajectory(atoms, positions_old, trajectory):
    """
    Function that updates the MD trajectory if the atom have shifted max(r1-r2)>0.2
    """
    # Check if the atoms have shifted more than the threshold
    positions_current, is_add_to_trajectory = check_coordinate_shift(atoms, positions_old)
    # Update trajectory if the atoms have shifted more
    if is_add_to_trajectory:
        temp = atoms.copy()
        trajectory.append(temp.copy())
    return positions_current, trajectory


def verlet_step(atoms, cell_atoms, dt, forces, lattice_force):
    '''
    Performing one Verlet step
    '''
    # Get positions, velocites and masses
    velocities = atoms.get_velocities()
    positions = atoms.get_positions()
    masses = get_masses(atoms)
    
    # Update the positions
    atoms.set_positions(positions + dt * velocities + 0.5 * dt * dt * (forces / masses))

    # if the system is periodic update the cell vectors
    if cell_atoms is not None:
        update_lattice_positions(atoms, cell_atoms, lattice_force, dt)
    
    # Calculate the new forces
    forces_new = atoms.get_forces()
    # update the velocites
    atoms.set_velocities(velocities + 0.5 * dt * ((forces + forces_new) / masses))
    
    # If system is periodic update the cell velocites
    if cell_atoms is not None:
        lattice_force_new = update_lattice_velocities(atoms, cell_atoms, lattice_force, dt)
    else:
        lattice_force_new = None


    return forces_new, lattice_force_new


def update_lattice_positions(atoms, cell_atoms, lattice_force, dt):
    '''
    Update of the lattice postions and moving the atoms accordingly
    '''
    # Get the postions, velocities and masses
    positions = atoms.get_positions()
    lattice = cell_atoms.positions
    cell_masses = get_cell_masses(cell_atoms)
    # Update the cell positions
    # reduced_positions = lat_opt.cart2frac(positions, lattice)
    cell_atoms.positions = cell_atoms.positions + dt * cell_atoms.velocities + 0.5 * dt * dt * (lattice_force / cell_masses)
    atoms.set_cell(cell_atoms.positions, scale_atoms=True, apply_constraint=False)
    # lattice = cell_atoms.positions
    # positions = lat_opt.frac2cart(reduced_positions, lattice)
    # atoms.set_positions(positions)


def update_lattice_velocities(atoms, cell_atoms, lattice_force, dt):
    '''
    Update the lattice velocities
    '''
    # Get stress tensor, lattice and cell masses
    stress_tensor = atoms.get_stress(voigt=False, apply_constraint=False)
    lattice = atoms.get_cell()
    cell_masses = get_cell_masses(cell_atoms)
    # Get new lattice derivatives
    lattice_force_new = lat_opt.lattice_derivative(stress_tensor, lattice)
    # update cell velocities
    cell_atoms.velocities = cell_atoms.velocities + 0.5 * dt * ((lattice_force + lattice_force_new) / cell_masses)
    return lattice_force_new


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
    return positions_current ,append_traj


def check(atoms, e_pot, i_steps, i_max, n_max, n_change, sign_old):
    '''
    Check if a new maximum is found or if 10000 steps are reached
    '''
    if i_steps > 10000:
        warning_msg = "MD did not overcome {:d} maxima in 10000 steps".format(n_max)
        warnings.warn(warning_msg, UserWarning)
        n_change = n_max
    else:
        e_pot_new = atoms.get_potential_energy()
        sign = int(np.sign(e_pot - e_pot_new))
        if sign_old != sign:
            sign_old = sign
            n_change += 1
            if n_change%2 == 0:
                i_max += 1
    return e_pot_new, sign_old, n_change, i_max


def write_log(atoms, e_pot, e_kin, e_tot, i_steps, dt, outpath, forces):
    '''
    Write each MD step into a file and print epot, ekin and etot. The file is overwritten each time the MD is
    started
    '''
    md_msg = "{:4d}      {:1.5f}      {:1.5f}       {:1.8f}       {:1.5f}    {:1.5f}\n".format(i_steps,
                                                                                            e_pot,
                                                                                            e_kin,
                                                                                            e_tot,
                                                                                            dt,
                                                                                            np.linalg.norm(forces))

    f = open(outpath+"MD_log.dat", "a")
    f.write(md_msg)
    f.close()
    write(outpath + "MD_trajectory.extxyz", atoms, append=True, parallel=False)


def adjust_dt(etot_max, etot_min, epot_max, epot_min, dt):
    """
    adjust the timestep according to energy conservation
    """
    _defcon = (etot_max - etot_min)#/(3 * self._nat)
    # print("DEBUGG:   ", (_defcon / (epot_max-epot_min)), etot_max, etot_min,  epot_max, epot_min, dt)
    if (_defcon / (epot_max-epot_min)) < 1e-2:
        dt *= 1.05
    else:
        dt *= 1.0/1.05
    return dt


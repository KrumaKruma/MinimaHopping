import numpy as np
import minimahopping.mh.lattice_operations as lat_opt
import warnings
from ase.io import write
import minimahopping.md.dbscan as dbscan
import minimahopping.logging.logger as logging



def md(atoms, calculator, outpath, fixed_cell_simulation=False,cell_atoms = None, dt = 0.001, n_max = 3, verbose = True, collect_md_file = None, dt_min = 0.0001, md_max_steps=10000, margin=0.3):
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
    
    if verbose:
        md_trajectory_file = open(outpath + "MD_trajectory.extxyz", "w")

        write(md_trajectory_file, atoms, parallel=False)
        
        md_log_file = open(outpath + "MD_log.dat", "w")
        msg = 'STEP      EPOT          EKIN          ETOT               DT\n'
        md_log_file.write(msg)
    else:
        md_trajectory_file = None
        md_log_file = None

    try:
        # Initialization of the MD. Get the energy and forces for the first MD step
        e_pot, forces, lattice_force = initialize(atoms, cell_atoms)
        # Run the MD until n_max minima have been visited
        etot_max, etot_min, e_pot_max, e_pot_min, trajectory, i_steps = run(atoms, 
                                                                            fixed_cell_simulation,
                                                                            cell_atoms, 
                                                                            dt, 
                                                                            forces, 
                                                                            lattice_force, 
                                                                            e_pot, 
                                                                            n_max, 
                                                                            verbose, 
                                                                            collect_md_file, 
                                                                            md_trajectory_file, 
                                                                            md_log_file, 
                                                                            md_max_steps,
                                                                            margin)

        # adjust the time step for the next MD
        new_dt = adjust_dt(etot_max, etot_min, e_pot_max, e_pot_min, dt, dt_min)
        # attach the last structure to the MD trajectory
        trajectory.append(atoms.copy())
        trajectory[-1].info['energy'] = atoms.get_potential_energy()
        trajectory[-1].info.pop('label', None)
    finally:
        if verbose:
            md_trajectory_file.close()
            md_log_file.close()
    return atoms.get_positions(), atoms.get_cell(), new_dt, trajectory, e_pot_max, i_steps



def initialize(atoms, cell_atoms):
    '''
    Initialization of the MD before the iterative part starts
    '''
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


def calc_etot_and_ekin(atoms, cell_atoms):
    """
    Calculation of the total and kinetic energy
    """
    # get the masses of the system (the mass is always 1 here)
    masses = get_masses(atoms=atoms)
    # get potential energy
    e_pot = atoms.get_potential_energy()
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


def run(atoms, fixed_cell_simulation, cell_atoms, dt, forces, lattice_force, e_pot, n_max, verbose, collect_md_file, md_trajectory_file, md_log_file, md_max_steps, margin):
    '''
    Running the MD over n_max maxima. If this is not reached after 10'000 steps the MD stops
    '''
    # Initialization of some variables
    i_steps = 0
    sign_old = -1
    n_change = 0
    i_max = 0
    e_pot_old = atoms.get_potential_energy()
    trajectory = [atoms.copy()]
    trajectory[0].info['energy'] = e_pot_old
    trajectory[0].info.pop('label', None)
    is_one_cluster = True
    epot_max = np.NINF
    epot_min = np.Inf
    etot_max = np.NINF
    etot_min = np.Inf

    initial_positions = atoms.get_positions()
    initial_lattice = atoms.get_cell()
    initial_velocities = atoms.get_velocities()
    positions_old = initial_positions
    lattice_old = initial_lattice

    e_pot_old, e_kin_old, e_tot_old = calc_etot_and_ekin(atoms, cell_atoms)
    
    while i_steps < md_max_steps:
        # perform velocity verlet step
        forces_new, lattice_force_new = verlet_step(atoms, fixed_cell_simulation,cell_atoms, dt, forces, lattice_force)
        i_steps += 1
        # update the trajectory
        positions_old, lattice_old, trajectory = update_trajectory(atoms, positions_old, lattice_old, trajectory, collect_md_file)
        # update the maximal/minimal potential energy
        epot_min_new, epot_max_new = update_epot_minmax(atoms, epot_min, epot_max)

        # check if a new minimum was found
        sign_new, n_change, i_max = check(atoms, fixed_cell_simulation,cell_atoms, lattice_force_new, i_max, n_change, sign_old)
        # calculate the kinetic and total energy
        e_pot, e_kin, e_tot = calc_etot_and_ekin(atoms, cell_atoms)


        _defcon = abs(e_tot - e_tot_old)#/(3 * self._nat)
        energy_conservation = _defcon / len(atoms) #/ abs(e_pot - e_pot_old)
        if energy_conservation > 3.0:
            warning_msg = "MD failed. Restart with smaller dt"
            logging.logger.warning(warning_msg)
            dt = dt/10.
            i_steps = 0
            atoms.set_positions(initial_positions)
            atoms.set_cell(initial_lattice)
            atoms.set_velocities(initial_velocities)


        # update the maximal/minimal total energy
        e_tot_min_new, e_tot_max_new = update_etot_minmax(e_tot, etot_min, etot_max)

        if sum(atoms.pbc) == 0:
            if i_steps%5 == 0:
                is_one_cluster = check_and_fix_fragmentation(atoms, margin)

        # Write a log file if verbosity is True
        if verbose:
            write_log(atoms, e_pot, e_kin, e_tot, i_steps, dt, forces_new, md_trajectory_file, md_log_file, energy_conservation)

        # Update energy and forces
        forces = forces_new
        lattice_force = lattice_force_new
        epot_max = epot_max_new
        epot_min = epot_min_new
        etot_max = e_tot_max_new
        etot_min = e_tot_min_new
        sign_old = sign_new
        e_tot_old = e_tot

        # Stop MD if n_max minim have been visited
        if i_max > n_max:
            if is_one_cluster:
                break
    
    if i_steps >= md_max_steps:
        warning_msg = "MD finished after maximal given steps: {:d}".format(i_steps)
        logging.logger.warning(warning_msg)


    return etot_max, etot_min, epot_max, epot_min, trajectory, i_steps


def check_and_fix_fragmentation(atoms, margin):
    # check if the cluster has not fragmented
    positions = atoms.get_positions()
    elements = atoms.get_atomic_numbers()
    is_one_cluster = dbscan.one_cluster(positions, elements, margin)
    if not is_one_cluster:
        warning_msg = "Cluster fragmented: fixing fragmentation"
        logging.logger.warning(warning_msg)
        velocities = atoms.get_velocities()
        masses = atoms.get_masses()
        velocities = dbscan.adjust_velocities(positions, velocities, elements, masses, margin)
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


def update_trajectory(atoms, positions_old, lattice_old, trajectory, collect_md_file):
    """
    Function that updates the MD trajectory if the atom have shifted max(r1-r2)>0.2
    """
    # Check if the atoms have shifted more than the threshold
    positions_current, lattice_current, is_add_to_trajectory = check_coordinate_shift(atoms, positions_old, lattice_old)
    # Update trajectory if the atoms have shifted more
    if is_add_to_trajectory:
        temp = atoms.copy()
        trajectory.append(temp.copy())
        trajectory[-1].info['energy'] = atoms.get_potential_energy()
        trajectory[-1].info.pop('label', None)
        if collect_md_file is not None:
            atoms.info['energy'] = atoms.get_potential_energy()
            atoms.info['stress'] = atoms.get_stress()
            write(collect_md_file, atoms, parallel=False)
            collect_md_file.flush()

    return positions_current, lattice_current, trajectory


def verlet_step(atoms, fixed_cell_simulation, cell_atoms, dt, forces, lattice_force):
    '''
    Performing one Verlet step
    '''
    # Get positions, velocites and masses
    velocities = atoms.get_velocities()
    positions = atoms.get_positions()
    masses = get_masses(atoms)
    
    # if the system is periodic update the cell vectors
    if True in atoms.pbc and not fixed_cell_simulation:
        assert cell_atoms is not None, "Cell atom class not defined"
        # transform lattice force so that atoms are not moved
        lattice_force_transformed = transform_deralat(atoms, forces, lattice_force)
        # update the positions of the lattice
        update_lattice_positions(atoms, cell_atoms, lattice_force_transformed, dt)

    # Update the positions
    atoms.set_positions(positions + dt * velocities + 0.5 * dt * dt * (forces / masses))
    
    # Calculate the new forces
    forces_new = atoms.get_forces()
    # update the velocites
    atoms.set_velocities(velocities + 0.5 * dt * ((forces + forces_new) / masses))
    
    # If system is periodic update the cell velocites
    if True in atoms.pbc and not fixed_cell_simulation:
        assert cell_atoms is not None, "Cell atom class not defined"
        stress_tensor = atoms.get_stress(voigt=False, apply_constraint=False)
        lattice_force_new = lat_opt.lattice_derivative(stress_tensor, cell_atoms.positions)
        lattice_force_new_transformed = transform_deralat(atoms, forces_new, lattice_force_new)
        update_lattice_velocities(cell_atoms, lattice_force_transformed, lattice_force_new_transformed, dt)
    else:
        lattice_force_new = None
        


    return forces_new, lattice_force_new


def transform_deralat(atoms, forces, deralat):
    index = len(np.where(atoms.pbc==True)[0])
    new_deralat = deralat.copy()
    reduced_positions = atoms.get_scaled_positions(wrap=False)
    sumsum = np.zeros((index,index))
    sumsum = np.dot(forces.T[:index,:], reduced_positions[:,:index])
    new_deralat[:index,:index] = deralat[:index,:index] - sumsum.T
    return new_deralat


def update_lattice_positions(atoms, cell_atoms, lattice_force, dt):
    '''
    Update of the lattice postions and moving the atoms accordingly
    '''
    # Get the postions, velocities and masses
    cell_masses = get_cell_masses(cell_atoms)
    # Update the cell positions
    cell_atoms.positions = cell_atoms.positions + dt * cell_atoms.velocities + 0.5 * dt * dt * (lattice_force / cell_masses)
    atoms.cell = cell_atoms.positions


def update_lattice_velocities(cell_atoms, lattice_force, lattice_force_new, dt):
    '''
    Update the lattice velocities
    '''
    # get the masses of the cell
    cell_masses = get_cell_masses(cell_atoms)
    # update cell velocities
    cell_atoms.velocities = cell_atoms.velocities + 0.5 * dt * ((lattice_force + lattice_force_new) / cell_masses)



def check_coordinate_shift(atoms, positions_old, lattice_old):
    """
    checks if maximal coordinate shift is larger than 0.1
    """
    # Get the maximal coordinate shift
    positions_cur = atoms.get_positions()
    lattice_cur = atoms.get_cell()
    pos_diff = np.max(np.abs(positions_cur-positions_old))
    lat_diff = 0
    if atoms.get_cell() is not None:
        lat_diff = np.max(np.abs(atoms.get_cell() - lattice_old))
    max_diff = max(pos_diff, lat_diff)
    # if maximal shift is larger than 0.1 -> write to trajectory and update current position
    if max_diff > 0.1:
        append_traj = True
        positions_current = positions_cur
        lattice_current = lattice_cur
    else:
        positions_current = positions_old
        lattice_current = lattice_old
        append_traj = False
    return positions_current, lattice_current, append_traj


def check(atoms, fixed_cell_simulation,cell_atoms, cell_forces, i_max, n_change, sign_old):
    '''
    Check if a new maximum is found
    '''
    sign = calculate_sign(atoms, fixed_cell_simulation,cell_atoms, cell_forces) 
    if sign_old != sign:
        sign_old = sign
        n_change += 1
        if n_change%2 == 0:
            i_max += 1
    return sign_old, n_change, i_max


def calculate_sign(atoms, fixed_cell_simulation,cell_atoms, cell_forces):
    '''
    Calculation whether the energy goes up or down with dot product of forces and veliocities.
    '''
    if True in atoms.pbc and not fixed_cell_simulation:
        f = np.concatenate([atoms.get_forces(),cell_forces]).flatten()
        v = np.concatenate([atoms.get_velocities(),cell_atoms.velocities]).flatten()
    else:
        f = atoms.get_forces().flatten()
        v = atoms.get_velocities().flatten()
    dot_product = np.dot(v,f)
    sign = int(np.sign(dot_product))
    return sign


def write_log(atoms, e_pot, e_kin, e_tot, i_steps, dt, forces, md_trajectory_file, md_log_file, energy_conservation):
    '''
    Write each MD step into a file and print epot, ekin and etot. The file is overwritten each time the MD is
    started
    '''
    md_msg = "{:4d}      {:1.5f}      {:1.5f}       {:1.8f}       {:1.5f}    {:1.5f}    {:1.5f}\n".format(i_steps,
                                                                                            e_pot,
                                                                                            e_kin,
                                                                                            e_tot,
                                                                                            dt,
                                                                                            np.linalg.norm(forces),
                                                                                            energy_conservation)

    md_log_file.write(md_msg)
    md_log_file.flush()
    write(md_trajectory_file, atoms, parallel=False)
    md_trajectory_file.flush()


def adjust_dt(etot_max, etot_min, epot_max, epot_min, dt, dt_min):
    """
    adjust the timestep according to energy conservation
    """
    _defcon = (etot_max - etot_min)#/(3 * self._nat)
    # print("DEBUGG:   ", (_defcon / (epot_max-epot_min)), etot_max, etot_min,  epot_max, epot_min, dt)
    if (_defcon / (epot_max-epot_min)) < 1e-2:
        dt *= 1.05
    else:
        if dt > dt_min:
            dt *= 1.0/1.05
    return dt


import numpy as np
import vcs_md.lattice_operations as lat_opt
import warnings
from ase.io import write
import logging



def md(atoms, calculator, outpath, cell_atoms, n_write=1,dt = 0.001, n_steps = 1000):
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

    # If verbosity is true make the output files
    md_trajectory_file = open(outpath + "MD_trajectory.extxyz", "w")

    write(md_trajectory_file, atoms, parallel=False)
    
    md_log_file = open(outpath + "MD_log.dat", "w")
    msg = 'STEP      EPOT          EKIN          ETOT               DT\n'
    md_log_file.write(msg)

    try:
        # Initialization of the MD. Get the energy and forces for the first MD step
        e_pot, forces, lattice_force = initialize(atoms, cell_atoms)
        # Run the MD until n_max minima have been visited
        etot_max, etot_min, e_pot_max, e_pot_min, trajectory, i_steps = run(atoms, cell_atoms, dt, forces, lattice_force, e_pot, n_steps, n_write, md_trajectory_file, md_log_file)

        # attach the last structure to the MD trajectory
        trajectory.append(atoms.copy())
        trajectory[-1].info['energy'] = atoms.get_potential_energy()
        trajectory[-1].info.pop('label', None)
    finally:
        md_trajectory_file.close()
        md_log_file.close()
    return atoms.get_positions(), atoms.get_cell(), trajectory, e_pot_max, i_steps



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


def run(atoms, cell_atoms, dt, forces, lattice_force, e_pot, n_steps, n_write,md_trajectory_file, md_log_file):
    '''
    Running the MD over n_max maxima. If this is not reached after 10'000 steps the MD stops
    '''
    # Initialization of some variables
    i_steps = 0
    e_pot_old = atoms.get_potential_energy()
    trajectory = [atoms.copy()]
    trajectory[0].info['energy'] = e_pot_old
    trajectory[0].info.pop('label', None)
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

    for i_steps in range(n_steps):
        # perform velocity verlet step
        forces_new, lattice_force_new = verlet_step(atoms, cell_atoms, dt, forces, lattice_force)
        # update the maximal/minimal potential energy
        epot_min_new, epot_max_new = update_epot_minmax(atoms, epot_min, epot_max)

        # calculate the kinetic and total energy
        e_pot, e_kin, e_tot = calc_etot_and_ekin(atoms, cell_atoms)


        _defcon = abs(e_tot - e_tot_old)#/(3 * self._nat)
        energy_conservation = _defcon / len(atoms) #/ abs(e_pot - e_pot_old)
        if energy_conservation > 3.0:
            warning_msg = "MD failed. Restart with smaller dt"
            logging.warning(warning_msg)
            dt = dt/10.
            i_steps = 0
            atoms.set_positions(initial_positions)
            atoms.set_cell(initial_lattice)
            atoms.set_velocities(initial_velocities)


        # update the maximal/minimal total energy
        e_tot_min_new, e_tot_max_new = update_etot_minmax(e_tot, etot_min, etot_max)

        # Write a log file if verbosity is True
        if i_steps%n_write == 0:
            write_log(atoms, e_pot, e_kin, e_tot, i_steps, dt, forces_new, md_trajectory_file, md_log_file, energy_conservation)

        # Update energy and forces
        forces = forces_new
        lattice_force = lattice_force_new
        epot_max = epot_max_new
        epot_min = epot_min_new
        etot_max = e_tot_max_new
        etot_min = e_tot_min_new
        e_tot_old = e_tot

    return etot_max, etot_min, epot_max, epot_min, trajectory, i_steps



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



def verlet_step(atoms, cell_atoms, dt, forces, lattice_force):
    '''
    Performing one Verlet step
    '''
    # Get positions, velocites and masses
    velocities = atoms.get_velocities()
    positions = atoms.get_positions()
    masses = get_masses(atoms)

    # if the system is periodic update the cell vectors
    if cell_atoms is not None:
        lattice_force_transformed = transform_deralat(atoms, forces, lattice_force)
        update_lattice_positions(atoms=atoms, cell_atoms=cell_atoms, lattice_force=lattice_force_transformed, dt=dt)
        #cell_masses = get_cell_masses(cell_atoms)
        #cell_atoms.positions = cell_atoms.positions + dt * cell_atoms.velocities + 0.5 * dt * dt * (lattice_force_transformed / cell_masses)
        #atoms.cell = cell_atoms.positions 

    
    # Update the positions
    atoms.set_positions(positions + dt * velocities + 0.5 * dt * dt * (forces / masses))

    
    # Calculate the new forces
    forces_new = atoms.get_forces()

    atoms.set_velocities(velocities + 0.5 * dt * ((forces + forces_new) / masses))

    if cell_atoms is not None:
        stress_tensor = atoms.get_stress(voigt=False, apply_constraint=False)
        lattice_force_new = lat_opt.lattice_derivative(stress_tensor, cell_atoms.positions)
        lattice_force_new_transformed = transform_deralat(atoms, forces_new, lattice_force_new)
        update_lattice_velocities(cell_atoms=cell_atoms, lattice_force=lattice_force_transformed, lattice_force_new=lattice_force_new_transformed, dt=dt)
    
    return forces_new, lattice_force_new

def transform_deralat(atoms, forces, deralat):
    nat = len(atoms)
    reduced_positions = atoms.get_scaled_positions(wrap=False)
    sumsum = np.zeros((3,3))
    for iat in range(nat):
        for i in range(3):
            for j in range(3):
                sumsum[i,j] = sumsum[i,j] + forces[iat,i] * reduced_positions[iat,j]

    new_deralat = deralat - sumsum.T
    return new_deralat


def update_lattice_positions(atoms, cell_atoms, lattice_force, dt):
    '''
    Update of the lattice postions and moving the atoms accordingly
    '''
    # Get the postions, velocities and masses
    cell_masses = get_cell_masses(cell_atoms)
    # Update the cell positions
    cell_atoms.positions = cell_atoms.positions + dt * cell_atoms.velocities + 0.5 * dt * dt * (lattice_force / cell_masses)
    # atoms.set_cell(cell_atoms.positions, scale_atoms=False, apply_constraint=False)
    atoms.cell = cell_atoms.positions #, scale_atoms=True, apply_constraint=False)


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




import numpy as np
import minimahopping.mh.lattice_operations as lat_opt
import warnings
from ase.io import write
import minimahopping.md.dbscan as dbscan
import minimahopping.logging.logger as logging
from minimahopping.mh.cell_atom import Cell_atom
import typing
import ase.atom



def md(atoms: ase.atom.Atom, 
       calculator: ase.calculators.calculator.Calculator, 
       outpath: str, fixed_cell_simulation: bool = False,
       cell_atoms: Cell_atom = None, 
       dt: float = 0.001, 
       n_max: int = 3, 
       verbose: bool = True, 
       collect_md_file: str = None, 
       dt_min: float = 0.0001, 
       md_max_steps: int = 10000, 
       margin: float = 0.3):
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
        collect_md_file: path if the MD is to be collected and None if MD is not collected
        dt_min: minimum timestep for MD
        md_max_steps: maximal number of timesteps for MD
        margin: margin for the fixing of fragmented clusters (2*(1+margin)*max(rcov))
    Return:
        atoms.get_positions(): positions after MD
        atoms.get_cell(): Lattice/Cell after MD
        new_dt: time step for the next MD
        trajectory: trajectory of the MD
        e_pot_max: maximal potentiel energy during MD run
        i_steps: number of steps to finish the MD
    """
    # Make a copy of the atoms object and attach a calculator
    atoms = atoms.copy()
    atoms.calc = calculator
    
    if verbose:
        # if verbosity is True open the log and trajectory file
        # this file is overwritten if a new MD is started
        md_trajectory_file = open(outpath + "MD_trajectory.extxyz", "w")
        write(md_trajectory_file,atoms, parallel=False)
        md_log_file = open(outpath + "MD_log.dat", "w")
        msg = 'STEP      EPOT          EKIN          ETOT               DT\n'
        md_log_file.write(msg)
    else:
        md_trajectory_file = None
        md_log_file = None

    try:
        # Initialization of the MD. Get the energy and forces for the first MD step
        e_pot, forces, lattice_force = initialize(atoms = atoms, cell_atoms = cell_atoms)
        # Run the MD until n_max minima have been visited
        etot_max, etot_min, e_pot_max, e_pot_min, trajectory, i_steps = run(atoms, 
                                                                            fixed_cell_simulation,
                                                                            cell_atoms, 
                                                                            dt, 
                                                                            forces, 
                                                                            lattice_force, 
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



def initialize(atoms: ase.atom.Atom, cell_atoms: Cell_atom):
    '''
    Initialization of the MD before the iterative part starts
    Input:
        atoms: ASE atoms object
        cell_atoms: cell object
    Output:
        potential_energy: initial potential energy
        forces: initial forces
        lattice_forces: initial lattice/cell forces
    '''
    # Get the initial energy and forces
    forces = atoms.get_forces()
    potential_energy = atoms.get_potential_energy()

    # Get the initial lattice derivatives if the structure is periodic
    lattice_force = 0.
    if cell_atoms is not None:
        stress_tensor = atoms.get_stress(voigt=False,apply_constraint=False)
        lattice = atoms.get_cell()
        lattice_force = lat_opt.lattice_derivative(stress_tensor, lattice)

    return potential_energy, forces, lattice_force


def calc_etot_and_ekin(atoms: ase.atom.Atom, cell_atoms: Cell_atom):
    """
    Calculation of the total and kinetic energy
    Input:
        atoms: ASE atoms object
        cell_atoms: Cell atom object
    Return:
        potential_energy: potential energy
        kinetic_energy: kinetic energy
        total_energy: total energy
    """
    # get the masses of the system (the mass is always 1 here)
    masses = get_masses(atoms=atoms)
    # get potential energy
    potential_energy = atoms.get_potential_energy()
    # calculate the kinetic energy
    kinetic_energy = 0.5 * np.sum(masses * atoms.get_velocities() * atoms.get_velocities())
    # calculate the kinetic energy contribution of the cell atoms if the system is periodic
    if cell_atoms is not None:
        cell_masses = get_cell_masses(cell_atoms=cell_atoms)
        kinetic_energy = kinetic_energy + 0.5 * np.sum(cell_masses * cell_atoms.velocities * cell_atoms.velocities)
    total_energy = kinetic_energy + potential_energy

    return potential_energy, kinetic_energy, total_energy


def update_etot_minmax(total_energy: float, minimal_total_energy: float, maximal_total_energy: float):
    """
    Function that updates the maximal and minimal total energy
    Input:
        total_energy: current total energy
        minimal_total_energy: current minimal total energy
        maximal_total_energy: current maximal total energy
    Return:
        minimal_total_energy: updated minimal total energy
        maximal_total_energy: updated maximal total energy
    """
    if total_energy > maximal_total_energy:
        maximal_total_energy = total_energy

    if total_energy < minimal_total_energy:
        minimal_total_energy = total_energy
    return minimal_total_energy, maximal_total_energy


def get_masses(atoms: ase.atom.Atom):
    """
    get the masses of the system in the shape (3,nat)
    Input: 
        atoms: ASE atoms object
    Output:
        masses: masses in a (nat, 3) array
    """
    masses = atoms.get_masses()[:, np.newaxis]/ atoms.get_masses()[:, np.newaxis]
    return masses


def get_cell_masses(cell_atoms: Cell_atom):
    """
    get the masses of the cell atoms in the shape (3,3)
    Input:
        cell_atoms: Cell atoms object
    Return:
        cell_masses: pseudo masses of the cell
    """
    cell_masses = cell_atoms.masses[:, np.newaxis]
    return cell_masses


def run(atoms: ase.atom.Atom, 
        fixed_cell_simulation: bool, 
        cell_atoms: Cell_atom, 
        dt: float, 
        forces: np.ndarray, 
        lattice_force: np.ndarray, 
        n_max: int, 
        verbose: bool, 
        collect_md_file: typing.IO, 
        md_trajectory_file: typing.IO, 
        md_log_file: typing.IO, 
        md_max_steps: int, 
        margin: float):
    '''
    Running the MD over n_max maxima. If this is not reached after 10'000 steps the MD stops
    Input:
        atoms: ASE atoms object
        fixed_cell_simulation: if true the simulation cell is not moved during the MD
        cell_atoms: object containing the postions and the velocity of the cell atoms
        dt: time step
        forces: initial forces for the MD
        lattice_forces: initial lattice/cell forces
        n_max: number of minima that are visited by the md
        verbose: if True and MD_log and the MD positions are written
        collect_md_file: file for the MD trajectory collection
        md_trajectory_file: debug file for the MD trajectory (only if verbose True)
        md_log_file: debug file with MD log (only if verbose True)
        md_max_steps: maximal number of timesteps for MD
        margin: margin for the fixing of fragmented clusters (2*(1+margin)*max(rcov))
    Return:
        etot_max: maximal total energy
        etot_min: minimal total energy
        epot_max: maximal potential energy
        epot_min: minimal potential energy
        trajectory: MD trajectory 
        i_steps: number of MD steps
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
    epot_max = - np.inf
    epot_min = np.inf
    etot_max = - np.inf
    etot_min = np.inf

    initial_positions = atoms.get_positions()
    initial_lattice = atoms.get_cell()
    initial_velocities = atoms.get_velocities()
    positions_old = initial_positions
    lattice_old = initial_lattice

    e_pot_old, e_kin_old, e_tot_old = calc_etot_and_ekin(atoms = atoms, cell_atoms = cell_atoms)

    # start of the MD loop 
    while i_steps < md_max_steps:
        # perform velocity verlet step
        forces_new, lattice_force_new = verlet_step(atoms = atoms,
                                                    fixed_cell_simulation = fixed_cell_simulation,
                                                    cell_atoms = cell_atoms, 
                                                    dt = dt, 
                                                    forces = forces, 
                                                    lattice_force = lattice_force)
        i_steps += 1
        # update the trajectory
        positions_old, lattice_old, trajectory = update_trajectory(atoms = atoms, 
                                                                   positions_old = positions_old, 
                                                                   lattice_old = lattice_old, 
                                                                   trajectory = trajectory, 
                                                                   collect_md_file = collect_md_file)
        # update the maximal/minimal potential energy
        epot_min_new, epot_max_new = update_epot_minmax(atoms = atoms, minimal_potential_energy = epot_min, maximal_potential_energy = epot_max)

        # check if a new minimum was found
        sign_new, n_change, i_max = check(atoms = atoms, 
                                          fixed_cell_simulation = fixed_cell_simulation,
                                          cell_atoms = cell_atoms, 
                                          cell_forces = lattice_force_new, 
                                          i_max = i_max, 
                                          n_change = n_change, 
                                          sign_old = sign_old)
        
        # calculate the kinetic and total energy
        e_pot, e_kin, e_tot = calc_etot_and_ekin(atoms = atoms, cell_atoms = cell_atoms)

        # check if MD is crashed
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
        e_tot_min_new, e_tot_max_new = update_etot_minmax(total_energy = e_tot, minimal_total_energy = etot_min, maximal_total_energy = etot_max)

        if sum(atoms.pbc) == 0:
            if i_steps%5 == 0:
                is_one_cluster = check_and_fix_fragmentation(atoms = atoms, margin = margin)

        # Write a log file if verbosity is True
        if verbose:
            write_log(atoms = atoms, 
                      e_pot = e_pot, 
                      e_kin = e_kin, 
                      e_tot = e_tot, 
                      i_steps = i_steps, 
                      dt = dt, 
                      forces = forces_new, 
                      md_trajectory_file = md_trajectory_file, 
                      md_log_file = md_log_file, 
                      energy_conservation = energy_conservation)

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


def check_and_fix_fragmentation(atoms: ase.atom.Atom, margin: float):
    """
    check if the cluster has not fragmented
    Input:
        atoms: ASE atoms object
        margin: margin for cluster fragmentation
    Return:
        is_one_cluster: if True the cluster is not fragmented
    """
    positions = atoms.get_positions()
    elements = atoms.get_atomic_numbers()
    is_one_cluster = dbscan.one_cluster(positions = positions, elements = elements, margin = margin)
    if not is_one_cluster:
        warning_msg = "Cluster fragmented: fixing fragmentation"
        logging.logger.warning(warning_msg)
        velocities = atoms.get_velocities()
        masses = atoms.get_masses()
        velocities = dbscan.adjust_velocities(positions = positions, velocities = velocities, elements = elements, masses = masses, margin = margin)
        atoms.set_velocities(velocities)
    return is_one_cluster



def update_epot_minmax(atoms: ase.atom.Atom, minimal_potential_energy: float, maximal_potential_energy: float):
    """
    Function that updates the maximal and minimal potential energy
    Input:
        atoms: ASE atoms object
        minimal_potential_energy: current minimal potential energy
        maximal_potential_energy: current maximal potential energy
    Return:
        minimal_potential_energy: updated minimal potential energy
        maximal_potential_energy: updated maximal potential energy
    """
    potential_energy = atoms.get_potential_energy()
    if potential_energy > maximal_potential_energy:
        maximal_potential_energy = potential_energy
    if potential_energy < minimal_potential_energy:
        minimal_potential_energy = potential_energy
    return minimal_potential_energy, maximal_potential_energy


def update_trajectory(atoms: ase.atom.Atom, positions_old: np.ndarray, lattice_old: np.ndarray, trajectory: list, collect_md_file:typing.IO):
    """
    Function that updates the MD trajectory if the atom have shifted max(r1-r2)>0.2
    Input:
        atoms: ASE atoms object
        positions_old: old positions from previous MD step
        lattice_old: old lattice/cell from previous MD step
        trajectory: trajectory list
        collect_md_file: file to collect the MD trajectories
    Output:
        positions_current: current positions
        lattice_current: current lattice
        trajectory: MD trajectory list
    """
    # Check if the atoms have shifted more than the threshold
    positions_current, lattice_current, is_add_to_trajectory = check_coordinate_shift(atoms = atoms, positions_old = positions_old, lattice_old = lattice_old)
    # Update trajectory if the atoms have shifted more
    if is_add_to_trajectory:
        temp = atoms.copy()
        trajectory.append(temp.copy())
        trajectory[-1].info['energy'] = atoms.get_potential_energy()
        trajectory[-1].info.pop('label', None)
        if collect_md_file is not None:
            atoms.info['energy'] = atoms.get_potential_energy()
            if True in atoms.pbc:
                atoms.info['stress'] = atoms.get_stress()
            write(filename = collect_md_file, images = atoms, parallel=False)
            collect_md_file.flush()

    return positions_current, lattice_current, trajectory


def verlet_step(atoms: ase.atom.Atom, fixed_cell_simulation: bool, cell_atoms: Cell_atom, dt: float, forces: np.ndarray, lattice_force: np.ndarray):
    '''
    Performing one Verlet step
    Input:
        atoms: ASE atoms object
        fixed_cell_simulation: if true the simulation cell is not moved during the MD
        cell_atoms: object containing the postions and the velocity of the cell atoms
        dt: time step
        forces: initial forces for the MD
        lattice_forces: initial lattice/cell forces
    Output:
        forces_new: forces after the MD step
        lattice_force_new: lattice/cell forces after the MD step
    '''
    # Get positions, velocites and masses
    velocities = atoms.get_velocities()
    positions = atoms.get_positions()
    masses = get_masses(atoms = atoms)
    
    # if the system is periodic update the cell vectors
    if True in atoms.pbc and not fixed_cell_simulation:
        assert cell_atoms is not None, "Cell atom class not defined"
        # transform lattice force so that atoms are not moved
        lattice_force_transformed = transform_deralat(atoms = atoms, forces = forces, deralat = lattice_force)
        # update the positions of the lattice
        update_lattice_positions(atoms = atoms, cell_atoms = cell_atoms, lattice_force = lattice_force_transformed, dt = dt)
    
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
        lattice_force_new_transformed = transform_deralat(atoms = atoms, forces = forces_new, deralat = lattice_force_new)
        update_lattice_velocities(cell_atoms = cell_atoms, lattice_force = lattice_force_transformed, lattice_force_new = lattice_force_new_transformed, dt = dt)
    else:
        lattice_force_new = None

    return forces_new, lattice_force_new


def transform_deralat(atoms: ase.atom.Atom, forces: np.ndarray, deralat: np.ndarray):
    """
    Transformation of the lattice derivative (see Gubler et al., J.Chem.Phys X, 17, 2023)
    Input:
        atoms: ASE atoms object
        forces: forces of the atoms
        deralat: lattice/cell derivative
    Output:
        new_deralat: transformed lattice derivative
    """
    index = len(np.where(atoms.pbc==True)[0])
    new_deralat = deralat.copy()
    reduced_positions = atoms.get_scaled_positions(wrap=False)
    sumsum = np.zeros((index,index))
    sumsum = np.dot(forces.T[:index,:], reduced_positions[:,:index])
    new_deralat[:index,:index] = deralat[:index,:index] - sumsum.T
    return new_deralat


def update_lattice_positions(atoms: ase.atom.Atom, cell_atoms: Cell_atom, lattice_force: np.ndarray, dt: float):
    '''
    Update of the lattice postions and moving the atoms accordingly
    Input:
        atoms: ASE atoms object
        cell_atoms: Cell atoms object
        lattice_force: force applite to the lattice
        dt: time step
    '''
    # Get the postions, velocities and masses
    cell_masses = get_cell_masses(cell_atoms = cell_atoms)
    # Update the cell positions
    cell_atoms.positions = cell_atoms.positions + dt * cell_atoms.velocities + 0.5 * dt * dt * (lattice_force / cell_masses)
    atoms.cell = cell_atoms.positions


def update_lattice_velocities(cell_atoms: Cell_atom, lattice_force: np.ndarray, lattice_force_new: np.ndarray, dt: float):
    '''
    Update the lattice velocities
    Input:
        cell_atoms: Cell atom object
        lattice_force: lattice/cell forces
        lattice_force_new: updated lattice/cell forces
        dt: time step
    '''
    # get the masses of the cell
    cell_masses = get_cell_masses(cell_atoms)
    # update cell velocities
    cell_atoms.velocities = cell_atoms.velocities + 0.5 * dt * ((lattice_force + lattice_force_new) / cell_masses)



def check_coordinate_shift(atoms: ase.atom.Atom, positions_old: np.ndarray, lattice_old: np.ndarray):
    """
    checks if maximal coordinate shift is larger than 0.1
    Input:
        atoms: ASE atoms object
        positions_old: old positions from previous MD step
        lattice_old: old lattice/cell from previous MD step
    Output:
        positions_current: updated positions if shift is greater than 0.1
        lattice_current: updated lattice/cell if shift is greater than 0.1
        append_trajectory: True if shift is greater than 0.1
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
    if max_diff > 0.01:
        append_traj = True
        positions_current = positions_cur
        lattice_current = lattice_cur
    else:
        positions_current = positions_old
        lattice_current = lattice_old
        append_traj = False
    return positions_current, lattice_current, append_traj


def check(atoms: ase.atom.Atom, fixed_cell_simulation: bool, cell_atoms: Cell_atom, cell_forces: np.ndarray, i_max: int, n_change: int, sign_old: int):
    '''
    Check if a new maximum is found
    Input:
        atoms: ASE atoms object
        fixed_cell_simulation: if True the cell is fixed
        cell_atoms: Cell atom object
        cell_forces: lattice/cell forces
        i_max: number of minima ocillations traversed
        n_change: number of sign changes (number of maxima and minima traversed)
        sign_old: preveous sign
    Return:
        sign_old: updated sign
        n_change: updated number of sign changes
        i_max: updated number of minima oscillations
    '''
    sign = calculate_sign(atoms = atoms, fixed_cell_simulation = fixed_cell_simulation, cell_atoms = cell_atoms, cell_forces = cell_forces) 
    if sign_old != sign:
        sign_old = sign
        n_change += 1
        if n_change%2 == 0:
            i_max += 1
    return sign_old, n_change, i_max


def calculate_sign(atoms: ase.atom.Atom, fixed_cell_simulation: bool, cell_atoms: Cell_atom, cell_forces: np.ndarray):
    '''
    Calculation whether the energy goes up or down with dot product of forces and veliocities.
    Input:
        atoms: ASE atoms object
        fixed_cell_simulations: if True cell is fixed and cell forces are not considered
        cell_atoms: Cell atom object
        cell_forces: Lattice/cell forces
    Return:
        sign: negative/positive slope
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


def write_log(atoms: ase.atom.Atom, 
              e_pot: float, 
              e_kin: float, 
              e_tot: float, 
              i_steps: int, 
              dt: float, 
              forces: np.ndarray, 
              md_trajectory_file: typing.IO, 
              md_log_file: typing.IO, 
              energy_conservation: float):
    '''
    Write each MD step into a file and print epot, ekin and etot. The file is overwritten each time the MD is
    started
    Input:
        atoms: ASE atoms object
        e_pot: potential energy
        e_kin: kinetic energy
        e_tot: total energy
        i_steps: current MD step
        dt: time step
        forces: current forces
        md_trajectory_file: trajectory output file
        md_log_file: log output file
        energy_conservation: energy conservation value
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
    write(filename = md_trajectory_file, images = atoms, parallel=False)
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


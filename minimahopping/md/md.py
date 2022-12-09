import numpy as np
from copy import deepcopy
import minimahopping.mh.lattice_operations as lat_opt
import warnings
from ase.io import read, write
import minimahopping.md.dbscan as dbscan

class MD():
    '''
    Velocity Verlet MD which visits n_max maxima for clusters and variable cell shape velocity Verlet MD for bulk
    systems
    '''
    def __init__(self, atoms, outpath, cell_atoms=None, dt=0.001, n_max=3, verbose=True):
        self._atoms = deepcopy(atoms)
        self._dt = dt
        self._n_max = n_max
        self._verbose = verbose
        self._outpath = outpath
        if cell_atoms is not None:
            self._cell_atoms = deepcopy(cell_atoms)
            self._nat = len(self._atoms) + 3
        else:
            self._cell_atoms = cell_atoms
            self._nat = len(self._atoms)



    def run(self):
        '''
        Running the MD over n_max maxima. If this is not reached after 10'000 steps the MD stops
        '''
        atoms = self._atoms
        positions_old = atoms.get_positions() # initial positions to calcluate the shift in the md to decide when to add to trajectory
        cell_atoms = self._cell_atoms
        e_pot_old, forces_old, lattice_force = self._initialize(atoms=atoms, cell_atoms=cell_atoms)
        is_one_cluster = True
        for i in range(10000):
            forces_new, lattice_force_new, positions_traj = self._verlet_step(atoms, cell_atoms, e_pot_old, forces_old, lattice_force, positions_old)
            self._i_steps += 1
            e_pot_new = self._check(atoms, e_pot_old)
            e_pot, e_kin, e_tot = self._calc_etot_and_ekin(atoms, cell_atoms, e_pot_new)
            
            if self._i_steps%5 == 0:
                positions = atoms.get_positions()
                elements = atoms.get_atomic_numbers()
                is_one_cluster = dbscan.one_cluster(positions, elements)
                if not is_one_cluster:
                    velocities = atoms.get_velocities()
                    masses = atoms.get_masses()
                    velocities = dbscan.adjust_velocities(positions, velocities, elements, masses)
                    atoms.set_velocities(velocities)



            if self._verbose:
                self._write(e_pot, e_kin, e_tot)

            forces_old = forces_new
            lattice_force = lattice_force_new
            e_pot_old = e_pot_new
            positions_old = positions_traj

            if self._i_max > self._n_max:
                if is_one_cluster:
                    break


        self._adjust_dt()
        temp = deepcopy(atoms)
        self._trajectory.append(temp.copy())
        if self._cell_atoms is not None:
            return atoms.get_positions(), atoms.get_cell(), self._dt, self._trajectory, self._epot_max
        else:
            return atoms.get_positions(), self._dt, self._trajectory, self._epot_max


    def _initialize(self, atoms, cell_atoms):
        '''
        Initialization of the MD before the iterative part starts
        '''
        self._trajectory = []
        if self._verbose:
            write(self._outpath + "MD.extxyz", self._atoms)
            f = open(self._outpath + "MD_log.dat", "w")
            msg = 'STEP      EPOT          EKIN          ETOT               DT\n'
            f.write(msg)
            f.close()
        self._masses = atoms.get_masses()[:, np.newaxis]/ self._atoms.get_masses()[:, np.newaxis]  # masses in the md are set to 1
        forces = atoms.get_forces()
        e_pot = atoms.get_potential_energy()
        self._sign_old = -1
        self._n_change = 0
        self._i_steps = 0
        lattice_force = 0.
        if cell_atoms is not None:
            self._cell_masses = cell_atoms.masses[:, np.newaxis]
            stress_tensor = atoms.get_stress(voigt=False,apply_constraint=False)
            lattice = atoms.get_cell()
            lattice_force = lat_opt.lattice_derivative(stress_tensor, lattice)
        self._etot_max = -1e10
        self._epot_max = -1e10
        self._etot_min = 1e10
        self._epot_min = 1e10
        self._i_max = 0
        self._calc_etot_and_ekin(atoms, cell_atoms, e_pot)

        return e_pot, forces, lattice_force


    def _verlet_step(self, atoms, cell_atoms,e_pot, forces, lattice_force, positions_old):
        '''
        Performing one Verlet step
        '''
        velocities = atoms.get_velocities()

        positions = atoms.get_positions()
        atoms.set_positions(positions + self._dt * velocities + 0.5 * self._dt * self._dt * (forces / self._masses))

        if self._cell_atoms is not None:
            self._update_lattice_positions(atoms, cell_atoms, lattice_force)

        forces_new = atoms.get_forces()
        atoms.set_velocities(velocities + 0.5 * self._dt * ((forces + forces_new) / self._masses))
        
        if self._cell_atoms is not None:
            lattice_force = self._update_lattice_velocities(atoms, cell_atoms, lattice_force)

        positions_current, is_add_to_trajectory = self._check_coordinate_shift(atoms, positions_old)
        if is_add_to_trajectory:
            temp = deepcopy(atoms)
            self._trajectory.append(temp.copy())

        _e_pot = atoms.get_potential_energy()
        if _e_pot > self._epot_max:
            self._epot_max = _e_pot
        if _e_pot < self._epot_min:
            self._epot_min = _e_pot

        return forces_new, lattice_force, positions_current


    def _check(self, atoms, e_pot):
        '''
        Check if a new maximum is found or if 10000 steps are reached
        '''
        if self._i_steps > 10000:
            warning_msg = "MD did not overcome {:d} maxima in 10000 steps".format(self._n_max)
            warnings.warn(warning_msg, UserWarning)
            self._n_change = self._n_max
        else:
            e_pot_new = atoms.get_potential_energy()
            sign = int(np.sign(e_pot - e_pot_new))
            if self._sign_old != sign:
                self._sign_old = sign
                self._n_change += 1
                if self._n_change%2 == 0:
                    self._i_max += 1
        return e_pot_new


    def _calc_etot_and_ekin(self, atoms, cell_atoms, e_pot):
        _e_kin = 0.5 * np.sum(self._masses * atoms.get_velocities() * atoms.get_velocities())
        if cell_atoms is not None:
            _e_kin = _e_kin + 0.5 * np.sum(self._cell_masses * cell_atoms.velocities * cell_atoms.velocities)
        e_kin = _e_kin
        e_tot = e_kin + e_pot

        if e_tot > self._etot_max:
            self._etot_max = e_tot

        if e_tot < self._etot_min:
            self._etot_min = e_tot

        return e_pot, e_kin, e_tot


    def _adjust_dt(self):
        _defcon = (self._etot_max - self._etot_min)#/(3 * self._nat)
        #print("DEBUGG:   ", (_defcon / (self._epot_max-self._epot_min)), self._etot_max, self._etot_min,  self._epot_max, self._epot_min, self._dt)
        if (_defcon / (self._epot_max-self._epot_min)) < 1e-2:
            self._dt *= 1.05
        else:
            self._dt *= 1./1.05


    def _write(self, e_pot, e_kin, e_tot):
        '''
        Write each MD step into a file and print epot, ekin and etot. The file is overwritten each time the MD is
        started
        '''
        _i = self._i_steps

        md_msg = "{:4d}      {:1.5f}      {:1.5f}       {:1.8f}       {:1.5f}\n".format(_i,
                                                                                             e_pot,
                                                                                             e_kin,
                                                                                             e_tot,
                                                                                             self._dt)

        f = open(self._outpath+"MD_log.dat", "a")
        f.write(md_msg)
        f.close()
        write(self._outpath + "MD.extxyz", self._atoms, append=True)


    def _check_coordinate_shift(self, atoms, positions_old):
        positions_cur = atoms.get_positions()
        pos_diff = np.abs(positions_cur-positions_old)
        max_diff = np.max(pos_diff)
        if max_diff > 0.01:
            append_traj = True
            self._atoms_old = deepcopy(self._atoms)
            positions_current = positions_cur
        else:
            positions_current = positions_old
            append_traj = False
        return positions_current ,append_traj


    def _update_lattice_positions(self, atoms, cell_atoms, lattice_force):
        '''
        Update of the lattice postions and moving the atoms accordingly
        '''
        positions = atoms.get_positions()
        lattice = cell_atoms.positions
        reduced_positions = lat_opt.cart2frac(positions, lattice)
        cell_atoms.positions = cell_atoms.positions + self._dt * cell_atoms.velocities + 0.5 * self._dt * self._dt * (lattice_force / self._cell_masses)
        atoms.set_cell(cell_atoms.positions, scale_atoms=False, apply_constraint=False)
        lattice = cell_atoms.positions
        positions = lat_opt.frac2cart(reduced_positions, lattice)
        atoms.set_positions(positions)


    def _update_lattice_velocities(self,atoms, cell_atoms, lattice_force):
        '''
        Update the lattice velocities
        '''
        stress_tensor = atoms.get_stress(voigt=False, apply_constraint=False)
        lattice = atoms.get_cell()
        lattice_force_new = lat_opt.lattice_derivative(stress_tensor, lattice)
        cell_atoms.velocities = cell_atoms.velocities + 0.5 * self._dt * ((lattice_force + lattice_force_new) / self._cell_masses)
        return lattice_force_new






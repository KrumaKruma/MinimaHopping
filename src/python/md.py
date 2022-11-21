import numpy as np
from copy import deepcopy
import lattice_operations as lat_opt
import warnings
from ase.io import read, write

#todo: adjust dt automatically
class MD():
    '''
    Velocity Verlet MD which visits n_max maxima for clusters and variable cell shape velocity Verlet MD for bulk
    systems
    '''
    def __init__(self, atoms, outpath, cell_atoms=None, dt=0.001, n_max=3, verbose=True):
        self._atoms = deepcopy(atoms)
        self._atoms_old = deepcopy(atoms)
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
        self._initialize()
        while self._n_change < self._n_max:
            self._verlet_step()
            self._i_steps += 1
            self._check()
            self._calc_etot_and_ekin()
            if self._verbose:
                self._write()
        self._adjust_dt()
        temp = deepcopy(self._atoms)
        self._trajectory.append(temp.copy())
        if self._cell_atoms is not None:
            return self._atoms.get_positions(), self._atoms.get_cell(), self._dt, self._trajectory, self._epot_max
        else:
            return self._atoms.get_positions(), self._dt, self._trajectory, self._epot_max


    def _initialize(self,):
        '''
        Initialization of the MD before the iterative part starts
        '''
        self._trajectory = []
        if self._verbose:
            write(self._outpath + "MD.extxyz", self._atoms)
            f = open(self._outpath + "MD_log.dat", "w")
            f.close()
        self._masses = self._atoms.get_masses()[:, np.newaxis]/ self._atoms.get_masses()[:, np.newaxis]  # for the moment no masses
        self._forces = self._atoms.get_forces()
        self._e_pot = self._atoms.get_potential_energy()
        self._sign_old = -1
        self._n_change = 0
        self._i_steps = 0
        if self._cell_atoms is not None:
            self._cell_masses = self._cell_atoms.masses[:, np.newaxis]
            _stress_tensor = self._atoms.get_stress(voigt=False,apply_constraint=False)
            _lattice = self._atoms.get_cell()
            self._lattice_force = lat_opt.lattice_derivative(_stress_tensor, _lattice)
        self._etot_max = -1e10
        self._epot_max = -1e10
        self._etot_min = 1e10
        self._calc_etot_and_ekin()
        self._target_e_kin = self._e_kin


    def _verlet_step(self):
        '''
        Performing one Verlet step
        '''
        _velocities = self._atoms.get_velocities()
        _positions = self._atoms.get_positions()
        self._atoms.set_positions(_positions + self._dt * _velocities + 0.5 * self._dt * self._dt * (self._forces / self._masses))

        if self._cell_atoms is not None:
            self._update_lattice_positions()

        _forces_new = self._atoms.get_forces()
        self._atoms.set_velocities(_velocities + 0.5 * self._dt * ((self._forces + _forces_new) / self._masses))
        self._forces = _forces_new

        if self._cell_atoms is not None:
            self._update_lattice_velocities()

        if self._check_coordinate_shift():
            temp = deepcopy(self._atoms)
            self._trajectory.append(temp.copy())

        _e_pot = self._atoms.get_potential_energy()
        if _e_pot > self._epot_max:
            self._epot_max = _e_pot


    def _check(self):
        '''
        Check if a new maximum is found or if 10000 steps are reached
        '''
        if self._i_steps > 10000:
            warning_msg = "MD did not overcome {:d} maxima in 10000 steps".format(self._n_max)
            warnings.warn(warning_msg, UserWarning)
            self._n_change = self._n_max
        else:
            _e_pot_new = self._atoms.get_potential_energy()
            sign = int(np.sign(self._e_pot - _e_pot_new))
            if self._sign_old != sign:
                self._sign_old = sign
                self._n_change += 1
            self._e_pot = _e_pot_new


    def _calc_etot_and_ekin(self):
        _e_kin = 0.5 * np.sum(self._masses * self._atoms.get_velocities() * self._atoms.get_velocities())
        if self._cell_atoms is not None:
            _e_kin = _e_kin + 0.5 * np.sum(self._cell_masses * self._cell_atoms.velocities * self._cell_atoms.velocities)
        self._e_kin = _e_kin
        self._e_tot = self._e_kin + self._e_pot

        if self._e_tot > self._etot_max:
            self._etot_max = self._e_tot

        if self._e_tot < self._etot_min:
            self._etot_min = self._e_tot


    def _adjust_dt(self):
        _defcon = (self._etot_max - self._etot_min)/(3 * self._nat)
        #print('devcon', _defcon, self._target_e_kin, _defcon/self._target_e_kin)
        if (_defcon/self._target_e_kin) < 1e-4:
            self._dt *= 1.05
        else:
            self._dt *= 1./1.05


    def _write(self):
        '''
        Write each MD step into a file and print epot, ekin and etot. The file is overwritten each time the MD is
        started
        '''
        _i = self._i_steps
        md_msg = "MD STEP:  {:d}   e_pot: {:1.5f}  e_kin:  {:1.5f}   e_tot:  {:1.5f}  dt:  {:1.5f}\n".format(_i,
                                                                                             self._e_pot,
                                                                                             self._e_kin,
                                                                                             self._e_tot,
                                                                                             self._dt)

        f = open(self._outpath+"MD_log.dat", "a")
        f.write(md_msg)
        f.close()
        write(self._outpath + "MD.extxyz", self._atoms, append=True)


    def _check_coordinate_shift(self,):
        positions_old = self._atoms_old.get_positions()
        positions_cur = self._atoms.get_positions()
        pos_diff = np.abs(positions_cur-positions_old)
        max_diff = np.max(pos_diff)
        if max_diff > 0.01:
            append_traj = True
            self._atoms_old = deepcopy(self._atoms)
        else:
            append_traj = False
        return append_traj


    def _update_lattice_positions(self):
        '''
        Update of the lattice postions and moving the atoms accordingly
        '''
        _positions = self._atoms.get_positions()
        _lattice = self._cell_atoms.positions
        _reduced_positions = lat_opt.cart2frac(_positions, _lattice)
        self._cell_atoms.positions = self._cell_atoms.positions + self._dt * self._cell_atoms.velocities + 0.5 * self._dt * self._dt * (self._lattice_force / self._cell_masses)
        self._atoms.set_cell(self._cell_atoms.positions, scale_atoms=False, apply_constraint=False)
        _lattice = self._cell_atoms.positions
        _positions = lat_opt.frac2cart(_reduced_positions, _lattice)
        self._atoms.set_positions(_positions)


    def _update_lattice_velocities(self,):
        '''
        Update the lattice velocities
        '''
        _stress_tensor = self._atoms.get_stress(voigt=False, apply_constraint=False)
        _lattice = self._atoms.get_cell()
        _lattice_force_new = lat_opt.lattice_derivative(_stress_tensor, _lattice)
        self._cell_atoms.velocities = self._cell_atoms.velocities + 0.5 * self._dt * ((self._lattice_force + _lattice_force_new) / self._cell_masses)
        self._lattice_force = _lattice_force_new






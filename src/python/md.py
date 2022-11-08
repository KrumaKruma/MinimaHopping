import numpy as np
from copy import deepcopy
import lattice_operations as lat_opt
import warnings
from ase.io import read, write


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
        else:
            self._cell_atoms = cell_atoms


    def run(self):
        '''
        Running the MD over n_max maxima. If this is not reached after 10'000 steps the MD stops
        '''
        self._initialize()
        while self._n_change < self._n_max:
            self._verlet_step()
            self._i_steps += 1
            self._check()
            if self._verbose:
                self._write()

        if self._cell_atoms is not None:
            return self._atoms.get_positions(), self._atoms.get_cell()
        else:
            return self._atoms.get_positions()


    def _initialize(self,):
        '''
        Initialization of the MD before the iterative part starts
        '''
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

    def _write(self):
        '''
        Write each MD step into a file and print epot, ekin and etot. The file is overwritten each time the MD is
        started
        '''
        _e_kin = 0.5 * np.sum(self._masses * self._atoms.get_velocities() * self._atoms.get_velocities())
        if self._cell_atoms is not None:
            _e_kin = _e_kin + 0.5 * np.sum(self._cell_masses * self._cell_atoms.velocities * self._cell_atoms.velocities)
        _e_pot = self._e_pot
        _i = self._i_steps
        md_msg = "MD STEP:  {:d}   e_pot: {:1.5f}  e_kin:  {:1.5f}   e_tot:  {:1.5f}".format(_i,
                                                                                             _e_pot,
                                                                                             _e_kin,
                                                                                             _e_pot + _e_kin)
        print(md_msg)
        write(self._outpath + "MD.extxyz", self._atoms, append=True)


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






import numpy as np
from copy import deepcopy
import minimahopping.mh.lattice_operations as lat_opt

import minimahopping.biomode.biomode as biomode 

class Softening():
    """
    Softening the velocites along the softest modes of the postitions and the lattice in case of periodic boundary
    conditions.
    Reference fortran implementation originally written by Hannes Huber
    """
    def __init__(self, atoms, calculator, cell_atoms = None):
        self._atoms = atoms.copy()
        self._atoms.calc = calculator
        if cell_atoms is not None:
            self._cell_atoms = deepcopy(cell_atoms)
            self._alpha_pos = 1e-3
            self._alpha_lat = 1e-3
        else:
            self._cell_atoms = cell_atoms
            self._alpha_pos = 1e-1

        self._eps_dd = 1e-2



    def run(self,  nsoft):
        '''
        Running the softening of the velocities iteratively
        '''
        self._initialize()
        for i in range(nsoft):
            self._update_velocities()
            if self._res < (self._curve * self._eps_dd * 0.5):
                break

        self._velocities /= self._norm_const
        # Restore initial positions in atoms object
        self._atoms.set_positions(self._pos_in)


        if self._cell_atoms is not None:
            self._cell_velocities /= self._norm_const
            return self._velocities, self._cell_velocities
        else:
            return self._velocities


    def _initialize(self):
        '''
        Initialization before the iterative part of the softening
        '''
        self._masses = self._atoms.get_masses()
        self._pos_in = self._atoms.get_positions()
        self._e_pot_in = self._atoms.get_potential_energy()
        self._velocities = self._atoms.get_velocities()

        if self._cell_atoms is not None:
            self._cell_masses = self._cell_atoms.masses
            self._cell_positions_in = self._cell_atoms.positions
            self._cell_velocities = self._cell_atoms.velocities

        self._norm_velocities()




    def _norm_velocities(self):
        '''
        Normalization of the velocities
        '''
        if self._cell_atoms is None:
            self._norm_const = self._eps_dd / np.sqrt(np.sum(self._velocities ** 2))
            self._velocities *= self._norm_const
        else:
            self._norm_const = self._eps_dd / np.sqrt(np.sum(self._velocities ** 2) + np.sum(self._cell_velocities ** 2))
            self._velocities *= self._norm_const
            self._cell_velocities *= self._norm_const


    def _update_velocities(self):
        '''
        Performing one softening steps of the velocities
        '''
        self._pos = self._pos_in + self._velocities
        self._atoms.set_positions(self._pos)

        if self._cell_atoms is not None:
            self._cell_positions = self._cell_positions_in + self._cell_velocities
            _reduced_positions = lat_opt.cart2frac(self._atoms.get_positions(), self._atoms.get_cell())
            self._atoms.set_cell(self._cell_positions, scale_atoms=False, apply_constraint=False)
            _positions = lat_opt.frac2cart(_reduced_positions, self._atoms.get_cell())
            self._atoms.set_positions(_positions)


        _e_pot = self._atoms.get_potential_energy()
        self._forces = self._atoms.get_forces()

        # atomnames = self._atoms.get_chemical_symbols()
        # lattice = self._atoms.get_cell()
        # forces_covalent, forces_rest = biomode.split_bond_forces(_positions, atomnames, lattice, self._forces)
        # self._forces = forces_rest
        # Only a parameter for check
        _fd2 = 2* (_e_pot - self._e_pot_in) / self._eps_dd ** 2




        _sdf = np.sum(self._velocities * self._forces)
        _sdd = np.sum(self._velocities * self._velocities)

        if self._cell_atoms is not None:
            _stress_tensor = self._atoms.get_stress(voigt=False,apply_constraint=False)
            _cell = self._atoms.get_cell()
            _deralat = lat_opt.lattice_derivative(_stress_tensor, _cell)
            _sdf += np.sum(self._cell_velocities * _deralat)
            _sdd += np.sum(self._cell_velocities * self._cell_velocities)


        self._curve = -_sdf/_sdd

        _tt = np.sum(self._forces * self._forces) # only for debugging reasons
        self._forces += self._curve * self._velocities
        self._res = np.sum(self._forces * self._forces)


        if self._cell_atoms is not None:
            _tt += np.sum(_deralat*_deralat)
            _deralat = _deralat + self._curve * self._cell_velocities
            self._res = self._res + np.sum(_deralat*_deralat)

        _tt = np.sqrt(_tt)
        self._res = np.sqrt(self._res)

        # print(_tt, self._res, self._curve , _fd2, _e_pot - self._e_pot_in)


        self._pos = self._pos + self._alpha_pos * self._forces
        self._velocities = self._pos - self._pos_in

        if self._cell_atoms is not None:
            #self._cell_positions = self._cell_positions + self._alpha_lat*_deralat
            self._cell_positions[0, 0] = self._cell_positions[0, 0] + self._alpha_lat * _deralat[0, 0]
            self._cell_positions[0, 1] = self._cell_positions[0, 1] + self._alpha_lat * _deralat[0, 1]
            self._cell_positions[0, 2] = self._cell_positions[0, 2] + self._alpha_lat * _deralat[0, 2]
            self._cell_positions[1, 0] = self._cell_positions[1, 0] + self._alpha_lat * _deralat[1, 0]
            self._cell_positions[1, 1] = self._cell_positions[1, 1] + self._alpha_lat * _deralat[1, 1]
            self._cell_positions[2, 0] = self._cell_positions[2, 0] + self._alpha_lat * _deralat[2, 0]

            self._cell_velocities = self._cell_positions - self._cell_positions_in
            self._velocities = self._elim_moment(self._velocities)
            self._cell_velocities = self._elim_torque(self._cell_velocities, self._cell_positions, self._cell_masses)
            _divisor = np.sqrt(np.sum(self._velocities ** 2) + np.sum(self._cell_velocities ** 2))
        else:
            self._velocities = self._elim_moment(self._velocities)
            self._velocities = self._elim_torque(self._velocities, self._pos, self._masses)
            _divisor = np.sqrt(np.sum(self._velocities ** 2))

        _sdd = self._eps_dd / _divisor

        self._velocities *= _sdd
        if self._cell_atoms is not None:
            self._cell_velocities *= _sdd


    def _elim_moment(self, velocities):
        """
        Elimination of the momentum in the velocities
        """
        # eliminiation of momentum
        _s = np.sum(velocities, axis=0) / velocities.shape[0]
        velocities -= _s
        return velocities


    def _elim_torque(self, velocities, positions, masses):
        """
        Elimination of the torque in the velocites
        """
        # elimination of torque
        # calculate center of mass and subtracti it from positions
        total_mass = np.sum(masses)
        masses_3d = np.vstack([masses] * 3).T
        weighted_positions = positions * masses_3d
        cm = np.sum(weighted_positions, axis=0)
        cm /= total_mass
        weighted_positions -= cm

        evaleria, teneria = self._moment_of_inertia(positions, masses)

        vrot = np.zeros((positions.shape[0], 3, 3))
        for iat, at in enumerate(positions):
            vrot[iat, :, 0] = np.cross(teneria[:, 0], at)
            vrot[iat, :, 1] = np.cross(teneria[:, 1], at)
            vrot[iat, :, 2] = np.cross(teneria[:, 2], at)

        velocities = velocities.flatten()
        vrot = vrot.reshape((positions.shape[0] * 3, 3), order="C")

        for i, vec in enumerate(vrot.T):
            vrot[:, i] = self._normalize(vec)

        weighted_positions += cm

        for i, eval in enumerate(evaleria):
            if abs(eval) > 1e-10:
                alpha = np.dot(vrot[:, i], velocities)
                velocities -= alpha * vrot[:, i]

        velocities = velocities.reshape((positions.shape[0], 3))

        # For debugging reasons this can be switched on to controle if torque is eliminated
        # get_torque(weighted_positions, velocities, masses)
        return velocities



    def _moment_of_inertia(self, positions, masses):
        '''
        Calcualtion of the eigenvalues and eigenvectors of the inertia tensor
        '''
        inertia_tensor = np.zeros((3, 3))
        for at, mass in zip(positions, masses):
            inertia_tensor[0, 0] += mass * (at[1] ** 2 + at[2] ** 2)
            inertia_tensor[1, 1] += mass * (at[0] ** 2 + at[2] ** 2)
            inertia_tensor[2, 2] += mass * (at[0] ** 2 + at[1] ** 2)
            inertia_tensor[0, 1] -= mass * (at[0] * at[1])
            inertia_tensor[0, 2] -= mass * (at[0] * at[2])
            inertia_tensor[1, 2] -= mass * (at[1] * at[2])

        inertia_tensor[1, 0] = inertia_tensor[0, 1]
        inertia_tensor[2, 0] = inertia_tensor[0, 2]
        inertia_tensor[2, 1] = inertia_tensor[1, 2]

        eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
        return eigenvalues, eigenvectors



    def _normalize(self,v):
        """
        Function that normalized a vector of arbitrary length
        """
        norm = np.linalg.norm(v)
        if norm == 0:
            return v
        return v / norm


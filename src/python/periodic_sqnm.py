#!/usr/bin/env python3

# The variable cell shape optimization method is based on the following 
# paper: https://arxiv.org/abs/2206.07339
# More details about the SQNM optimization method are available here:
# https://comphys.unibas.ch/publications/Schaefer2015.pdf
# Author of this document: Moritz Gubler 

# Copyright (C) 2022 Moritz Gubler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import sqnm
import sys


class periodic_sqnm:
    """
    Implementation of the vc-sqnm method. More informations about the algorithm can be found here: https://arxiv.org/abs/2206.07339
    """

    def __init__(self, nat, init_lat, initial_step_size, nhist_max, lattice_weigth, alpha_min, eps_subsp):
        """
        Construct a periodic optimizer object that can be used for variable cell shape optimization.
        Parameters
        ----------
        nat: int
            Number of atoms
        init_lat: 3*3 numpy matrix 
            Matrix containing initial lattice vectors stored columnwise.
        initial_step_size: double
            initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
            * 2.5. If a system only contains weaker bonds a value up to 5.0 may speed up the convergence.
        nhist_max: int
            Maximal number of steps that will be stored in the history list. Use a value between 3 and 20. Must be <= than 3*nat + 9.
        lattice_weight: double
            weight / size of the supercell that is used to transform lattice derivatives. Use a value between 1 and 2. Default is 2.
        alpha_min: double
            Lower limit on the step size. 1.e-2 is the default.
        eps_subsp: double 
            Lower limit on linear dependencies of basis vectors in history list. Default 1.e-4.
        """

        self.nat = nat
        self.ndim = 3 * nat + 9
        self.lattice_weight = lattice_weigth
        self.initial_lat = init_lat
        self.initial_lat_inverse = np.linalg.inv(init_lat)
        self.lattice_transformer = np.diag(1 / np.linalg.norm(self.initial_lat, axis=0)) * self.lattice_weight * np.sqrt(nat)
        self.lattice_transformer_inv = np.linalg.inv(self.lattice_transformer)
        self.optimizer = sqnm.SQNM(self.ndim, nhist_max, initial_step_size, eps_subsp, alpha_min)
        self.fluct = 0.0
        self.a_inv = np.zeros((3,3))
        self.q = np.zeros(3 * nat)
        self.df_dq = np.zeros(3 * nat)
        self.a_tilde = np.zeros((3,3))
        self.df_da_tilde = np.zeros((3, 3))
        self.q_and_lat = np.zeros( 3*nat + 9)
        self.dq_and_dlat = np.zeros( 3*nat + 9)
        self.dd = np.zeros( 3*nat + 9)

    def optimizer_step(self, pos, alat, epot, forces, deralat):
        """
        Calculates new atomic coordinates that are closer to the local minimum. Variable cell shape optimization.
        This function should be used the following way:
        1. calculate energies, forces and stress tensor at positions r and lattice vectors a, b, c.
        2. call the step function to update positions r and lattice vectors.
        3. repeat.
        Parameters
        ----------
        pos: numpy matrix, dimension(3, nat)
            Input: atomic coordinates, dimension(3, nat). 
            Output: improved coordinates that are calculated based on forces from this and previous iterations.
        alat: Numpy 3*3 matrix
            Matrix containing lattice vectors stored columnwise.
        epot: double
            Potential energy of current geometry.
        forces: numpy matrix, dimension(3, nat)
            Forces of current geometry
        deralat: numpy matrix, dimension(3, nat)
            Derivative of energy with repect to lattice vectors. Is connected to the stress tensor by:
            dE / dA = -det(A) * stress A^{-1}^T
            See equation 9 of the vc-sqnm paper. https://arxiv.org/abs/2206.07339
        """

        # check for noise in forces using eq. 23 of vc-sqnm paper
        fnoise = np.linalg.norm(np.sum(forces, axis=1)) / np.sqrt(3 * self.nat)
        if self.fluct == 0.0:
            self.fluct = fnoise
        else:
            self.fluct = .8 * self.fluct + .2 * fnoise
        if self.fluct > 0.2 * np.max( np.abs(forces) ):
            print("""Warning: noise in forces is larger than 0.2 times the largest force component. 
            Convergence is not guaranteed.""", file=sys.stderr)

        self.a_inv = np.linalg.inv(alat)

        self.q = ((self.initial_lat @ self.a_inv) @ pos).reshape(3 * self.nat)
        self.df_dq = (- (alat @ self.initial_lat_inverse) @ forces).reshape(3 * self.nat)

        self.a_tilde = (alat @ self.lattice_transformer).reshape(9)
        self.df_da_tilde = (- deralat @ self.lattice_transformer_inv).reshape(9)

        self.q_and_lat = np.concatenate((self.q, self.a_tilde))
        self.dq_and_dlat = np.concatenate((self.df_dq, self.df_da_tilde))

        self.dd = self.optimizer.sqnm_step(self.q_and_lat, epot, self.dq_and_dlat)

        self.q_and_lat = self.q_and_lat + self.dd

        alat = (self.q_and_lat[(3*self.nat):].reshape(3, 3)) @ self.lattice_transformer_inv
        pos = (alat @ self.initial_lat_inverse) @ (self.q_and_lat[:(3 * self.nat)].reshape(3, self.nat))

        return pos, alat

    def lower_bound(self):
        """ Returns an estimate of a lower bound for the local minumum.
        The estimate is only accurate when the optimization is converged.
        """

        return self.optimizer.lower_bound()

# the rest of this file can be used for testing only

def _energyandforces(nat, pos, alat):
    import bazant
    import random
    epot, forces, deralat = bazant.energyandforces_bazant(alat, pos, nat)
    return epot, forces, deralat

def _rand_vec(nat, sigma):
    import random
    x = np.zeros((3, nat))
    for i in range(nat):
        for j in range(3):
            x[j, i] = random.gauss(0.0, sigma)
    return x

def _tests():
    from ase import io
    import sys
    import time
    b2a = 0.52917721067

    filename = sys.argv[1]

    at = io.read(filename)
    pos = at.get_positions().T / b2a
    lat = at.get_cell().T / b2a
    nat = at.get_global_number_of_atoms()
    lattice_weight = 2.0
    nhist_max = 10

    # if alpha is negative, initial step size will be estimated using eq. 24 and 25 of the
    # of the vc-sqnm paper: https://arxiv.org/abs/2206.07339
    alpha = -0.1
    print('initial step size', alpha)

    opt = periodic_sqnm(nat, lat, alpha, nhist_max, lattice_weight, 1e-2, 1e-3)

    for i in range(30):
        epot, forces, deralat = _energyandforces(nat, pos, lat)
        pos, lat = opt.optimizer_step(pos, lat, epot, forces, deralat)
        print('Epot: %12.10f, Force and lattice derivative norm: %.2E' % (epot, max(np.max(np.linalg.norm(forces, axis=0)), np.linalg.norm(deralat))) )
    print('The current energy is: ', epot)
    print('The estimated lower bound of the ground state is:', opt.lower_bound())
    print('The estimated energy error is:', epot - opt.lower_bound())


if __name__ == "__main__":
    _tests()

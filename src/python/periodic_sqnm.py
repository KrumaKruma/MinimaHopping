#!/usr/bin/env python3
from ast import Lambda
import numpy as np
import sqnm
import sys


class periodic_sqnm:

    def __init__(self, nat, init_lat, initial_step_size, nhist_max, lattice_weigth, alpha_min, eps_subsp):
        self.nat = nat
        self.ndim = 3 * nat + 9
        self.lattice_weight = lattice_weigth
        self.initial_lat = init_lat
        self.initial_lat_inverse = np.linalg.inv(init_lat)
        self.lattice_transformer = np.diag(1 / np.linalg.norm(self.initial_lat, axis=0)) * self.lattice_weight * np.sqrt(nat)
        self.lattice_transformer_inv = np.linalg.inv(self.lattice_transformer)
        self.optimizer = sqnm.SQNM(self.ndim, nhist_max, initial_step_size, eps_subsp, alpha_min)
        self.fluct = 0.0

    def optimizer_step(self, pos, alat, epot, forces, deralat):
        # check for noise in forces using eq. 23 of vc-sqnm paper
        fnoise = np.linalg.norm(np.sum(forces, axis=1)) / np.sqrt(3 * self.nat)
        if self.fluct == 0.0:
            self.fluct = fnoise
        else:
            self.fluct = .8 * self.fluct + .2 * fnoise
        if self.fluct > 0.2 * np.max( np.abs(forces) ):
            print("""Warning: noise in forces is larger than 0.2 times the largest force component. Convergence is not guaranteed.""", file=sys.stderr)

        a_inv = np.linalg.inv(alat)

        q = ((self.initial_lat @ a_inv) @ pos).reshape(3 * self.nat)
        df_dq = (- (alat @ self.initial_lat_inverse) @ forces).reshape(3 * self.nat)

        a_tilde = (alat @ self.lattice_transformer).reshape(9)
        df_da_tilde = (- deralat @ self.lattice_transformer).reshape(9)

        q_and_lat = np.concatenate((q, a_tilde))
        dq_and_dlat = np.concatenate((df_dq, df_da_tilde))

        dd = self.optimizer.sqnm_step(q_and_lat, epot, dq_and_dlat)

        q_and_lat = q_and_lat + dd

        q = q_and_lat[:(3 * self.nat)].reshape(3, self.nat)
        a_tilde = q_and_lat[(3*self.nat):].reshape(3, 3)

        alat = a_tilde @ self.lattice_transformer_inv
        pos = (alat @ self.initial_lat_inverse) @ q

        return pos, alat

    def lower_bound(self):
        return self.optimizer.lower_bound()


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

    # calculate optimal stepsize
    beta = 0.1
    e0, f0, d0 = _energyandforces(nat, pos, lat)
    p1 = pos + beta * f0
    e1, f1, d0 = _energyandforces(nat, p1, lat)
    gtg = np.linalg.norm(f0)**2    
    lmax = ( 2* (e1 - e0 + beta * gtg) / ( gtg * beta**2 ))
    lmax1 = np.linalg.norm(f1 - f0) / (beta * np.linalg.norm(f0))
    alpha = 1 / max(lmax, lmax1)
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


#!/usr/bin/env python3
import random
import numpy as np
import sqnm


class free_sqnm:

    def __init__(self, nat, initial_step_size, nhist_max, alpha_min, eps_subsp):
        self.nat = nat
        self.ndim = 3 * nat
        if self.ndim < nhist_max:
            print("Number of subspace dimensions bigger than number of dimensions.")
            print("Number of subspace dimensions will be reduced")
            nhist_max = self.ndim
        self.optimizer = sqnm.SQNM(self.ndim, nhist_max, initial_step_size, eps_subsp, alpha_min)
        self.fluct = 0.0

    def optimizer_step(self, pos, epot, forces):
        # check for noise in forces using eq. 23 of vc-sqnm paper
        fnoise = np.linalg.norm(np.sum(forces, axis=1)) / np.sqrt(3 * self.nat)
        if self.fluct == 0.0:
            self.fluct = fnoise
        else:
            self.fluct = .8 * self.fluct + .2 * fnoise
        if self.fluct > 0.2 * np.max( np.abs(forces) ):
            print("""Warning: noise in forces is larger than 0.2 times the largest force component. Convergence is not guaranteed.""", file=sys.stderr)
        pos = pos.reshape(3 * self.nat)
        pos = pos + self.optimizer.sqnm_step(pos, epot, -forces.reshape(3 * self.nat))
        pos = pos.reshape((3, self.nat))
        return pos

    def lower_bound(self):
        return self.optimizer.lower_bound()


def _energyandforces(nat, pos, alat):
    import bazant
    epot, forces, deralat = bazant.energyandforces_bazant(alat, pos, nat)
    return epot, forces, deralat

def _tests():
    import bazant
    from ase import io
    import sys
    import time
    b2a = Bohr_Ang = 0.52917721067

    filename = sys.argv[1]

    at = io.read(filename)
    pos = at.get_positions().T / b2a
    lat = at.get_cell().T / b2a
    nat = at.get_global_number_of_atoms()
    alpha = 2


    opt = free_sqnm(nat, alpha, 10, 1e-2, 1e-4)

    for i in range(30):
        epot, forces, deralat = _energyandforces(nat, pos, lat)
        print(epot, np.linalg.norm(forces), np.linalg.norm(deralat))
        pos = opt.optimizer_step(pos, epot, forces)

    print('The current energy is: ', epot)
    print('The estimated lower bound of the ground state is:', opt.lower_bound())
    print('The estimated energy error is:', epot - opt.lower_bound())

if __name__ == "__main__":
    _tests()

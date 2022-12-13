import numpy as np
from copy import deepcopy
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from minimahopping.md.md import MD
from minimahopping.opt.optim import Opt
from minimahopping.mh.minimum import Minimum
from minimahopping.mh.cell_atom import Cell_atom

class adjust_fp():
    def __init__(self, 
                    fmax,
                    iterations=10, 
                    temperature=500, 
                    dt=0.01, 
                    md_min=1, 
                    ns_orb=1, 
                    np_orb=1, 
                    width_cutoff=5.5, 
                    maxnatsphere = 100,
                    exclude=[]):


        self.n_steps = iterations
        self._temperature = temperature
        self._dt = dt
        self._mdmin = md_min
        self._fmax = fmax
        self._ns_orb = ns_orb
        self._np_orb = np_orb
        self._width_cutoff = width_cutoff
        self._verbose = True
        self._exclude = exclude
        self._maxnatsphere = maxnatsphere
        self.minimalist = []
        self.fpd_lst = []
        self.energy_lst = []
        
        self._outpath = "./"

    def run(self,atoms):

        for i in range(self.n_steps):
            atom = deepcopy(atoms)
            MaxwellBoltzmannDistribution(atom, temperature_K=self._temperature)
            if True in atom.pbc:
                mass = .75 * np.sum(atoms.get_masses()) / 10.
                cell_atoms = Cell_atom(mass=mass, positions=atom.get_cell())
                cell_atoms.set_velocities_boltzmann(temperature=self._temperature)
                md = MD(atoms=atom, outpath=self._outpath, cell_atoms=self._cell_atoms, dt=self._dt, n_max=self._mdmin, verbose=self._verbose)
                _positions, _cell , dt, _md_trajectory, _epot_max = md.run()
                atom.set_positions(_positions)
                atom.set_cell(_cell)
                opt = Opt(atoms=atom, outpath=self._outpath, max_froce_threshold=self._fmax, verbose=self._verbose)
                _positions, _lattice, _noise, _opt_trajectory = opt.run()
                atom.set_positions(_positions)
                atom.set_cell(_lattice)
            
            else:
                md = MD(atoms=atom, outpath=self._outpath, cell_atoms=None, dt=self._dt, n_max=self._mdmin, verbose=self._verbose)
                _positions , dt, _md_trajectory, _epot_max = md.run()
                atom.set_positions(_positions)
                opt = Opt(atoms=atom, outpath=self._outpath, max_froce_threshold=self._fmax, verbose=self._verbose)
                _positions, self._noise, _opt_trajectory = opt.run()
                atom.set_positions(_positions)

            mini = Minimum(atom,
                        s = self._ns_orb,
                        p = self._np_orb, 
                        width_cutoff = self._width_cutoff,
                        maxnatsphere = self._maxnatsphere,
                        epot = atom.get_potential_energy(),
                        T=self._temperature,
                        ediff=0)
            
            self.minimalist.append(mini)

        for i in range(len(self.minimalist)):
            mini1 = self.minimalist[i]
            for j in range(i,len(self.minimalist), 1):
                mini2 = self.minimalist[j]
                self.fpd_lst.append(mini1.__equals__(mini2))
                self.energy_lst.append(mini1.__compareto__(mini2))
        
        fpds = np.array(self.fpd_lst)
        f_mean = np.mean(fpds)
        f_std = np.std(fpds)
        f_maximum = np.max(fpds)

        energies = np.array(self.energy_lst)
        e_mean = np.mean(energies)
        e_std = np.mean(energies)
        e_maximum = np.max(energies)
        
        outdict = { 'fp': {
                        'mean' : f_mean,
                        'std' : f_std,
                        'max' : f_maximum
                    },
                    'energy': {
                        'mean' : e_mean,
                        'std' : e_std,
                        'max' : e_maximum
                    }
        }

        return outdict



# def main():

#     atoms = Icosahedron('Na', 2, latticeconstant=None)
#     write("test.extxyz", atoms)
#     calculator = EAM(potential="Na_v2.eam.fs")
#     atoms.calc = calculator


#     fnrm =  0.000005
#     adjust = adjust_fp(atoms, fnrm,)
#     fp_max, fp_mean, fp_std = adjust.run()
#     msg = 'Maximal fingerprint distance between the same local minima:\n' + str(fp_max)
#     msg += '\n Mean fingerprint distance between the same local minima:\n' + str(fp_mean)
#     msg += '\n Standard deviaton of the fingerprint distances:\n' + str(fp_std)
#     print(msg)




# if __name__  == '__main__':
#     main()

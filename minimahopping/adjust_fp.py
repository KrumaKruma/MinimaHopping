import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import minimahopping.md.md as md
import minimahopping.opt.optim as opt
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
        self._verbose = False
        self._exclude = exclude
        self._maxnatsphere = maxnatsphere
        self.minimalist = []
        self.fpd_lst = []
        self.energy_lst = []
        
        self._outpath = "./"

    def run(self, atoms):
        msg = "START MD AND OPTIMIZATION FOR {:d} CYCLES: ".format(self.n_steps)
        print(msg)
        calculator = atoms.calc
        for i in range(self.n_steps):
            msg = "  START CYCLE {:d}".format(i)
            print(msg)
            atom = atoms.copy()
            MaxwellBoltzmannDistribution(atom, temperature_K=self._temperature)
            if True in atom.pbc:
                mass = .75 * np.sum(atoms.get_masses()) / 10.
                cell_atoms = Cell_atom(mass=mass, positions=atom.get_cell())
                cell_atoms.set_velocities_boltzmann(temperature=self._temperature)

                positions, lattice , dt, _md_trajectory, _epot_max, number_of_md_steps = md.md(atoms = atom, 
                                                                                                        calculator = calculator,
                                                                                                        outpath = self._outpath, 
                                                                                                        cell_atoms = self._cell_atoms,
                                                                                                        dt = self._temperature, 
                                                                                                        n_max = self._mdmin,
                                                                                                        verbose = self._verbose)

                atom.set_positions(positions)
                atom.set_cell(lattice)

                positions, lattice, self._noise, _opt_trajectory, number_of_opt_steps = opt.optimization(atoms=atom, 
                                                                        calculator=calculator, 
                                                                        max_force_threshold=self._fmax, 
                                                                        outpath=self._outpath, 
                                                                        verbose=self._verbose)

                atom.set_positions(positions)
                atom.set_cell(lattice)
            
            else:
                print(atom.get_velocities())
                positions, dt, _md_trajectory, _epot_max, number_of_md_steps = md.md(atoms = atom, 
                                                                                                        calculator = calculator,
                                                                                                        outpath = self._outpath, 
                                                                                                        cell_atoms = None,
                                                                                                        dt = self._dt, 
                                                                                                        n_max = self._mdmin,
                                                                                                        verbose = True)
                                                                                                        
                atom.set_positions(positions)

                positions, lattice, self._noise, _opt_trajectory, number_of_opt_steps = opt.optimization(atoms=atom, 
                                        calculator=calculator, 
                                        max_force_threshold=self._fmax, 
                                        outpath=self._outpath, 
                                        verbose=self._verbose)

                atom.set_positions(positions)

            atom.calc = calculator

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

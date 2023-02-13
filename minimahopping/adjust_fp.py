import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import minimahopping.md.md as md
import minimahopping.opt.optim as opt
from minimahopping.mh.minimum import Minimum
from minimahopping.mh.cell_atom import Cell_atom
import minimahopping.mh.lattice_operations as lattice_operations
import minimahopping.mh.parameters
import ase.atom


class adjust_fp():
    def __init__(self,initial_configuration : ase.atom.Atom, iterations : int,**kwargs):

        self.initial_configuration = initial_configuration

        initalParameters = minimahopping.mh.parameters.minimaHoppingParameters(**kwargs)

        self.parameters = initalParameters
        self.parameters._dt = self.parameters.dt0
        self.parameters._T = self.parameters.T0

        self.n_steps = iterations
        self._verbose = False
        self.minimalist = []
        self.fpd_lst = []
        self.energy_lst = []
        
        self._outpath = "./"

    def run(self):
        atoms = self.initial_configuration
        msg = "START MD AND OPTIMIZATION FOR {:d} CYCLES: ".format(self.n_steps)
        print(msg)
        calculator = atoms.calc
        for i in range(self.n_steps):
            msg = "  START CYCLE {:d}".format(i)
            print(msg)
            atom = atoms.copy()
            atom.set_masses(self.initial_configuration.get_masses())
            MaxwellBoltzmannDistribution(atom, temperature_K=self.parameters._T)

            # check that periodic boundaries are the same in all directions (no mixed boundary conditions)
            _pbc = list(set(atom.pbc))
            assert len(_pbc) == 1, "mixed boundary conditions"

            # if periodic boundary conditions create cell atom object
            if True in _pbc:
                # calculate mass for cell atoms
                # Formula if for the MD real masses are used
                # mass = .75 * np.sum(atoms.get_masses()) / 10. # Formula if for the MD real masses are used
                # Formula if mass 1 is used in the MD
                mass = .75 * np.sum(len(atom)) / 10.
                # set position and mass of cell atoms
                cell_atoms = Cell_atom(mass=mass, positions=atom.get_cell())
                # set velocities of the cell atoms
                cell_atoms.set_velocities_boltzmann(temperature=self.parameters._T)
            else:
                # IMPORTANT: cell atoms has to be None for soften/md/geopt if no pbc
                cell_atoms = None
            positions, lattice, dt, _md_trajectory, epot_max_md, number_of_md_steps = md.md(atoms = atom, 
                                                                                                calculator = calculator,
                                                                                                outpath = self._outpath, 
                                                                                                cell_atoms = cell_atoms,
                                                                                                dt = self.parameters._dt, 
                                                                                                n_max = self.parameters.mdmin,
                                                                                                verbose = True, #self._verbose,
                                                                                                collect_md_file = None)


            atom.set_positions(positions)
            # If pbc set new lattice and reshape cell
            if True in _pbc:
                atom.set_cell(lattice)
                lattice_operations.reshape_cell(atom, self.parameters.symprec)


            positions, lattice, self._noise, _opt_trajectory, number_of_opt_steps, epot_max_geopt = opt.optimization(atoms=atom, 
                                                                                                    calculator=calculator, 
                                                                                                    max_force_threshold=self.parameters.fmax, 
                                                                                                    outpath=self._outpath, 
                                                                                                    verbose=self._verbose)


            # Set optimized positions
            atom.set_positions(positions)
            # If Pbc set optimized lattice 
            if True in _pbc:
                atom.set_cell(lattice)




            atom.calc = calculator

            mini = Minimum(atom,
                        s = self.parameters.n_S_orbitals,
                        p = self.parameters.n_P_orbitals, 
                        width_cutoff = self.parameters.width_cutoff,
                        maxnatsphere = self.parameters.maxnatsphere,
                        epot = atom.get_potential_energy(),
                        T=self.parameters._T,
                        ediff=0)
            
            self.minimalist.append(mini)

        for i in range(len(self.minimalist)):
            mini1 = self.minimalist[i]
            for j in range(i,len(self.minimalist), 1):
                mini2 = self.minimalist[j]
                self.fpd_lst.append(mini1.fingerprint_distance(mini2))
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

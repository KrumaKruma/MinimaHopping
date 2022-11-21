import numpy as np
import warnings
from ase.io import read, write
import periodic_sqnm
import free_or_fixed_cell_sqnm
import lattice_operations as lat_opt
from copy import deepcopy



class Opt():
    '''
    SQNM optimization for clusters and variable cell shape SQNM for bulk systems
    '''
    def __init__(self, atoms, max_froce_threshold,outpath,initial_step_size=None, nhist_max=10, lattice_weight=2, alpha_min=1e-3, eps_subsop=1e-3, verbose=True):
        self._atoms = deepcopy(atoms)
        self._max_force_threshold = max_froce_threshold
        self._initial_step_size = initial_step_size
        self._nhist_max = nhist_max
        self._lattice_weight = lattice_weight
        self._alpha_min = alpha_min
        self._eps_subsop = eps_subsop
        self._verbose = verbose
        self._outpath = outpath

        self._nat = self._atoms.get_positions().shape[0]
        self._i_step = 0

        self._atoms_old = deepcopy(atoms)
        self._trajectory = []



    def run(self):
        '''
        Performing either variable cell shape or cluster optimization
        '''
        _pbc = list(set(self._atoms.pbc))
        assert len(_pbc) == 1, "mixed boundary conditions"
        if self._verbose:
            write(self._outpath + "OPT.extxyz", self._atoms)
            f = open(self._outpath + "OPT_log.dat", "w")
            f.close()
        if True in _pbc:
            self._vcs_geom_opt()
            temp = deepcopy(self._atoms)
            self._trajectory.append(temp.copy())
            return self._atoms.get_positions(), self._atoms.get_cell(), self._optim.lower_bound(), self._trajectory
        else:
            self._geom_opt()
            temp = deepcopy(self._atoms)
            self._trajectory.append(temp.copy())
            return self._atoms.get_positions(), self._optim.lower_bound(), self._trajectory




    def _vcs_geom_opt(self,):
        '''
        variable cell shape geometry optimization
        '''
        if self._initial_step_size is None:
            self._get_init_step()

        _init_lat = self._atoms.get_cell().T
        self._optim = periodic_sqnm.periodic_sqnm(self._nat, _init_lat, self._initial_step_size, self._nhist_max, self._lattice_weight, self._alpha_min, self._eps_subsop)
        self._max_force_comp = 100
        while self._max_force_comp > self._max_force_threshold:
            self._vcs_optim_step()
            self._i_step += 1
            if self._verbose:
                self._write()
            if self._check_coordinate_shift():
                temp = deepcopy(self._atoms)
                self._trajectory.append(temp.copy())
            self._check()


    def _geom_opt(self,):
        '''
        Cluster or fixed cell geometry optimization
        '''

        if self._initial_step_size is None:
            self._get_init_step()

        self._optim = free_or_fixed_cell_sqnm.free_sqnm(nat=self._nat, initial_step_size=self._initial_step_size, nhist_max=self._nhist_max,alpha_min=self._alpha_min, eps_subsp=self._eps_subsop)
        self._max_force_comp = 100
        while self._max_force_comp > self._max_force_threshold:
            self._optim_step()
            self._i_step += 1
            if self._verbose:
                self._write()
            if self._check_coordinate_shift():
                temp = deepcopy(self._atoms)
                self._trajectory.append(temp.copy())
            self._check()



    def _optim_step(self,):
        '''
        Cluster or fixed cell geometry optimization step
        '''
        _energy = self._atoms.get_potential_energy()
        _forces = self._atoms.get_forces()

        self._max_force_comp = np.max(_forces)

        _pos = self._atoms.get_positions()

        _pos_new = self._optim.optimizer_step(_pos.T, _energy, _forces.T)

        self._atoms.set_positions(_pos_new.T)


    def _vcs_optim_step(self):
        '''
        variable cell shape geometry optimization step
        '''
        _energy = self._atoms.get_potential_energy()
        _forces = self._atoms.get_forces()
        _stress_tensor = self._atoms.get_stress(voigt=False, apply_constraint=False)
        _lattice = self._atoms.get_cell()
        _deralat = lat_opt.lattice_derivative(_stress_tensor, _lattice)

        _max_force_comp = np.max(_forces)
        _max_deralat_comp = np.max(_deralat)
        self._max_force_comp = np.maximum(_max_force_comp, _max_deralat_comp)

        _pos = self._atoms.get_positions().T
        _lat = self._atoms.get_cell().T

        _pos_new, _lat_new = self._optim.optimizer_step(_pos, _lat, _energy, _forces.T, _deralat)

        self._atoms.set_positions(_pos_new.T)
        self._atoms.set_cell(_lat_new.T)


    def _get_init_step(self,):
        _f0 = self._atoms.get_forces()
        _alpha = 0.001
        _x1 = self._atoms.get_positions() + _alpha*_f0
        _atom = deepcopy(self._atoms)
        _atom.set_positions(_x1)
        _f1 = _atom.get_forces()
        self._initial_step_size = 1. / (np.linalg.norm(_f1-_f0) / (_alpha*np.linalg.norm(_f0)))


    def _check(self):
        '''
        Check if the geometry optimization has reached the limit of 10000 optimization steps.
        '''
        if self._i_step > 10000:
            self._max_force_comp = -1
            warning_msg = "Geometry did not converge in {:d} optimizations steps".format(self._i_step)
            warnings.warn(warning_msg, UserWarning)


    def _write(self):
        '''
        If verbose is True each optimization step is written to a file and energy and the max force component is
        printed
        '''
        _energy = self._atoms.get_potential_energy()
        opt_msg = "OPT Step: {:d}   energy: {:1.8f}  max_force_comp:  {:1.5e}   gain ratio:   {:1.5e}   stepsize:   {:1.5e}\n".format(self._i_step,
                                                                                                                                    _energy,
                                                                                                                                    self._max_force_comp,
                                                                                                                                    self._optim.optimizer.gainratio,
                                                                                                                                    self._optim.optimizer.alpha)
        f = open(self._outpath+"OPT_log.dat", "a")
        f.write(opt_msg)
        f.close()
        write(self._outpath + "OPT.extxyz", self._atoms, append=True)

    def _check_coordinate_shift(self, ):
        positions_old = self._atoms_old.get_positions()
        positions_cur = self._atoms.get_positions()
        pos_diff = np.abs(positions_cur - positions_old)
        max_diff = np.max(pos_diff)
        if max_diff > 0.01:
            append_traj = True
            self._atoms_old = deepcopy(self._atoms)
        else:
            append_traj = False
        return append_traj



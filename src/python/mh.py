import numpy as np
import warnings
import os
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import lattice_operations as lat_opt
from OverlapMatrixFingerprint import OverlapMatrixFingerprint as OMFP
from copy import deepcopy
from soften import Softening
from md import MD
from optim import Opt
from minimum import Minimum
from cell_atom import Cell_atom
from database import Database
import time
import json
import file_handling 

"""
MH Software written by Marco Krummenacher (marco.krummenacher@unibas.ch)
Parts of the software were originally developped (some in Fortran) from other people:
  -- VCSMD: Martin Sommer-Joergenson
  -- VCS Softening: Hannes Huber
  -- VCS optimizer: Moritz Gubler
  -- OMFP in python: Jonas Finkler
"""


class Minimahopping:
    def __init__(self, atoms,
                        T0 = 5000,
                        beta_decrease = 1./1.1,
                        beta_increase = 1.1,
                        Ediff0 = .01,
                        alpha_a = 0.95,
                        alpha_r = 1.05,
                        n_soft = 20,
                        ns_orb = 1,
                        np_orb = 1,
                        width_cutoff = 3.5,
                        maxnatsphere = 100,
                        exclude = [],
                        dt = 0.1,
                        mdmin = 2,
                        fmax = 0.00005,
                        enhanced_feedback = False,
                        energy_threshold = 0.00005,
                        n_poslow = 30,
                        minima_threshold = 5e-4,
                        verbose = True,
                        new_start = False,
                        run_time = 'infinit'):
        """Initialize with an ASE atoms object and keyword arguments."""
        self._atoms = atoms
        self._T0 = T0 
        self._beta_decrease = beta_decrease
        self._beta_increase = beta_increase
        self._Ediff0 = Ediff0 
        self._alpha_a = alpha_a
        self._alpha_r = alpha_r
        self._n_soft = n_soft
        self._ns_orb = ns_orb
        self._np_orb = np_orb
        self._width_cutoff = width_cutoff
        self._maxnatsphere = maxnatsphere
        self._exclude = exclude
        self._dt = dt
        self._mdmin = mdmin
        self._fmax = fmax
        self._enhanced_feedback = enhanced_feedback
        self._energy_threshold = energy_threshold
        self._n_poslow = n_poslow
        self._minima_threshold = minima_threshold
        self._verbose = verbose
        self._new_start = new_start
        self._run_time = run_time

        self._temperature = self._T0
        self._Ediff = self._Ediff0

        self._counter = 0

        self.restart_dict = {"dt" : self._dt,
                             "T" : self._temperature, 
                             "Ediff" : self._Ediff,}

    def __call__(self, totalsteps = None):

        # Check if all the files are there for a restart.
        self._outpath = 'output/' 
        self.restart_path = self._outpath + "restart/"
        self._minima_path = 'minima/'
        self._is_restart = file_handling.restart(self._outpath, self.restart_path, self._minima_path)
        # Check if new start is desired
        if self._new_start:
            self._is_restart = False

        # initialize database
        with Database(self._energy_threshold, self._minima_threshold, self._is_restart, self.restart_path) as self.data:
            # Start up minimahopping 
            atoms = deepcopy(self._atoms)
            struct = self._startup(atoms,)  # gets an atoms object and a minimum object is returned.
            struct_cur = struct.__deepcopy__()

            # Start hopping loop
            while (self._counter <= totalsteps):
                msg = "START HOPPING STEP NR.  {:d}".format(self._counter)
                print(msg)
                print("  Start escape loop")
                print("  ---------------------------------------------------------------")
                struct_inter, n_visits, _epot_max, _md_trajectory, _opt_trajectory = self._escape(struct,) 
                print("  ---------------------------------------------------------------")
                print("  New minimum found!")

                # write the lowest n minima to files
                self.data._write_poslow(self._n_poslow ,self._minima_path)

                # write output
                self._hoplog(struct_inter)
                status = 'Inter'
                self._history_log(struct_inter, status, n_visits)

                # check if new proposed structure is accepted or rejected
                struct_cur, is_accepted =  self._acc_rej_step(struct_cur, struct_inter)

                # if structure is not accepted store the first rejected structre as currently lowest rejected structure
                if not is_accepted :
                    status = "Rejected"
                else:
                    status = "Accepted"
                
                
                log_msg = "  New minimum has been found {:d} time(s)".format(n_visits)
                print(log_msg)
                # adjust the temperature according to the number of visits
                self._adj_temperature(struct, n_visits)

                # write to history output file
                self._history_log(struct_inter, status)

                # if not accepted the first rejected is the current lowest rejected
                if not is_accepted:
                    lowest_rej = struct_inter.__deepcopy__()
                
                # if not accepted start reject loop until a minimum has been accepted
                while not is_accepted:
                    msg = "START HOPPING STEP NR.  {:d}".format(self._counter)
                    print(msg)
                    print("  Start escape loop")
                    print("  ---------------------------------------------------------------")
                    struct_inter, n_visits, _epot_max, _md_trajectory, _opt_trajectory = self._escape(struct,) 
                    print("  ---------------------------------------------------------------")
                    print("  New minimum found!")
                    # write the lowest n minima to files
                    self.data._write_poslow(self._n_poslow ,self._minima_path)

                    # write output
                    self._hoplog(struct_inter)
                    status = 'Inter'
                    self._history_log(struct_inter, status, n_visits)
                    struct_cur, is_accepted =  self._acc_rej_step(struct_cur, struct_inter)

                    if not is_accepted :
                        status = "Rejected"
                    else:
                        status = "Accepted"

                    epot_prev = lowest_rej.e_pot
                    epot_new = struct_inter.e_pot
                    if epot_new < epot_prev:
                        # if new found minimum is lower in energy that's the lowest rejected minimum 
                        lowest_rej = struct_inter.__deepcopy__()
                    struct_inter = lowest_rej.__deepcopy__()

                    log_msg = "  New minimum has been found {:d} time(s)".format(n_visits)
                    print(log_msg)
                    # adjust the temperature according to the number of visits
                    self._adj_temperature(struct, n_visits)

                    # write to history output file
                    self._history_log(struct_inter, status)
                    self._counter += 1


                struct = struct_cur.__deepcopy__()
                self._counter += 1
                print("DONE")
                print("=================================================================")

                _elapsed_time = time.time() - self._time_in

                if self._run_time is not "infinit":
                    if _elapsed_time > self._run_time_sec:
                        msg = 'Simulation stopped because the given time is over\n'
                        msg += 'Run terminated after {:d} steps'.format(self._counter)
                        print(msg)
                        print("=================================================================")
                        return

            _elapsed_time = time.time() - self._time_in
            day = _elapsed_time // (24 * 3600)
            _elapsed_time = _elapsed_time % (24 * 3600)
            hour = _elapsed_time // 3600
            _elapsed_time %= 3600
            minutes = _elapsed_time // 60
            _elapsed_time %= 60
            seconds = _elapsed_time
            msg = 'Run terminated after {:d} steps in {:d}D {:d}H {:d}M {:d}S'.format(totalsteps,
                                                                                    int(day),
                                                                                    int(hour),
                                                                                    int(minutes),
                                                                                    int(seconds))
            print(msg)


        
        

    def _startup(self, atoms):
        print("=================================================================")
        print("MINIMAHOPPING SETUP START")

        # Convert given time to seconds
        if self._run_time is not "infinit":
            self._run_time_sec = self._get_sec()
        self._time_in = time.time()

        # Initialization of global counters
        
        self._i_step = 0
        self._n_min = 1
        self._n_unique = 0
        self._n_notunique = 0
        self._n_same = 0

        # Check if restart
        if self._is_restart:
            msg = '  Restart MH run'
            print(msg)

            # Get calculator from input structure
            calc = atoms.calc

            # Read current structure
            filename = self.restart_path + "poscur.extxyz"
            atoms = read(filename)
            atoms.calc = calc

            label = self.data.nstructs + 1
            struct_cur = Minimum(atoms,
                        n_visit=1,
                        s = self._ns_orb,
                        p = self._np_orb, 
                        width_cutoff = self._width_cutoff,
                        maxnatsphere = self._maxnatsphere,
                        epot = atoms.get_potential_energy(),
                        T=self._temperature,
                        ediff=self._Ediff,
                        label=label)

            # Read parameters and set them
            filename = self.restart_path + "params.json"
            f = open(filename)
            self.restart_dict = json.load(f)
            f.close()
            self._dt = self.restart_dict["dt"]
            self._Ediff = self.restart_dict["Ediff"]
            self._temperature = self.restart_dict["T"]

        else:

            msg = '  New MH run is started'
            print(msg)

            _positions, _lattice = self._restart_opt(atoms)
            atoms.set_positions(_positions)
            atoms.set_cell(_lattice)
            write(self._outpath + "acc.extxyz", atoms, append=True)
            write(self.restart_path + "poscur.extxyz", atoms)

            # add input structure to database after optimization

            struct_cur = Minimum(atoms,
                        n_visit=1,
                        s = self._ns_orb,
                        p = self._np_orb, 
                        width_cutoff = self._width_cutoff,
                        maxnatsphere = self._maxnatsphere,
                        epot = atoms.get_potential_energy(),
                        T=self._temperature,
                        ediff=self._Ediff,
                        label=0)
        
        
        status = 'Initial'
        self._history_log(struct_cur, status, n_visits=1)

        self.data.addElement(struct_cur)

        print("DONE")
        print("=================================================================")
        return struct_cur



    def _restart_opt(self, atoms,):
        _atoms = deepcopy(atoms)
        opt = Opt(atoms=_atoms, outpath=self._outpath,max_froce_threshold=self._fmax, verbose=self._verbose)
        if True in atoms.pbc:
            _positions, _lattice, _noise, _opt_trajectory = opt.run()
        else:
            _positions, _noise, _opt_trajectory = opt.run()
            _lattice = np.zeros((3,3))
        return _positions, _lattice

    def _get_sec(self,):
        """Get seconds from time."""
        nd, d = self._run_time.split('-')
        h, m, s = d.split(':')
        return int(nd) * 86400 + int(h) * 3600 + int(m) * 60 + int(s)


    def _escape(self, struct):
        """
        Escape loop to find a new minimum
        """
        _escape = 0.0
        _escape_energy = 0.0
        atoms = deepcopy(struct.atoms)
        _beta_s = 1.05
        _i_steps = 0
        while _escape < self._minima_threshold or _escape_energy < self._energy_threshold:
            # if the loop not escaped (no new minimum found) rise temperature
            if _i_steps > 0:
                self._n_same += 1
                status = 'Same'
                self._history_log(struct_prop, status)
                self._temperature *= _beta_s
                log_msg = "    Same minimum found with fpd {:1.2e} {:d} time(s). Increase temperature to {:1.5f}".format(_escape, self._n_same, self._temperature)
                print(log_msg)

            # set the temperature according to Boltzmann distribution
            MaxwellBoltzmannDistribution(atoms, temperature_K=self._temperature)

            # in case of periodic system do variable cell shape md and optimization
            if True in atoms.pbc:

                # initialize the cell vectors as atoms
                _mass = .75 * np.sum(atoms.get_masses()) / 10.
                self._cell_atoms = Cell_atom(mass=_mass, positions=atoms.get_cell())
                self._cell_atoms.set_velocities_boltzmann(temperature=self._temperature)

                # softening of the velocities
                softening = Softening(atoms, self._cell_atoms)
                _velocities, _cell_velocities = softening.run(self._n_soft)
                atoms.set_velocities(_velocities)
                self._cell_atoms.velocities = _cell_velocities

                print("    VCS MD Start")

                md = MD(atoms=atoms, outpath=self._outpath, cell_atoms=self._cell_atoms, dt=self._dt, n_max=self._mdmin, verbose=self._verbose)
                _positions, _cell , self._dt, _md_trajectory, _epot_max = md.run()

                log_msg = "    VCS MD finished after {:d} steps visiting {:d} maxima. New dt is {:1.5f}".format(md._i_steps, self._mdmin, self._dt)
                print(log_msg)
                atoms.set_positions(_positions)
                atoms.set_cell(_cell)

                atoms = lat_opt.reshape_cell2(atoms, 6)

                print("    VCS OPT start")

                opt = Opt(atoms=atoms, outpath=self._outpath, max_froce_threshold=self._fmax, verbose=self._verbose)
                _positions, _lattice, self._noise, _opt_trajectory = opt.run()
                atoms.set_positions(_positions)
                atoms.set_cell(_lattice)

                log_msg = "    VCS OPT finished after {:d} steps         {:d}".format(opt._i_step, len(_opt_trajectory))
                print(log_msg)

            # in case of a non-periodic system do md and optimization
            else:
                #start softening
                softening = Softening(atoms)
                _velocities = softening.run(self._n_soft)
                atoms.set_velocities(_velocities)

                print("    MD Start")

                md = MD(atoms=atoms, outpath=self._outpath, cell_atoms=None, dt=self._dt, n_max=self._mdmin, verbose=self._verbose)
                _positions , self._dt, _md_trajectory, _epot_max = md.run()
                atoms.set_positions(_positions)

                log_msg = "    MD finished after {:d} steps visiting {:d} maxima. New dt is {:1.5f}".format(md._i_steps, self._mdmin, self._dt)
                print(log_msg)
                print("    OPT start")
                opt = Opt(atoms=atoms, outpath=self._outpath, max_froce_threshold=self._fmax, verbose=self._verbose)
                _positions, self._noise, _opt_trajectory = opt.run()
                atoms.set_positions(_positions)
                log_msg = "    OPT finished after {:d} steps".format(opt._i_step, len(_opt_trajectory))
                print(log_msg)

            # check if the energy threshold is below the optimization noise
            self._check_energy_threshold()

            struct_prop = Minimum(atoms,
                        n_visit=1,
                        s = self._ns_orb,
                        p = self._np_orb, 
                        width_cutoff = self._width_cutoff,
                        maxnatsphere = self._maxnatsphere,
                        epot = atoms.get_potential_energy(),
                        T=self._temperature,
                        ediff=self._Ediff,
                        label=0)

            # check if proposed structure is the same to the initial structure
            _escape_energy = struct.__compareto__(struct_prop)
            _escape = struct.__equals__(struct_prop)

            write(self._outpath + 'locm.extxyz', atoms, append=True)

            _i_steps += 1
            self._n_min += 1

            # update and save restart dict
            self.restart_dict["dt"] = self._dt
            self.restart_dict["T"] = self._temperature
            f = open(self.restart_path+"params.json", "w")
            json.dump(self.restart_dict,f)
            f.close()
        

        log_msg = "    New minimum found with fpd {:1.2e} after looping {:d} time(s)".format(_escape, _i_steps)
        print(log_msg)

        # add new minimum to database
        
        n_visits = self.data.addElement(struct = struct_prop)

        return struct_prop,n_visits,_epot_max, _md_trajectory, _opt_trajectory


    def _hoplog(self, struct):
        atoms = struct.atoms
        self._i_step += 1
        log_msg = "  Epot:  {:1.5f}   E_diff:  {:1.5f}    Temp:   {:1.5f} ".format(atoms.get_potential_energy(),
                                                                                             self._Ediff,
                                                                                             self._temperature)
        print(log_msg)



    def _acc_rej_step(self, struct_cur, struct):

        _e_pot_cur = struct_cur.e_pot
        _e_pot = struct.e_pot

        _ediff_in = self._Ediff

        if _e_pot - _e_pot_cur < self._Ediff:
            self._Ediff *= self._alpha_a
            struct_cur = struct.__deepcopy__()
            is_accepted = True
            ediff_acc = _e_pot - _e_pot_cur
        else:
            self._Ediff *= self._alpha_r
            is_accepted = False
            ediff_rej = _e_pot - _e_pot_cur


        if is_accepted:
            log_msg = "  Minimum was accepted:  Enew - Ecur = {:1.5f} < {:1.5f} = Ediff".format(ediff_acc,
                                                                                                _ediff_in)
            print(log_msg)
            write(self._outpath + "acc.extxyz", struct_cur.atoms, append=True)
            write(self.restart_path + "poscur.extxyz", struct_cur.atoms)
        else:
            log_msg = "  Minimum was rejected:  Enew - Ecur = {:1.5f} > {:1.5f} = Ediff".format(ediff_rej,
                                                                                                _ediff_in)
            print(log_msg)

        # update restart dict
        self.restart_dict["Ediff"] = self._Ediff
        f = open(self.restart_path+"params.json", "w")
        json.dump(self.restart_dict,f)
        f.close()

        return struct_cur, is_accepted



    def _adj_temperature(self,struct, n_visits):
        if n_visits > 1:
            self._n_notunique += 1
            if self._enhanced_feedback:
                self._temperature = self._temperature * self._beta_increase * (1. + 1. * np.log(float(n_visits)))
            else:
                self._temperature = self._temperature * self._beta_increase
        else:
            self._n_unique += 1
            self._temperature = self._temperature * self._beta_decrease
            atoms = struct.atoms
            write(self._outpath + "min.extxyz", atoms,  append=True)


    def _history_log(self, struct, status, n_visits = 0):
        atoms = struct.atoms
        _notunique_frac = float(self._n_notunique)/float(self._n_min)
        _same_frac = float(self._n_same)/float(self._n_min)
        _unique_frac = 1. - (_notunique_frac+_same_frac)

        history_msg = "{:1.9f}  {:d}  {:1.5f}  {:1.5f}  {:1.2f}  {:1.2f}  {:1.2f} {:s} \n".format(atoms.get_potential_energy(),
                                                                        n_visits,
                                                                        self._temperature,
                                                                        self._Ediff,
                                                                        _same_frac,
                                                                        _notunique_frac,
                                                                        _unique_frac,
                                                                        status)

        history_file = open(self._outpath + 'history.dat', 'a')
        history_file.write(history_msg)
        history_file.close()


    def _check_energy_threshold(self):
        if self._energy_threshold < self._noise:
            _warning_msg = 'Energy threshold is below the noise level'
            warnings.warn(_warning_msg, UserWarning)



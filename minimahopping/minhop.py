import numpy as np
import warnings
from ase.io import read, write
import ase.atom
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import minimahopping.mh.lattice_operations as lat_opt
from copy import deepcopy
from minimahopping.md.soften import Softening
from minimahopping.md.md import MD
from minimahopping.opt.optim import Opt
from minimahopping.mh.minimum import Minimum
from minimahopping.mh.cell_atom import Cell_atom
from minimahopping.mh.database import Database
import time
import json
import minimahopping.mh.file_handling as file_handling 
import minimahopping.graph.graph as graph

"""
MH Software written by Marco Krummenacher (marco.krummenacher@unibas.ch)
Parts of the software were originally developped (some in Fortran) from other people:
  -- VCSMD: Martin Sommer-Joergenson
  -- VCS Softening: Hannes Huber
  -- VCS optimizer: Moritz Gubler
  -- OMFP in python: Jonas Finkler
"""


class Minimahopping:
    def __init__(self, initial_configuration : ase.atom.Atom,
                        T0 = 1000,
                        beta_decrease = 1./1.05,
                        beta_increase = 1.05,
                        Ediff0 = .01,
                        alpha_accept = 0.95,
                        alpha_reject = 1.05,
                        n_soft = 20,
                        n_S_orbitals = 1, 
                        n_P_orbitals = 1, 
                        width_cutoff = 3.5, 
                        maxnatsphere = 100, 
                        exclude = [],
                        dt = 0.01,
                        mdmin = 2,
                        fmax = 0.001, 
                        enhanced_feedback = False,
                        energy_threshold = 0.001, #5 the noise
                        output_n_lowest_minima = 30,
                        fingerprint_threshold = 1e-3,
                        verbose_output = True,
                        new_start = False,
                        run_time = 'infinite', 
                        use_intermediate_mechanism = False,
                        restart_interval_minutes = 2,
                        overwriteParametersOnRestart = False):
        """Initialize with an ASE atoms object and keyword arguments."""

        self.initial_configuration = initial_configuration

        # time of last restart
        self._last_restart = time.time()

        self._outpath = 'output/' 
        self.restart_path = self._outpath + "restart/"
        self._minima_path = 'minima/'

        self.isRestart = self.checkIfRestart(new_start)

        if self.isRestart and not overwriteParametersOnRestart:
            self._read_changing_parameters("params.json")

        else:
            self.parameter_dictionary = {
                "beta_decrease" : beta_decrease,
                "beta_increase" : beta_increase,
                "alpha_accept" : alpha_accept,
                "alpha_reject" : alpha_reject,
                "n_softening_steps" : n_soft,
                "n_S_orbitals" : n_S_orbitals,
                "n_P_orbitals" : n_P_orbitals,
                "width_cutoff": width_cutoff,
                "max_atoms_in_cutoff_sphere" : maxnatsphere,
                "dt" : dt,
                "mdmin" : mdmin,
                "fmax" : fmax,
                "enhanced_feedback" : enhanced_feedback,
                "energy_threshold" : energy_threshold,
                "output_n_lowest_minima": output_n_lowest_minima,
                "fingerprint_threshold" : fingerprint_threshold,
                "verbose_output" : verbose_output,
                "T" : T0, 
                "energy_difference_to_accept" : Ediff0,
                "exclude" : exclude,
                "run_time" : run_time,
                "use_intermediate_mechanism" : use_intermediate_mechanism,
                "restart_interval_minutes": restart_interval_minutes,
            }

        # Initialization of global counters        
        self._n_min = 1
        self._n_notunique = 0
        self._n_same = 0


    def __call__(self, totalsteps = None):
        counter = 0

        # initialize database
        with Database(self.parameter_dictionary["energy_threshold"], self.parameter_dictionary["fingerprint_threshold"], self.isRestart, self.restart_path) as self.data, graph.MinimaHoppingGraph('output/graph.dat', 'output/trajectory.dat', self.isRestart) as self.mh_graph:
            # Start up minimahopping 
            atoms = deepcopy(self.initial_configuration)
            current_minimum = self._startup(atoms,)  # gets an atoms object and a minimum object is returned.
            #struct_cur = struct.__deepcopy__()

            # Start hopping loop
            while (counter <= totalsteps):
                is_accepted = False
                is_first_accept_iteration = True
                
                # if not accepted start reject loop until a minimum has been accepted
                msg = "START HOPPING STEP NR.  {:d}".format(counter)
                print(msg)
                while not is_accepted:
                    print("  Start escape loop")
                    print("  ---------------------------------------------------------------")
                    escaped_minimum, n_visits, epot_max, md_trajectory, opt_trajectory = self._escape(current_minimum) 
                    print("  ---------------------------------------------------------------")
                    print("  New minimum found!")
                    # write the lowest n minima to files
                    self.data._write_poslow(self.parameter_dictionary["output_n_lowest_minima"] ,self._minima_path)

                    # add new structure as edge to graph
                    print('current_minimum.label, escaped_minimum.label', current_minimum.label, escaped_minimum.label)
                    self.mh_graph.addStructure(current_minimum.label, escaped_minimum.label, md_trajectory + opt_trajectory, current_minimum.e_pot, escaped_minimum.e_pot, epot_max)

                    # write output
                    self._hoplog(escaped_minimum)

                    # check if intermediate or escaped minimum should be considered for accepting.
                    if is_first_accept_iteration or not self.parameter_dictionary["use_intermediate_mechanism"]: # after first iteration, the escaped minimum must be considered for accecpting
                        intermediate_minimum = escaped_minimum.__deepcopy__()
                        intermediate_minimum_is_escaped_minimum = True
                        is_first_accept_iteration = False
                    else: # After several iterations, the check if the escaped or the intermediate minimum has a lower energy
                        intermediate_minimum_is_escaped_minimum = False
                        if escaped_minimum.e_pot < intermediate_minimum.e_pot:
                            intermediate_minimum = escaped_minimum.__deepcopy__()
                            intermediate_minimum_is_escaped_minimum = True
                    # accept minimum if necessary
                    current_minimum, is_accepted =  self._acc_rej_step(current_minimum, intermediate_minimum)

                    # write log messages
                    if is_accepted:
                        if intermediate_minimum_is_escaped_minimum:
                            status = "Accepted minimum after escaping"
                            self._history_log(escaped_minimum, status, n_visits)
                        else:
                            status = 'Accepted intermediate minimum'
                            self._history_log(intermediate_minimum, status, n_visits)
                    else:
                        if intermediate_minimum_is_escaped_minimum:
                            status = "Rejected -> new intermediate minimum"
                        else:
                            status = "Rejected"
                        self._history_log(escaped_minimum, status, n_visits)

                    log_msg = "  New minimum has been found {:d} time(s)".format(n_visits)
                    print(log_msg)

                    # adjust the temperature according to the number of visits
                    self._adj_temperature(escaped_minimum, n_visits)
                    temp_atoms_towrite = escaped_minimum.atoms.copy()
                    temp_atoms_towrite.info = {}
                    temp_atoms_towrite.set_momenta(None)
                    temp_atoms_towrite.info['energy'] = escaped_minimum.e_pot
                    temp_atoms_towrite.info['label'] = escaped_minimum.label
                    write(self._outpath + "min.extxyz", temp_atoms_towrite,  append=True)
                    self._write_restart_if_necessary()

                counter += 1
                print("DONE")
                print("=================================================================")

                _elapsed_time = time.time() - self._time_in

                if self.parameter_dictionary["run_time"] != "infinite":
                    if _elapsed_time > self._run_time_sec:
                        msg = 'Simulation stopped because the given time is over\n'
                        msg += 'Run terminated after {:d} steps'.format(counter)
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


    def _write_restart(self):
        self.mh_graph.write_restart_files()
        self.data.write_restart_files()
        self._last_restart = time.time()
        # TODO write poscur and params.json


    def _write_restart_if_necessary(self):
        minutes_since_restart = (self._last_restart - time.time()) / 60.0
        if minutes_since_restart > self.parameter_dictionary["restart_interval_minutes"]:
            self._write_restart()


    def _startup(self, atoms):
        print("=================================================================")
        print("MINIMAHOPPING SETUP START")


        if not isinstance(atoms, list):
            atoms = [deepcopy(atoms)]

        for atom in atoms:
            calc = atom.calc
            _positions, _lattice = self._restart_opt(atom,)
            atom.set_positions(_positions)
            atom.set_cell(_lattice)

            struct = Minimum(atom,
                        n_visit=1,
                        s = self.parameter_dictionary['n_S_orbitals'],
                        p = self.parameter_dictionary["n_P_orbitals"], 
                        width_cutoff = self.parameter_dictionary["width_cutoff"],
                        maxnatsphere = self.parameter_dictionary["max_atoms_in_cutoff_sphere"],
                        epot = atom.get_potential_energy(),
                        T=self.parameter_dictionary['T'],
                        ediff=self.parameter_dictionary["energy_difference_to_accept"],
                        label=0,
                        exclude= self.parameter_dictionary["exclude"])
            
            self.data.addElement(struct)

        # Convert given time to seconds
        if self.parameter_dictionary["run_time"] != "infinite":
            self._run_time_sec = self._get_sec()
        self._time_in = time.time()

        # Check if restart
        if self.isRestart:
            print('  Restart MH run')

            # Read current structure
            filename = self.restart_path + "poscur.extxyz"
            atoms = read(filename)
            atoms.calc = calc
            
            struct_cur = Minimum(atoms,
                        n_visit=1,
                        s = self.parameter_dictionary['n_S_orbitals'],
                        p = self.parameter_dictionary["n_P_orbitals"], 
                        width_cutoff = self.parameter_dictionary["width_cutoff"],
                        maxnatsphere = self.parameter_dictionary["max_atoms_in_cutoff_sphere"],
                        epot = atoms.get_potential_energy(),
                        T=self.parameter_dictionary['T'],
                        ediff=self.parameter_dictionary["energy_difference_to_accept"],
                        label=-1,
                        exclude= self.parameter_dictionary["exclude"])
            index = self.data.get_element(struct_cur)
            struct_cur.set_label(self.data.unique_minima_sorted[index].label)
        else:

            msg = '  New MH run is started'
            print(msg)

            # add input structure to database after optimization
            struct_cur = self.data.unique_minima_sorted[0].__copy__()
            struct_cur.atoms.calc = calc
            atoms = struct_cur.atoms
            write(self._outpath + "acc.extxyz", atoms, append=True)
            write(self.restart_path + "poscur.extxyz", atoms)


        self.data.addElement(struct_cur)
        status = 'Initial'
        self._history_log(struct_cur, status, n_visits=struct_cur.n_visit)

        print("DONE")
        print("=================================================================")
        return struct_cur


    def _read_changing_parameters(self, param_dict_name):
        # Read parameters and set them
        filename = self.restart_path + param_dict_name
        f = open(filename)
        self.parameter_dictionary = json.load(f)
        f.close()


    def _restart_opt(self, atoms,):
        _atoms = deepcopy(atoms)
        opt = Opt(atoms=_atoms, outpath=self._outpath,max_froce_threshold=self.parameter_dictionary["fmax"], verbose=self.parameter_dictionary["verbose_output"])
        if True in atoms.pbc:
            _positions, _lattice, _noise, _opt_trajectory = opt.run()
        else:
            _positions, _noise, _opt_trajectory = opt.run()
            _lattice = np.zeros((3,3))
        return _positions, _lattice

    def _get_sec(self,):
        """Get seconds from time."""
        nd, d = self.parameter_dictionary["run_time"].split('-')
        h, m, s = d.split(':')
        return int(nd) * 86400 + int(h) * 3600 + int(m) * 60 + int(s)


    def _escape(self, struct):
        """
        Escape loop to find a new minimum
        """
        self._n_same = 0
        _escape = 0.0
        _escape_energy = 0.0
        atoms = deepcopy(struct.atoms)
        _beta_s = 1.05
        _i_steps = 0
        while _escape < self.parameter_dictionary["fingerprint_threshold"] or _escape_energy < self.parameter_dictionary["energy_threshold"]:
            # if the loop not escaped (no new minimum found) rise temperature
            if _i_steps > 0:
                self._n_same += 1
                status = 'Same'
                self._history_log(struct_prop, status)
                self.parameter_dictionary['T'] *= _beta_s
                log_msg = "    Same minimum found with fpd {:1.2e} {:d} time(s). Increase temperature to {:1.5f}".format(_escape, self._n_same, self.parameter_dictionary['T'])
                print(log_msg)

            # set the temperature according to Boltzmann distribution
            MaxwellBoltzmannDistribution(atoms, temperature_K=self.parameter_dictionary['T'])

            # in case of periodic system do variable cell shape md and optimization
            if True in atoms.pbc:

                # initialize the cell vectors as atoms
                _mass = .75 * np.sum(atoms.get_masses()) / 10.
                self._cell_atoms = Cell_atom(mass=_mass, positions=atoms.get_cell())
                self._cell_atoms.set_velocities_boltzmann(temperature=self.parameter_dictionary['T'])

                # softening of the velocities
                softening = Softening(atoms, self._cell_atoms)
                _velocities, _cell_velocities = softening.run(self.parameter_dictionary['n_softening_steps'])
                atoms.set_velocities(_velocities)
                self._cell_atoms.velocities = _cell_velocities

                print("    VCS MD Start")

                md = MD(atoms=atoms, outpath=self._outpath, cell_atoms=self._cell_atoms, dt=self.parameter_dictionary['dt'], n_max=self.parameter_dictionary["mdmin"], verbose=self.parameter_dictionary["verbose_output"])
                _positions, _cell , self.parameter_dictionary['dt'], _md_trajectory, _epot_max = md.run()

                log_msg = "    VCS MD finished after {:d} steps visiting {:d} maxima. New dt is {:1.5f}".format(md._i_steps, self.parameter_dictionary["mdmin"], self.parameter_dictionary['dt'])
                print(log_msg)
                atoms.set_positions(_positions)
                atoms.set_cell(_cell)

                atoms = lat_opt.reshape_cell2(atoms, 6)

                print("    VCS OPT start")

                opt = Opt(atoms=atoms, outpath=self._outpath, max_froce_threshold=self.parameter_dictionary["fmax"], verbose=self.parameter_dictionary["verbose_output"])
                _positions, _lattice, self._noise, _opt_trajectory = opt.run()
                atoms.set_positions(_positions)
                atoms.set_cell(_lattice)

                log_msg = "    VCS OPT finished after {:d} steps         {:d}".format(opt._i_step, len(_opt_trajectory))
                print(log_msg)

            # in case of a non-periodic system do md and optimization
            else:
                #start softening
                softening = Softening(atoms)
                _velocities = softening.run(self.parameter_dictionary['n_softening_steps'])
                atoms.set_velocities(_velocities)

                print("    MD Start")

                md = MD(atoms=atoms, outpath=self._outpath, cell_atoms=None, dt=self.parameter_dictionary['dt'], n_max=self.parameter_dictionary["mdmin"], verbose=self.parameter_dictionary["verbose_output"])
                _positions , self.parameter_dictionary['dt'], _md_trajectory, _epot_max = md.run()
                atoms.set_positions(_positions)
                log_msg = "    MD finished after {:d} steps visiting {:d} maxima. New dt is {:1.5f}".format(md._i_steps, self.parameter_dictionary["mdmin"], self.parameter_dictionary['dt'])
                print(log_msg)
                print("    OPT start")
                opt = Opt(atoms=atoms, outpath=self._outpath, max_froce_threshold=self.parameter_dictionary["fmax"], verbose=self.parameter_dictionary["verbose_output"])
                _positions, self._noise, _opt_trajectory = opt.run()
                atoms.set_positions(_positions)
                log_msg = "    OPT finished after {:d} steps".format(opt._i_step, len(_opt_trajectory))
                print(log_msg)
            # check if the energy threshold is below the optimization noise
            self._check_energy_threshold()

            struct_prop = Minimum(atoms,
                        n_visit=1,
                        s = self.parameter_dictionary['n_S_orbitals'],
                        p = self.parameter_dictionary["n_P_orbitals"], 
                        width_cutoff = self.parameter_dictionary["width_cutoff"],
                        maxnatsphere = self.parameter_dictionary["max_atoms_in_cutoff_sphere"],
                        epot = atoms.get_potential_energy(),
                        T=self.parameter_dictionary['T'],
                        ediff=self.parameter_dictionary["energy_difference_to_accept"],
                        label=0,
                        exclude= self.parameter_dictionary["exclude"])

            # check if proposed structure is the same to the initial structure
            _escape_energy = struct.__compareto__(struct_prop)
            _escape = struct.__equals__(struct_prop)

            write(self._outpath + 'locm.extxyz', atoms, append=True)

            _i_steps += 1
            self._n_min += 1

            # update and save restart dict
            self.parameter_dictionary["dt"] = self.parameter_dictionary['dt']
            self.parameter_dictionary["T"] = self.parameter_dictionary['T']
            f = open(self.restart_path+"params.json", "w")
            json.dump(self.parameter_dictionary,f)
            f.close()
        

        log_msg = "    New minimum found with fpd {:1.2e} after looping {:d} time(s)".format(_escape, _i_steps)
        print(log_msg)

        # add new minimum to database
        
        self.data.addElement(struct = struct_prop)

        return struct_prop, struct_prop.n_visit, _epot_max, _md_trajectory, _opt_trajectory


    def _hoplog(self, struct):
        atoms = struct.atoms
        log_msg = "  Epot:  {:1.5f}   E_diff:  {:1.5f}    Temp:   {:1.5f} ".format(atoms.get_potential_energy(),
                                                                                             self.parameter_dictionary["energy_difference_to_accept"],
                                                                                             self.parameter_dictionary['T'])
        print(log_msg)



    def _acc_rej_step(self, struct_cur, struct):

        _e_pot_cur = struct_cur.e_pot
        _e_pot = struct.e_pot

        _ediff_in = self.parameter_dictionary["energy_difference_to_accept"]

        if _e_pot - _e_pot_cur < self.parameter_dictionary["energy_difference_to_accept"]:
            self.parameter_dictionary["energy_difference_to_accept"] *= self.parameter_dictionary['alpha_accept']
            struct_cur = struct.__deepcopy__()
            is_accepted = True
            ediff_acc = _e_pot - _e_pot_cur
        else:
            self.parameter_dictionary["energy_difference_to_accept"] *= self.parameter_dictionary['alpha_reject']
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
        self.parameter_dictionary["Ediff"] = self.parameter_dictionary["energy_difference_to_accept"]
        f = open(self.restart_path+"params.json", "w")
        json.dump(self.parameter_dictionary,f)
        f.close()

        return struct_cur, is_accepted



    def _adj_temperature(self,struct, n_visits):
        if n_visits > 1:
            self._n_notunique += 1
            if self.parameter_dictionary["enhanced_feedback"]:
                self.parameter_dictionary['T'] = self.parameter_dictionary['T'] * self.parameter_dictionary['beta_increase'] * (1. + 1. * np.log(float(n_visits)))
            else:
                self.parameter_dictionary['T'] = self.parameter_dictionary['T'] * self.parameter_dictionary['beta_increase']


    def _history_log(self, struct, status, n_visits = 0):
        atoms = struct.atoms
        _notunique_frac = float(self._n_notunique)/float(self._n_min)
        _same_frac = float(self._n_same)/float(self._n_min)
        _unique_frac = 1. - (_notunique_frac+_same_frac)

        history_msg = "{:1.9f}  {:d}  {:1.5f}  {:1.5f}  {:1.2f}  {:1.2f}  {:1.2f} {:s} \n".format(atoms.get_potential_energy(),
                                                                        n_visits,
                                                                        self.parameter_dictionary['T'],
                                                                        self.parameter_dictionary["energy_difference_to_accept"],
                                                                        _same_frac,
                                                                        _notunique_frac,
                                                                        _unique_frac,
                                                                        status)

        history_file = open(self._outpath + 'history.dat', 'a')
        history_file.write(history_msg)
        history_file.close()


    def _check_energy_threshold(self):
        if self.parameter_dictionary["energy_threshold"] < self._noise:
            _warning_msg = 'Energy threshold is below the noise level'
            warnings.warn(_warning_msg, UserWarning)

    def checkIfRestart(self, isNewStart):
        isRestart = file_handling.restart(self._outpath, self.restart_path, self._minima_path)
        # Check if new start is desired
        if isNewStart:
            isRestart = False
        return isRestart

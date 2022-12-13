import numpy as np
import warnings
from ase.io import read
import ase.atom
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import minimahopping.mh.lattice_operations as lat_opt
from copy import deepcopy
from minimahopping.md.soften import Softening
from minimahopping.md.md import MD
from minimahopping.opt.optim import Opt
from minimahopping.mh.minimum import Minimum
from minimahopping.mh.cell_atom import Cell_atom
import time
import json
import minimahopping.mh.file_handling as file_handling 
import minimahopping.MPI_database.mpi_database_master as MPI_server

"""
MH Software written by Marco Krummenacher (marco.krummenacher@unibas.ch)
Parts of the software were originally developped (some in Fortran) from other people:
  -- VCSMD: Martin Sommer-Joergenson
  -- VCS Softening: Hannes Huber
  -- VCS optimizer: Moritz Gubler
  -- OMFP in python: Jonas Finkler
"""


def importer(name, root_package=False, relative_globals=None, level=0):
    """ We only import modules, functions can be looked up on the module.
    Usage: 

    from foo.bar import baz
    >>> baz = importer('foo.bar.baz')

    import foo.bar.baz
    >>> foo = importer('foo.bar.baz', root_package=True)
    >>> foo.bar.baz

    from .. import baz (level = number of dots)
    >>> baz = importer('baz', relative_globals=globals(), level=2)
    """
    return __import__(name, locals=None, # locals has no use
                      globals=relative_globals, 
                      fromlist=[] if root_package else [None],
                      level=level)

class Minimahopping:
    def __init__(self, initial_configuration : ase.atom.Atom,
                        T0 = 1000,
                        beta_decrease = 1./1.02,
                        beta_increase = 1.02,
                        Ediff0 = .5,
                        alpha_accept = 1/1.02,
                        alpha_reject = 1.02,
                        n_soft = 50,
                        n_S_orbitals = 1, 
                        n_P_orbitals = 1, 
                        width_cutoff = 3.5, 
                        maxnatsphere = 100, 
                        exclude = [],
                        dt = 0.01,
                        mdmin = 5,
                        fmax = 0.001,
                        enhanced_feedback = False,
                        energy_threshold = 0.0001, #5 the noise
                        output_n_lowest_minima = 20,
                        fingerprint_threshold = 1e-3,
                        verbose_output = True,
                        new_start = False,
                        run_time = 'infinite', 
                        use_intermediate_mechanism = True,
                        overwriteParametersOnRestart = False,
                        write_graph_output = True,
                        use_MPI = False):
        """Initialize with an ASE atoms object and keyword arguments."""

        self.initial_configuration = initial_configuration

        try:
            from mpi4py import MPI
            self.mpiRank = MPI.COMM_WORLD.Get_rank()
            self.mpiSize = MPI.COMM_WORLD.Get_size()
        except ImportError:
            if use_MPI:
                print("use_MPI was set to true but mpi4py is not installed. Install mpi4py and start again.")
                quit()
            self.mpiRank = 0
            self.mpiSize = 1
        
        if self.mpiSize == 1 or not use_MPI:
            self.isMaster = False
            self.database = importer("minimahopping.mh.database")
            self._outpath = 'output/' 
            self.restart_path = self._outpath + "restart/"
            self._minima_path = 'minima/'
        else:
            if self.mpiRank == 0:
                self.isMaster = True
                self._outpath = 'output/master/'
                self.restart_path = self._outpath + "restart/"
                self._minima_path = 'minima/'
                
            else:
                self.isMaster = False
                self.database = importer("minimahopping.MPI_database.mpi_database_worker")
                self._outpath = 'output/worker_' + str(self.mpiRank) + '/' 
                self.restart_path = self._outpath + "restart/"
                self._minima_path = 'minima/worker_' + str(self.mpiRank) + '/'
                


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
                "write_graph_output" : write_graph_output,
                "use_MPI": use_MPI
            }

        if self.isMaster:
            f = open(self.restart_path+"params.json", "w")
            json.dump(self.parameter_dictionary,f)
            f.close()

        # Initialization of global counters        
        self._n_min = 1
        self._n_notunique = 0
        self._n_same = 0


    def __call__(self, totalsteps = None):
        counter = 0

        if self.parameter_dictionary['use_MPI'] and self.mpiRank == 0:
            print('starting mpi server on rank', self.mpiRank)
            MPI_server.MPI_database_server_loop(self.parameter_dictionary['energy_threshold'], self.parameter_dictionary['fingerprint_threshold']
                , self.parameter_dictionary['output_n_lowest_minima'], self.isRestart, self.restart_path, self._minima_path
                , self.parameter_dictionary['write_graph_output'])
            print('this shoul never be called')
            quit()
        else:
            time.sleep(1)
            print('Worker starting work ', self.mpiRank)
        # initialize database
        with self.database.Database(self.parameter_dictionary["energy_threshold"], self.parameter_dictionary["fingerprint_threshold"]\
                , self.parameter_dictionary["output_n_lowest_minima"], self.isRestart, self.restart_path, self._minima_path\
                , self.parameter_dictionary["write_graph_output"])\
                as self.data:
            # Start up minimahopping 
            atoms = deepcopy(self.initial_configuration)
            current_minimum = self._startup(atoms,)  # gets an atoms object and a minimum object is returned.
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
                    escaped_minimum, epot_max, md_trajectory, opt_trajectory = self._escape(current_minimum) 
                    print("  ---------------------------------------------------------------")
                    print("  New minimum found!")

                    self.data.addElementandConnectGraph(current_minimum, escaped_minimum, md_trajectory + opt_trajectory, epot_max)

                    # write output
                    self._hoplog(escaped_minimum)

                    # check if intermediate or escaped minimum should be considered for accepting.
                    if is_first_accept_iteration or not self.parameter_dictionary["use_intermediate_mechanism"]: 
                        # after first iteration, the escaped minimum must be considered for accecpting
                        intermediate_minimum = escaped_minimum.__deepcopy__()
                        intermediate_minimum_is_escaped_minimum = True
                        is_first_accept_iteration = False
                    else: # After several iterations, the check if the escaped or the intermediate minimum has a lower energy
                        intermediate_minimum_is_escaped_minimum = False
                        if escaped_minimum.e_pot < intermediate_minimum.e_pot:
                            intermediate_minimum = escaped_minimum.__deepcopy__()
                            intermediate_minimum_is_escaped_minimum = True

                    # accept minimum if necessary
                    is_accepted =  self._accept_reject_step(current_minimum, intermediate_minimum)

                    # write log messages
                    if is_accepted:
                        current_minimum = intermediate_minimum.__deepcopy__()
                        if intermediate_minimum_is_escaped_minimum:
                            status = "Accepted minimum after escaping"
                            self._history_log(escaped_minimum, status, escaped_minimum.n_visit)
                        else:
                            status = 'Accepted intermediate minimum'
                            self._history_log(intermediate_minimum, status, escaped_minimum.n_visit)
                    else:
                        if intermediate_minimum_is_escaped_minimum:
                            status = "Rejected -> new intermediate minimum"
                        else:
                            status = "Rejected"
                        self._history_log(escaped_minimum, status, escaped_minimum.n_visit)

                    log_msg = "  New minimum has been found {:d} time(s)".format(escaped_minimum.n_visit)
                    print(log_msg)

                    # adjust the temperature according to the number of visits
                    self._adj_temperature(escaped_minimum.n_visit)
                    counter += 1
                    self._write_restart(escaped_minimum, intermediate_minimum, is_accepted)

                print("DONE")
                print("=================================================================")

                if self.parameter_dictionary["run_time"] != "infinite":
                    if time.time() - self._time_in > self._run_time_sec:
                        msg = 'Simulation stopped because the given time is over\n'
                        msg += 'Run terminated after {:d} steps'.format(counter)
                        print(msg)
                        print("=================================================================")
                        return
                
        self.print_elapsed_time(totalsteps)
        return


    def _write_restart(self, escaped_minimum: Minimum, intermediate_minimum: Minimum, isAccepted: bool):
        f = open(self.restart_path+"params.json", "w")
        json.dump(self.parameter_dictionary,f)
        f.close()
        if isAccepted:
            intermediate_minimum.write(self.restart_path + "poscur.extxyz", append=False)
            intermediate_minimum.write(self._outpath + "accepted_minima.extxyz", append=True)
        escaped_minimum.write(self._outpath + "all_minima.extxyz", append=True)
        if escaped_minimum.n_visit == 1:
            escaped_minimum.write(self._outpath + "all_minima_no_duplicates.extxyz", append=True)
        

    def _startup(self, atoms):
        print("=================================================================")
        print("MINIMAHOPPING SETUP START")

        if not isinstance(atoms, list):
            atoms = [deepcopy(atoms)]

        for atom in atoms:
            calc = atom.calc
            self.calc = calc
            _positions, _lattice = self._restart_opt(atom,)
            atom.set_positions(_positions)
            atom.set_cell(_lattice)

            struct = Minimum(atom,
                        n_visit=1,
                        s = self._ns_orb,
                        p = self._np_orb, 
                        width_cutoff = self._width_cutoff,
                        maxnatsphere = self._maxnatsphere,
                        epot = atom.get_potential_energy(),
                        T=self._temperature,
                        ediff=self._Ediff,
                        label=0)
            
            self.data.addElement(struct)
        # Convert given time to seconds
        if self.parameter_dictionary["run_time"] != "infinite":
            self._run_time_sec = self._get_sec()
        self._time_in = time.time()

        if not isinstance(atoms, list):
            atoms = [deepcopy(atoms)]
        calc = atoms[0].calc

        # Check if this is a fresh start
        if not self.isRestart:
            print('  New MH run is started')
            for atom in atoms:
                _positions, _lattice = self._restart_opt(atom,)
                atom.set_positions(_positions)
                atom.set_cell(_lattice)
                struct = Minimum(atom,
                            s = self.parameter_dictionary['n_S_orbitals'],
                            p = self.parameter_dictionary["n_P_orbitals"], 
                            width_cutoff = self.parameter_dictionary["width_cutoff"],
                            maxnatsphere = self.parameter_dictionary["max_atoms_in_cutoff_sphere"],
                            epot = atom.get_potential_energy(),
                            T=self.parameter_dictionary['T'],
                            ediff=self.parameter_dictionary["energy_difference_to_accept"],
                            exclude= self.parameter_dictionary["exclude"])
                self.data.addElement(struct)
            # add input structure to database after optimization
            struct_cur = self.data.get_element(0)
            struct_cur.atoms.calc = calc
            self._write_restart(struct_cur, struct_cur, True)
        else:
            print('  Restart MH run')

            # Read current structure
            filename = self.restart_path + "poscur.extxyz"
            atoms = read(filename, parallel= False)
            atoms.calc = calc
            
            struct_cur = Minimum(atoms,
                        s = self.parameter_dictionary['n_S_orbitals'],
                        p = self.parameter_dictionary["n_P_orbitals"], 
                        width_cutoff = self.parameter_dictionary["width_cutoff"],
                        maxnatsphere = self.parameter_dictionary["max_atoms_in_cutoff_sphere"],
                        epot = atoms.get_potential_energy(),
                        T=self.parameter_dictionary['T'],
                        ediff=self.parameter_dictionary["energy_difference_to_accept"],
                        exclude= self.parameter_dictionary["exclude"])
        
            database_index = self.data.get_element_index(struct_cur)
            if database_index == -1:
                print("restart structure not in database, quitting")
                quit()
            struct_cur = self.data.get_element(database_index)
            struct_cur.atoms = atoms

        status = 'Initial'
        self._history_log(struct_cur, status, n_visits=struct_cur.n_visit)
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
        atoms = struct.atoms.copy()
        atoms.calc = self.calc

        _i_steps = 0
        while _escape < self.parameter_dictionary["fingerprint_threshold"] or _escape_energy < self.parameter_dictionary["energy_threshold"]:
            atoms = deepcopy(struct.atoms)
            # if the loop not escaped (no new minimum found) rise temperature
            if _i_steps > 0:
                self._n_same += 1
                status = 'Same'
                self._history_log(proposed_structure, status)
                self.parameter_dictionary['T'] *= self.parameter_dictionary['beta_increase']
                log_msg = "    Same minimum found with fpd {:1.2e} {:d} time(s). Increase temperature to {:1.5f}".format(_escape, self._n_same, self.parameter_dictionary['T'])
                print(log_msg)

            # set the temperature according to Boltzmann distribution
            MaxwellBoltzmannDistribution(atoms, temperature_K=self.parameter_dictionary['T'], communicator='serial')

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

                md = MD(atoms=atoms, calculator=self.calc, outpath=self._outpath, cell_atoms=self._cell_atoms, dt=self.parameter_dictionary['dt'], n_max=self.parameter_dictionary["mdmin"], verbose=self.parameter_dictionary["verbose_output"])
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

            proposed_structure = Minimum(atoms,
                        s = self.parameter_dictionary['n_S_orbitals'],
                        p = self.parameter_dictionary["n_P_orbitals"], 
                        width_cutoff = self.parameter_dictionary["width_cutoff"],
                        maxnatsphere = self.parameter_dictionary["max_atoms_in_cutoff_sphere"],
                        epot = atoms.get_potential_energy(),
                        T=self.parameter_dictionary['T'],
                        ediff=self.parameter_dictionary["energy_difference_to_accept"],
                        exclude= self.parameter_dictionary["exclude"])

            # check if proposed structure is the same to the initial structure
            _escape_energy = struct.__compareto__(proposed_structure)
            _escape = struct.__equals__(proposed_structure)

            _i_steps += 1
            self._n_min += 1
        

        log_msg = "    New minimum found with fpd {:1.2e} after looping {:d} time(s)".format(_escape, _i_steps)
        print(log_msg)

        return proposed_structure, _epot_max, _md_trajectory, _opt_trajectory


    def _hoplog(self, struct):
        atoms = struct.atoms
        log_msg = "  Epot:  {:1.5f}   E_diff:  {:1.5f}    Temp:   {:1.5f} ".format(atoms.get_potential_energy(),
                                                                                             self.parameter_dictionary["energy_difference_to_accept"],
                                                                                             self.parameter_dictionary['T'])
        print(log_msg)



    def _accept_reject_step(self, struct_cur: Minimum, struct: Minimum):
        _e_pot_cur = struct_cur.e_pot
        _e_pot = struct.e_pot

        _ediff_in = self.parameter_dictionary["energy_difference_to_accept"]

        if _e_pot - _e_pot_cur < self.parameter_dictionary["energy_difference_to_accept"]:
            self.parameter_dictionary["energy_difference_to_accept"] *= self.parameter_dictionary['alpha_accept']
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
        else:
            log_msg = "  Minimum was rejected:  Enew - Ecur = {:1.5f} > {:1.5f} = Ediff".format(ediff_rej,
                                                                                                _ediff_in)
            print(log_msg)

        self.parameter_dictionary["energy_difference_to_accept"] = self.parameter_dictionary["energy_difference_to_accept"]
        return is_accepted


    def _adj_temperature(self, n_visits):
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
        isRestart = file_handling.restart(self._outpath, self.restart_path, self._minima_path, self.isMaster)
        # Check if new start is desired
        if isNewStart:
            isRestart = False
        return isRestart


    def print_elapsed_time(self, totalsteps):
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
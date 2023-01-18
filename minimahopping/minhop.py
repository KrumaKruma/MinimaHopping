import numpy as np
import warnings
from ase.io import read
import ase.atom
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import minimahopping.mh.lattice_operations as lat_opt
import minimahopping.md.soften as softening
import minimahopping.md.md as md
import minimahopping.opt.optim as opt
from minimahopping.mh.minimum import Minimum
from minimahopping.mh.cell_atom import Cell_atom
import time
import json
import minimahopping.mh.file_handling as file_handling 
import minimahopping.MPI_database.mpi_database_master as MPI_server
import os
from mpi4py import MPI
from minimahopping.MPI_database import mpi_messages
from ase import Atoms

"""
MH Software written by Marco Krummenacher (marco.krummenacher@unibas.ch), Moritz Gubler and Jonas Finkler
Parts of the software were originally developped (some in Fortran) from other people:
  -- VCSMD: Martin Sommer-Joergenson
  -- VCS Softening: Hannes Huber
  -- VCS optimizer: Moritz Gubler
  -- OMFP in python: Jonas Finkler
"""

# TODO: test mh with periodic bounrdary conditions and bazant


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
    """A minima hopping object. Must be accessed through a context manager

    Attributes:
        data            database object that handles all local minima.
        mpiRank         Rank of mpi
        mpiSize         size of MPI_COMM_WORLD
        _minima_path    path where database stores all minima
        isMaster        boolean that stores whether process is master or not
        database        library that contains databse implementation
    """

    data = None
    mpiRank = MPI.COMM_WORLD.Get_rank()
    mpiSize = MPI.COMM_WORLD.Get_size()
    _minima_path = 'minima/'
    isMaster = None
    database = None

    def __init__(self, initial_configuration : ase.atom.Atom,
                        T0 = 1000,
                        beta_decrease = 1./1.02,
                        beta_increase = 1.02,
                        Ediff0 = .1,
                        alpha_accept = 1/1.02,
                        alpha_reject = 1.02,
                        n_soft = 20,
                        soften_positions = 1e-2,
                        soften_lattice = 1e-3,
                        soften_biomode = False,
                        n_S_orbitals = 1, 
                        n_P_orbitals = 1, 
                        width_cutoff = 3.5, 
                        maxnatsphere = 100, 
                        exclude = [],
                        dt = 0.01,
                        mdmin = 2,
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
                        use_MPI = False,
                        totalWorkers = None
                        ):
        """Initialize with an ASE atoms object and keyword arguments."""

        self.initial_configuration = initial_configuration

        
        if self.mpiSize == 1 or not use_MPI: # no mpi should be used
            self.isMaster = False
            self.database = importer("minimahopping.mh.database")
            self._outpath = 'output/' 
            self.restart_path = self._outpath + "restart/"
            # self._minima_path = 'minima/'
            if use_MPI:
                print('UseMPI is true but only one rank is present. simulation will be stopped.')
                quit()
        else: # mpi parallelized simulation
            if totalWorkers is None:
                print("""In an mpi simulation, the total number of workers parameter must be set to the correct number""")
                MPI.COMM_WORLD.Abort()
            
            if self.mpiRank == 0:
                self.isMaster = True
                self._outpath = 'output/master/'
                self.restart_path = self._outpath + "restart/"
                # self._minima_path = 'minima/'
                if not os.path.exists("output"):
                    os.mkdir('output')
                if not os.path.exists(self._minima_path):
                    os.mkdir(self._minima_path)
                
            comm = MPI.COMM_WORLD
            comm.Barrier()
            if self.mpiRank > 0:
                self.isMaster = False
                self.database = importer("minimahopping.MPI_database.mpi_database_worker")
                self._outpath = 'output/worker_' + str(self.mpiRank) + '/' 
                self.restart_path = self._outpath + "restart/"
                # self._minima_path = 'minima/worker_' + str(self.mpiRank) + '/'


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
                "alpha_soften_positions" : soften_positions,
                "alpha_soften_lattice" : soften_lattice,
                "soften_biomode" : soften_biomode,
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
                "use_MPI": use_MPI,
                "totalWorkers": totalWorkers
            }

        if self.isMaster:
            f = open(self.restart_path+"params.json", "w")
            json.dump(self.parameter_dictionary,f)
            f.close()

        # Initialization of global counters        
        self._n_min = 1
        self._n_notunique = 0
        self._n_same = 0


    def __enter__(self):
        if self.mpiRank > 0 or not self.parameter_dictionary["use_MPI"]: 
            self.data = self.database.Database(self.parameter_dictionary["energy_threshold"], self.parameter_dictionary["fingerprint_threshold"]\
                    , self.parameter_dictionary["output_n_lowest_minima"], self.isRestart, self.restart_path, self._minima_path\
                    , self.parameter_dictionary["write_graph_output"])
            self.data.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        if self.mpiRank > 0: 
            self.data.__exit__(exc_type, exc_value, exc_traceback)
        if self.parameter_dictionary['use_MPI'] and not self.isMaster: # slave threads must send exit signals to master.
            MPI.COMM_WORLD.send((mpi_messages.clientWorkDone, None), 0)
            print('client set work done signal to server')


    def __call__(self, totalsteps = None):
        counter = 0
        if totalsteps is None:
            print('provide integer for the total number of minima hopping steps')
            quit()

        if self.isMaster:
            print('starting mpi server on rank', self.mpiRank)
            MPI_server.MPI_database_server_loop(self.parameter_dictionary['energy_threshold'], self.parameter_dictionary['fingerprint_threshold']
                , self.parameter_dictionary['output_n_lowest_minima'], self.isRestart, self.restart_path, self._minima_path
                , self.parameter_dictionary['write_graph_output'], totalWorkers=self.parameter_dictionary["totalWorkers"], maxTimeHours=self._get_sec() / 3600)
            print('All clients have left and the server will shut down as well.')
            print("sending mpi_abort to comm_world to make sure that all clients stop working")
            # time.sleep(5)
            # MPI.COMM_WORLD.Abort()
            # quit()
        else:
            # time.sleep(1)
            print('Worker starting work ', self.mpiRank)

        if self.data is None:
            print("The minimahopping class must be accessed through a context manager.")
            quit()

        # Start up minimahopping 
        structure_list, calculator = self._initialize_structures(self.initial_configuration)
        self.calculator = calculator
        current_minimum = self._startup(structure_list)  # gets an atoms object and a minimum object is returned.

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

                n_visit, label, continueSimulation = self.data.addElementandConnectGraph(current_minimum, escaped_minimum, md_trajectory + opt_trajectory, epot_max)

                # write output
                self._hoplog(escaped_minimum)

                # check if intermediate or escaped minimum should be considered for accepting.
                if is_first_accept_iteration or not self.parameter_dictionary["use_intermediate_mechanism"]: 
                    # after first iteration, the escaped minimum must be considered for accecpting
                    intermediate_minimum = escaped_minimum.__copy__()
                    intermediate_minimum_is_escaped_minimum = True
                    is_first_accept_iteration = False
                else: # After several iterations, the check if the escaped or the intermediate minimum has a lower energy
                    intermediate_minimum_is_escaped_minimum = False
                    if escaped_minimum.e_pot < intermediate_minimum.e_pot:
                        intermediate_minimum = escaped_minimum.__copy__()
                        intermediate_minimum_is_escaped_minimum = True

                # accept minimum if necessary
                is_accepted =  self._accept_reject_step(current_minimum, intermediate_minimum)

                # write log messages
                if is_accepted:
                    current_minimum = intermediate_minimum.__copy__()
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
                if not continueSimulation:
                    print("Client got shut down signal from server and is shutting down.")
                    return

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
        # if self.parameter_dictionary['use_MPI']:
        #     print("An MPI worker is not allowed to leave the above loop because the server might freeze. Sending MPI_abort to comm_world")
        #     from mpi4py import MPI
        #     MPI.COMM_WORLD.Abort()
        return


    def _initialize_structures(self, atoms_in):
        """ 
        initializes the sturcture
        Input: 
            atoms_in: Initial structure (ASE atoms object with and attachted calculator)
        Return:
            atoms_list: List of ASE atoms objects (without calculator)
            calculator: ASE calculator
        """

        atoms_list = []
        if not isinstance(atoms_in, list):
            atoms_out, calculator = self._split_atoms_and_calculator(atoms_in)
            atoms_list.append(atoms_out.copy())
        else:
            for atom in atoms_in:
                atoms_out, calculator = self._split_atoms_and_calculator(atom)
                atoms_list.append(atoms_out.copy())
        return atoms_list, calculator

    
    def _split_atoms_and_calculator(self, atoms_in):
        """
        Extract calculator and only nessecairy information from atoms object
        Input:
            atoms_in: ASE atoms object with an attached calculator
        Return:
            atoms_out: ASE atoms object only with nessecairy information
            calulator: ASE calculator
        """
        calculator = atoms_in.calc
        positions = atoms_in.get_positions()
        cell = atoms_in.get_cell()
        pbc = atoms_in.pbc
        elements = atoms_in.get_chemical_symbols()
        atoms_out = Atoms(symbols=elements, positions=positions, pbc=pbc, cell=cell)
        return atoms_out, calculator


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
        """
        Startup of the minimahopping algorithm
        Input:
            atoms: list of ASE atoms objects
        Return:
            struct_cur: current structure (ASE atoms object)
        """

        # Input is a list of ASE atoms objects
        print("=================================================================")
        print("MINIMAHOPPING SETUP START")

        # Convert given time to seconds
        if self.parameter_dictionary["run_time"] != "infinite":
            self._run_time_sec = self._get_sec()
        self._time_in = time.time()

        # Check if this is a fresh start
        if not self.isRestart:
            print('  New MH run is started')
            for atom in atoms:
                print(len(atom))
                _positions, _lattice = self._restart_opt(atom,)
                atom.set_positions(_positions)
                atom.set_cell(_lattice)
                atom.calc = self.calculator
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
            self._write_restart(struct_cur, struct_cur, True)
        else:
            print('  Restart MH run')

            # Read current structure
            filename = self.restart_path + "poscur.extxyz"
            atoms = read(filename, parallel= False)
            
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
        """
        Optimization wrapper for the startup
        """
        positions, lattice, noise, trajectory, number_of_steps = opt.optimization(atoms=atoms, 
                                                                        calculator=self.calculator, 
                                                                        max_force_threshold=self.parameter_dictionary["fmax"], 
                                                                        outpath=self._outpath, 
                                                                        verbose=self.parameter_dictionary["verbose_output"])
        return positions, lattice


    def _get_sec(self,):
        """
        Get seconds from time.
        """
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
        _i_steps = 0

        is_escape = True

        while is_escape:
            atoms = struct.atoms.copy()
            atoms.calc = self.calculator
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

            # check that periodic boundaries are the same in all directions (no mixed boundary conditions)
            _pbc = list(set(atoms.pbc))
            assert len(_pbc) == 1, "mixed boundary conditions"

            # if periodic boundary conditions create cell atom object
            if True in _pbc:
                print("    VARIABLE CELL SHAPE SOFTENING, MD AND OPTIMIZATION ARE PERFORMED")
                # calculate mass for cell atoms
                # Formula if for the MD real masses are used
                # mass = .75 * np.sum(atoms.get_masses()) / 10. # Formula if for the MD real masses are used
                # Formula if mass 1 is used in the MD
                mass = .75 * np.sum(len(atoms)) / 10.
                # set position and mass of cell atoms
                cell_atoms = Cell_atom(mass=_mass, positions=atoms.get_cell())
                # set velocities of the cell atoms
                cell_atoms.set_velocities_boltzmann(temperature=self.parameter_dictionary['T'])
            else:
                # IMPORTANT: cell atoms has to be None for soften/md/geopt if no pbc
                cell_atoms = None

            # softening of the velocities
            velocities, cell_velocities = softening.soften(atoms=atoms, 
                            calculator=self.calculator, 
                            nsoft=self.parameter_dictionary['n_softening_steps'],
                            alpha_pos = self.parameter_dictionary['alpha_soften_positions'], 
                            cell_atoms = cell_atoms,
                            alpha_lat =  self.parameter_dictionary['alpha_soften_lattice'])
            
            # set softened velocities
            atoms.set_velocities(velocities)

            # set cell velocities if pbc
            if True in _pbc:
                cell_atoms.velocities = cell_velocities
   
            # Perfom MD run
            print("    MD Start")
            positions, lattice, self.parameter_dictionary['dt'], _md_trajectory, _epot_max, number_of_md_steps = md.md(atoms = atoms, 
                                                                                                        calculator = self.calculator,
                                                                                                        outpath = self._outpath, 
                                                                                                        cell_atoms = cell_atoms,
                                                                                                        dt = self.parameter_dictionary['dt'], 
                                                                                                        n_max = self.parameter_dictionary["mdmin"],
                                                                                                        verbose = self.parameter_dictionary["verbose_output"])

            log_msg = "    MD finished after {:d} steps visiting {:d} maxima. New dt is {:1.5f}".format(number_of_md_steps, self.parameter_dictionary["mdmin"], self.parameter_dictionary['dt'])

            print(log_msg)
            # Set new positions after the MD
            atoms.set_positions(positions)
            # If pbc set new lattice and reshape cell
            if True in _pbc:
                atoms.set_cell(lattice)
                atoms = lat_opt.reshape_cell2(atoms, 6)

            print("    OPT start")
            positions, lattice, self._noise, _opt_trajectory, number_of_opt_steps = opt.optimization(atoms=atoms, 
                                                                    calculator=self.calculator, 
                                                                    max_force_threshold=self.parameter_dictionary["fmax"], 
                                                                    outpath=self._outpath, 
                                                                    verbose=self.parameter_dictionary["verbose_output"])
            # Set optimized positions
            atoms.set_positions(positions)
            # If Pbc set optimized lattice 
            if True in _pbc:
                atoms.set_cell(lattice)

            log_msg = "    OPT finished after {:d} steps.".format(number_of_opt_steps)
            print(log_msg)
            # TODO: testing cluster and VCS

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

            if  _escape > self.parameter_dictionary["fingerprint_threshold"]:
                is_escape = False
            elif _escape_energy > self.parameter_dictionary["energy_threshold"]:
                is_escape = False

        log_msg = "    New minimum found with fpd {:1.2e} after looping {:d} time(s)".format(_escape, _i_steps)
        print(log_msg)

        return proposed_structure, _epot_max, _md_trajectory, _opt_trajectory


    def _hoplog(self, struct):
        """
        Print log information of the minimahopping
        """
        atoms = struct.atoms
        atoms.calc = self.calculator
        log_msg = "  Epot:  {:1.5f}   E_diff:  {:1.5f}    Temp:   {:1.5f} ".format(atoms.get_potential_energy(),
                                                                                             self.parameter_dictionary["energy_difference_to_accept"],
                                                                                             self.parameter_dictionary['T'])
        print(log_msg)



    def _accept_reject_step(self, struct_cur: Minimum, struct: Minimum):
        """
        Accept/Reject step in the algorithm and adjustment of E_diff. 
        Input:
            sturct_cur: Minimum object of the old minimum
            struct: Minimum object of the newly found minimum
        Return:
            is_accepted: bool, True if struct is accepted
        """
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
        """
        Adjust the temperature depending if minimum was found previously
        Input:
            n_visits: number of times the minimum has been visited
        """
        if n_visits > 1:
            self._n_notunique += 1
            if self.parameter_dictionary["enhanced_feedback"]:
                self.parameter_dictionary['T'] = self.parameter_dictionary['T'] * self.parameter_dictionary['beta_increase'] * (1. + 1. * np.log(float(n_visits)))
            else:
                self.parameter_dictionary['T'] = self.parameter_dictionary['T'] * self.parameter_dictionary['beta_increase']
        else:
            self.parameter_dictionary['T'] = self.parameter_dictionary['T'] * self.parameter_dictionary['beta_decrease']


    def _history_log(self, struct, status, n_visits = 0):
        """
        Writing log message in the history file
        Input:
            struct: Minimum object
            status: str (Accept, Reject, Intermediate, etc.)
            n_visits: number of time the minimum has been visited
        """
        atoms = struct.atoms
        atoms.calc = self.calculator
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
        """
        Function that checks if the energy_threshold is below the noise
        """
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
        """
        Function that prints the elapsed time in the form DD-HH-MM-SS
        """
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
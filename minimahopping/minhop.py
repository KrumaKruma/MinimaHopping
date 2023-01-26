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
import minimahopping.mh.parameters

import minimahopping.mh.database
import minimahopping.MPI_database.mpi_database_worker

"""
MH Software written by Marco Krummenacher (marco.krummenacher@unibas.ch), Moritz Gubler and Jonas Finkler
Parts of the software were originally developped (some in Fortran) from other people:
  -- VCSMD: Martin Sommer-Joergenson
  -- VCS Softening: Hannes Huber
  -- VCS optimizer: Moritz Gubler
  -- OMFP in python: Jonas Finkler
"""

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
    isWorker = None


    def __init__(self, initial_configuration : ase.atom.Atom, **kwargs):
        """Initialize with an ASE atoms object and keyword arguments."""

        self.initial_configuration = initial_configuration

        initalParameters = minimahopping.mh.parameters.minimaHoppingParameters(**kwargs)

        self.createPathsAndSetMPIVariables(initalParameters.use_MPI, initalParameters.totalWorkers)

        self.isRestart = self.checkIfRestart(initalParameters.new_start)

        if self.isRestart:
            filename = self.restart_path + "params.json"
            f = open(filename)
            parameter_dictionary = json.load(f)
            f.close()
            self.parameters = minimahopping.mh.parameters.minimaHoppingParameters(**parameter_dictionary)
        else:
            self.parameters = initalParameters
            self.parameters._dt = self.parameters.dt0
            self.parameters._T = self.parameters.T0
            self.parameters._eDiff = self.parameters.Ediff0

        if self.isMaster:
            f = open(self.restart_path+"params.json", "w")
            json.dump(self.parameters.to_dict(),f)
            f.close()

        # Initialization of global counters        
        self._n_min = 1
        self._n_notunique = 0
        self._n_same = 0


    def __enter__(self):

        # master does not need setting up here
        if self.isMaster:
            return self

        if not self.parameters.use_MPI:
            self.data = minimahopping.mh.database.Database(self.parameters.energy_threshold, self.parameters.fingerprint_threshold\
                    , self.parameters.output_n_lowest_minima, self.isRestart, self.restart_path, self._minima_path\
                    , self.parameters.write_graph_output)
        elif self.isWorker:
            self.data = minimahopping.MPI_database.mpi_database_worker.Database(self.parameters.energy_threshold, self.parameters.fingerprint_threshold\
                    , self.parameters.output_n_lowest_minima, self.isRestart, self.restart_path, self._minima_path\
                    , self.parameters.write_graph_output)

        self.data.__enter__()

        if self.parameters.collect_md_data:
            self.collect_md_file = open(self._outpath + "MD_collection.extxyz", "a")
        else:
            self.collect_md_file = None
        # open history file
        self.history_file = open(self._outpath + 'history.dat', 'a')
        return self


    def __exit__(self, exc_type, exc_value, exc_traceback):
        # master does not need to close anything here
        if self.isMaster:
            return
        self.data.__exit__(exc_type, exc_value, exc_traceback)
        if self.isWorker: # slave threads must send exit signals to master.
            MPI.COMM_WORLD.send((mpi_messages.clientWorkDone, None), 0)
            print('client set work done signal to server', flush=True)
        # close MD collection file if MD is collected for ML
        if self.parameters.collect_md_data:
            self.collect_md_file.close()
        # close histroy file
        self.history_file.close()


    def __call__(self, totalsteps = None):
        counter = 0
        if totalsteps is None:
            print('provide integer for the total number of minima hopping steps', flush=True)
            quit()

        if self.isMaster:
            print('starting mpi server on rank', self.mpiRank, flush=True)
            MPI_server.MPI_database_server_loop(self.parameters.energy_threshold, self.parameters.fingerprint_threshold
                , self.parameters.output_n_lowest_minima, self.isRestart, self.restart_path, self._minima_path
                , self.parameters.write_graph_output, totalWorkers=self.parameters.totalWorkers, maxTimeHours=self._get_sec() / 3600)
            print('All clients have left and the server will shut down as well.', flush=True)
            # print("sending mpi_abort to comm_world to make sure that all clients stop working", flush=True)
            # time.sleep(5)
            # MPI.COMM_WORLD.Abort()
            # quit()
            return
        else:
            # time.sleep(1)
            print('Worker starting work ', self.mpiRank, flush=True)

        if self.data is None:
            print("The minimahopping class must be accessed through a context manager.", flush=True)
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
            print(msg, flush=True)
            while not is_accepted:
                print("  Start escape loop", flush=True)
                print("  ---------------------------------------------------------------", flush=True)
                escaped_minimum, epot_max, md_trajectory, opt_trajectory = self._escape(current_minimum) 
                print("  ---------------------------------------------------------------",flush=True)
                print("  New minimum found!", flush=True)

                n_visit, label, continueSimulation = self.data.addElementandConnectGraph(current_minimum, escaped_minimum, md_trajectory + opt_trajectory, epot_max)

                # write output
                self._hoplog(escaped_minimum)

                # check if intermediate or escaped minimum should be considered for accepting.
                if is_first_accept_iteration or not self.parameters.use_intermediate_mechanism: 
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
                print(log_msg, flush=True)

                # adjust the temperature according to the number of visits
                self._adj_temperature(escaped_minimum.n_visit)
                counter += 1
                self._write_restart(escaped_minimum, intermediate_minimum, is_accepted)
                if not continueSimulation:
                    print("Client got shut down signal from server and is shutting down.", flush=True)
                    return

            print("DONE", flush=True)
            print("=================================================================", flush=True)

            if self.parameters.run_time != "infinite":
                if time.time() - self._time_in > self._run_time_sec:
                    msg = 'Simulation stopped because the given time is over\n'
                    msg += 'Run terminated after {:d} steps'.format(counter)
                    print(msg, flush=True)
                    print("=================================================================", flush=True)
                    return
        
        self.print_elapsed_time(totalsteps)
        # if self.parameters.use_MPI:
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
        json.dump(self.parameters.to_dict(),f)
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
        print("=================================================================", flush=True)
        print("MINIMAHOPPING SETUP START", flush=True)

        # Convert given time to seconds
        if self.parameters.run_time != "infinite":
            self._run_time_sec = self._get_sec()
        self._time_in = time.time()

        # Check if this is a fresh start
        if not self.isRestart:
            print('  New MH run is started', flush=True)
            for atom in atoms:
                _positions, _lattice = self._restart_opt(atom,)
                atom.set_positions(_positions)
                atom.set_cell(_lattice)
                atom.calc = self.calculator
                struct = Minimum(atom,
                            s = self.parameters.n_S_orbitals,
                            p = self.parameters.n_P_orbitals, 
                            width_cutoff = self.parameters.width_cutoff,
                            maxnatsphere = self.parameters.maxnatsphere,
                            epot = atom.get_potential_energy(),
                            T=self.parameters._T,
                            ediff=self.parameters._eDiff,
                            exclude= self.parameters.exclude)
                self.data.addElement(struct)
            # add input structure to database after optimization
            struct_cur = self.data.get_element(0)
            self._write_restart(struct_cur, struct_cur, True)
        else:
            print('  Restart MH run', flush=True)

            # Read current structure
            filename = self.restart_path + "poscur.extxyz"
            atoms = read(filename, parallel= False)
            
            struct_cur = Minimum(atoms,
                        s = self.parameters.n_S_orbitals,
                        p = self.parameters.n_P_orbitals, 
                        width_cutoff = self.parameters.width_cutoff,
                        maxnatsphere = self.parameters.maxnatsphere,
                        epot = atoms.get_potential_energy(),
                        T=self.parameters._T,
                        ediff=self.parameters._eDiff,
                        exclude= self.parameters.exclude)
        
            database_index = self.data.get_element_index(struct_cur)
            if database_index == -1:
                print("restart structure not in database, quitting", flush=True)
                quit()
            struct_cur = self.data.get_element(database_index)
            struct_cur.atoms = atoms

        status = 'Initial'
        self._history_log(struct_cur, status, n_visits=struct_cur.n_visit)

        return struct_cur


    def _restart_opt(self, atoms,):
        """
        Optimization wrapper for the startup
        """
        positions, lattice, noise, trajectory, number_of_steps = opt.optimization(atoms=atoms, 
                                                                        calculator=self.calculator, 
                                                                        max_force_threshold=self.parameters.fmax, 
                                                                        outpath=self._outpath, 
                                                                        verbose=self.parameters.verbose_output)
        return positions, lattice


    def _get_sec(self,):
        """
        Get seconds from time.
        """
        if self.parameters.run_time is 'infinite':
            return np.inf
        nd, d = self.parameters.run_time.split('-')
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
            try:
                atoms.calc.recalculateBasis(atoms)
            except:
                pass
            # if the loop not escaped (no new minimum found) rise temperature
            if _i_steps > 0:
                self._n_same += 1
                status = 'Same'
                self._history_log(proposed_structure, status)
                self.parameters._T *= self.parameters.beta_increase
                log_msg = "    Same minimum found with fpd {:1.2e} {:d} time(s). Increase temperature to {:1.5f}".format(_escape, self._n_same, self.parameters._T)
                print(log_msg, flush=True)

            # set the temperature according to Boltzmann distribution
            MaxwellBoltzmannDistribution(atoms, temperature_K=self.parameters._T, communicator='serial')

            # check that periodic boundaries are the same in all directions (no mixed boundary conditions)
            _pbc = list(set(atoms.pbc))
            assert len(_pbc) == 1, "mixed boundary conditions"

            # if periodic boundary conditions create cell atom object
            if True in _pbc:
                print("    VARIABLE CELL SHAPE SOFTENING, MD AND OPTIMIZATION ARE PERFORMED", flush=True)
                # calculate mass for cell atoms
                # Formula if for the MD real masses are used
                # mass = .75 * np.sum(atoms.get_masses()) / 10. # Formula if for the MD real masses are used
                # Formula if mass 1 is used in the MD
                mass = .75 * np.sum(len(atoms)) / 10.
                # set position and mass of cell atoms
                cell_atoms = Cell_atom(mass=mass, positions=atoms.get_cell())
                # set velocities of the cell atoms
                cell_atoms.set_velocities_boltzmann(temperature=self.parameters._T)
            else:
                # IMPORTANT: cell atoms has to be None for soften/md/geopt if no pbc
                cell_atoms = None

            # softening of the velocities
            velocities, cell_velocities = softening.soften(atoms=atoms, 
                            calculator=self.calculator, 
                            nsoft = self.parameters.n_soft,
                            alpha_pos = self.parameters.soften_positions, 
                            cell_atoms = cell_atoms,
                            alpha_lat =  self.parameters.soften_lattice)
            
            # set softened velocities
            atoms.set_velocities(velocities)

            # set cell velocities if pbc
            if True in _pbc:
                cell_atoms.velocities = cell_velocities
   
            # Perfom MD run
            print("    MD Start", flush=True)
            positions, lattice, self.parameters._dt, _md_trajectory, _epot_max, number_of_md_steps = md.md(atoms = atoms, 
                                                                                                        calculator = self.calculator,
                                                                                                        outpath = self._outpath, 
                                                                                                        cell_atoms = cell_atoms,
                                                                                                        dt = self.parameters._dt, 
                                                                                                        n_max = self.parameters.mdmin,
                                                                                                        verbose = self.parameters.verbose_output,
                                                                                                        collect_md_file = self.collect_md_file)

            log_msg = "    MD finished after {:d} steps visiting {:d} maxima. New dt is {:1.5f}".format(number_of_md_steps, self.parameters.mdmin, self.parameters._dt)

            print(log_msg, flush=True)
            # Set new positions after the MD
            atoms.set_positions(positions)
            # If pbc set new lattice and reshape cell
            if True in _pbc:
                atoms.set_cell(lattice)
                lat_opt.reshape_cell2(atoms, 6)

            try:
                atoms.calc.recalculateBasis(atoms)
            except:
                pass

            print("    OPT start", flush=True)
            positions, lattice, self._noise, _opt_trajectory, number_of_opt_steps = opt.optimization(atoms=atoms, 
                                                                    calculator=self.calculator, 
                                                                    max_force_threshold=self.parameters.fmax, 
                                                                    outpath=self._outpath, 
                                                                    verbose=self.parameters.verbose_output)
            # Set optimized positions
            atoms.set_positions(positions)
            # If Pbc set optimized lattice 
            if True in _pbc:
                atoms.set_cell(lattice)

            log_msg = "    OPT finished after {:d} steps.".format(number_of_opt_steps)
            print(log_msg, flush=True)

            # check if the energy threshold is below the optimization noise
            self._check_energy_threshold()

            proposed_structure = Minimum(atoms,
                        s = self.parameters.n_S_orbitals,
                        p = self.parameters.n_P_orbitals, 
                        width_cutoff = self.parameters.width_cutoff,
                        maxnatsphere = self.parameters.maxnatsphere,
                        epot = atoms.get_potential_energy(),
                        T=self.parameters._T,
                        ediff=self.parameters._eDiff,
                        exclude= self.parameters.exclude)

            # check if proposed structure is the same to the initial structure
            _escape_energy = struct.__compareto__(proposed_structure)
            _escape = struct.__equals__(proposed_structure)

            _i_steps += 1
            self._n_min += 1

            if  _escape > self.parameters.fingerprint_threshold:
                is_escape = False
            elif _escape_energy > self.parameters._eDiff:
                is_escape = False

        log_msg = "    New minimum found with fpd {:1.2e} after looping {:d} time(s)".format(_escape, _i_steps)
        print(log_msg, flush=True)

        return proposed_structure, _epot_max, _md_trajectory, _opt_trajectory


    def _hoplog(self, struct):
        """
        Print log information of the minimahopping
        """
        atoms = struct.atoms
        atoms.calc = self.calculator
        log_msg = "  Epot:  {:1.5f}   E_diff:  {:1.5f}    Temp:   {:1.5f} ".format(atoms.get_potential_energy(),
                                                                                             self.parameters._eDiff,
                                                                                             self.parameters._T)
        print(log_msg, flush=True)


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

        _ediff_in = self.parameters._eDiff

        if _e_pot - _e_pot_cur < self.parameters._eDiff:
            self.parameters._eDiff *= self.parameters.alpha_accept
            is_accepted = True
            ediff_acc = _e_pot - _e_pot_cur
        else:
            self.parameters._eDiff *= self.parameters.alpha_reject
            is_accepted = False
            ediff_rej = _e_pot - _e_pot_cur

        if is_accepted:
            log_msg = "  Minimum was accepted:  Enew - Ecur = {:1.5f} < {:1.5f} = Ediff".format(ediff_acc,
                                                                                                _ediff_in)
            print(log_msg, flush=True)
        else:
            log_msg = "  Minimum was rejected:  Enew - Ecur = {:1.5f} > {:1.5f} = Ediff".format(ediff_rej,
                                                                                                _ediff_in)
            print(log_msg, flush=True)

        return is_accepted


    def _adj_temperature(self, n_visits):
        """
        Adjust the temperature depending if minimum was found previously
        Input:
            n_visits: number of times the minimum has been visited
        """
        if n_visits > 1:
            self._n_notunique += 1
            if self.parameters.enhanced_feedback:
                self.parameters._T = self.parameters._T * self.parameters.beta_increase * (1. + 1. * np.log(float(n_visits)))
            else:
                self.parameters._T = self.parameters._T * self.parameters.beta_increase
        else:
            self.parameters._T = self.parameters._T * self.parameters.beta_decrease


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
                                                                        self.parameters._T,
                                                                        self.parameters._eDiff,
                                                                        _same_frac,
                                                                        _notunique_frac,
                                                                        _unique_frac,
                                                                        status)

        # history_file = open(self._outpath + 'history.dat', 'a')
        self.history_file.write(history_msg)
        self.history_file.flush()
       

    def _check_energy_threshold(self):
        """
        Function that checks if the energy_threshold is below the noise
        """
        if self.parameters._eDiff < self._noise:
            _warning_msg = 'Energy threshold is below the noise level'
            warnings.warn(_warning_msg, UserWarning)


    def checkIfRestart(self, isNewStart):
        isRestart = file_handling.restart(self._outpath, self.restart_path, self._minima_path, self.isMaster)
        # Check if new start is desired
        if isNewStart:
            isRestart = False
        return isRestart


    def createPathsAndSetMPIVariables(self, use_MPI, totalWorkers):
        if self.mpiSize == 1 or not use_MPI: # no mpi should be used
            self.isMaster = False
            self.isWorker = False
            self._outpath = 'output/' 
            self.restart_path = self._outpath + "restart/"
            # self._minima_path = 'minima/'
            if use_MPI:
                print('UseMPI is true but only one rank is present. simulation will be stopped.', flush=True)
                quit()
        else: # mpi parallelized simulation
            if totalWorkers == 1:
                print("""In an mpi simulation, the total number of workers parameter must be set to the correct number""", flush=True)
                MPI.COMM_WORLD.Abort()
            
            if self.mpiRank == 0:
                self.isMaster = True
                self.isWorker = False
                self._outpath = 'output/master/'
                self.restart_path = self._outpath + "restart/"
                # self._minima_path = 'minima/'
                if not os.path.exists("output"):
                    os.mkdir('output')
                if not os.path.exists(self._minima_path):
                    os.mkdir(self._minima_path)
                
            if not os.path.exists('output'):
                time.sleep(1)
                if not os.path.exists('output'):
                    time.sleep(4)
                    if not os.path.exists('output'):
                        print("Failed to create an output directory. Aborting")
                        quit()

            if self.mpiRank > 0:
                self.isMaster = False
                self.isWorker = True
                self._outpath = 'output/worker_' + str(self.mpiRank) + '/' 
                self.restart_path = self._outpath + "restart/"
                # self._minima_path = 'minima/worker_' + str(self.mpiRank) + '/'


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
        print(msg, flush=True)
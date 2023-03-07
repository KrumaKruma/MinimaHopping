import numpy as np
import warnings
from ase.io import read
import ase.atom
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import minimahopping.mh.lattice_operations as lattice_operations
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
import signal
import sys

import logging

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
    logger = None


    def __init__(self, initial_configuration : ase.atom.Atom, **kwargs):
        """Initialize with an ASE atoms object and keyword arguments."""

        self.initial_configuration = initial_configuration

        initalParameters = minimahopping.mh.parameters.minimaHoppingParameters(**kwargs)

        self.createPathsAndSetMPIVariables(initalParameters.use_MPI, initalParameters.logLevel)

        self.isRestart = self.checkIfRestart(initalParameters.new_start)

        if self.isRestart:
            filename = self.restart_path + "params.json"
            f = open(filename)
            parameter_dictionary = json.load(f)
            f.close()
            compinedDictionary = {**parameter_dictionary, **kwargs}

            for i in initalParameters.getFixedParameterList():
                if compinedDictionary[i] != parameter_dictionary[i]:
                    logging.error("Parameter %s was changed between restart. Changing %s between restarts is not safe. Therefore, the previous value will be used"%(i, i))
                    compinedDictionary[i] = parameter_dictionary[i]

            self.parameters = minimahopping.mh.parameters.minimaHoppingParameters(**compinedDictionary)
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

        if self.parameters.compare_energies:
            self.parameters.calculate_fingerprint = False
        else:
            self.parameters.calculate_fingerprint = True


    def __enter__(self):

        # set up sigtermcatcher
        signal.signal(signal.SIGTERM, self.sigTermCatcher)

        # master does not need setting up here
        if self.isMaster:
            return self

        if not self.parameters.use_MPI:
            self.data = minimahopping.mh.database.Database(self.parameters.energy_threshold, self.parameters.fingerprint_threshold\
                    , self.parameters.output_n_lowest_minima, self.isRestart, self.restart_path, self._minima_path\
                    , self.parameters.write_graph_output, self.parameters.compare_energies)
        elif self.isWorker:
            self.data = minimahopping.MPI_database.mpi_database_worker.Database(self.parameters.energy_threshold, self.parameters.fingerprint_threshold\
                    , self.parameters.output_n_lowest_minima, self.isRestart, self.restart_path, self._minima_path\
                    , self.parameters.write_graph_output, self.parameters.compare_energies)

        self.data.__enter__()

        if self.parameters.collect_md_data:
            self.collect_md_file = open(self._outpath + "MD_collection.extxyz", "a")
        else:
            self.collect_md_file = None
        # open history file
        self.history_file = open(self._outpath + 'history.dat', 'a')
        self.history_file.write("Energy[eV], number of visits, label, temperature [K], delta acceptance, ratio of accepted minima, ratio of rejected minima, ratio of unsucessful escapes (Same), status")
        return self


    def __exit__(self, exc_type, exc_value, exc_traceback):
        # master does not need to close anything here
        if self.isMaster:
            return
        self.data.__exit__(exc_type, exc_value, exc_traceback)
        if self.isWorker: # slave threads must send exit signals to master.
            MPI.COMM_WORLD.send((mpi_messages.clientWorkDone, None), 0)
            logging.info('client set work done signal to server')
        # close MD collection file if MD is collected for ML
        if self.parameters.collect_md_data:
            self.collect_md_file.close()
        # close histroy file
        self.history_file.close()


    def __call__(self, totalsteps = np.inf):
        counter = 0

        if self.isMaster:
            logging.info('starting mpi server on rank %i'%self.mpiRank)
            MPI_server.MPI_database_server_loop(self.parameters.energy_threshold, self.parameters.fingerprint_threshold
                , self.parameters.output_n_lowest_minima, self.isRestart, self.restart_path, self._minima_path
                , self.parameters.write_graph_output, maxTimeHours=self._get_sec() / 3600, compare_energies=self.parameters.compare_energies)
            logging.info('All clients have left and the server will shut down as well.')
            return
        else:
            # time.sleep(1)
            logging.info('Worker %i starting work '%self.mpiRank)

        if self.data is None:
            logging.info("The minimahopping class must be accessed through a context manager.")
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

            logging.info("START HOPPING STEP NR.  {:d}".format(counter))
            while not is_accepted:
                logging.info("  Start escape loop")
                logging.info("  ---------------------------------------------------------------")
                escaped_minimum, epot_max, md_trajectory, opt_trajectory = self._escape(current_minimum) 
                logging.info("  ---------------------------------------------------------------")
                logging.info("  Succesfully escaped minimum")
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
                    self.parameters._n_accepted += 1
                    current_minimum = intermediate_minimum.__copy__()
                    if intermediate_minimum_is_escaped_minimum:
                        status = "Accepted minimum after escaping"
                        self._history_log(escaped_minimum, status)
                    else:
                        status = 'Accepted intermediate minimum'
                        self._history_log(intermediate_minimum, status)
                else:
                    self.parameters._n_rejected += 1
                    if intermediate_minimum_is_escaped_minimum and self.parameters.use_intermediate_mechanism:
                        status = "Rejected -> new intermediate minimum"
                    else:
                        status = "Rejected"
                    self._history_log(escaped_minimum, status)

                logging.info("  New minimum has been found {:d} time(s)".format(escaped_minimum.n_visit))

                # adjust the temperature according to the number of visits
                self._adj_temperature(escaped_minimum.n_visit)
                counter += 1
                self._write_restart(escaped_minimum, intermediate_minimum, is_accepted)
                if not continueSimulation:
                    logging.info("Client got shut down signal from server and is shutting down.")
                    return

            logging.info("DONE")
            logging.info("=================================================================")

            if self.parameters.run_time != "infinite":
                if time.time() - self._time_in > self._run_time_sec:
                    msg = 'Simulation stopped because the given time is over\n'
                    msg += 'Run terminated after {:d} steps'.format(counter)
                    logging.info(msg)
                    logging.info("=================================================================")
                    return
        
        self.print_elapsed_time(totalsteps)
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
            # set the temperature according to Boltzmann distribution
            atoms_out.set_masses(np.ones(len(atoms_out)))
            atoms_list.append(atoms_out.copy())
        else:
            for atom in atoms_in:
                atoms_out, calculator = self._split_atoms_and_calculator(atom)
                atoms_out.set_masses(np.ones(len(atoms_out)))
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
        logging.info("=================================================================")
        logging.info("MINIMAHOPPING SETUP START")

        # Convert given time to seconds
        if self.parameters.run_time != "infinite":
            self._run_time_sec = self._get_sec()
        self._time_in = time.time()

        # Check if this is a fresh start
        if not self.isRestart:
            logging.info('  New MH run is started')
            for atom in atoms:
                try:
                    self.calc.recalculateBasis(atom)
                except:
                    pass
                atom.calc = self.calculator
                _positions, _lattice = self._restart_opt(atom,)
                atom.set_positions(_positions)
                atom.set_cell(_lattice)
                struct = Minimum(atom,
                            s = self.parameters.n_S_orbitals,
                            p = self.parameters.n_P_orbitals, 
                            width_cutoff = self.parameters.width_cutoff,
                            maxnatsphere = self.parameters.maxnatsphere,
                            epot = atom.get_potential_energy(),
                            T=self.parameters._T,
                            ediff=self.parameters._eDiff,
                            exclude= self.parameters.exclude,
                            calculate_fingerprint=self.parameters.calculate_fingerprint)
                logging.debug("Before initial database request in startup")
                n_visit, label, continueSimulation = self.data.addElement(struct)
                self.parameters._n_accepted += 1
                logging.debug("After initial database request in startup")
                if not continueSimulation:
                    logging.info("received shutdown signal after adding first element.")
                    quit()
            # add input structure to database after optimization
            struct_cur = self.data.get_element(0)
            self._write_restart(struct_cur, struct_cur, True)
            try:
                self.calc.recalculateBasis(atom)
            except:
                pass
        else:
            logging.info('  Restart MH run')

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
                        exclude= self.parameters.exclude,
                        calculate_fingerprint=self.parameters.calculate_fingerprint)
        
            database_index = self.data.get_element_index(struct_cur)
            if database_index == -1:
                logging.critical("restart structure not in database, quitting")
                quit()
            struct_cur = self.data.get_element(database_index)
            struct_cur.atoms = atoms

        status = 'Initial'
        self._history_log(struct_cur, status)

        return struct_cur


    def _restart_opt(self, atoms,):
        """
        Optimization wrapper for the startup
        """
        positions, lattice, noise, trajectory, number_of_steps, epot_max = opt.optimization(atoms=atoms, 
                                                                        calculator=self.calculator, 
                                                                        max_force_threshold=self.parameters.fmax, 
                                                                        outpath=self._outpath, 
                                                                        verbose=self.parameters.verbose_output)
        return positions, lattice


    def _get_sec(self,):
        """
        Get seconds from time.
        """
        if self.parameters.run_time == 'infinite':
            return np.inf
        nd, d = self.parameters.run_time.split('-')
        h, m, s = d.split(':')
        return int(nd) * 86400 + int(h) * 3600 + int(m) * 60 + int(s)


    def _escape(self, struct: Minimum):
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
                self._history_log(struct, status)
                self.parameters._T *= self.parameters.beta_increase
                log_msg = "    Same minimum found with fpd {:1.2e} {:d} time(s). Increase temperature to {:1.5f}".format(_escape, self._n_same, self.parameters._T)
                logging.info(log_msg)

            MaxwellBoltzmannDistribution(atoms, temperature_K=self.parameters._T, communicator='serial')

            # check that periodic boundaries are the same in all directions (no mixed boundary conditions)
            _pbc = list(set(atoms.pbc))
            assert len(_pbc) == 1, "mixed boundary conditions"

            # if periodic boundary conditions create cell atom object
            if True in _pbc:
                logging.info("    VARIABLE CELL SHAPE SOFTENING, MD AND OPTIMIZATION ARE PERFORMED")
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
            logging.info("    MD Start")
            positions, lattice, self.parameters._dt, _md_trajectory, epot_max_md, number_of_md_steps = md.md(atoms = atoms, 
                                                                                                        calculator = self.calculator,
                                                                                                        outpath = self._outpath, 
                                                                                                        cell_atoms = cell_atoms,
                                                                                                        dt = self.parameters._dt, 
                                                                                                        n_max = self.parameters.mdmin,
                                                                                                        verbose = self.parameters.verbose_output,
                                                                                                        collect_md_file = self.collect_md_file)

            log_msg = "    MD finished after {:d} steps visiting {:d} maxima. New dt is {:1.5f}".format(number_of_md_steps, self.parameters.mdmin, self.parameters._dt)

            logging.info(log_msg)
            # Set new positions after the MD
            atoms.set_positions(positions)
            # If pbc set new lattice and reshape cell
            if True in _pbc:
                atoms.set_cell(lattice)
                lattice_operations.reshape_cell(atoms, self.parameters.symprec)

            try:
                atoms.calc.recalculateBasis(atoms)
            except:
                pass

            logging.info("    OPT start")
            positions, lattice, self._noise, _opt_trajectory, number_of_opt_steps, epot_max_geopt = opt.optimization(atoms=atoms, 
                                                                    calculator=self.calculator, 
                                                                    max_force_threshold=self.parameters.fmax, 
                                                                    outpath=self._outpath, 
                                                                    verbose=self.parameters.verbose_output)

            if epot_max_geopt > epot_max_md:
                _epot_max = epot_max_geopt
                msg = "maximal potential energy in geometry optimization"
                logging.warning(msg)
            else:
                _epot_max = epot_max_md

            # Set optimized positions
            atoms.set_positions(positions)
            # If Pbc set optimized lattice 
            if True in _pbc:
                atoms.set_cell(lattice)

            log_msg = "    OPT finished after {:d} steps.".format(number_of_opt_steps)
            logging.info(log_msg)

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
                        exclude= self.parameters.exclude,
                        calculate_fingerprint=self.parameters.calculate_fingerprint)

            # check if proposed structure is the same to the initial structure
            _escape_energy = struct.__compareto__(proposed_structure)
            if not self.parameters.compare_energies:
                _escape = struct.fingerprint_distance(proposed_structure)

            _i_steps += 1
            self._n_min += 1

            if not self.parameters.compare_energies:
                if  _escape > self.parameters.fingerprint_threshold:
                    is_escape = False
            elif _escape_energy > self.parameters._eDiff:
                is_escape = False
            else: # not escaped, same minimum found
                self.parameters._n_same += 1

        if self.parameters.compare_energies:
            log_msg = "    New minimum found with energy difference {:1.2e} after looping {:d} time(s)".format(_escape_energy, _i_steps)
        else:
            log_msg = "    New minimum found with fpd {:1.2e} after looping {:d} time(s)".format(_escape, _i_steps)
        logging.info(log_msg)

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
        logging.info(log_msg)


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

        n_visits = struct.n_visit

        if _e_pot - _e_pot_cur < self.parameters._eDiff:
            self.parameters._eDiff *= self.parameters.alpha_accept
            is_accepted = True
            ediff_acc = _e_pot - _e_pot_cur
        else:
            if self.parameters.enhanced_feedback:
                self.parameters._eDiff = self.parameters._eDiff * self.parameters.alpha_reject * (1. + 0.2 * np.log(n_visits))
            else:
                self.parameters._eDiff *= self.parameters.alpha_reject

            is_accepted = False
            ediff_rej = _e_pot - _e_pot_cur

        if is_accepted:
            log_msg = "  Minimum was accepted:  Enew - Ecur = {:1.5f} < {:1.5f} = Ediff".format(ediff_acc,
                                                                                                _ediff_in)
            logging.info(log_msg)
        else:
            log_msg = "  Minimum was rejected:  Enew - Ecur = {:1.5f} > {:1.5f} = Ediff".format(ediff_rej,
                                                                                                _ediff_in)
            logging.info(log_msg)

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
                self.parameters._T = self.parameters._T * self.parameters.beta_increase * (1. + 0.2 * np.log(float(n_visits)))
            else:
                self.parameters._T = self.parameters._T * self.parameters.beta_increase
        else:
            self.parameters._T = self.parameters._T * self.parameters.beta_decrease


    def _history_log(self, struct: Minimum, status):
        """
        Writing log message in the history file
        Input:
            struct: Minimum object
            status: str (Accept, Reject, Intermediate, etc.)
        """
        # _notunique_frac = float(self._n_notunique) / float(self._n_min)
        # _same_frac = float(self._n_same) / float(self._n_min)
        # _unique_frac = 1.0 - (_notunique_frac+_same_frac)

        totalMinima = self.parameters._n_accepted + self.parameters._n_rejected + self.parameters._n_same
        accept_ratio = self.parameters._n_accepted / totalMinima
        rejected_ratio = self.parameters._n_rejected / totalMinima
        same_ratio = self.parameters._n_same / totalMinima

        n_visits = struct.n_visit
        if n_visits is None:
            n_visits = 'None'

        label = struct.label
        if label is None:
            label = 'None'

        history_msg = "%15.8f %s %s %.2f %8.4f %.2f %.2f %.2f %s \n"%(
                                            struct.e_pot,
                                            n_visits,
                                            label,
                                            self.parameters._T,
                                            self.parameters._eDiff,
                                            accept_ratio,
                                            rejected_ratio,
                                            same_ratio,
                                            status
        )

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


    def createPathsAndSetMPIVariables(self, use_MPI, logLevel):
        if self.mpiSize > 1 and not use_MPI:
            print("Detected multiple MPI Processes but use_MPI parameter was set to false. Is this on purpose?")
        if self.mpiSize == 1 or not use_MPI: # no mpi should be used
            logging.basicConfig(level=logLevel, 
                format='%(message)s'    
            )
            self.isMaster = False
            self.isWorker = False
            self._outpath = 'output/'
            self.restart_path = self._outpath + "restart/"
            if not os.path.exists('output'):
                os.mkdir(self._outpath)
                os.mkdir(self.restart_path)
            if not os.path.exists(self._minima_path):
                os.mkdir(self._minima_path)
            if use_MPI:
                logging.error('UseMPI is true but only one rank is present. simulation will be stopped.')
                quit()
            
        else: # mpi parallelized simulation

            if self.mpiRank == 0:
                self.isMaster = True
                self.isWorker = False
                self._outpath = 'output/master/'
                self.restart_path = self._outpath + "restart/"
                # self._minima_path = 'minima/'
                if not os.path.exists("output"):
                    os.mkdir('output')
                    os.mkdir(self._outpath)
                    os.mkdir(self.restart_path)
                if not os.path.exists(self._minima_path):
                    os.mkdir(self._minima_path)
                
            if not os.path.exists('output'):
                time.sleep(1)
                if not os.path.exists('output'):
                    time.sleep(4)
                    if not os.path.exists('output'):
                        print("Failed to create an output directory. Aborting", flush=True)
                        quit()
            if not os.path.exists(self._minima_path):
                time.sleep(1)
                if not os.path.exists(self._minima_path):
                    time.sleep(4)
                    if not os.path.exists(self._minima_path):
                        print("Failed to create an output directory. Aborting", flush=True)
                        quit()

            if self.mpiRank > 0:
                self.isMaster = False
                self.isWorker = True
                self._outpath = 'output/worker_' + str(self.mpiRank) + '/' 
                self.restart_path = self._outpath + "restart/"
                if not os.path.exists(self._outpath):
                    os.mkdir(self._outpath)
                if not os.path.exists(self.restart_path):
                    os.mkdir(self.restart_path)
                MPI.COMM_WORLD.send((mpi_messages.loginRequestFromClient, 1), dest=0)
            
            logging.basicConfig(filename=self._outpath + 'minimahopping.log', 
                level=logLevel, 
                format='%(message)s'    
            )


    def sigTermCatcher(self, *args):
        if self.isMaster:
            logging.info("Received sigterm on master. Closing files and send mpi abort to comm_world")
        else:
            logging.info("Received sigterm. I will close files and exit.")
        # if self.isMaster:
        #     MPI.COMM_WORLD.Abort()
        sys.exit()

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
        logging.info(msg)
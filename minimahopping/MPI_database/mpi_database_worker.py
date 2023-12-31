import minimahopping.mh.minimum as minimum
from mpi4py import MPI
import minimahopping.MPI_database.mpi_messages as message
import minimahopping.logging.logger as logging


class Database():
    def __init__(self, energy_threshold: float, minima_threshold: float, output_n_lowest_minima: int
                 , is_restart: bool = False, outpath: str = './'
                 , minima_path: str = "lowest_minima/", write_graph_output: bool = True):
        self.energy_threshold = energy_threshold
        self.minima_threshold = minima_threshold
        self.output_n_lowest_minima = output_n_lowest_minima
        self.is_restart = is_restart
        self.outpath = outpath
        self.minima_path = minima_path
        self.write_graph_output = write_graph_output

        self.comm_world = MPI.COMM_WORLD

    def __enter__(self):
        return self

    def __exit__(self,exc_type, exc_value, exc_traceback):
        pass

    def addElement(self,struct: minimum.Minimum):
        minimum_copy = struct.__copy__()
        logging.logger.debug("Sending minimum to master: %s from rank %i"%(str(message.addelement), self.comm_world.Get_rank()))
        self.comm_world.send((message.addelement, minimum_copy), dest=0)
        logging.logger.debug("Waiting for answer from add element of master")
        n_visit, label, continueSimulation = self.comm_world.recv(source=0)
        logging.logger.debug("Received answer from add element of master, n_vistit: %i, label: %i, continue: %r"%(n_visit, label, continueSimulation))
        struct.n_visit = n_visit
        struct.label = label
        return n_visit, label, continueSimulation

    def addElementandConnectGraph(self, currentMinimum: minimum.Minimum, escapedMinimum: minimum.Minimum, trajectory: list, epot_max: float):
        self.comm_world.send((message.addElementandConnectGraph, [currentMinimum.__copy__(), escapedMinimum.__copy__(), trajectory, epot_max]), dest=0)
        n_visit, label, continueSimulation = self.comm_world.recv(source=0)
        escapedMinimum.n_visit = n_visit
        escapedMinimum.label = label
        return n_visit, label, continueSimulation

    def get_element(self, index: int):
        self.comm_world.send((message.get_element, index), dest=0)
        element = self.comm_world.recv(source=0)
        return element

    def get_element_index(self, structure: minimum.Minimum):
        self.comm_world.send((message.get_element_index, structure.__copy__()), dest = 0)
        return self.comm_world.recv(source = 0)



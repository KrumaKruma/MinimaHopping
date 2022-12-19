import minimahopping.mh.minimum as minimum
from mpi4py import MPI
import minimahopping.MPI_database.mpi_messages as message


class Database():
    def __init__(self,energy_threshold, minima_threshold, output_n_lowest_minima, is_restart = False, outpath='./', minima_path= "lowest_minima/", write_graph_output = True):
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
        self.comm_world.send((message.addelement, minimum_copy), dest=0)
        n_visit, label = self.comm_world.recv(source=0)
        struct.n_visit = n_visit
        struct.label = label
        return n_visit, label

    def addElementandConnectGraph(self, currentMinimum: minimum.Minimum, escapedMinimum: minimum.Minimum, trajectory, epot_max):
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



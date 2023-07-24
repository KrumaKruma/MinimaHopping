import networkx as nx
import shelve
import pickle
import os
import matplotlib.pyplot as plt
from ase import Atoms
import copy

graphDotName = 'graph.dot'

class MinimaHoppingGraph:
    
    def __init__(self, graphFileName, trajectoryDatabaseName, is_restart) -> None:
        self.graphFileName = graphFileName
        self.trajectoryDatabaseName = trajectoryDatabaseName
        self.is_restart = is_restart
        self.trajectoryDict = None
        self.graph = None
        self.lucky_counter = 0

    def __enter__(self):
        return self.read_from_disk()

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.write_to_disk()

    def read_from_disk(self):
        if self.is_restart:
            graph_pickle = open(self.graphFileName, "rb")
            self.graph = pickle.load(graph_pickle)
            graph_pickle.close()
        else:
            self.graph = nx.DiGraph()
            try:
                os.remove(self.trajectoryDatabaseName)
            except FileNotFoundError:
                pass
        self.trajectoryDict = shelve.open(self.trajectoryDatabaseName)
        return self

    def write_to_disk(self):
        """
        Writes the graph to the disk and closes the trajectory data shelve.
        """
        self.trajectoryDict.close()
        self.write_restart_files()
        nx.drawing.nx_pydot.write_dot(self.graph, graphDotName)
    

    def write_restart_files(self):
        """
        Writes the graph to the disk and updates the trajectory data shelve.
        """
        # self.trajectoryDict.sync()
        graph_pickle = open (self.graphFileName, "wb")
        pickle.dump(self.graph, graph_pickle)
        graph_pickle.close()


    def addStructure(self, initialStuctureLabel, structureLabel, trajectory, e_old, e_new, e_max):
        if initialStuctureLabel == structureLabel:
            return
        weight_old = e_max - e_old
        weight_new = e_max - e_new
        if weight_new < 0:
            weight_new = float('inf')
            return
        if weight_old < 0:
            weight_old = float('inf')
            return
        mt = copy.deepcopy(trajectory)
        reverse_trajectory = copy.deepcopy(trajectory)
        reverse_trajectory.reverse()

        restart_file_update_necessary = True

        if not self.graph.has_node(structureLabel):
            if self.graph.size() == 0:
                self.graph.add_node(initialStuctureLabel, energy = e_old)
            self.graph.add_node(structureLabel, energy = e_new)

        if self.graph.has_edge(initialStuctureLabel, structureLabel):
            if self.graph[initialStuctureLabel][structureLabel]['weight'] > weight_old:
                self.graph.remove_edge(initialStuctureLabel, structureLabel)
                self.graph.remove_edge(structureLabel, initialStuctureLabel)
                self._add_edge(initialStuctureLabel, structureLabel, weight=weight_old, trajectory=mt)
                self._add_edge(structureLabel, initialStuctureLabel, weight=weight_new, trajectory=reverse_trajectory)
            else:
                restart_file_update_necessary = False
        else: 
            self._add_edge(initialStuctureLabel, structureLabel, weight=weight_old, trajectory=mt)
            self._add_edge(structureLabel, initialStuctureLabel, weight=weight_new, trajectory=reverse_trajectory)
        if restart_file_update_necessary:
            self.write_restart_files()

    def get_lowest_energy(self):
        emin = float("inf")
        ind = -1
        for n in self.graph.nodes:
            e = float(self.graph.nodes[n]['energy'])
            if emin > float(e):
                emin = e
                ind = n
        self.graph.nodes[ind]['color'] = 'red'
        return emin, ind

    def shift_energy_to_zero(self):
        emin, ind = self.get_lowest_energy()
        for n in self.graph.nodes:
            self.graph.nodes[n]['shifted_energy'] = float(self.graph.nodes[n]['energy']) - emin

    def remove_leaves(self, number_of_iterations=1):
        graph_copy = copy.deepcopy(self.graph)
        nx.set_node_attributes(graph_copy, 0.5, 'width')
        nx.set_node_attributes(graph_copy, 0.5, 'height')
        for i in range(number_of_iterations):
            remove = [node for node, degree in graph_copy.degree() if degree <= 2]
            for i in remove:
                for v in graph_copy.neighbors(i):
                    # print(i, graph_copy.edges(i), v)
                    graph_copy.nodes[v]['width'] = graph_copy.nodes[v]['width'] + 0.05 / graph_copy.nodes[v]['width']
                    graph_copy.nodes[v]['height'] = graph_copy.nodes[v]['height'] + 0.05 / graph_copy.nodes[v]['height']
            graph_copy.remove_nodes_from(remove)
        return graph_copy
    
    def _add_edge(self, initialStuctureLabel, structureLabel, weight, trajectory):
        self.graph.add_edge(initialStuctureLabel, structureLabel, weight=weight)
        self.trajectoryDict.sync()
        self.trajectoryDict[self._getEdgeString(initialStuctureLabel, structureLabel)] = trajectory
        self.trajectoryDict.sync()

    def _getEdgeString(self, initialStuctureLabel, structureLabel):
        edgeString = str(initialStuctureLabel) + ',' + str(structureLabel)
        return edgeString

    def shortestPath(self, a, b):
        return nx.shortest_path(self.graph, a, b, weight="weight")

    def getTrajectoryList(self, a, b):
        path = self.shortestPath(a, b)
        TList = []
        for i in range(1, len(path)):
            TList = TList + copy.deepcopy(self.trajectoryDict[self._getEdgeString(path[i - 1], path[i])])
        return TList

    def draw(self):
        nx.draw(self.graph, with_labels=True, font_weight='bold')
        plt.show()
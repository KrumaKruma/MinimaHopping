import networkx as nx
import shelve
import pickle
import os
import matplotlib.pyplot as plt
from ase import Atoms
import copy


class MinimaHoppingGraph:
    
    def __init__(self, graphFileName, trajectoryDatabaseName, restart) -> None:
        self.graphFileName = graphFileName
        self.trajectoryDatabaseName = trajectoryDatabaseName
        self.restart = restart
        self.trajectoryDict = None
        self.graph = None

    def __enter__(self):
        if self.restart:
            graph_pickle = open(self.graphFileName, "rb")
            self.graph = pickle.load(graph_pickle)
            graph_pickle.close()
        else:
            self.graph = nx.DiGraph()
            try:
                os.remove(self.trajectoryDatabaseName)
            except FileNotFoundError:
                #print('File not found')
                pass
        self.trajectoryDict = shelve.open(self.trajectoryDatabaseName,  writeback=True)
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.trajectoryDict.close()
        graph_pickle = open (self.graphFileName, "wb")
        pickle.dump(self.graph, graph_pickle)
        graph_pickle.close()


    def addStructure(self, initialStuctureLabel, structureLabel, trajectory, e_old, e_new, e_max):
        if initialStuctureLabel == structureLabel:
            return
        weight_old = e_max - e_old
        weight_new = e_max - e_new
        mt = copy.deepcopy(trajectory)
        reverse_trajectory = copy.deepcopy(trajectory)
        reverse_trajectory.reverse()

        if not self.graph.has_node(structureLabel):
            if self.graph.size() == 0:
                self.graph.add_node(initialStuctureLabel, energy = e_old)
            else:
                pass
            self.graph.add_node(structureLabel, energy = e_new)

        if self.graph.has_edge(initialStuctureLabel, structureLabel):
            if self.graph[initialStuctureLabel][structureLabel]['weight'] > weight_old:
                self.graph.remove_edge(initialStuctureLabel, structureLabel)
                self.graph.remove_edge(structureLabel, initialStuctureLabel)
                self._add_edge(initialStuctureLabel, structureLabel, weight=weight_old, trajectory=mt)
                self._add_edge(structureLabel, initialStuctureLabel, weight=weight_new, trajectory=reverse_trajectory)
        else: 
            self._add_edge(initialStuctureLabel, structureLabel, weight=weight_old, trajectory=mt)
            self._add_edge(structureLabel, initialStuctureLabel, weight=weight_new, trajectory=reverse_trajectory)


    def _add_edge(self, initialStuctureLabel, structureLabel, weight, trajectory):
        self.graph.add_edge(initialStuctureLabel, structureLabel, weight=weight)
        self.trajectoryDict.sync()
        self.trajectoryDict[self._getEdgeString(initialStuctureLabel, structureLabel)] = trajectory
        self.trajectoryDict.sync()

    def _getEdgeString(self, initialStuctureLabel, structureLabel):
        #minimum = min(structureLabel, initialStuctureLabel)
        #maximum = max(structureLabel, initialStuctureLabel)
        edgeString = str(initialStuctureLabel) + ',' + str(structureLabel)
        return edgeString

    def shortestPath(self, a, b):
        return nx.shortest_path(self.graph, a, b, weight="weight")

    def getTrajectoryList(self, a, b):
        path = self.shortestPath(a, b)
        TList = []
        #print(dict(self.trajectoryDict))
        for i in range(1, len(path)):
            self.trajectoryDict[self._getEdgeString(path[i], path[i-1])]
            #temp_trajectory = self.trajectoryDict[self._getEdgeString(path[i-1], path[i])]
            #print(temp_trajectory)
            TList = TList + copy.deepcopy(self.trajectoryDict[self._getEdgeString(path[i - 1], path[i])])
            #print(i-1, i)
        return TList

    def draw(self):
        nx.draw(self.graph, with_labels=True, font_weight='bold')
        plt.show()
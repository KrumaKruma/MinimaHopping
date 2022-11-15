import networkx as nx
import shelve
import pickle
import os
import matplotlib.pyplot as plt
from ase import Atoms


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
            self.graph = nx.Graph()
            try:
                os.remove(self.trajectoryDatabaseName)
            except FileNotFoundError:
                print('File not found')
        self.trajectoryDict = shelve.open(self.trajectoryDatabaseName)
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.trajectoryDict.close()
        graph_pickle = open (self.graphFileName, "wb")
        pickle.dump(self.graph, graph_pickle)
        graph_pickle.close()


    def addStructure(self, structureLabel, initialStuctureLabel, weight, trajectory):
        if self.graph.has_edge(structureLabel, initialStuctureLabel):
            if self.graph[structureLabel][initialStuctureLabel]['weight'] > weight:
                self.graph.remove_edge(structureLabel, initialStuctureLabel)
                self._add_edge(structureLabel, initialStuctureLabel, weight=weight, trajectory=trajectory)
        else:
            self._add_edge(structureLabel, initialStuctureLabel, weight=weight, trajectory=trajectory)


    def _add_edge(self, structureLabel, initialStuctureLabel, weight, trajectory):
        self.graph.add_edge(structureLabel, initialStuctureLabel, weight=weight)
        self.trajectoryDict[self._getEdgeString(initialStuctureLabel, structureLabel)] = trajectory

    def _getEdgeString(self, structureLabel, initialStuctureLabel):
        minimum = min(structureLabel, initialStuctureLabel)
        maximum = max(structureLabel, initialStuctureLabel)
        edgeString = str(minimum) + ',' + str(maximum)
        return edgeString

    def shortestPath(self, a, b):
        return nx.shortest_path(self.graph, a, b, weight=lambda x, y, z: z['weight'] + 11.37)

    def getTrajectoryList(self, a, b):
        path = self.shortestPath(a, b)
        TList = []
        #print(dict(self.trajectoryDict))
        for i in range(1, len(path)):
            self.trajectoryDict[self._getEdgeString(path[i], path[i-1])]
            temp_trajectory = self.trajectoryDict[self._getEdgeString(path[i], path[i-1])]
            #print(temp_trajectory)
            TList = TList + self.trajectoryDict[self._getEdgeString(path[i], path[i-1])] 
        return TList

    def draw(self):
        nx.draw(self.graph, with_labels=True, font_weight='bold')
        plt.show()
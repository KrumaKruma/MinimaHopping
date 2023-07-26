import networkx as nx
import numpy as np
import shelve
import pickle
import os
from ase import Atoms
import copy
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
import matplotlib
import matplotlib.pyplot as plt

graphDotName = 'output/graph.dot'

class MinimaHoppingGraph:
    
    def __init__(self, graphFileName: str, trajectoryDatabaseName: str, is_restart: bool) -> None:
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


    def addStructure(self, initialStuctureLabel: int, structureLabel: int, trajectory: list, e_old: float, e_new: float, e_max: float):
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
                self.graph.add_node(initialStuctureLabel, energy = e_old
                                    , removed_leaves = 0, width=0.5, height = 0.5, num_atoms=len(trajectory[1]))
            self.graph.add_node(structureLabel, energy = e_new
                                , removed_leaves = 0, width=0.5, height = 0.5, num_atoms=len(trajectory[1]))

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
        return get_lowest_energy_static(self.graph)

    def shift_energy_to_zero(self):
        shift_energy_to_zero_static(self.graph)

    def remove_leaves(self, number_of_iterations: int = 1):
        return remove_leaves_static(self.graph)
    
    def _add_edge(self, initialStuctureLabel: int, structureLabel: int, weight: float, trajectory: list):
        self.graph.add_edge(initialStuctureLabel, structureLabel, weight=weight)
        self.trajectoryDict.sync()
        self.trajectoryDict[self._getEdgeString(initialStuctureLabel, structureLabel)] = trajectory
        self.trajectoryDict.sync()

    def _getEdgeString(self, initialStuctureLabel: int, structureLabel: int):
        edgeString = str(initialStuctureLabel) + ',' + str(structureLabel)
        return edgeString

    def shortestPath(self, a: int, b: int):
        return nx.shortest_path(self.graph, a, b, weight="weight")

    def getTrajectoryList(self, a: int, b: int):
        path = self.shortestPath(a, b)
        TList = []
        for i in range(1, len(path)):
            TList = TList + copy.deepcopy(self.trajectoryDict[self._getEdgeString(path[i - 1], path[i])])
        return TList
    
    def getTrajectoryListFromPath(self, path: list):
        TList = []
        for i in range(1, len(path)):
            TList = TList + copy.deepcopy(self.trajectoryDict[self._getEdgeString(path[i - 1], path[i])])
        return TList

    def draw(self):
        nx.draw(self.graph, with_labels=True, font_weight='bold')
        plt.show()

def shift_energy_to_zero_static(graph: nx.DiGraph):
    emin, ind = get_lowest_energy_static(graph)
    for n in graph.nodes:
        graph.nodes[n]['shifted_energy'] = float(graph.nodes[n]['energy']) - emin

def get_lowest_energy_static(graph: nx.DiGraph):
    emin = float("inf")
    ind = -1
    for n in graph.nodes:
        e = float(graph.nodes[n]['energy'])
        if emin > float(e):
            emin = e
            ind = n
    # graph.nodes[ind]['color'] = 'red'
    return emin, ind

def color_graph(graph: nx.DiGraph):
    cmap = matplotlib.cm.get_cmap('Reds')
    shift_energy_to_zero_static(graph)
    nx.set_node_attributes(graph, 'blue', 'fillcolor')
    nx.set_node_attributes(graph, 'filled', 'style')
    for n in graph.nodes():
        intensity = graph.nodes()[n]['shifted_energy'] / graph.nodes()[n]['num_atoms']
        intensity = intensity * 10
        graph.nodes()[n]['fillcolor'] = matplotlib.colors.rgb2hex(cmap( intensity ))

    # todo: create a colorbar pdf
    cm = 1/2.54
    gradient = np.linspace(0, 1, 101, endpoint=True)
    gradient = np.vstack((gradient, gradient))
    fig, axs = plt.subplots(figsize=( 6.8 * cm, 2.8* cm))
    # fig, axs = plt.subplots()
    # fig.set_size_inches(6*cm, 5*cm)
    axs.set_xlabel('meV per Atom', fontsize=10)
    plt.yticks([])
    plt.xticks([0, 20, 40, 60, 80, 100])
    axs.imshow(gradient, aspect=7.0, cmap=cmap)
    # plt.show()
    fig.savefig('coloourbar.pdf')

def remove_leaves_static(graph, number_of_iterations: int = 1):
    graph_copy = copy.deepcopy(graph)
    for i in range(number_of_iterations):
        remove = [node for node, degree in graph_copy.degree() if degree <= 2]
        for i in remove:
            for v in graph_copy.neighbors(i):
                # print(i, graph_copy.edges(i), v)
                graph_copy.nodes[v]['width'] = graph_copy.nodes[v]['width'] + 0.05 / graph_copy.nodes[v]['width']
                graph_copy.nodes[v]['height'] = graph_copy.nodes[v]['height'] + 0.05 / graph_copy.nodes[v]['height']
                graph_copy.nodes[v]['removed_leaves'] += 1
        graph_copy.remove_nodes_from(remove)
    return graph_copy

def contract(g_directed: nx.DiGraph, contraction_iterations: int = 1):
    """
    Contract chains of neighbouring vertices with degree 2 into a single edge.

    Arguments:
    ----------
    g -- networkx.Graph or networkx.DiGraph instance

    contraction_iterations -- int

    Returns:
    --------
    h -- networkx.Graph or networkx.DiGraph instance
        the contracted graph

    """

    # create subgraph of all nodes with degree 2
    print('number of nodes in iteration 0: %i'%(len(g_directed)))
    g = remove_leaves_static(g_directed)
    g = g.to_undirected()
    for i in range(contraction_iterations):
        print('number of nodes in iteration %i: %i'%(i+1, len(g)))
        for e in g.edges():
            g[e[0]][e[1]]['weight'] = 1.0
        is_chain = [node for node, degree in g.degree() if degree == 2]
        chains = g.subgraph(is_chain)

        # contract connected components (which should be chains of variable length) into single node
        # components = list(nx.components.connected_component_subgraphs(chains))
        components = [chains.subgraph(c) for c in nx.components.connected_components(chains)]

        hyper_edges = []
        for component in components:
            end_points = [node for node, degree in component.degree()
                          if degree < 2]
            candidates = set(
                [neighbor for node in end_points for neighbor in g.neighbors(node)])
            connectors = candidates - set(list(component.nodes()))
            hyper_edge = list(connectors)
            weights = [component.get_edge_data(
                *edge)['weight'] for edge in component.edges()]
            if len(hyper_edge) >= 2: # the if clause is necessary because the graph may contain loops
                hyper_edges.append((hyper_edge, np.sum(weights)))

        # initialise new graph with all other nodes
        not_chain = [node for node in g.nodes() if not node in is_chain]
        h = g.subgraph(not_chain).copy()
        for hyper_edge, weight in hyper_edges:
            # print(*hyper_edge, weight)
            h.add_edge(*hyper_edge, weight=weight)
        g = remove_leaves_static(h)
    return h

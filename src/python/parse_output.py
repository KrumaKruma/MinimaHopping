#!/usr/bin/env python3
import networkx as nx
import matplotlib.pyplot as plt
from ase.io import write
import pydot



import graph

with graph.MinimaHoppingGraph('graph.dat', 'trajectories.dat', True) as g:
    g.draw()
    l = g.shortestPath(0, 14)
    tl = g.getTrajectoryList(0, 14)

    print(l)

    f = open('good_trajectory.extxyz', 'w')
    f.close()

    write('good_trajectory.extxyz', tl, append = True)
    
    graph = nx.drawing.nx_pydot.to_pydot(g.graph)
    graph.write_png('output.pdf')
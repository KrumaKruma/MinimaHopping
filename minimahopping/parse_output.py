#!/usr/bin/env python3
import networkx as nx
import matplotlib.pyplot as plt
from ase.io import write
import pydot
import pygraphviz


from minimahopping.graph import graph


def draw_pygraphviz(g, filename, layout='fdp'):
    g.graph_attr['concentrate'] = 'true'
    # g.graph_attr['beautify'] = 'true'
    g.layout(layout)
    g.draw(filename)
    # g.write('file.dot')
    

g = graph.MinimaHoppingGraph('output/graph.dat', 'output/trajectory.dat', True)
g.read_from_disk()
n = 28
# l = g.shortestPath(0, n)
# tl = g.getTrajectoryList(0, n)

# print(l)

# f = open('output/good_trajectory.extxyz', 'w')
# f.close()

# write('output/good_trajectory.extxyz', tl, append = True)

emin, ind = g.get_lowest_energy()

g.shift_energy_to_zero()

stripped_graph = g.remove_leaves()

draw_pygraphviz(nx.nx_agraph.to_agraph(stripped_graph), 'no_leaves.pdf')
draw_pygraphviz(nx.nx_agraph.to_agraph(g.graph), 'with_leaves.pdf')

# graph.write('file.dot')


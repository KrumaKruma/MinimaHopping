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
    

g = graph.MinimaHoppingGraph('output/master/restart/graph.dat', 'output/master/restart/trajectory.dat', True)
g.read_from_disk()
n1 = 289
n2 = 432
print('size', g.graph.size())
l = g.shortestPath(n1, n2)
tl = g.getTrajectoryList(n1, n2)

print(l)

f = open('output/good_trajectory.extxyz', 'w')
f.close()

write('output/good_trajectory.extxyz', tl, append = True)

emin, ind = g.get_lowest_energy()

g.shift_energy_to_zero()

stripped_graph = g.remove_leaves()

draw_pygraphviz(nx.nx_agraph.to_agraph(stripped_graph), 'no_leaves.pdf', layout='fdp')
draw_pygraphviz(nx.nx_agraph.to_agraph(g.graph), 'with_leaves.pdf', layout='fdp')

# graph.write('file.dot')


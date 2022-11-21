#!/usr/bin/env python3
import networkx as nx
import matplotlib.pyplot as plt
from ase.io import write
import pydot
import pygraphviz


import graph

with graph.MinimaHoppingGraph('output/graph.dat', 'output/trajectory.dat', True) as g:
    # g.draw()
    n = 40
    l = g.shortestPath(0, n)
    tl = g.getTrajectoryList(0, n)

    print(l)

    f = open('output/good_trajectory.extxyz', 'w')
    f.close()

    write('output/good_trajectory.extxyz', tl, append = True)


    graph = nx.nx_agraph.to_agraph(g.graph)

    graph.nodes()

    emin = 10000.0
    ind = -1

    for n in graph.nodes():
        try:
            e = float(n.attr['energy'])
            #print(e)
            if e == '':
                e = 0.0
            if emin > float(e):
                emin = float(n.attr['energy'])
                ind = n.get_name()
        except:
            pass

    graph.get_node(ind).attr['color'] = 'red'
    graph.get_node(ind).attr['root'] = 'true'


    # e0 = emin + .8

    # for n in graph.nodes():
    #     graph.get_node(n.get_name()).attr['width'] = '.5'
    #     graph.get_node(n.get_name()).attr['heigt'] = '.5'
    #     try:
    #         e = float(n.attr['energy'])
    #         if e > e0:
    #             ind1 = n.get_name()
    #             graph.remove_node(ind1)
    #     except:
    #         pass


    graph.graph_attr['concentrate'] = 'true'
    #graph.graph_attr['K'] = 0.1
    #graph.graph_attr['repulsiveforce'] = 2.0
    #graph.graph_attr['beautify'] = 'true'
    graph.layout('fdp')
    graph.draw('networkx_graph.pdf')
    graph.write('file.dot')


    #graph = nx.drawing.nx_pydot.to_pydot(g.graph)
    #graph.set_strict(True)
    #graph.concentrate = True
    #print(graph.get_strict(True))
    #graph.write_png('output.pdf')
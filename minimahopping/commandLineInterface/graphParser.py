import networkx as nx
import argparse
from ase.io import read, write
from minimahopping.graph import graph as mh_graph
import os
import pickle

def main():
    parser = argparse.ArgumentParser(description ="""Utility that parses graph output, visualizes the graph and writes a
                                     .dot and a binary .dat file of the graph with all modifications.
                                     The graph_binary.dat file can be used in the graphName command line option if this script.""")
    subparsers = parser.add_subparsers(help="""Command execute.
                                       Use the -h option in combination with the command
                                       itself for a detailed command documentation.""",
                                       dest='command', required=False)
    parser.add_argument('-g', '--graphName', dest ='graphName',
                action ='store', help ='input filename.', required=False,
                default='output/master/restart/graph.dat', type=argparse.FileType('r'))
    parser.add_argument('-t', '--trajectoryName', dest ='trajectoryName',
                action ='store', help ='trajectory filename.', required=False,
                default='output/master/restart/trajectory.dat', type=argparse.FileType('r'))
    parser.add_argument('--removeLeaves', help="""Removes leaves (nodes with only one edge)
                        from the graph and writes a dot file of the new graph.
                        This is done, before any other operation is applied to the graph"""
                        , action='store_true'
                        , required=False)
    parser.add_argument('--shift2zero', help="Shifts the energies of all nodes such that the lowest energy is zero."
                        , action='store_true'
                        , required=False)

    shortestPathParser = subparsers.add_parser("shortestPath", help = "Calculates the shortest path between two nodes.")
    shortestPathParser.add_argument('-n1', type=int, required=True, help='Label of first minima', action='store', dest='n1')
    shortestPathParser.add_argument('-n2', type=int, required=True, help='Label of second minima', action='store', dest='n1')

    plotParser = subparsers.add_parser('plotGraph'
                    , help="Creates a pdf represantation of the graph. Requires graphviz and pygraphviz installed on your system.")
    plotParser.add_argument('--layout', help="""pygraphviz layout of graph. 
                            Must be one of: dot, neato, fdp, sfdp, circo,
                            twopi, nop, nop2,osage, patchwork"""
                            , type=str, required=False, default='fdp', action='store', dest='layout')

    listParser = subparsers.add_parser("listPath", help="""Returns a trajectory (.extxyz) file
                                       that connects the nodes that are specified in remaining arguments. If there is no connection between
                                       any of the nodes specified, this operation fails.
                                       
                                       Example usage:
                                       graphParsery listPath 5 3 7
                                       creates an extxyz file that contains the edge 5->3, 3->7""")
    listParser.add_argument('edges', type=int, nargs=argparse.REMAINDER, help="""
                            List of edges.""")

    args = parser.parse_args()

    g = mh_graph.MinimaHoppingGraph(args.graphName.name, args.trajectoryName.name, True)
    g.read_from_disk()
    graph = g.graph
    
    # apply all opations:
    if args.shift2zero:
        g.shift_energy_to_zero()
    if args.removeLeaves:
        graph = g.remove_leaves()

    # write dot file with all operations applied in text and binary form.
    nx.drawing.nx_pydot.write_dot(graph, 'graph.dot')
    with open('graph_binary.dat', 'wb') as graph_pickle:
        pickle.dump(graph, graph_pickle)


    # execute commands
    if args.command == 'shortestPath':
        shortestPath(g, args.n1, args.n2)
    elif args.command == 'plotGraph':
        pygraphviz_graph = nx.nx_agraph.to_agraph(graph)
        pygraphviz_graph.graph_attr['concentrate'] = 'true'
        pygraphviz_graph.layout(args.layout)
        pygraphviz_graph.draw('graph.pdf')
        stripped_graph = g.remove_leaves()
        pygraphviz_graph = nx.nx_agraph.to_agraph(stripped_graph)
        pygraphviz_graph.graph_attr['concentrate'] = 'true'
        pygraphviz_graph.layout(args.layout)
        pygraphviz_graph.draw('stripped_graph.pdf')
    elif args.command == 'listPath':
        if len(args.edges) < 2:
            print("At least two nodes required to make the trajectory. Aborting")
            quit()
        trajectory_list = g.getTrajectoryListFromPath(args.edges)
        listString = ''
        for i in args.edges:
            listString += "_%i"%(i)
        filename = "connection" + listString + '.extxyz'
        if os.path.exists(filename):
            os.remove(filename)
        write(filename, trajectory_list, append=True)
        

def shortestPath(g: mh_graph.MinimaHoppingGraph, n1: int, n2: int):
    trajectory_list = g.getTrajectoryList(n1, n2)
    filename = "connection_%i_%i.extxyz"%(n1, n2)
    if os.path.exists(filename):
        os.remove(filename)
    write(filename, trajectory_list, append=True)

if __name__ =='__main__':
    main()
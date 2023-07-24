import networkx as nx
import argparse
from ase.io import read, write
from minimahopping.graph import graph

def main():
    parser = argparse.ArgumentParser(description ='Utility that parses graph output, visualizes the graph and writes .dot file')
    subparsers = parser.add_subparsers(help="Command that database will execute", dest='command', required=False)
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

    args = parser.parse_args()

    g = graph.MinimaHoppingGraph(args.graphName, args.trajectroyName, True)
    g.read_from_disk()
    graph = g.graph
    
    # apply all opations:
    if args.shift2zero:
        g.shift_energy_to_zero()
    if args.removeLeaves:
        graph = g.remove_leaves()

    # write dot file with all operations applied
    nx.drawing.nx_pydot.write_dot(graph, 'graph.dot')

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

def shortestPath(g: graph.MinimaHoppingGraph, n1: int, n2: int):
    tl = g.getTrajectoryList(n1, n2)
    filename = "connection_%i_%i.extxyz"%(n1, n2)


if __name__ =='__main__':
    main()
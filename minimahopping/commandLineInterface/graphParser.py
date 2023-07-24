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
    parser.add_argument('--removeLeaves', help="Removes leaves (nodes with only one edge) from the graph", action='store_true')
    

    shortestPathParser = subparsers.add_parser("shortestPath", help = "Calculates the shortest path between two nodes.")
    shortestPathParser.add_argument('-n1', type=int, required=True, help='Label of first minima', action='store', dest='n1')
    shortestPathParser.add_argument('-n2', type=int, required=True, help='Label of second minima', action='store', dest='n1')

    args = parser.parse_args()

    g = graph.MinimaHoppingGraph(args.graphName, args.trajectroyName, True)
    g.read_from_disk()

    # pygraphviz_graph = nx.nx_agraph.to_agraph(g.graph)
    # pygraphviz_graph.graph_attr['concentrate'] = 'true'
    # pygraphviz_graph.layout(layout)
    # pygraphviz_graph.draw(filename)


    if args.command == 'shortestPath':
        shortestPath(g, args.n1, args.n2)


def shortestPath(g: graph.MinimaHoppingGraph, n1: int, n2: int):
    tl = g.getTrajectoryList(n1, n2)
    filename = "connection_%i_%i.extxyz"%(n1, n2)


if __name__ =='__main__':
    main()
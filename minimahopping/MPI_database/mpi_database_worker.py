import minimahopping.mh.minimum as minimum
import minimahopping.graph.graph

class Database():
    def __init__(self,energy_threshold, minima_threshold, output_n_lowest_minima, is_restart = False, outpath='./', minima_path= "lowest_minima/", write_graph_output = True):
        self.energy_threshold = energy_threshold
        self.minima_threshold = minima_threshold
        self.output_n_lowest_minima = output_n_lowest_minima
        self.is_restart = is_restart
        self.outpath = outpath
        self.minima_path = minima_path
        self.write_graph_output = write_graph_output

    def __enter__(self):
        return self

    def __exit__(self,exc_type, exc_value, exc_traceback):
        pass

    def addElement(self,struct: minimum.Minimum):
        pass

    def addElementandConnectGraph(self, currentMinimum: minimum.Minimum, escapedMinimum: minimum.Minimum, trajectory, epot_max):
        pass

    def get_element(self, struct: minimum.Minimum):
        index = 0
        return index



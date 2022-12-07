import bisect
import shelve
import minimahopping.mh.minimum as minimum
import minimahopping.graph.graph

class Database():
    def __init__(self,energy_threshold, minima_threshold, output_n_lowest_minima, is_restart = False, outpath='./', minima_path= "lowest_minima/"):
        self.unique_minima_sorted = []
        self.nstructs = 0

        self.energy_threshold = energy_threshold
        self.minima_threshold = minima_threshold
        self.output_n_lowest_minima = output_n_lowest_minima
        self.is_restart = is_restart
        self.outpath = outpath
        self.minima_path = minima_path

        self.minima_shelve = None


        self.graphFilename = self.outpath + "graph.dat"
        self.graphTrajectoryFilename = self.outpath + 'trajectory.dat'
        self.graph = minimahopping.graph.graph.MinimaHoppingGraph(self.graphFilename, self.graphTrajectoryFilename, self.is_restart)


    def __enter__(self):
        self.read_restart_files()
        self.graph.read_from_disk()
        return self


    def __exit__(self,exc_type, exc_value, exc_traceback):
        self.graph.write_to_disk()
        self.minima_shelve.close()
        

    def read_restart_files(self):        
        filename = self.outpath + "minima.pickle.shelve"
        self.minima_shelve = shelve.open(filename)
        if self.is_restart:
            # print('asdf', list(dict(self.minima_shelve).values()).sort())
            self.unique_minima_sorted = list(dict(self.minima_shelve).values())
            self.unique_minima_sorted.sort()
            self.nstructs = len(self.unique_minima_sorted)


    def addElement(self,struct: minimum.Minimum):
        index = self.get_element(struct=struct)
        already_found = self.contains(index=index)

        if already_found:
            self.unique_minima_sorted[index].n_visit += 1
            struct.n_visit = self.unique_minima_sorted[index].n_visit
            struct.label = self.unique_minima_sorted[index].label
            self.minima_shelve[str(self.unique_minima_sorted[index].label)] = self.unique_minima_sorted[index]
        else:
            label = self.nstructs
            struct.set_label(label)
            self.nstructs += 1
            struct.n_visit = 1
            struct1 = struct.__copy__()
            struct1.atoms.set_momenta(None)
            struct1.atoms.info['energy'] = struct.e_pot
            struct1.atoms.info['label'] = label
            index = bisect.bisect_left(self.unique_minima_sorted, struct1)
            self.unique_minima_sorted.insert(index, struct1)

            if index < self.output_n_lowest_minima:
                self._write_poslow(self.output_n_lowest_minima, self.minima_path)

            self.minima_shelve[str(label)] = struct1
        return

    def addElementandConnectGraph(self, currentMinimum: minimum.Minimum, escapedMinimum: minimum.Minimum, trajectory, epot_max):
        self.addElement(escapedMinimum)
        self.graph.addStructure(currentMinimum.label, escapedMinimum.label, trajectory, currentMinimum.e_pot, escapedMinimum.e_pot, epot_max)


    def get_element(self, struct: minimum.Minimum):
        
        indices = self.get_index_energyrange(struct)

        min_dist = 10e10
        index = -1

        for i_compare in indices:
            s = self.unique_minima_sorted[i_compare]
            energy_difference = struct.__compareto__(s)
            fp_dist = struct.__equals__(s)
            if fp_dist < self.minima_threshold:
                if fp_dist < min_dist:
                    min_dist = energy_difference
                    index = i_compare

        return index

    def get_index_energyrange(self, struct):
        if self.nstructs == 0:
            return []

        i_start = bisect.bisect(self.unique_minima_sorted, struct)

        if i_start >= len(self.unique_minima_sorted):
            i_start -= 1

        indices = []
        
        # #backward
        energy_difference = 0
        for i_compare in range(i_start-1, -1, -1):
            if energy_difference > self.energy_threshold:
                break
            else:
                s = self.unique_minima_sorted[i_compare]
                energy_difference = struct.__compareto__(s)
                #if fp_dist < self.minima_threshold:
                if energy_difference < self.energy_threshold:
                    indices.append(i_compare)

        energy_difference = 0
        #forward
        i_compare = i_start
        for i_compare in range(i_start, self.nstructs, 1):
            if energy_difference > self.energy_threshold:
                break
            else:
                s = self.unique_minima_sorted[i_compare]
                energy_difference = struct.__compareto__(s)
                if energy_difference < self.energy_threshold:
                    indices.append(i_compare)

        return indices

    def contains(self, index):
        if index < 0:
            already_found = False
        else:
            already_found = True
        return already_found

    def _write_poslow(self, n_poslow, path):
        i_poslow = 0
        for s in self.unique_minima_sorted:
            filename = 'min'+str(i_poslow).zfill(6)
            if True in s.atoms.pbc:
                filename += '.ascii'
            else:
                filename += '.xyz'
            filename = path + filename
            s.write(filename, append=False)
            i_poslow += 1
            if i_poslow-1 > n_poslow:
                break

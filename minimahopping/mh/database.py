import bisect
import shelve
import minimahopping.mh.minimum as minimum
import minimahopping.graph.graph
import logging
import time

class Database():
    def __init__(self,energy_threshold, minima_threshold, output_n_lowest_minima, is_restart = False, outpath='./', minima_path= "lowest_minima/", write_graph_output = True, compare_energies = False):
        self.unique_minima_sorted = []
        self.nstructs = 0

        self.energy_threshold = energy_threshold
        self.minima_threshold = minima_threshold
        self.output_n_lowest_minima = output_n_lowest_minima
        self.is_restart = is_restart
        self.outpath = outpath
        self.minima_path = minima_path
        self.write_graph_output = write_graph_output
        self.compare_energies = compare_energies

        self.minima_shelve = None

        self.noDublicatesFileName = self.minima_path + "all_minima_no_duplicates.extxyz"
        self.allMinimaFilename = self.minima_path + "all_minima.extxyz"

        if logging.root.level <= logging.DEBUG:
            self.verbosity = True
        else:
            self.verbosity = False

        if self.write_graph_output:
            self.graphFilename = self.outpath + "graph.dat"
            self.graphTrajectoryFilename = self.outpath + 'trajectory.dat'
            self.graph = minimahopping.graph.graph.MinimaHoppingGraph(self.graphFilename, self.graphTrajectoryFilename, self.is_restart)


    def __enter__(self):
        self.read_restart_files()
        self.noDublicatesFile = open(self.noDublicatesFileName, mode='a')
        self.allMinimaFile = open(self.allMinimaFilename, mode='a')
        if self.write_graph_output:
            self.graph.read_from_disk()
        return self


    def __exit__(self,exc_type, exc_value, exc_traceback):
        self.minima_shelve.close()
        self.noDublicatesFile.close()
        self.allMinimaFile.close()
        if self.write_graph_output:
            self.graph.trajectoryDict.close()


    def read_restart_files(self):   
        filename = self.outpath + "minima.pickle.shelve.dat"
        self.minima_shelve = shelve.open(filename)
        if self.is_restart:
            # print('asdf', list(dict(self.minima_shelve).values()).sort())
            self.unique_minima_sorted = list(dict(self.minima_shelve).values())
            self.unique_minima_sorted.sort()
            self.nstructs = len(self.unique_minima_sorted)


    def addElement(self,struct: minimum.Minimum):
        if self.verbosity:
            t1 = time.time()
        index = self.get_element_index(struct=struct)
        if self.verbosity:
            t2 = time.time()
            finding_time = t2 - t1
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
            struct.write(self.noDublicatesFile, append=True)

            if index < self.output_n_lowest_minima:
                self._write_poslow(self.output_n_lowest_minima, self.minima_path)

            self.minima_shelve[str(label)] = struct1
        

        struct.write(self.allMinimaFile, append=True)
        if self.verbosity:
            t1 = time.time()
            db_time = t1 - t2
            if already_found:
                timingMessage = "Database search time: %.4f, adjusting minima shelve time: %.3f"%(finding_time, db_time)
            else:
                timingMessage = "Database search time: %.4f, adding minima shelve time: %.3f"%(finding_time, db_time)
            logging.info(timingMessage)

        return self.unique_minima_sorted[index].n_visit, self.unique_minima_sorted[index].label, True

    def addElementandConnectGraph(self, currentMinimum: minimum.Minimum, escapedMinimum: minimum.Minimum, trajectory, epot_max):
        n_vistit, label, _ = self.addElement(escapedMinimum)
        if self.write_graph_output:
            self.graph.addStructure(currentMinimum.label, escapedMinimum.label, trajectory, currentMinimum.e_pot, escapedMinimum.e_pot, epot_max)
        return n_vistit, label, True # last return determines if worker should continue. Since this class is not used with mpi True must be returned

    def get_element(self, index: int):
        return self.unique_minima_sorted[index].__copy__()

    def get_element_index(self, struct: minimum.Minimum):
        
        index = -1
        if self.compare_energies:
            i_start = bisect.bisect(self.unique_minima_sorted, struct)
            # if i_start is at the end of the list compare only to one lower
            if i_start >= len(self.unique_minima_sorted) and len(self.unique_minima_sorted) != 0:
                i_compare = i_start - 1
                print("DD1P:    ", i_compare)
                s = self.unique_minima_sorted[i_compare]
                energy_difference = struct.__compareto__(s)
                if energy_difference < self.energy_threshold:
                    index = i_compare
            # elif i_start is at the beginning only compare to the first structure
            elif i_start == 0 and len(self.unique_minima_sorted) != 0:
                i_compare = i_start
                print("DD2P:    ", i_compare)
                s = self.unique_minima_sorted[i_compare]
                energy_difference = struct.__compareto__(s)
                if energy_difference < self.energy_threshold:
                    index = i_compare
            # else compare one structure up and one down and the index itself
            elif len(self.unique_minima_sorted) != 0:
                i_compare = i_start - 1
                print("DD3P:    ", i_compare)
                s = self.unique_minima_sorted[i_compare]
                energy_difference_down = struct.__compareto__(s)
                if energy_difference_down < self.energy_threshold:
                    index = i_compare

                i_compare = i_start + 1
                print("DD4P:    ", i_compare)
                s = self.unique_minima_sorted[i_compare]
                energy_difference_up = struct.__compareto__(s)
                if energy_difference_up < self.energy_threshold:
                    index = i_compare

                i_compare = i_start
                print("DD5P:    ", i_compare)
                s = self.unique_minima_sorted[i_compare]
                energy_difference_itself = struct.__compareto__(s)
                if energy_difference_itself < self.energy_threshold:
                    index = i_compare
            else:
                index = -1
                    
        else:
            indices = self.get_index_energyrange(struct)
            for i_compare in indices:
                s = self.unique_minima_sorted[i_compare]      
                fp_dist = struct.fingerprint_distance(s)
                if fp_dist < self.minima_threshold:
                    index = i_compare
                    break

        return index

    def get_index_energyrange(self, struct):
        if self.nstructs == 0:
            return []

        i_start = bisect.bisect(self.unique_minima_sorted, struct)

        if i_start >= len(self.unique_minima_sorted):
            i_start -= 1

        indices = []
        break_backward = False
        break_forward = False
        i = 0

        while i < self.nstructs:

            # backward
            i_compare = i_start - i
            if i_compare >= 0:
                s = self.unique_minima_sorted[i_compare]
                energy_difference = struct.__compareto__(s)
                if energy_difference < self.energy_threshold:
                    indices.append(i_compare)
                else:
                    if i != 0:
                        break_backward = True
            else:
                break_backward = True

            i += 1
            # forward
            i_compare = i_start + i
            if i_compare < self.nstructs:
                s = self.unique_minima_sorted[i_compare]
                energy_difference = struct.__compareto__(s)
                if energy_difference < self.energy_threshold:
                    indices.append(i_compare)
                else:
                    break_forward = True
            else:
                break_forward = True
            # check if break
            if break_backward and break_forward:
                  break

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

import bisect
from ase.io import read, write
import pickle
import os

class Database():
    def __init__(self,energy_threshold, minima_threshold, is_restart = False, outpath='./'):
        self.unique_minima_sorted = []
        self.nstructs = 0

        self.energy_threshold = energy_threshold
        self.minima_threshold = minima_threshold
        self.is_restart = is_restart
        self.outpath = outpath


    def __enter__(self):
        if self.is_restart:
            filename = self.outpath + "data.pickle"
            listpickle = open(filename, "rb")
            self.unique_minima_sorted = pickle.load(listpickle)
            listpickle.close()
        return self

    def __exit__(self,exc_type, exc_value, exc_traceback):
        filename = self.outpath + "data.pickle"
        listpickle = open(filename, "wb")
        pickle.dump(self.unique_minima_sorted, listpickle)
        listpickle.close()
        return
        



    def addElement(self, struct,):
        struct = struct.__copy__()
        index = self.get_element(struct=struct)
        already_found = self.contains(index=index)

        if already_found:
            self.unique_minima_sorted[index].n_visit += 1
            n_visits = self.unique_minima_sorted[index].n_visit
        else:
            self.nstructs += 1
            label = self.nstructs
            struct.set_label(label)
            bisect.insort(self.unique_minima_sorted, struct)
            n_visits = 1

        return n_visits


    def get_element(self, struct):
        
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
            write(filename,s.atoms)
            i_poslow += 1
            if i_poslow-1 > n_poslow:
                break

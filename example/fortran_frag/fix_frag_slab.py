import numpy as np
from ase.io import read, write
from ase.data import covalent_radii, atomic_numbers
from ase.atoms import Atoms
from ase.atom import Atom
# from ase.neighborlist import NeighborList, build_neighbor_list, NewPrimitiveNeighborList
try:
    from numba import njit
except ImportError:
    # todo: raise warning
    print("numba not installed")
    def njit(f):
        return f

import sys
import time

def fix_frag_free(initial_structure: Atoms, debug = False, threshold = 1.3):

    # get all covalent radii from ASE
    rcovs = np.empty(len(initial_structure))
    i = 0
    for atom in initial_structure:
        rcovs[i] = covalent_radii[atomic_numbers[atom.symbol]]
        # print(i, atom.symbol, rcovs[i])
        i+=1

    nat = len(initial_structure)
    rxyz = initial_structure.positions
    t1 = time .time()
    fix_frag_numba(nat, rxyz, rcovs, threshold)
    t2 = time.time()
    print(t2 - t1)

# @njit()
def fix_frag_numba(nat, rxyz, rcovs, threshold = 1.3):
    cutoffs = threshold * rcovs


 
def main():
    import time

    filename = sys.argv[1]

    initial_structure = read(filename)

    pbcfail = False
    if not initial_structure.pbc[0]:
        pbcfail = True
    if not initial_structure.pbc[1]:
        pbcfail = True
    if initial_structure.pbc[2]:
        pbcfail = True
    if pbcfail:
        print("atoms.pbc must be [True, True, False]")
        quit()


    t1 = time.time()
    fix_frag_free(initial_structure)
    t2 = time.time()
    print("ela", t2 - t1)


    initial_structure = read(filename)
    t1 = time.time()
    fix_frag_free(initial_structure)
    t2 = time.time()
    print("ela (compiled)", t2 - t1)


    write('fixed.xyz', initial_structure)


if __name__ == "__main__":
    main()
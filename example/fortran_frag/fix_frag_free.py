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

@njit()
def fix_frag_numba(nat, rxyz, rcovs, threshold = 1.3):
    cutoffs = threshold * rcovs
    belong = np.empty(nat, dtype=np.bool)
    loggrow = False
    logexit = True
    vec= np.empty(3)

    icount = 0
    scale = 0.0

    belong[:] = False
    belong[0] = True
    while(True):
        # print(icount, '\n')
        # icount += 1
        # if icount > 10:
        #     return
        loggrow =False
        for iiat in range(nat):
            if belong[iiat]:
                for iat in range(iiat):
                    d2 = np.sum((rxyz[iat, :] - rxyz[iiat, :])**2)
                    if d2 <= (cutoffs[iiat] + cutoffs[iat]) * (cutoffs[iiat] + cutoffs[iat]):
                        if not belong[iat]:
                            loggrow = True
                        belong[iat] = True
                for iat in range(iiat + 1, nat):
                    d2 = np.sum((rxyz[iat, :] - rxyz[iiat, :])**2)
                    if d2 <= (cutoffs[iiat] + cutoffs[iat]) * (cutoffs[iiat] + cutoffs[iat]):
                        if not belong[iat]:
                            loggrow = True
                        belong[iat] = True

        # print("loggrow", loggrow)
        # print("updatebel", belong)
        if loggrow:
            continue
        # print("belong", belong)

        logexit = True
        for iat in range(nat):
            if not belong[iat]:
                logexit = False
        if logexit:
            return
        
        # All the atoms in the first fragment have been found
        # Now find shortest distance between atoms in the two fragments

        dmin = np.inf
        for iat in range(nat):
            if belong[iat]:
                for jat in range(nat):
                    if not belong[jat]:
                        d2 = np.sum((rxyz[iat, :] - rxyz [jat,:])**2)
                        if d2 <= dmin:
                            dmin = d2
                            iiat = iat
                            jjat =jat

        # print("dmin", np.sqrt(dmin))
        
        vec = rxyz[iiat, : ] - rxyz[jjat, :]
        veclength = np.linalg.norm(vec)
        scale = (veclength - rcovs[iiat] - rcovs[jjat]) / veclength
        vec = vec * scale
        for iat in range(nat):
            if not belong[iat]:
                rxyz[iat, :] = rxyz[iat, :] + vec
        
        # belong[:] = False
        # belong[0] = True

 
def main():
    import time

    filename = sys.argv[1]

    initial_structure = read(filename)

    if True in initial_structure.pbc: # check if one dimension is periodic
        print("this program cannot do pbc")
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
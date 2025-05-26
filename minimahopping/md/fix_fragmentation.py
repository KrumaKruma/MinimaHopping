
import numpy as np
from ase.data import covalent_radii, atomic_numbers
from ase.atoms import Atoms
import minimahopping.logging.logger as logging

# from ase.neighborlist import NeighborList, build_neighbor_list, NewPrimitiveNeighborList
try:
    from numba import njit
except ImportError:
    # todo: raise warning
    print("numba not installed")
    def njit(f):
        return f


def fix_frag_free(initial_structure: Atoms, threshold = 1.3):
    # extract arrays to run on numba

    # get all covalent radii from ASE
    rcovs = np.empty(len(initial_structure))
    i = 0
    for atom in initial_structure:
        rcovs[i] = covalent_radii[atomic_numbers[atom.symbol]]
        # print(i, atom.symbol, rcovs[i])
        i+=1

    nat = len(initial_structure)
    rxyz = initial_structure.positions
    did_fix = fix_frag_free_numba(nat, rxyz, rcovs, threshold)

    return did_fix

@njit()
def fix_frag_free_numba(nat, rxyz, rcovs, threshold = 1.3):
    cutoffs = threshold * rcovs
    belong = np.empty(nat, dtype=np.bool)
    loggrow = False
    logexit = True
    vec= np.empty(3)

    icount = 0
    scale = 0.0

    fix_stuff = False

    belong[:] = False
    belong[0] = True
    while(True):
        # print(icount, '\n')
        icount += 1
        if icount > 100:
            print("This should really not happen in fix fragmentation free. aborting")
            quit()
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
            return fix_stuff
        
        # All the atoms in the first fragment have been found
        # Now find shortest distance between atoms in the two fragments

        # at this point cluster is for sure fragmented:
        fix_stuff = True

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

def fix_frag_slab(initial_structure: Atoms, threshold = 1.3):

    # get all covalent radii from ASE
    rcovs = np.empty(len(initial_structure))
    i = 0
    for atom in initial_structure:
        rcovs[i] = covalent_radii[atomic_numbers[atom.symbol]]
        # print(i, atom.symbol, rcovs[i])
        i+=1

    nat = len(initial_structure)
    rxyz = initial_structure.positions
    return fix_frag_slab_numba(nat, rxyz, rcovs, threshold)



def fix_frag_slab_numba(nat, rxyz, rcovs, threshold=1.3, constrains_freeze = []):
    """
    Fix fragmented slabs by translating fragments in z-direction to bring them closer together.

    Parameters:
    - nat: int, number of atoms
    - rxyz: np.ndarray, shape (nat, 3), atomic coordinates
    - rcovs: np.ndarray, shape (nat,), covalent radii
    - threshold: float, distance factor for bonding
    - constrains_freeze: list, indices of atoms to freeze

    Returns:
    - None (in-place modification of rxyz)
    """
    cutoffs = threshold * rcovs
    max_cut = np.max(cutoffs)
    min_cut = max_cut / threshold
    did_fix = False

    # Sort atoms by z to check for gaps
    zvalues = rxyz[:, 2].copy()
    sorted_indices = np.argsort(zvalues)
    z_sorted = zvalues[sorted_indices]

    zdiff = np.diff(z_sorted)

    # If no big gap, return
    if np.all(zdiff <= max_cut):
        return

    # Fragmented: identify fragments
    fragments = np.zeros(nat, dtype=np.int32)
    current_frag = 0
    fragments[sorted_indices[0]] = current_frag

    for i in range(1, nat):
        idx_prev = sorted_indices[i - 1]
        idx_curr = sorted_indices[i]
        if abs(rxyz[idx_curr, 2] - rxyz[idx_prev, 2]) <= max_cut:
            fragments[idx_curr] = current_frag
        else:
            current_frag += 1
            fragments[idx_curr] = current_frag

    # list of fragments that contain at least one atom that should be frozen
    frozen_frags = []
    for i_constr in constrains_freeze:
        frozen_frags.append(fragments[sorted_indices[i_constr]])


    # if frozen_frags is not empty
    if len(frozen_frags) > 0:
        first_frag = frozen_frags[0]
        # covert frozen_frags to set
        frozen_frags = set(frozen_frags)
        # if more than one fragment frozen, abort.
        if len(frozen_frags) > 1:
            print("More than one fragment frozen. Aborting.")
            quit()

        if first_frag != 0:
            print("When using slabs, particles with the lowest z values can be frozen.")
            quit()

    # Move fragments above the first to close gaps
    for frag_id in range(1, current_frag + 1):
        did_fix = True
        frag_mask = fragments == frag_id
        prev_mask = fragments < frag_id

        # Minimum z of current fragment
        z_min = np.min(rxyz[frag_mask, 2])
        # Maximum z of previous fragment(s)
        z_max_prev = np.max(rxyz[prev_mask, 2])

        # Shift current fragment to close gap
        dz = z_max_prev + min_cut - z_min
        rxyz[frag_mask, 2] += dz
    
    return did_fix

if __name__ == "__main__":
    import time
    import sys
    import time
    from ase.io import read, write

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
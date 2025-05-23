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


def fix_frag_numba(nat, rxyz, rcovs, threshold=1.3, constrains_freeze = []):
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

    # covert frozen_frags to set
    frozen_frags = set(frozen_frags)
    # if more than one fragment frozen, abort.
    if len(frozen_frags) > 1:
        print("More than one fragment frozen. Aborting.")
        quit()

    if frozen_frags[0] != 0:
        print("When using slabs, particles with the lowest z values can be frozen.")
        quit()

    # Move fragments above the first to close gaps
    for frag_id in range(1, current_frag + 1):
        frag_mask = fragments == frag_id
        prev_mask = fragments < frag_id

        # Minimum z of current fragment
        z_min = np.min(rxyz[frag_mask, 2])
        # Maximum z of previous fragment(s)
        z_max_prev = np.max(rxyz[prev_mask, 2])

        # Shift current fragment to close gap
        dz = z_max_prev + min_cut - z_min
        rxyz[frag_mask, 2] += dz
        

    



 
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


    write('fixed_slab.xyz', initial_structure)


if __name__ == "__main__":
    main()
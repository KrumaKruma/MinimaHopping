from pymatgen.core import Composition
from ase.data import covalent_radii
from ase.neighborlist import NeighborList
from ase import Atoms
import numpy as np
from ase.symbols import symbols2numbers


def get_random_packed(composition: str | Composition, scale: int):
    if isinstance(composition, str):
        composition = Composition(composition)

    elements = sum([[str(el)] * int(composition.element_composition[el]) for el in composition], [])

    radii = covalent_radii[symbols2numbers(elements)] * scale
    atomic_volumes = 3 * (4/3) * np.pi * radii ** 3
    cell_vol = np.sum(atomic_volumes) 
        

    cell = np.eye(3) * cell_vol ** (1 / 3) 
    nat = len(elements)

    pos = np.random.rand(len(elements), 3) @ cell
    ats = Atoms(elements, cell=cell, pbc=True, positions=pos)
    skin = 0.0
    nl = NeighborList(radii, self_interaction=False, bothways=True, skin=skin)

    for it in range(500):
        nl.update(ats)
        pos = ats.get_positions()
        dx = np.zeros(pos.shape)
        dsum = 0
        for i in range(nat):
            indices, offsets = nl.get_neighbors(i)
            rs = pos[indices, :] + offsets @ cell - pos[i, :]
            # ds is the overlap
            ds = np.linalg.norm(rs, axis=1) - (radii[indices] + radii[i])
            if np.any(ds > 1.e-10):
                print('Assertion failed: ds <= 0')
                print(np.max(ds))
                quit()
            # sum overlaps 
            dsum += np.sum(ds) 
            ds -= skin
            # move atoms away from each other by overlap amount
            dx[i,:] = np.sum(rs / np.linalg.norm(rs, axis=1)[:, None] * ds[:, None], axis=0)
            
        # print(it, dsum, np.linalg.norm(dx))
        ats.set_positions(pos + dx)
        if dsum >= -1.e-5:
            break
        # print(i, ds)
    else: 
        raise RuntimeError('Cell packing not converged')

    ats.wrap()
    return ats


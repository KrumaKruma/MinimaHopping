import numpy as np

try:
    from numba import njit
except ImportError:
    # todo: raise warning
    def njit(f):
        return f

# todo: here we could use the ase neighorblist (probably faster)
@njit
def findNeighbors_jit(ats, els, rcut, lat=None):
    nat = ats.shape[0]
    # center atom is first in list
    neiats = [[ats[i, :]] for i in range(nat)]
    neiels = [[els[i]] for i in range(nat)]
    if lat is None:
        for i in range(nat):
            for j in range(i + 1, nat):
                d = np.linalg.norm(ats[i, :] - ats[j, :])
                if d < rcut:
                    neiats[i].append(ats[j, :])
                    neiats[j].append(ats[i, :])
                    neiels[i].append(els[j])
                    neiels[j].append(els[i])
    else:
        nc = ncells(lat, rcut)
        for i in range(nat):
            for j in range(nat):
                for ix in range(-nc[0], nc[0] + 1):
                    for iy in range(-nc[1], nc[1] + 1):
                        for iz in range(-nc[2], nc[2] + 1):
                            if i != j or ix != 0 or iy != 0 or iz != 0:
                                dlat = lat[0, :] * ix + lat[1, :] * iy + lat[2, :] * iz
                                d = np.linalg.norm(ats[i, :] - (ats[j, :] + dlat))
                                if d < rcut:
                                    neiats[i].append(ats[j, :] + dlat)
                                    neiels[i].append(els[j])

    #neiats = [np.array(x) for x in neiats]
    return neiats, neiels

# wrapper
def findNeighbors(ats, els, rcut, lat=None):
    neiats, neiels = findNeighbors_jit(ats, els, rcut, lat)
    neiats = [np.array(x) for x in neiats] # because numba does not like this line
    return neiats, neiels

@njit
def ncells(lat, rcut):
    n = []
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        c = np.cross(lat[j, :], lat[k, :])
        c = c / np.linalg.norm(c)
        cc = int(np.ceil(rcut / np.abs(np.dot(lat[i, :], c))))
        n.append(cc)
    return n

import numpy as np

try:
    from numba import njit
except ImportError:
    # todo: raise warning
    def njit(f):
        return f

# Vectorized version of the OM.
# This implementation changes the order of the orbitals (but returns the applied permutation)
# Numba does not make this much faster

@njit
def sphericalHarmonicsOverlap_vectorized(pa, pb, ra, rb, na, nb, ia, ib): # p->position, r->radius, n->spdf_idx, i->xyz_idx
    # coefficients for spherical harmonics
    # if you need f orbitals: copy the coefficients from here: http://openmopac.net/manual/real_spherical_harmonics.html
    # coefficients are in this functions for numba optimization of compile time constants
    coeffs = [
        [[(np.sqrt(1. / (4. * np.pi)), (0, 0, 0))]],                                                     # s
        [[(np.sqrt(3. / (4. * np.pi)), (1, 0, 0))],                                                      # px
         [(np.sqrt(3. / (4. * np.pi)), (0, 1, 0))],                                                      # py
         [(np.sqrt(3. / (4. * np.pi)), (0, 0, 1))]],                                                     # pz
        [[(-np.sqrt(5. / (16. * np.pi)), (2, 0, 0)), (-np.sqrt(5. / (16. * np.pi)), (0, 2, 0)),
          (2. * np.sqrt(5. / (16. * np.pi)), (0, 0, 2))],                                                # dz2
         [(np.sqrt(15. / (4. * np.pi)), (1, 0, 1))],                                                     # dxz
         [(np.sqrt(15. / (4. * np.pi)), (0, 1, 1))],                                                     # dyz
         [(np.sqrt(15. / (4. * np.pi)), (1, 1, 0))],                                                     # dxy
         [(np.sqrt(15. / (16. * np.pi)), (2, 0, 0)), (-np.sqrt(15. / (16. * np.pi)), (0, 2, 0))]]        # dx2-y2
    ]
    norb = ra.size
    o = np.zeros((norb,))
    for ci, lmni in coeffs[na][ia]:
        for cj, lmnj in coeffs[nb][ib]:
            # compute Gaussian exponent from rcov
            ai = 1. / (2. * ra ** 2)
            aj = 1. / (2. * rb ** 2)
            o += overlap_vectorized(
                ai, lmni, pa,
                aj, lmnj, pb) * ci * cj
    return o

def buildOverlapMatrix_vectorized(orbpos, orbrad, orbname, orbidx):
    norb = len(orbname)

    # sort the orbital parameters indo lists for each orbital type
    orbsort_pos = {}
    orbsort_rad = {}
    permsort = {}
    for i in range(norb):
        key = (orbname[i], orbidx[i])
        if not key in orbsort_pos:
            orbsort_pos[key] = [orbpos[i, :]]
            orbsort_rad[key] = [orbrad[i]]
            permsort[key] = [i]
        else:
            orbsort_pos[key].append(orbpos[i, :])
            orbsort_rad[key].append(orbrad[i])
            permsort[key].append(i)
    keys = list(orbsort_rad.keys())

    # make a list of all types and keep track of the permutation that we applied
    perm = []
    nkeys = [0]
    for k in keys:
        nkeys.append(nkeys[-1] + len(orbsort_pos[k]))
        orbsort_pos[k] = np.array(orbsort_pos[k])
        orbsort_rad[k] = np.array(orbsort_rad[k])
        perm.extend(permsort[k])

    O = np.zeros((norb, norb))
    for ia, ka in enumerate(keys):
        na = len(orbsort_rad[ka])
        for iib, kb in enumerate(keys[ia:]):   # only half is needed
            # calculate a part of the OM
            ib = iib + ia  # ia and ib are needed to put the part in the right place of O
            nb = len(orbsort_rad[kb])
            # make a list of all na*nb overlap integrals
            ii, jj = np.mgrid[:na, :nb]
            # get the correct locations and radii
            posa = orbsort_pos[ka][ii, :]
            posb = orbsort_pos[kb][jj, :]
            rada = orbsort_rad[ka][ii]
            radb = orbsort_rad[kb][jj]

            # this makes it slower
            # Maybe the whole array reordering is too expensive
            # Just don't use it now. Maybe it can be done better
            # It's a little unsafe any ways
            if False:  #ia == ib or isSymmetric(ka[0], kb[0], ka[1], kb[1]):  # the parts on the diagonal are symmetric (Other parts too. Need a clever way to check which ones.)
                oo = np.zeros((na, nb))
                tridx = np.triu_indices(na)  # indices for upper triangular part of matrix
                oo[tridx] = sphericalHarmonicsOverlap_vectorized(
                    posa[tridx],
                    posb[tridx],
                    rada[tridx],
                    radb[tridx],
                    ka[0],
                    kb[0],
                    ka[1],
                    kb[1])
                oo = oo + oo.transpose() - np.diag(np.diag(oo))  # make the full matrix from half
            else:
                # compute all matrix elements in one go and then reshape to the matrix form
                oo = sphericalHarmonicsOverlap_vectorized(
                    posa.reshape((na * nb, 3)),
                    posb.reshape((na * nb, 3)),
                    rada.reshape((na * nb,)),
                    radb.reshape((na * nb,)),
                    ka[0],
                    kb[0],
                    ka[1],
                    kb[1]).reshape((na, nb))
            O[nkeys[ia]:nkeys[ia+1], nkeys[ib]:nkeys[ib+1]] = oo
            if ia != ib:  # its a symmetric matrix
                O[nkeys[ib]:nkeys[ib+1], nkeys[ia]:nkeys[ia+1]] = oo.transpose()

    return O, np.array(perm)  # perm is the permutation that was applied to the orbitals

# is the overlap matrix of these orbitals symmetric?
def isSymmetric(namea, nameb, idxa, idxb):
    if namea == nameb:  # Careful: This makes assumptions (for ex. order of px and py is the same ...)
        return True
    return False

# taken from: https://github.com/jjgoings/McMurchie-Davidson/blob/master/mmd/integrals/reference.py
# todo: respect the license!
@njit
def E_vectorized(i, j, t, Qx, a, b):
    ''' Recursive definition of Hermite Gaussian coefficients.
        Returns a float.
        a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
        i,j: orbital angular momentum number on Gaussian 'a' and 'b'
        t: number nodes in Hermite (depends on type of integral,
           e.g. always zero for overlap integrals)
        Qx: distance between origins of Gaussian 'a' and 'b'
    '''

    if (t < 0) or (t > (i + j)):
        # out of bounds for t
        return np.zeros(a.shape)
    p = a + b
    q = a * b / p
    if i == j == t == 0:
        # base case
        return np.exp(-q * Qx * Qx)  # K_AB
    else:
        if i > j:
            # decrement index i
            return (1. / (2. * p)) * E_vectorized(i - 1, j, t - 1, Qx, a, b) - \
                   (q * Qx / a) * E_vectorized(i - 1, j, t, Qx, a, b) + \
                   (t + 1) * E_vectorized(i - 1, j, t + 1, Qx, a, b)
        else:
            # decrement index j
            return (1. / (2. * p)) * E_vectorized(i, j - 1, t - 1, Qx, a, b) + \
                   (q * Qx / b) * E_vectorized(i, j - 1, t, Qx, a, b) + \
                   (t + 1) * E_vectorized(i, j - 1, t + 1, Qx, a, b)

# taken from: https://github.com/jjgoings/McMurchie-Davidson/blob/master/mmd/integrals/reference.py
# todo: respect the license!

@njit
def overlap_vectorized(a, lmn1, A, b, lmn2, B):
    ''' Evaluates overlap integral between two Gaussians
        Returns a float.
        a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
              for Gaussian 'a'
        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
        A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
        B:    list containing origin of Gaussian 'b'
    '''
    l1, m1, n1 = lmn1  # shell angular momentum on Gaussian 'a'
    l2, m2, n2 = lmn2  # shell angular momentum on Gaussian 'b'
    S1 = E_vectorized(l1, l2, 0, A[:, 0] - B[:, 0], a, b)  # X
    S2 = E_vectorized(m1, m2, 0, A[:, 1] - B[:, 1], a, b)  # Y
    S3 = E_vectorized(n1, n2, 0, A[:, 2] - B[:, 2], a, b)  # Z
    return S1 * S2 * S3 * np.power(np.pi / (a + b), 1.5)

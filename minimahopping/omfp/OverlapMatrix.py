
import numpy as np

try:
    from numba import njit
except ImportError:
    # todo: raise warning
    def njit(f):
        return f


@njit
def sphericalHarmonicsOverlap(pa, pb, ra, rb, na, nb, ia, ib): # p->position, r->radius, n->spdf_idx, i->xyz_idx
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
    o = 0.
    for ci, lmni in coeffs[na][ia]:
        for cj, lmnj in coeffs[nb][ib]:
            # compute Gaussian exponent from rcov
            ai = 1. / (2. * ra ** 2)
            aj = 1. / (2. * rb ** 2)
            o += overlap(
                ai, lmni, pa,
                aj, lmnj, pb) * ci * cj
    return o

@njit
def buildOverlapMatrix(orbpos, orbrad, orbname, orbidx):
    norb = len(orbname)
    O = np.empty((norb, norb))
    for i in range(norb):
        for j in range(i, norb):
            O[i,j] = sphericalHarmonicsOverlap(orbpos[i,:], orbpos[j,:], orbrad[i], orbrad[j], orbname[i], orbname[j], orbidx[i], orbidx[j])
            O[j,i] = O[i,j]

    return O


# taken from: https://github.com/jjgoings/McMurchie-Davidson/blob/master/mmd/integrals/reference.py
# todo: respect the license!
@njit
def E(i, j, t, Qx, a, b):
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
        return 0.0
    p = a + b
    q = a * b / p
    if i == j == t == 0:
        # base case
        return np.exp(-q * Qx * Qx)  # K_AB
    else:
        if i > j:
            # decrement index i
            return (1. / (2. * p)) * E(i - 1, j, t - 1, Qx, a, b) - \
                   (q * Qx / a) * E(i - 1, j, t, Qx, a, b) + \
                   (t + 1) * E(i - 1, j, t + 1, Qx, a, b)
        else:
            # decrement index j
            return (1. / (2. * p)) * E(i, j - 1, t - 1, Qx, a, b) + \
                   (q * Qx / b) * E(i, j - 1, t, Qx, a, b) + \
                   (t + 1) * E(i, j - 1, t + 1, Qx, a, b)

# taken from: https://github.com/jjgoings/McMurchie-Davidson/blob/master/mmd/integrals/reference.py
# todo: respect the license!
@njit
def overlap(a, lmn1, A, b, lmn2, B):
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
    S1 = E(l1, l2, 0, A[0] - B[0], a, b)  # X
    S2 = E(m1, m2, 0, A[1] - B[1], a, b)  # Y
    S3 = E(n1, n2, 0, A[2] - B[2], a, b)  # Z
    return S1 * S2 * S3 * np.power(np.pi / (a + b), 1.5)
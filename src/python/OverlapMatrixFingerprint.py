import numpy as np
import warnings

class OverlapMatrixFingerprint:
    def __init__(self, lmn, rcut=-1, nex_cutoff=2, fplen=-1):
        self.lmn = lmn
        self.rcut = rcut
        self.nex_cutoff = nex_cutoff
        self.fplen = fplen

    # For compatibility with Stefans fortran version
    @staticmethod
    def stefansOMFP(s, p, width_cutoff, maxnatsphere=-1):
        lmn = {}
        for el in OverlapMatrixFingerprint.rcovs:
            iel = OverlapMatrixFingerprint.elementSymbolToNumber[el]
            #1: [(rH, (0, 0, 0)), (rH, (1, 0, 0)), (rH, (0, 1, 0)), (rH, (0, 0, 1))],
            lmn[iel] = []
            for ss in range(s):
                #cs(i)=sqrt(2.d0)**(i-1)
                lmn[iel].append((OverlapMatrixFingerprint.getRcov(iel) * np.sqrt(2.)**(float(ss) - ((float(s) + 1.)/2.)), 's'))
            for pp in range(p):
                lmn[iel].append((OverlapMatrixFingerprint.getRcov(iel) * np.sqrt(2.)**(float(pp) - ((float(p) + 1.)/2.)), 'p'))
        nex_cutoff = 2
        rcut = np.sqrt(2 * nex_cutoff) * width_cutoff
        return OverlapMatrixFingerprint(lmn, rcut=rcut, nex_cutoff=nex_cutoff, fplen=maxnatsphere*(s+3*p))

    # Counts how many orbitals are present for one element
    def getNOrbs(self, els):
        n = 0
        for el in els:
            for l in self.lmn[el]:
                if type(l[1]) is str:
                    if l[1] == 's':
                        n += 1
                    elif l[1] == 'p':
                        n += 3
                    elif l[1] == 'd':
                        n += 5
                else:
                    n += 1
        return n

    # makes a list of all orbitals for a given element
    # also entangles the two options for giving the orbitals
    def listOrbitals(self, el):
        orbs = []
        origin = np.zeros((3,))
        for l in self.lmn[el]:
            if type(l[1]) is str:
                if l[1] == 's':
                    orbs.append([(l[0], np.sqrt(1./(4. * np.pi)), (0, 0, 0))])
                elif l[1] == 'p':
                    orbs.append([
                        (l[0], np.sqrt(3. / (4. * np.pi)), (1, 0, 0)),
                    ])
                    orbs.append([
                        (l[0], np.sqrt(3. / (4. * np.pi)), (0, 1, 0)),
                    ])
                    orbs.append([
                        (l[0], np.sqrt(3. / (4. * np.pi)), (0, 0, 1)),
                    ])
                elif l[1] == 'd':
                    orbs.append([
                        (l[0], -np.sqrt(5. / (16. * np.pi)), (2, 0, 0)),
                        (l[0], -np.sqrt(5. / (16. * np.pi)), (0, 2, 0)),
                        (l[0], 2. * np.sqrt(5. / (16. * np.pi)), (0, 0, 2)),
                    ])
                    orbs.append([
                        (l[0], np.sqrt(15. / (4. * np.pi)), (1, 0, 1)),
                    ])
                    orbs.append([
                        (l[0], np.sqrt(15. / (4. * np.pi)), (0, 1, 1)),
                    ])
                    orbs.append([
                        (l[0], np.sqrt(15. / (4. * np.pi)), (1, 1, 0)),
                    ])
                    orbs.append([
                        (l[0], np.sqrt(15. / (16. * np.pi)), (2, 0, 0)),
                        (l[0], -np.sqrt(15. / (16. * np.pi)), (0, 2, 0)),
                    ])
                    # if you need f orbitals: copy the coefficients from here: http://openmopac.net/manual/real_spherical_harmonics.html
            else:
                (rci, oi) = l
                ai = 1. / (2 * rci ** 2)
                norm = 1.
                #norm = 1. / np.sqrt(self.overlap(
                #        ai, oi, origin,
                #        ai, oi, origin))
                orbs.append([(rci, norm, oi)])
        return orbs

    # Someone should clean everything up
    def overlapMatrixSpHar(self, ats, els, fcuts=None):
        origin = np.zeros((3,))
        nat = ats.shape[0]
        norbs = self.getNOrbs(els)

        # normalize, such that diagonal elements are 1.
        norms = np.zeros((norbs,))
        c = -1
        for i in range(nat):
            for ii, compsi in enumerate(self.listOrbitals(els[i])):
                c += 1
                # Since multiple Gauss functions per orbital, we need double loop to get the cross terms for normalization
                for rci, normi, oi in compsi:
                    ai = 1. / (2 * rci ** 2)
                    for rcj, normj, oj in compsi:
                        norms[c] += self.overlap(
                            ai, oi, origin,
                            ai, oj, origin) * (normi * normj)
                norms[c] = 1. / np.sqrt(norms[c]) #  diag = 1
                if fcuts is not None:
                    norms[c] *= fcuts[i] # cutoff function also goes into normalization

        # todo: only need to do half the matrix
        O = np.zeros((norbs, norbs))
        ci = -1
        # loop over atoms
        for i in range(nat):
            # for each atom, loop over orbitals
            for ii, compsi in enumerate(self.listOrbitals(els[i])):
                ci += 1
                cj = -1
                # now the same again
                for j in range(nat):
                    for jj, compsj in enumerate(self.listOrbitals(els[j])):
                        cj += 1
                        O[ci, cj] += self.orbOverlap(ats[i,:], ats[j,:], compsi, compsj) * norms[ci] * norms[cj]
        return O

    def orbOverlap(self, posA, posB, orbA, orbB):
        o = 0
        for rci, normi, oi in orbA:
            for rcj, normj, oj in orbB:
                ai = 1. / (2 * rci ** 2)
                aj = 1. / (2 * rcj ** 2)
                o += self.overlap(
                    ai, oi, posA,
                    aj, oj, posB) * normi * normj
        return o

    # old code
    def overlapMatrix(self, ats, els, fcuts=None):
        nat = ats.shape[0]
        norbs = self.getNOrbs(els)
        #for el in els:
        #    norbs += len(self.lmn[el])

        # normalize, such that diagonal elements are 1.
        norms = np.zeros((norbs,))
        c = -1
        for i in range(nat):
            for ii, (rci, oi) in enumerate(self.lmn[els[i]]):
                c += 1
                ai = 1. / (2 * rci ** 2)
                norms[c] = 1. / np.sqrt(self.overlap(
                    ai, oi, ats[i, :],
                    ai, oi, ats[i, :]))
                if fcuts is not None:
                    norms[c] *= fcuts[i]

        # todo: only need to do half the matrix
        O = np.zeros((norbs, norbs))
        ci = -1
        for i in range(nat):
            for ii, (rci, oi) in enumerate(self.lmn[els[i]]):
                ci += 1
                cj = -1
                for j in range(nat):
                    for jj, (rcj, oj) in enumerate(self.lmn[els[j]]):
                        cj += 1
                        ai = 1. / (2 * rci ** 2)
                        aj = 1. / (2 * rcj ** 2)
                        O[ci, cj] = self.overlap(
                            ai, oi, ats[i, :],
                            aj, oj, ats[j, :]) * norms[ci] * norms[cj]
        return O

    @staticmethod
    def diag(O):
        evals = np.linalg.eigvalsh(O)
        evals = np.sort(evals)[::-1]
        return evals

    def globalFingerprint(self, ats, els):
        O = self.overlapMatrixSpHar(ats, els)
        np.set_printoptions(threshold=np.inf, linewidth=np.inf)
        return self.diag(O)

    def fingerprint(self, ats, els, lat=None):
        nat = ats.shape[0]
        neiats, neiels = self.findNeighbors(ats, els, lat)
        fps = []
        for iat in range(nat):
            nneis = len(neiats[iat])
            fcuts = [self.fcut(self.rcut, np.linalg.norm(neiats[iat][0,:] - neiats[iat][i,:])) for i in range(nneis)]
            O = self.overlapMatrixSpHar(neiats[iat], neiels[iat], fcuts)
            fps.append(self.adjustFPlen(self.diag(O)))
        return fps

    def adjustFPlen(self, fp):
        if self.fplen < 0:
            return fp
        if fp.size < self.fplen:
            f = np.zeros((self.fplen,))
            f[:fp.size] = fp
            return f
        else:
            warning_msg = "Fingerprint of length {:d} is truncated to length {:d}".format(fp.shape[0], self.fplen)
            warnings.warn(warning_msg, FutureWarning)
            return fp[:self.fplen]

    @staticmethod
    def getRcov(el):
        return OverlapMatrixFingerprint.rcovs[OverlapMatrixFingerprint.elementSymbols[el]]

    @staticmethod
    def fcut(rcut, d, nex_cutoff=2):
        if d > rcut:
            return 0.
        wcut = rcut / np.sqrt(2. * nex_cutoff)
        fcut = 1. / (2. * nex_cutoff * wcut**2)
        return (1 - d**2 * fcut)**(nex_cutoff - 1) * (1. - d**2 * fcut)

    def findNeighbors(self, ats, els, lat=None):
        nat = ats.shape[0]
        # center atom is first in list
        neiats = [[ats[i,:]] for i in range(nat)]
        neiels = [[els[i]] for i in range(nat)]
        if lat is None:
            for i in range(nat):
                for j in range(i+1, nat):
                    d = np.linalg.norm(ats[i, :] - ats[j, :])
                    if d < self.rcut:
                        neiats[i].append(ats[j,:])
                        neiats[j].append(ats[i,:])
                        neiels[i].append(els[j])
                        neiels[j].append(els[i])
        else:
            nc = self.ncells(lat)
            for i in range(nat):
                for j in range(nat):
                    for ix in range(-nc[0], nc[0]+1):
                        for iy in range(-nc[1], nc[1]+1):
                            for iz in range(-nc[2], nc[2]+1):
                                if i!=j or ix != 0 or iy != 0 or iz != 0:
                                    dlat = lat[0,:] * ix + lat[1,:] * iy + lat[2,:] * iz
                                    d = np.linalg.norm(ats[i, :] - (ats[j, :] + dlat))
                                    if d < self.rcut:
                                        neiats[i].append(ats[j, :] + dlat)
                                        neiels[i].append(els[j])

        neiats = [np.array(x) for x in neiats]
        return neiats, neiels

    def ncells(self, lat):
        n = []
        for i in range(3):
            j = (i+1) % 3
            k = (i+2) % 3
            c = np.cross(lat[j, :], lat[k, :])
            c = c / np.linalg.norm(c)
            cc = int(np.ceil(self.rcut / np.abs(np.dot(lat[i,:], c))))
            n.append(cc)
        return n

    rcovs = {'H': 0.699198669167688, 'He': 0.6047123625234059, 'Li': 2.532233018066762, 'Be': 1.7007535195970789, 'B': 1.5495754289662274,
             'C': 1.4550891223219453, 'N': 1.4172945996642325, 'O': 1.3795000770065196, 'F': 1.3417055543488066, 'Ne': 1.3039110316910938,
             'Na': 2.9101782446438906, 'Mg': 2.4566439727513365, 'Al': 2.229876836805059, 'Si': 2.097596007503064, 'P': 2.003109700858782,
             'S': 1.9275206555433562, 'Cl': 1.870828871556787, 'Ar': 1.833034348899074, 'K': 3.703863220455861, 'Ca': 3.2881234712210192,
             'Sc': 2.721205631355326, 'Ti': 2.570027540724475, 'V': 2.3621576661070542, 'Cr': 2.399952188764767, 'Mn': 2.626719324711044,
             'Fe': 2.3621576661070542, 'Co': 2.3810549274359105, 'Ni': 2.2865686207916283, 'Cu': 2.6078220633821876, 'Zn': 2.4755412340801928,
             'Ga': 2.3810549274359105, 'Ge': 2.3054658821204845, 'As': 2.2487740981339153, 'Se': 2.192082314147346, 'Br': 2.154287791489633,
             'Kr': 2.078698746174208, 'Rb': 3.987322140388707, 'Sr': 3.628274175140435, 'Y': 3.061356335274742, 'Zr': 2.796794676670752,
             'Nb': 2.5889248020533313, 'Mo': 2.740102892684183, 'Tc': 2.9479727673016036, 'Ru': 2.3810549274359105, 'Rh': 2.5511302793956188,
             'Pd': 2.4755412340801928, 'Ag': 2.8912809833150344, 'Cd': 2.796794676670752, 'In': 2.721205631355326, 'Sn': 2.664513847368757,
             'Sb': 2.6078220633821876, 'Te': 2.5511302793956188, 'I': 2.5133357567379058, 'Xe': 2.4566439727513365, 'Cs': 4.251883798992697,
             'Ba': 3.741657743113574, 'Lu': 3.023561812617029, 'Hf': 2.834589199328465, 'Ta': 2.6078220633821876, 'W': 2.759000154013039,
             'Re': 3.004664551288173, 'Os': 2.4188494500936235, 'Ir': 2.5889248020533313, 'Pt': 2.4188494500936235, 'Au': 2.721205631355326,
             'Hg': 2.8156919379996084, 'Tl': 2.796794676670752, 'Pb': 2.7778974153418954, 'Bi': 2.759000154013039, 'Rn': 2.740102892684183, }

    elementSymbolToNumber = {'H': 1,'He': 2,'Li': 3,'Be': 4,'B': 5,'C': 6,'N': 7,'O': 8,'F': 9,'Ne': 10,
                            'Na': 11,'Mg': 12,'Al': 13,'Si': 14,'P': 15,'S': 16,'Cl': 17,'Ar': 18,'K': 19,'Ca': 20,
                            'Sc': 21,'Ti': 22,'V': 23,'Cr': 24,'Mn': 25,'Fe': 26,'Co': 27,'Ni': 28,'Cu': 29,'Zn': 30,
                            'Ga': 31,'Ge': 32,'As': 33,'Se': 34,'Br': 35,'Kr': 36,'Rb': 37,'Sr': 38,'Y': 39,'Zr': 40,
                            'Nb': 41,'Mo': 42,'Tc': 43,'Ru': 44,'Rh': 45,'Pd': 46,'Ag': 47,'Cd': 48,'In': 49,'Sn': 50,
                            'Sb': 51,'Te': 52,'I': 53,'Xe': 54,'Cs': 55,'Ba': 56,'La': 57,'Ce': 58,'Pr': 59,'Nd': 60,
                            'Pm': 61,'Sm': 62,'Eu': 63,'Gd': 64,'Tb': 65,'Dy': 66,'Ho': 67,'Er': 68,'Tm': 69,'Yb': 70,
                            'Lu': 71,'Hf': 72,'Ta': 73,'W': 74,'Re': 75,'Os': 76,'Ir': 77,'Pt': 78,'Au': 79,'Hg': 80,
                            'Tl': 81,'Pb': 82,'Bi': 83,'Po': 84,'At': 85,'Rn': 86,'Fr': 87,'Ra': 88,'Ac': 89,'Th': 90,
                            'Pa': 91,'U': 92,'Np': 93,'Pu': 94,'Am': 95,'Cm': 96,'Bk': 97,'Cf': 98,'Es': 99,'Fm': 100,
                            'Md': 101,'No': 102,'Lr': 103,'Rf': 104,'Db': 105,'Sg': 106,'Bh': 107,'Hs': 108,'Mt': 109,'Ds': 110,
                            'Rg': 111,'Cn': 112,'Nh': 113,'Fl': 114,'Mc': 115,'Lv': 116,'Ts': 117,'Og': 118}

    elementSymbols = [' ', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S',
          'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
          'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
          'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
          'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
          'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
          'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
          'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

    # taken from: https://github.com/jjgoings/McMurchie-Davidson/blob/master/mmd/integrals/reference.py
    # todo: respect the license!
    @staticmethod
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
        p = a + b
        q = a * b / p
        if (t < 0) or (t > (i + j)):
            # out of bounds for t
            return 0.0
        elif i == j == t == 0:
            # base case
            return np.exp(-q * Qx * Qx)  # K_AB
        elif j == 0:
            # decrement index i
            return (1 / (2 * p)) * OverlapMatrixFingerprint.E(i - 1, j, t - 1, Qx, a, b) - \
                   (q * Qx / a) * OverlapMatrixFingerprint.E(i - 1, j, t, Qx, a, b) + \
                   (t + 1) * OverlapMatrixFingerprint.E(i - 1, j, t + 1, Qx, a, b)
        else:
            # decrement index j
            return (1 / (2 * p)) * OverlapMatrixFingerprint.E(i, j - 1, t - 1, Qx, a, b) + \
                   (q * Qx / b) * OverlapMatrixFingerprint.E(i, j - 1, t, Qx, a, b) + \
                   (t + 1) * OverlapMatrixFingerprint.E(i, j - 1, t + 1, Qx, a, b)

    # taken from: https://github.com/jjgoings/McMurchie-Davidson/blob/master/mmd/integrals/reference.py
    # todo: respect the license!
    @staticmethod
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
        S1 = OverlapMatrixFingerprint.E(l1, l2, 0, A[0] - B[0], a, b)  # X
        S2 = OverlapMatrixFingerprint.E(m1, m2, 0, A[1] - B[1], a, b)  # Y
        S3 = OverlapMatrixFingerprint.E(n1, n2, 0, A[2] - B[2], a, b)  # Z
        return S1 * S2 * S3 * np.power(np.pi / (a + b), 1.5)

import numpy as np
import warnings
from minimahopping.omfp.OverlapMatrix import buildOverlapMatrix
from minimahopping.omfp.OverlapMatrixVectorized import buildOverlapMatrix_vectorized
from minimahopping.omfp.Neighborlist import findNeighbors
import minimahopping.mh.periodictable as periodictable 
import minimahopping.logging.logger as logging

class OverlapMatrixFingerprint:
    def __init__(self, lmn, rcut=-1, nex_cutoff=2):
        self.lmn = lmn # of the form {element_name: [(r_c, 's'), (r_c, 'p'), ...]}
        self.rcut = rcut
        self.nex_cutoff = nex_cutoff

    # For compatibility with Stefans fortran version
    @staticmethod
    def stefansOMFP(s, p, width_cutoff):
        lmn = {}
        for el in periodictable.get_rcov_dict():
            iel = periodictable.get_elementSymbolToNumber()[el]
            #1: [(rH, (0, 0, 0)), (rH, (1, 0, 0)), (rH, (0, 1, 0)), (rH, (0, 0, 1))],
            lmn[iel] = []
            for ss in range(s):
                #cs(i)=sqrt(2.d0)**(i-1)
                #lmn[iel].append((OverlapMatrixFingerprint.getRcov(iel) * np.sqrt(2.)**((ss + 1) - ((s + 1.) / 2.)), 's'))
                lmn[iel].append((periodictable.getRcov_n(iel) * np.sqrt(2.) ** (ss), 's'))
            for pp in range(p):
                #lmn[iel].append((OverlapMatrixFingerprint.getRcov(iel) * np.sqrt(2.)**((pp + 1) - ((p + 1.) / 2.)), 'p'))
                lmn[iel].append((periodictable.getRcov_n(iel) * np.sqrt(2.) ** ((pp)), 'p'))
        nex_cutoff = 2
        rcut = np.sqrt(2 * nex_cutoff) * width_cutoff
        return OverlapMatrixFingerprint(lmn, rcut=rcut, nex_cutoff=nex_cutoff)

    orbitalIndex = {'s': 0, 'p': 1, 'd': 2, 'f': 3}

    def getNOrbs(self, els):
        n = 0
        for el in els:
            for l in self.lmn[el]:
                n += self.orbitalIndex[l[1]] * 2 + 1
        return n

    def overlapMatrixSpHar(self, ats, els, fcuts=None):
        nat = ats.shape[0]
        norbs = self.getNOrbs(els)

        # assemble a list of all orbitals
        orbpos = np.zeros((norbs, 3))  # positions
        orbnames = np.zeros((norbs,), dtype=int)  # spdf
        orbidx = np.zeros((norbs, ), dtype=int)  # lmn
        orbrad = np.zeros((norbs,))  # r_cov
        orbcuts = np.ones((norbs,))  # cutoff function (default is 1 -> no cutoff)

        c = -1
        for iat in range(nat):  # all atoms
            for rc, spd in self.lmn[els[iat]]:  # s, p, d, ...
                oidx = self.orbitalIndex[spd]
                for ico in range(oidx * 2 + 1):  # s, px, py, pz, ...
                    c += 1
                    orbpos[c,:] = ats[iat,:]
                    orbnames[c] = oidx
                    orbidx[c] = ico
                    orbrad[c] = rc
                    if fcuts is not None:
                        orbcuts[c] = fcuts[iat]

        #O = buildOverlapMatrix(orbpos, orbrad, orbnames, orbidx)  # almost the whole runtime is spent here! numba is almost 20x faster
        # This one is even faster (even without numba)
        O, perm = buildOverlapMatrix_vectorized(orbpos, orbrad, orbnames, orbidx)  # but it permutes the oder of the orbitals
        orbcuts = orbcuts[perm]  # we just apply the same permutation to the cutoff too
        normalization = np.reciprocal(np.sqrt(np.diagonal(O))) * orbcuts # compute normalization, such that diagonal is 1 and apply cutoff function.
        O = O * np.outer(normalization, normalization)

        return O

    @staticmethod
    def diag(O):
        evals = np.linalg.eigvalsh(O)
        evals = np.sort(evals)[::-1]
        return evals

    def globalFingerprint(self, ats, els):
        O = self.overlapMatrixSpHar(ats, els)
        return self.diag(O)

    def fingerprint(self, ats, els, lat=None):
        nat = ats.shape[0]
        neiats, neiels = findNeighbors(ats, np.array(els), self.rcut, lat)
        fps = []
        maxlen = 0
        for iat in range(nat):
            #nneis = len(neiats[iat])
            #fcuts = [self.fcut(self.rcut, np.linalg.norm(neiats[iat][0,:] - neiats[iat][i,:])) for i in range(nneis)]
            fcuts = self.fcut(self.rcut, np.linalg.norm(neiats[iat][:,:] - neiats[iat][0,:], axis=1))  # the above in numpy notation
            O = self.overlapMatrixSpHar(neiats[iat], neiels[iat], fcuts)
            fps.append(self.diag(O))
            maxlen=max(maxlen, len(neiels[iat]))
        logging.logger.info('    Maximal number of atoms in the sphere for fingerprint: %d'%maxlen)
        return fps

    def adjustFPlen(fps, fplen):
        return_fps = []
        if fplen < 0:
            return fps
        for fp in fps:
            if fp.size <= fplen:
                f = np.zeros((fplen,))
                f[:fp.size] = fp
                return_fps.append(f)
            else:
                warning_msg = "Fingerprint of length {:d} is truncated to length {:d}".format(fp.shape[0], fplen)
                warnings.warn(warning_msg, FutureWarning)
                return_fps.append(fp[:fplen])
        return return_fps

    @staticmethod
    def getRcov(el):
        return OverlapMatrixFingerprint.rcovs[OverlapMatrixFingerprint.elementSymbols[el]]

    @staticmethod
    def fcut(rcut, d, nex_cutoff=2):
        wcut = rcut / np.sqrt(2. * nex_cutoff)
        fcut = 1. / (2. * nex_cutoff * wcut**2)
        fcut = (1 - d**2 * fcut)**(nex_cutoff - 1) * (1. - d**2 * fcut)
        return np.where(d <= rcut, fcut, 0.)


    # rcovs = {'H': 0.699198669167688, 'He': 0.6047123625234059, 'Li': 2.532233018066762, 'Be': 1.7007535195970789, 'B': 1.5495754289662274,
    #          'C': 1.4550891223219453, 'N': 1.4172945996642325, 'O': 1.3795000770065196, 'F': 1.3417055543488066, 'Ne': 1.3039110316910938,
    #          'Na': 2.9101782446438906, 'Mg': 2.4566439727513365, 'Al': 2.229876836805059, 'Si': 2.097596007503064, 'P': 2.003109700858782,
    #          'S': 1.9275206555433562, 'Cl': 1.870828871556787, 'Ar': 1.833034348899074, 'K': 3.703863220455861, 'Ca': 3.2881234712210192,
    #          'Sc': 2.721205631355326, 'Ti': 2.570027540724475, 'V': 2.3621576661070542, 'Cr': 2.399952188764767, 'Mn': 2.626719324711044,
    #          'Fe': 2.3621576661070542, 'Co': 2.3810549274359105, 'Ni': 2.2865686207916283, 'Cu': 2.6078220633821876, 'Zn': 2.4755412340801928,
    #          'Ga': 2.3810549274359105, 'Ge': 2.3054658821204845, 'As': 2.2487740981339153, 'Se': 2.192082314147346, 'Br': 2.154287791489633,
    #          'Kr': 2.078698746174208, 'Rb': 3.987322140388707, 'Sr': 3.628274175140435, 'Y': 3.061356335274742, 'Zr': 2.796794676670752,
    #          'Nb': 2.5889248020533313, 'Mo': 2.740102892684183, 'Tc': 2.9479727673016036, 'Ru': 2.3810549274359105, 'Rh': 2.5511302793956188,
    #          'Pd': 2.4755412340801928, 'Ag': 2.8912809833150344, 'Cd': 2.796794676670752, 'In': 2.721205631355326, 'Sn': 2.664513847368757,
    #          'Sb': 2.6078220633821876, 'Te': 2.5511302793956188, 'I': 2.5133357567379058, 'Xe': 2.4566439727513365, 'Cs': 4.251883798992697,
    #          'Ba': 3.741657743113574, 'Lu': 3.023561812617029, 'Hf': 2.834589199328465, 'Ta': 2.6078220633821876, 'W': 2.759000154013039,
    #          'Re': 3.004664551288173, 'Os': 2.4188494500936235, 'Ir': 2.5889248020533313, 'Pt': 2.4188494500936235, 'Au': 2.721205631355326,
    #          'Hg': 2.8156919379996084, 'Tl': 2.796794676670752, 'Pb': 2.7778974153418954, 'Bi': 2.759000154013039, 'Rn': 2.740102892684183, }

    # elementSymbolToNumber = {'H': 1,'He': 2,'Li': 3,'Be': 4,'B': 5,'C': 6,'N': 7,'O': 8,'F': 9,'Ne': 10,
    #                         'Na': 11,'Mg': 12,'Al': 13,'Si': 14,'P': 15,'S': 16,'Cl': 17,'Ar': 18,'K': 19,'Ca': 20,
    #                         'Sc': 21,'Ti': 22,'V': 23,'Cr': 24,'Mn': 25,'Fe': 26,'Co': 27,'Ni': 28,'Cu': 29,'Zn': 30,
    #                         'Ga': 31,'Ge': 32,'As': 33,'Se': 34,'Br': 35,'Kr': 36,'Rb': 37,'Sr': 38,'Y': 39,'Zr': 40,
    #                         'Nb': 41,'Mo': 42,'Tc': 43,'Ru': 44,'Rh': 45,'Pd': 46,'Ag': 47,'Cd': 48,'In': 49,'Sn': 50,
    #                         'Sb': 51,'Te': 52,'I': 53,'Xe': 54,'Cs': 55,'Ba': 56,'La': 57,'Ce': 58,'Pr': 59,'Nd': 60,
    #                         'Pm': 61,'Sm': 62,'Eu': 63,'Gd': 64,'Tb': 65,'Dy': 66,'Ho': 67,'Er': 68,'Tm': 69,'Yb': 70,
    #                         'Lu': 71,'Hf': 72,'Ta': 73,'W': 74,'Re': 75,'Os': 76,'Ir': 77,'Pt': 78,'Au': 79,'Hg': 80,
    #                         'Tl': 81,'Pb': 82,'Bi': 83,'Po': 84,'At': 85,'Rn': 86,'Fr': 87,'Ra': 88,'Ac': 89,'Th': 90,
    #                         'Pa': 91,'U': 92,'Np': 93,'Pu': 94,'Am': 95,'Cm': 96,'Bk': 97,'Cf': 98,'Es': 99,'Fm': 100,
    #                         'Md': 101,'No': 102,'Lr': 103,'Rf': 104,'Db': 105,'Sg': 106,'Bh': 107,'Hs': 108,'Mt': 109,'Ds': 110,
    #                         'Rg': 111,'Cn': 112,'Nh': 113,'Fl': 114,'Mc': 115,'Lv': 116,'Ts': 117,'Og': 118}

    # elementSymbols = [' ', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S',
    #       'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
    #       'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    #       'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    #       'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    #       'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
    #       'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
    #       'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']


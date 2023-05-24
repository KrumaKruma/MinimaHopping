import numpy as np
from scipy import optimize
from copy import deepcopy
from minimahopping.omfp.OverlapMatrixFingerprint import OverlapMatrixFingerprint as OMFP
from ase.io import write
import logging
import time
try:
    from numba import njit
except ImportError:
    # todo: raise warning
    def njit(f):
        return f

@njit
def _costmatrix(desc1, desc2):
    """
    Cost matrix of the local fingerprints for the hungarian algorithm
    """
    # assert desc1.shape[0] == desc2.shape[0], "descriptor has not the same length"

    costmat = np.zeros((desc1.shape[0], desc2.shape[0]))

    for i, vec1 in enumerate(desc1):
        for j, vec2 in enumerate(desc2):
            costmat[i, j] = np.linalg.norm(vec1 - vec2)
    return costmat

class Minimum():
    """ 
    Minimum class for managing the database of the minima hopping. 
    """
    def __init__(self, atoms, epot, s, p, width_cutoff, T, ediff, n_visit = None, label = None, exclude=[], fingerprint = None):
        self.atoms = atoms.copy()
        self.atoms.wrap()
        self.e_pot = epot
        if fingerprint is None:
            self.fp = epot #use epot as criteria for different structures
            # self.fp = self._get_OMFP(s=s, p=p, width_cutoff=width_cutoff, exclude=exclude)
        else:
            self.fp = fingerprint.copy()
        self.maxNatInEnv = 0
        # for f in self.fp:
            # self.maxNatInEnv = max(self.maxNatInEnv, f.size)
        self.maxNatInEnv = 1 #max(self.maxNatInEnv, f.size)
        self.temperature = T
        self.ediff = ediff
        self.n_visit = n_visit
        self.label = label
        self.s = s
        self.p = p
        self.width_cutoff = width_cutoff
        self.exclude = exclude

    def set_label(self, label):
        self.label = label

    def __lt__(self, other):
        return self.e_pot < other.e_pot

    def __gt__(self, other):
        return self.e_pot > other.e_pot

    def __copy__(self,):
        return Minimum(self.atoms.copy(), self.e_pot, self.s, self.p, self.width_cutoff,
            self.temperature, self.ediff, self.n_visit, self.label, self.exclude, fingerprint= self.fp)

    def __compareto__(self, other):
        return abs(self.e_pot - other.e_pot)

    def fingerprint_distance(self, other):
        
        #Use difference in epot as distinguishing structures:
        return abs(self.e_pot - other.e_pot)
        
        
        """
         Calcualtes the fingerprint distance of 2 structures with local environment descriptors using the hungarian algorithm
         if a local environment descriptor is used. Else the distance is calculated using l2-norm.
         """

        n1 = len(self.fp)
        n2 = len(other.fp)

        assert n1 == n2, "Number of particles for vector 1 is and vector 2 is different"
        # assert n_dim1 < 3, "Dimension of vector 1 is larger that 2"
        # assert n_dim2 < 3, "Dimension of vector 2 is larger that 2"

        # if n_dim1 == 1 and n_dim2 == 1:
        #     fp_dist = np.linalg.norm(self.fp - other.fp) / len(self.atoms)
        # else:
        maxNatInEnv = max(self.maxNatInEnv, other.maxNatInEnv)
        fp1 = np.array(OMFP.adjustFPlen(self.fp, maxNatInEnv))
        fp2 = np.array(OMFP.adjustFPlen(other.fp, maxNatInEnv))

        costmat = _costmatrix(fp1, fp2)
        ans_pos = optimize.linear_sum_assignment(costmat)
        # use this formula for euclidian fingerprint distance
        # fp_dist = np.linalg.norm( self.fp[ans_pos[0], :] - other.fp[ans_pos[1], :]) / len(self.atoms)
        fp_dist = np.max( np.abs(fp1[ans_pos[0], :] - fp2[ans_pos[1], :]) )

        return fp_dist

    def write(self, filename, append = False, info_dict: dict = {}):
        temp_atoms = self.atoms.copy()
        temp_atoms.info = {}
        temp_atoms.set_momenta(None)
        temp_atoms.info['energy'] = self.e_pot
        temp_atoms.info['label'] = self.label
        temp_atoms.info['n_visit'] = self.n_visit
        temp_atoms.info.update(info_dict)
        write(filename, temp_atoms, append=append, parallel=False)


    def _get_OMFP(self,s=1, p=0, width_cutoff=1.5, exclude=[]):
            """
            Calculation of the Overlapmatrix fingerprint. For peridoic systems a local environment fingerprint is calculated
            and a hungarian algorithm has to be used for the fingerprint distance. For non-periodic systems a global fingerprint
            is calculated and a simple l2-norm is sufficient for as a distance measure.

            If you use that function please reference:

            @article{sadeghi2013metrics,
            title={Metrics for measuring distances in configuration spaces},
            author={Sadeghi, Ali and Ghasemi, S Alireza and Schaefer, Bastian and Mohr, Stephan and Lill, Markus A and Goedecker, Stefan},
            journal={The Journal of chemical physics},
            volume={139},
            number={18},
            pages={184118},
            year={2013},
            publisher={American Institute of Physics}
            }

            and

            @article{zhu2016fingerprint,
            title={A fingerprint based metric for measuring similarities of crystalline structures},
            author={Zhu, Li and Amsler, Maximilian and Fuhrer, Tobias and Schaefer, Bastian and Faraji, Somayeh and Rostami, Samare and Ghasemi, S Alireza and Sadeghi, Ali and Grauzinyte, Migle and Wolverton, Chris and others},
            journal={The Journal of chemical physics},
            volume={144},
            number={3},
            pages={034203},
            year={2016},
            publisher={AIP Publishing LLC}
            }
            """

            _pbc = list(set(self.atoms.pbc))
            assert len(_pbc) == 1, "mixed boundary conditions"
            _ang2bohr = 1.8897161646320724

            _symbols = self.atoms.get_chemical_symbols()
            _positions = self.atoms.get_positions()
            _elements = self.atoms.get_atomic_numbers()
            _selected_postions = []
            _selected_elem = []

            for symb,elem, pos in zip(_symbols, _elements,_positions):
                if symb not in exclude:
                    _selected_postions.append(pos)
                    _selected_elem.append(elem)
            _selected_postions = np.array(_selected_postions)


            _selected_positions = _selected_postions*_ang2bohr
            _omfpCalculator = OMFP.stefansOMFP(s=s, p=p, width_cutoff=width_cutoff)
            if True in _pbc:
                _lattice = self.atoms.get_cell()*_ang2bohr
                _omfp = _omfpCalculator.fingerprint(_selected_positions, _selected_elem, lat=_lattice)
            else:
                _omfp = _omfpCalculator.fingerprint(_selected_positions, _selected_elem)
            return _omfp


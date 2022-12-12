import numpy as np
from scipy import optimize
from copy import deepcopy
from minimahopping.omfp.OverlapMatrixFingerprint import OverlapMatrixFingerprint as OMFP
from ase.io import write

class Minimum():
    """ 
    Minimum class for managing the database of the minima hopping. 
    """
    def __init__(self, atoms, epot, s, p, width_cutoff, maxnatsphere, T, ediff, n_visit = None, label = None, exclude=[]):
        self.atoms = deepcopy(atoms)
        self.e_pot = epot
        self.fp = self._get_OMFP(s=s, p=p, width_cutoff=width_cutoff, maxnatsphere=maxnatsphere, exclude=exclude)
        self.temperature = T
        self.ediff = ediff
        self.n_visit = n_visit
        self.label = label
        self.s = s
        self.p = p
        self.width_cutoff = width_cutoff
        self.maxnatsphere = maxnatsphere
        self.exclude = exclude

    def set_label(self, label):
        self.label = label

    def __lt__(self, other):
        return self.e_pot < other.e_pot

    def __gt__(self, other):
        return self.e_pot > other.e_pot

    def __copy__(self,):
        return Minimum(self.atoms.copy(), self.e_pot, self.s, self.p, self.width_cutoff, self.maxnatsphere, self.temperature, self.ediff, self.n_visit, self.label, self.exclude)

    def __deepcopy__(self,):
        return Minimum(self.atoms, self.e_pot ,self.s, self.p, self.width_cutoff, self.maxnatsphere, self.temperature, self.ediff, self.n_visit, self.label, self.exclude)

    def __compareto__(self, other):
        return abs(self.e_pot - other.e_pot)

    def __equals__(self, other):
        """
         Calcualtes the fingerprint distance of 2 structures with local environment descriptors using the hungarian algorithm
         if a local environment descriptor is used. Else the distance is calculated using l2-norm.
         """

        n_dim1 = len(self.fp.shape)
        n_dim2 = len(other.fp.shape)

        assert n_dim1 == n_dim2, "Dimension of vector 1 is and vector 2 is different"
        assert n_dim1 < 3, "Dimension of vector 1 is larger that 2"
        assert n_dim2 < 3, "Dimension of vector 2 is larger that 2"

        if n_dim1 == 1 and n_dim2 == 1:
            fp_dist = np.linalg.norm(self.fp - other.fp)
        else:
            costmat = self._costmatrix(self.fp, other.fp)
            ans_pos = optimize.linear_sum_assignment(costmat)
            fp_dist = 0.
            for index1, index2 in zip(ans_pos[0], ans_pos[1]):
                fp_dist += np.dot((self.fp[index1, :] - other.fp[index2, :]), (self.fp[index1, :] - other.fp[index2, :]))
            fp_dist = np.sqrt(fp_dist)

        fp_dist /= len(self.atoms)
        return fp_dist

    def write(self, filename :str, append = False, info_dict: dict = {}):
        temp_atoms = self.atoms.copy()
        temp_atoms.info = {}
        temp_atoms.set_momenta(None)
        temp_atoms.info['energy'] = self.e_pot
        temp_atoms.info['label'] = self.label
        temp_atoms.info = temp_atoms.info | info_dict
        write(filename, temp_atoms, append=append, parallel=False)


    def _costmatrix(self, desc1, desc2):
        """
        Cost matrix of the local fingerprints for the hungarian algorithm
        """
        assert desc1.shape[0] == desc2.shape[0], "descriptor has not the same length"

        costmat = np.zeros((desc1.shape[0], desc2.shape[0]))

        for i, vec1 in enumerate(desc1):
            for j, vec2 in enumerate(desc2):
                costmat[i, j] = np.linalg.norm(vec1 - vec2)

        return costmat


    def _get_OMFP(self,s=1, p=0, width_cutoff=1.5, maxnatsphere=100, exclude=[]):
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


            if True in _pbc:
                _selected_positions = _selected_postions*_ang2bohr
                _lattice = self.atoms.get_cell()*_ang2bohr
                _omfpCalculator = OMFP.stefansOMFP(s=s, p=p, width_cutoff=width_cutoff, maxnatsphere=maxnatsphere)
                _omfp = _omfpCalculator.fingerprint(_selected_positions, _selected_elem, lat=_lattice)
                _omfp = np.array(_omfp)

            else:
                _selected_positions = _selected_postions*_ang2bohr
                _elements = self.atoms.get_atomic_numbers()
                _width_cutoff = 1000000
                _maxnatsphere = len(self.atoms)
                _omfpCalculator = OMFP.stefansOMFP(s=s, p=p, width_cutoff=_width_cutoff, maxnatsphere=_maxnatsphere)
                _omfp = _omfpCalculator.globalFingerprint(_selected_positions, _selected_elem)
                _omfp = np.array(_omfp)

            return _omfp


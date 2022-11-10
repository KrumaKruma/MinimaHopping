import numpy as np
import scipy
from ase.io import read, write
from ase.calculators.lj import LennardJones
from ase.calculators.eam import EAM
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from md import MD
from optim import Opt
from OverlapMatrixFingerprint import OverlapMatrixFingerprint as OMFP
from cell_atom import Cell_atom
from copy import deepcopy
from bazant_calc import BazantCalculator




class adjust_fp():
    def __init__(self, atoms, fmax,iterations=100, temperature=500, dt=0.1, md_min=1, s=1, p=1, width_cutoff=3.5, exclude=[]):
        self._atoms = deepcopy(atoms)
        self.n_steps = iterations
        self._temperature = temperature
        self._dt = dt
        self._mdmin = md_min
        self._fmax = fmax
        self._s = s
        self._p = p
        self._width_cutoff = width_cutoff
        self._verbose = True
        self._exclude = exclude


        self.md_structures = []
        self.opt_structures = []
        self.fps = []

    def run(self,):
        self._run_md()
        self._opt()
        self._fingerprints()
        self._get_distances()
        return self._max_dist, self._mean_dist, self._std_dist

    def _run_md(self,):
        for i in range(self.n_steps):
            MaxwellBoltzmannDistribution(self._atoms, temperature_K=self._temperature)
            if True in self._atoms.pbc:
                self._vcsmd()
                self.md_structures.append(deepcopy(self._atoms))
            else:
                self._md()
                self.md_structures.append(deepcopy(self._atoms))


    def _vcsmd(self,):
        _mass = .75 * np.sum(self._atoms.get_masses()) / 10.
        self._cell_atoms = Cell_atom(mass=_mass, positions=self._atoms.get_cell())
        self._cell_atoms.set_velocities_boltzmann(temperature=self._temperature)
        md = MD(atoms=self._atoms,outpath='./', cell_atoms=self._cell_atoms, dt=self._dt, n_max=self._mdmin, verbose=self._verbose)
        _positions, _cell, self._dtt = md.run()
        self._atoms.set_positions(_positions)
        self._atoms.set_cell(_cell)


    def _md(self,):
        md = MD(atoms=self._atoms, outpath='./',cell_atoms=None, dt=self._dt, n_max=self._mdmin, verbose=self._verbose)
        _positions = md.run()
        self._atoms.set_positions(_positions)


    def _opt(self,):
        if True in self._atoms.pbc:
            self._vcs_opt()
        else:
            self._c_opt()


    def _vcs_opt(self,):
        for _atom in self.md_structures:
            opt = Opt(atoms=_atom, outpath='./', max_froce_threshold=self._fmax, verbose=self._verbose)
            _positions, _lattice, self._noise = opt.run()
            _atom.set_positions(_positions)
            _atom.set_cell(_lattice)
            self.opt_structures.append(deepcopy(_atom))

    def _c_opt(self,):
        for _atom in self.md_structures:
            opt = Opt(atoms=_atom,outpath='./', max_froce_threshold=self._fmax, verbose=self._verbose)
            _positions, self._noise = opt.run()
            _atom.set_positions(_positions)
            self.opt_structures.append(deepcopy(_atom))


    def _fingerprints(self,):
        for _atom in self.opt_structures:
            fp = self._get_OMFP(_atom, s=self._s, p=self._p, width_cutoff=self._width_cutoff, exclude=self._exclude)
            self.fps.append(fp)


    def _get_distances(self,):
        _distances = []
        for i,fp1 in enumerate(self.fps):
            for j, fp2 in enumerate(self.fps[i+1:]):
                fp_dist = self.fp_distance(fp1, fp2)
                _distances.append(fp_dist)

        self._max_dist = max(_distances)
        self._mean_dist = np.mean(_distances)
        self._std_dist = np.std(_distances)

    def _get_OMFP(self, _atoms, s=1, p=1, width_cutoff=1.5, maxnatsphere=100, exclude=[]):
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

        Input:
            s: int
                number of s orbitals for which the fingerprint is calculated
            p: int
                number of p orbitals for which the fingerprint is calculated
            width_cutoff: float
                cutoff for searching neighbouring atoms
            maxnatsphere:
                maximum of the neighboring atoms which can be in the sphere
        Return:
            omfp: np array
                numpy array which contains the fingerprint
        """

        _pbc = list(set(self._atoms.pbc))
        assert len(_pbc) == 1, "mixed boundary conditions"
        _ang2bohr = 1.8897161646320724

        _symbols = _atoms.get_chemical_symbols()
        _positions = _atoms.get_positions()
        _elements = _atoms.get_atomic_numbers()
        _selected_postions = []
        _selected_elem = []

        for symb, elem, pos in zip(_symbols, _elements, _positions):
            if symb not in exclude:
                _selected_postions.append(pos)
                _selected_elem.append(elem)
        _selected_postions = np.array(_selected_postions)

        if True in _pbc:
            _selected_positions = _selected_postions * _ang2bohr
            _lattice = _atoms.get_cell() * _ang2bohr
            _omfpCalculator = OMFP.stefansOMFP(s=s, p=p, width_cutoff=width_cutoff, maxnatsphere=maxnatsphere)
            _omfp = _omfpCalculator.fingerprint(_selected_positions, _selected_elem, lat=_lattice)
            _omfp = np.array(_omfp)

        else:
            _selected_positions = _selected_postions #* _ang2bohr
            _elements = _atoms.get_atomic_numbers()
            _width_cutoff = 1000000
            _maxnatsphere = len(_atoms)
            _omfpCalculator = OMFP.stefansOMFP(s=s, p=p, width_cutoff=_width_cutoff, maxnatsphere=_maxnatsphere)
            _omfp = _omfpCalculator.globalFingerprint(_selected_positions, _selected_elem)
            _omfp = np.array(_omfp)

        return _omfp

    def fp_distance(self, desc1, desc2):
        """
        Calcualtes the fingerprint distance of 2 structures with local environment descriptors using the hungarian algorithm
        if a local environment descriptor is used. Else the distance is calculated using l2-norm.
        desc1: np array
            numpy array containing local environments of structure 1
        desc2: np array
            numpy array containing local environments of structure 2
        Return:
            Global fingerprint distance between structure 1 and structure 2
        """

        n_dim1 = len(desc1.shape)
        n_dim2 = len(desc2.shape)

        assert n_dim1 == n_dim2, "Dimension of vector 1 is and vector 2 is different"
        assert n_dim1 < 3, "Dimension of vector 1 is larger that 2"
        assert n_dim2 < 3, "Dimension of vector 2 is larger that 2"

        if n_dim1 == 1 and n_dim2 == 1:
            fp_dist = np.linalg.norm(desc1 - desc2)
        else:
            costmat = self._costmatrix(desc1, desc2)
            ans_pos = scipy.optimize.linear_sum_assignment(costmat)
            fp_dist = 0.
            for index1, index2 in zip(ans_pos[0], ans_pos[1]):
                fp_dist += np.dot((desc1[index1, :] - desc2[index2, :]), (desc1[index1, :] - desc2[index2, :]))
            fp_dist = np.sqrt(fp_dist)

        return fp_dist

    def _costmatrix(self, desc1, desc2):
        """
        Cost matrix of the local fingerprints for the hungarian algorithm
        desc1: np array
            numpy array containing local fingerprints of structure 1
        desc2: np array
            numpy array containing local fingerprints of structure 2
        Return:
            cost matrix of with the distances of the local fingerprints
        """
        assert desc1.shape[0] == desc2.shape[0], "descriptor has not the same length"

        costmat = np.zeros((desc1.shape[0], desc2.shape[0]))

        for i, vec1 in enumerate(desc1):
            for j, vec2 in enumerate(desc2):
                costmat[i, j] = np.linalg.norm(vec1 - vec2)

        return costmat


def main():
    filename = "../../data/SiC_in.extxyz"
    atoms = read(filename)
    # calculator = LennardJones()
    # calculator.parameters.epsilon = 1.0
    # calculator.parameters.sigma = 1.0
    # calculator.parameters.rc = 6.0
    #calculator = EAM(potential="Na_v2.eam.fs")
    calculator = BazantCalculator()
    atoms.calc = calculator


    fnrm =  0.000005
    adjust = adjust_fp(atoms, fnrm,)
    fp_max, fp_mean, fp_std = adjust.run()
    msg = 'Maximal fingerprint distance between the same local minima:\n' + str(fp_max)
    msg += '\n Mean fingerprint distance between the same local minima:\n' + str(fp_mean)
    msg += '\n Standard deviaton of the fingerprint distances:\n' + str(fp_std)
    print(msg)




if __name__  == '__main__':
    main()




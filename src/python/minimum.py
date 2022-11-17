import numpy as np
from scipy import optimize

class Minimum():
    def __init__(self, atoms, epot, n_visit, fingerprint, T, ediff, acc_rej, label):
        self.atoms = atoms.copy()
        self.e_pot = epot
        self.fp = fingerprint
        self.temperature = T
        self.ediff = ediff
        self.acc_rej = acc_rej
        self.n_visit = n_visit
        self.label = label


    def set_label(self, label):
        self.label = label

    def __lt__(self, other):
        return self.e_pot < other.e_pot

    def __gt__(self, other):
        return self.e_pot > other.e_pot


    def __compareto__(self, other):
        return abs(self.e_pot - other.e_pot)

    def __equals__(self, other):
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




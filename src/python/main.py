import numpy as np
from ase.io import read, write
from ase.calculators.lj import LennardJones
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import Atoms
from ase import units
#from nequip.ase import NequIPCalculator
from ase.calculators.espresso import Espresso
from flare.utils.parameter_helper import ParameterHelper
from flare.gp import GaussianProcess
from ase import units
import time
import scipy
#import bazant
import periodic_sqnm
import free_or_fixed_cell_sqnm
import warnings
from OverlapMatrixFingerprint import OverlapMatrixFingerprint as OMFP
#from dscribe.descriptors import SOAP
from ase_flare import mh_otf


"""
MH Software written by Marco Krummenacher (marco.krummenacher@unibas.ch)
Parts of the software were originally developped (some in Fortran) from other people:
  -- VCSMD: Martin Sommer-Joergenson
  -- VCS Softening: Hannes Huber
  -- VCS optimizer: Moritz Gubler
  -- OMFP in python: Jonas Finkler
"""


class cell_atom:
    def __init__(self, positions, mass=None, velocities=None):
        self.positions = positions
        self.masses = np.array([mass, mass, mass])
        self.velocities = velocities

    def set_velocities_boltzmann(self, temperature):
        """
        Set the velocity of the cell atoms for the MD part accorting to a temperature and the boltzmann distribution
        Input:
            temperature: float
                Temperature of the cell atoms
        Return:
            velocities of the cell atoms
        """
        xi = np.random.standard_normal((len(self.masses), 3))
        temp = units.kB * temperature
        self.velocities = xi * np.sqrt(self.masses * temp)[:, np.newaxis]
        self.velocities /= self.masses

        return None



def costmatrix(desc1, desc2):
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
            costmat[i,j] = np.linalg.norm(vec1-vec2)

    return costmat


# def get_SOAP(atoms, rcut = 5, nmax = 8, lmax = 6, crossover=True, sigma=1.):
#     """
#     atoms: ASE atoms object
#         ASE atoms object containing the atomic structure
#     rcut: float
#         cut off radius of the SOAP fingerprint
#     nmax: int
#         number of radial basis functions
#     lmax: int
#         maximum degree of spherical harmonics
#     crossover: bool
#         Determines if crossover of atomic types should be included in the power spectrum. If enabled, the power
#         spectrum is calculated over all unique species combinations Z and Z'. If disabled, the power spectrum does not
#         contain cross-species information and is only run over each unique species Z. Turned on by default to
#         correspond to the original definition.
#     sigma: float
#         standard deviation of the gaussian basis function
#     Return: np array
#         numpy array containing local environment SOAP descriptor
#
#     Ref:
#         "Willatt, M., Musil, F., & Ceriotti, M. (2018).
#         Feature optimization for atomistic machine learning
#         yields a data-driven construction of the periodic
#         table of the elements.  Phys. Chem. Chem. Phys., 20,
#         29661-29668.
#
#         "Caro, M. (2019). Optimizing many-body atomic
#         descriptors for enhanced computational performance of
#         machine learning based interatomic potentials.  Phys.
#         Rev. B, 100, 024112."
#
#     """
#
#
#     element_list = list(set(atoms.get_chemical_symbols()))
#
#     pbc = atoms.pbc
#     pbc_list = list(set(pbc))
#     assert len(pbc_list) == 1, "Mixed periodic boundary conditions"
#     pbc = pbc_list[0]
#
#     soap_desc = SOAP(species=element_list, rcut=rcut, nmax=nmax, lmax=lmax, crossover=crossover, periodic=pbc, sigma=sigma)
#     desc = soap_desc.create(atoms)
#
#     return desc


def get_OMFP(atoms, s=1, p=0, width_cutoff=3.5, maxnatsphere=100):
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
        atoms: ASE atoms object
            ASE atoms object with conatins either a bulk structure or a non-periodic structure
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


    _pbc = list(set(atoms.pbc))
    assert len(_pbc) == 1, "mixed boundary conditions"

    if True in _pbc:
        positions = atoms.get_positions()
        lattice = atoms.get_cell()
        elements = atoms.get_atomic_numbers()
        omfpCalculator = OMFP.stefansOMFP(s=s, p=p, width_cutoff=width_cutoff, maxnatsphere=maxnatsphere)
        omfp = omfpCalculator.fingerprint(positions, elements, lat=lattice)
        omfp = np.array(omfp)

    else:
        positions = atoms.get_positions()
        elements = atoms.get_atomic_numbers()
        width_cutoff = 1000000
        maxnatsphere = len(atoms)
        omfpCalculator = OMFP.stefansOMFP(s=s, p=p, width_cutoff=width_cutoff, maxnatsphere=maxnatsphere)
        omfp = omfpCalculator.globalFingerprint(positions, elements)
        omfp = np.array(omfp)
    print(len(omfp.shape))
    return omfp


def fp_distance(desc1, desc2):
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
        costmat = costmatrix(desc1,desc2)
        ans_pos = scipy.optimize.linear_sum_assignment(costmat)
        fp_dist = 0.
        for index1, index2 in zip(ans_pos[0], ans_pos[1]):
            fp_dist += np.dot((desc1[index1,:]-desc2[index2,:]), (desc1[index1,:]-desc2[index2,:]))
        fp_dist = np.sqrt(fp_dist)


    return fp_dist









def get_area2(lattice):
    """
    Calculates the area of the lattice
    Input:
        lattice: np array
            numpy array containing the lattice vectors
    Return:
        area: float
            area of the lattice
    :return:
    """

    area = 0.
    area += np.linalg.norm(np.cross(lattice[0, :], lattice[1, :]))
    area += np.linalg.norm(np.cross(lattice[0, :], lattice[2, :]))
    area += np.linalg.norm(np.cross(lattice[1, :], lattice[2, :]))

    return area




def alat2ascii(pos, alat_in):
    """
    function that rotates the cell so that the lattice is in the ascii format
    input:
        pos: np array
            numpy array which contains the atom positions
        alat_in: np array
            numpy array containing the lattice parameters
    Retrun:
        pos_out: np array
            numpy array containing the positions after the rotation
        alat_out: np array
            numpy array containing the lattice vecotr after the rotation
    """

    alat_out = np.zeros((3, 3))

    r1 = np.linalg.norm(alat_in[:, 0])
    r2 = np.linalg.norm(alat_in[:, 1])
    r3 = np.linalg.norm(alat_in[:, 2])

    alat_out[0, 0] = r1
    alat_out[0, 1] = np.dot(alat_in[:, 0], alat_in[:, 1]) / r1
    alat_out[0, 2] = np.dot(alat_in[:, 0], alat_in[:, 2]) / r1
    alat_out[1, 1] = np.sqrt(r2 ** 2 - alat_out[0, 1] ** 2)
    alat_out[1, 2] = (np.dot(alat_in[:, 1], alat_in[:, 2]) - alat_out[0, 1] * alat_out[0, 2]) / alat_out[1, 1]
    alat_out[2, 2] = np.sqrt(r3 ** 2 - alat_out[0, 2] ** 2 - alat_out[1, 2] ** 2)

    inv_alat_out = np.linalg.inv(alat_in)
    t = np.matmul(alat_out, inv_alat_out)
    pos_out = np.matmul(t, pos)
    return pos_out, alat_out


def reshape_cell_ascii(lattice_in, imax):
    """
    Reshapes the cell which is in ascii format
    Input:
        lattice_in: numpy array
            numpy array containing the lattice vectors
        imax: int
            maximum iterations for trying new cell combinations
    Return:
        lattice_min: np array
            numpy array containing the best lattice found
        success: bool
            True if a new cell is found, False otherwise
    """
    area_min = get_area2(lattice_in)
    lattice = np.zeros((3, 3))
    lattice_min = np.zeros((3, 3))
    lattice[:, 0] = lattice_in[:, 0].copy()
    success = False
    lattice_min[:,:] = lattice[:,:]
    for iba in range(-imax, imax, 1):
        lattice[:, 1] = lattice_in[:, 1] + iba * lattice_in[:, 0]
        for ica in range(-imax, imax, 1):
            for icb in range(-imax, imax, 1):
                lattice[:, 2] = lattice_in[:, 2] + ica * lattice_in[:, 0] + icb * lattice_in[:, 1]
                area = get_area2(lattice)
                if area < area_min:
                    area_min = area
                    lattice_min[:, :] = lattice[:, :].copy()
                    success = True
                    #print(area)
                    #print(lattice_min)

    # quit()
    return lattice_min, success


def reshape_cell2(atoms, imax):
    """
    Function that reshapes the cell so that the cell is as cubic as possible
    Input:
        atoms: ASE atoms object
            atoms object containing the lattice parameters
        imax: int
            maximum of the lattice expansion to try
    Return:
        atoms object with the changed cell
    """
    lattice_in = atoms.get_cell().T
    positions_in = atoms.get_positions().T
    nat = atoms.get_global_number_of_atoms()

    positions, lattice = alat2ascii(positions_in, lattice_in)
    lattice_temp, success = reshape_cell_ascii(lattice, imax)
    if success:
        lattice[:, :] = lattice_temp[:, :]

    permutation_matrix = np.zeros((3, 3))
    permutation_matrix[0, 1] = 1.
    permutation_matrix[1, 2] = 1.
    permutation_matrix[2, 0] = 1.
    inverse_permutation_matrix = np.linalg.inv(permutation_matrix)
    lattice_temp = np.matmul(lattice, permutation_matrix)
    lattice_temp, success = reshape_cell_ascii(lattice_temp, imax)
    lattice_temp = np.matmul(lattice_temp, inverse_permutation_matrix)

    if not success:
        positions[:, :] = positions_in[:, :]
    else:
        atoms.set_cell(lattice_temp.T, scale_atoms=False, apply_constraint=False)
        atoms.set_positions(positions.T)
        positions = atoms.get_positions(wrap=True).T
        positions_in[:, :] = positions[:, :]
        lattice[:, :] = lattice_temp[:, :]

    # positions_temp = positions
    permutation_matrix = np.zeros((3, 3))
    permutation_matrix[0, 2] = 1.
    permutation_matrix[1, 0] = 1.
    permutation_matrix[2, 1] = 1.
    inverse_permutation_matrix = np.linalg.inv(permutation_matrix)
    lattice_temp = np.matmul(lattice, permutation_matrix)
    lattice_temp, success = reshape_cell_ascii(lattice_temp, imax)
    lattice_temp = np.matmul(lattice_temp, inverse_permutation_matrix)

    if not success:
        positions[:, :] = positions_in[:, :]
    else:
        atoms.set_cell(lattice_temp.T, scale_atoms=False, apply_constraint=False)
        atoms.set_positions(positions.T)
        positions = atoms.get_positions(wrap=True).T
        lattice[:, :] = lattice_temp[:, :]


    atoms.set_positions(positions.T)
    atoms.set_cell(lattice.T, scale_atoms=False, apply_constraint=False)

    return None


def energyandforces(atoms):
    positions = atoms.get_positions()
    cell = atoms.get_cell(complete=False).T
    nat = positions.shape[0]
    e_pot, forces, deralat, stress_tensor = bazant.energyandforces_bazant(cell, positions.T, nat)

    return e_pot, forces.T, deralat, stress_tensor


def frac2cart(atoms, reduced_positions):
    """
    Conversion of the fractional coordinates to cartesian coordinates
    Input:
        atoms: ASE atoms object
            atoms object including the cell vectors
        reduced_positions: np array
            numpy array containing the reduced positions
    Return:
        position: np array
            numpy array containing the cartesian coordinates
    """
    cell = atoms.get_cell()
    positions = np.zeros(reduced_positions.shape)

    for i, at in enumerate(reduced_positions):
        positions[i, :] = np.matmul(cell, at)

    return positions


def cart2frac(atoms):
    """
    Conversiton of the cartesian coordinates to fractional coordinates
    Input:
        atoms: ASE atoms object
            atoms object containing the positions and cell vectors
    Return:
        reduced_positions: np array
            numpy array containing the reduced coordinates
    """

    positions = atoms.get_positions()
    cell = atoms.get_cell()
    inv_cell = np.linalg.inv(cell)
    reduced_positions = np.zeros(positions.shape)
    for i, at in enumerate(positions):
        reduced_positions[i, :] = np.matmul(inv_cell, at)

    return reduced_positions


def lattice_derivative(atoms):
    """
    Calculation of the lattice derivative from the stress tensor. This function cannot be used or has to be changed
    if the stress tensor is not included in the calculator used
    Input:
        atoms: ASE atoms object
            atoms object containing a cell vector and a calculator which has a stress tensor implemented
    Return:
        deralat: np array
            numpy array containing the lattice derivatives
    """
    # BAZANT TEST
    # ______________________________________________________________________________________
    #e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
    stress_tensor = atoms.get_stress(voigt=False,apply_constraint=False)
    # ______________________________________________________________________________________
    cell = atoms.get_cell(complete=False)

    inv_cell = np.linalg.inv(cell)
    prefact = np.linalg.det(cell)
    deralat = - prefact * np.matmul(stress_tensor, inv_cell)
    return deralat

def lattice_derivative_stress(cell,stress):
    """

    """
    # BAZANT TEST
    # ______________________________________________________________________________________
    #e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
    stress_tensor = stress
    # ______________________________________________________________________________________

    inv_cell = np.linalg.inv(cell)
    prefact = np.linalg.det(cell)
    deralat = - prefact * np.matmul(stress_tensor, inv_cell)
    return deralat


def get_moment(velocities):
    # eliminiation of momentum
    s = np.sum(velocities, axis=0) / velocities.shape[0]
    print(s)
    return None


def get_torque(positions, velocities, masses):
    total_mass = np.sum(masses)
    masses_3d = np.vstack([masses] * 3).T
    weighted_positions = positions * masses_3d
    cm = np.sum(weighted_positions, axis=0)
    cm /= total_mass

    tv = 0
    for at, v in zip(positions, velocities):
        tv += np.cross(at, v)
    print(tv)

    return None


def normalize(v):
    """
    Function that normalized a vector of arbitrary length
    Input:
        v: np array
            one dimensional numpy array
    Return:
        Normalized vector as a one dimensional numpy array
    """
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def moment_of_inertia(masses, positions):
    """
    Calculation of the eigenvalues and eigenvector of the moment of inertia tensor
    Input:
        masses: np array
            numpy array containing all the atom masses
        positions: np array
            numpy array containing the atom positions
    Return:
        eigenvalues: np array
            vector containing the eigenvalues of the inertia tensor
        eigenvectors: np array
            matrix containing the eigenvectors of the inertia tensor
    """

    inertia_tensor = np.zeros((3, 3))
    for at, mass in zip(positions, masses):
        inertia_tensor[0, 0] += mass * (at[1] ** 2 + at[2] ** 2)
        inertia_tensor[1, 1] += mass * (at[0] ** 2 + at[2] ** 2)
        inertia_tensor[2, 2] += mass * (at[0] ** 2 + at[1] ** 2)
        inertia_tensor[0, 1] -= mass * (at[0] * at[1])
        inertia_tensor[0, 2] -= mass * (at[0] * at[2])
        inertia_tensor[1, 2] -= mass * (at[1] * at[2])

    inertia_tensor[1, 0] = inertia_tensor[0, 1]
    inertia_tensor[2, 0] = inertia_tensor[0, 2]
    inertia_tensor[2, 1] = inertia_tensor[1, 2]

    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
    return eigenvalues, eigenvectors


def elim_moment(velocities):
    """
    Elimination of the momentum in the velocities
    Input:
        velocities: np array
            numpy array containing the velocities
    Return:
        velocities: np array
            numpy array containing the velocities without momentum
    """

    # eliminiation of momentum
    s = np.sum(velocities, axis=0) / velocities.shape[0]
    velocities -= s
    return velocities


def elim_torque(positions, velocities, masses):
    """
    Elimination of the torque in the velocites
    Input:
        positions: np array
            numpy array containing the atom positions
        velocities: np array
            numpy array containing the atom velocities
        masses: np array
            numpy array containing the atom masses
    Return:
        velocities: np array
            numpy array containing the velocities without torque
    """
    # elimination of torque
    # calculate center of mass and subtracti it from positions
    total_mass = np.sum(masses)
    masses_3d = np.vstack([masses] * 3).T
    weighted_positions = positions * masses_3d
    cm = np.sum(weighted_positions, axis=0)
    cm /= total_mass
    weighted_positions -= cm

    evaleria, teneria = moment_of_inertia(masses, positions)

    vrot = np.zeros((positions.shape[0], 3, 3))
    for iat, at in enumerate(positions):
        vrot[iat, :, 0] = np.cross(teneria[:, 0], at)
        vrot[iat, :, 1] = np.cross(teneria[:, 1], at)
        vrot[iat, :, 2] = np.cross(teneria[:, 2], at)

    velocities = velocities.flatten()
    vrot = vrot.reshape((positions.shape[0] * 3, 3), order="C")

    for i, vec in enumerate(vrot.T):
        vrot[:, i] = normalize(vec)

    weighted_positions += cm

    for i, eval in enumerate(evaleria):
        if abs(eval) > 1e-10:
            alpha = np.dot(vrot[:, i], velocities)
            velocities -= alpha * vrot[:, i]

    velocities = velocities.reshape((positions.shape[0], 3))

    # For debugging reasons this can be switched on to controle if torque is eliminated
    # get_torque(weighted_positions, velocities, masses)
    return velocities


def vcs_soften(atoms, cell_atoms, nsoft):
    """
    Softening the velocites along the softest modes of the postitions and the lattice

    Input:
        atoms: ASE atoms object
            contains the unoptimized atoms object
        cell_atoms: class
            contains unit cell atoms postitions, velocities and masses
        nsoft: int
            number of softening steps
    Return:
        atoms object with softened velocites
    """

    # Softening constants
    eps_dd = 1e-2
    alpha_pos = 1e-3
    alpha_lat = 1e-3
    cell_masses = cell_atoms.masses  # for the moment no masses
    masses = atoms.get_masses()
    # Get initial energy
    # BAZANT TEST
    # ___________________________________________________________________________________________________________
    e_pot_in   = atoms.get_potential_energy()
    #e_pot_in, forces, deralat, stress_tensor = energyandforces(atoms)
    # ___________________________________________________________________________________________________________
    positions_in = atoms.get_positions()
    cell_positions_in = cell_atoms.positions

    # Normalize initial guess
    velocities = atoms.get_velocities()
    cell_velocities = cell_atoms.velocities
    norm_const = eps_dd / (np.sqrt(np.sum(velocities ** 2) + np.sum(cell_velocities ** 2)))
    velocities *= norm_const
    cell_velocities *= norm_const

    # Softening cycle
    for it in range(nsoft):
        w_positions = positions_in + velocities
        w_cell_positions = cell_positions_in + cell_velocities

        atoms.set_positions(w_positions)
        reduced_positions = cart2frac(atoms)
        atoms.set_cell(w_cell_positions, scale_atoms=False, apply_constraint=False)
        positions = frac2cart(atoms, reduced_positions)
        atoms.set_positions(positions)

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        forces = atoms.get_forces()
        e_pot = atoms.get_potential_energy()
        # e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________

        deralat = lattice_derivative(atoms)

        # fd2 is only a control variable
        fd2 = 2. * (e_pot - e_pot_in) / (eps_dd ** 2)

        sdf = np.sum(velocities * forces) + np.sum(cell_velocities * deralat)
        sdd = np.sum(velocities * velocities) + np.sum(cell_velocities * cell_velocities)

        curve = -sdf / sdd
        if it == 0:
            curve0 = curve

        tt = np.sqrt(np.sum(forces * forces) + np.sum(deralat * deralat))
        forces += curve * velocities
        deralat += curve * cell_velocities
        res = np.sqrt(np.sum(forces * forces) + np.sum(deralat * deralat))

        # print("SOFTEN:   ", it, tt, res, curve/fd2, e_pot - e_pot_in)

        w_positions = w_positions + alpha_pos * forces
        w_cell_positions = w_cell_positions + alpha_lat * deralat
        velocities = w_positions - positions_in
        cell_velocities = w_cell_positions - cell_positions_in

        velocities = elim_moment(velocities)
        cell_velocities = elim_torque(w_cell_positions, cell_velocities, cell_masses)

        sdd = eps_dd / np.sqrt(np.sum(velocities ** 2) + np.sum(cell_velocities ** 2))
        if res < (curve * eps_dd * 0.5):
            break
        velocities *= sdd
        cell_velocities *= sdd

    velocities /= norm_const
    cell_velocities /= norm_const
    # Restore initial positions in atoms object
    atoms.set_positions(positions_in)
    atoms.set_cell(cell_positions_in)
    cell_atoms.positions = cell_positions_in

    return velocities, cell_velocities


def soften(atoms, nsoft):
    """
    Softening the velocites along the softest modes
    Input:
        atoms: ASE atoms object
            atoms object containing positions, velocities and a calculator
        nsoft: int
            number of softening iterations
    Return:
        ASE atoms object with softened velocities
    """

    # Softening constants
    eps_dd = 1e-2
    alpha = 1e-3
    # Get initial energy
    # BAZANT TEST
    # ___________________________________________________________________________________________________________
    e_pot_in = atoms.get_potential_energy()
    #e_pot_in, forces, deralat, stress_tensor = energyandforces(atoms)
    # ___________________________________________________________________________________________________________
    positions_in = atoms.get_positions()
    masses = atoms.get_masses()

    # Normalize initial guess
    velocities = atoms.get_velocities()
    norm_const = eps_dd / np.sqrt(np.sum(velocities ** 2))
    velocities *= norm_const

    # Softening cycle
    for it in range(nsoft):
        w_positions = positions_in + velocities
        atoms.set_positions(w_positions)

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        e_pot = atoms.get_potential_energy()
        forces = atoms.get_forces()
        #e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________

        # fd2 is only a control variable
        fd2 = 2. * (e_pot - e_pot_in) / eps_dd ** 2

        sdf = np.sum(velocities * forces)
        sdd = np.sum(velocities * velocities)
        curve = -sdf / sdd
        if it == 0:
            curve0 = curve

        tt = np.sqrt(np.sum(forces * forces))
        forces += curve * velocities
        res = np.sqrt(np.sum(forces * forces))

        # Print statement for debugging reasons
        # print(it, tt, res, curve/ fd2,e_pot - e_pot_in)

        w_positions = w_positions + alpha * forces
        velocities = w_positions - positions_in

        velocities = elim_moment(velocities)
        velocities = elim_torque(w_positions, velocities, masses)

        sdd = eps_dd / np.sqrt(np.sum(velocities ** 2))
        if res < (curve * eps_dd * 0.5):
            break
        velocities *= sdd
    velocities /= norm_const
    # Restore initial positions in atoms object
    atoms.set_positions(positions_in)

    return velocities


def md(atoms, dt, n_max=4, verbose=True):
    """
    velocity Verlet MD which visits n_max maxima

    Input:
        atoms: ASE atoms object
            contains the unoptimized atoms object
        cell_atoms: class
            contains unit cell atoms postitions, velocities and masses
        dt: float
            velocity Verlet timestep
        n_max: int
            number of maxima visited before MD stops
        verbose: bool
            if True the iterationstep, kinetic energy, potential energy and total enery is printed. Furthermore
            in each iteration a ascii file of the current structure is wirtten

    Return:
        atoms object containing all the information after the MD
    """
    # MD which visits at least three max
    masses = atoms.get_masses()[:, np.newaxis] / atoms.get_masses()[:, np.newaxis]  # for the moment no masses

    # BAZANT TEST
    # ___________________________________________________________________________________________________________
    epot_old = atoms.get_potential_energy()
    forces = atoms.get_forces()
    # epot_old, forces, deralat, stress_tensor = energyandforces(atoms)
    # ___________________________________________________________________________________________________________

    sign_old = -1
    n_change = 0
    i = 0
    while n_change < n_max:
        velocities = atoms.get_velocities()
        positions = atoms.get_positions()
        epot_old = atoms.get_potential_energy()
        # Update postions
        atoms.set_positions(positions + dt * velocities + 0.5 * dt * dt * (forces / (masses)))

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # epot_old  = atoms.get_potential_energy()
        new_forces = atoms.get_forces()
        # epot_old, new_forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________
        atoms.set_velocities(velocities + 0.5 * dt * ((forces + new_forces) / masses))

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        epot = atoms.get_potential_energy()
        # epot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________

        sign = int(np.sign(epot_old - epot))
        if sign_old != sign:
            sign_old = sign
            n_change += 1

        if verbose:
            e_kin = 0.5 * np.sum(masses * atoms.get_velocities() * atoms.get_velocities())
            md_msg = "MD STEP:  {:d}   e_pot: {:1.5f}  e_kin:  {:1.5f}   e_tot:  {:1.5f}".format(i, epot, e_kin,
                                                                                                 epot + e_kin)
            print(md_msg)
            filename = "MD{0:05d}.xyz".format(i)
            write(filename, atoms)

        if i > 10000:
            warning_msg = "Geometry did not converge in {:d} optimizations steps".format(i)
            warnings.warn(warning_msg, FutureWarning)
            break

        i += 1

    return None


def vcsmd(atoms, OFT, cell_atoms, dt, n_max=2, verbose=True):
    """
    Variable cell shape velocity Verlet MD which visits n_max maxima

    Input:
        atoms: ASE atoms object
            contains the unoptimized atoms object
        cell_atoms: class
            contains unit cell atoms postitions, velocities and masses
        dt: float
            velocity Verlet timestep
        n_max: int
            number of maxima visited before MD stops
        verbose: bool
            if True the iterationstep, kinetic energy, potential energy and total enery is printed. Furthermore
            in each iteration a ascii file of the current structure is wirtten

    Return:
        atoms object containing all the information after the MD

    """
    # MD which visits at least three max
    # Initializations
    masses = atoms.get_masses()[:, np.newaxis] / atoms.get_masses()[:, np.newaxis]  # for the moment no masses
    cell_masses = cell_atoms.masses[:, np.newaxis]

    sign_old = -1
    n_change = 0
    i = 0

    # BAZANT TEST
    # ___________________________________________________________________________________________________________
    # forces = atoms.get_forces()
    # epot_old  = atoms.get_potential_energy()
    epot_old, forces, stress, is_dft = OFT.energyandforces(atoms)
    #e_pot_cur, forces, deralat, stress_tensor = energyandforces(atoms)
    # ___________________________________________________________________________________________________________

    cell_forces = lattice_derivative_stress(atoms.get_cell(complete=False), stress)

    #while n_change < n_max:
    for kj in range(1000):
        velocities = atoms.get_velocities()
        positions = atoms.get_positions()

        # Update postions
        atoms.set_positions(positions + dt * velocities + 0.5 * dt * dt * (forces / (masses)))

        # Update lattice so that fractional coordinates remain invariant
        reduced_postitions = cart2frac(atoms)
        cell_positions = cell_atoms.positions
        cell_velocities = cell_atoms.velocities
        cell_atoms.positions = cell_positions + dt * cell_velocities + 0.5 * dt * dt * (cell_forces / cell_masses)
        atoms.set_cell(cell_atoms.positions, scale_atoms=False, apply_constraint=False)
        positions = frac2cart(atoms, reduced_postitions)
        atoms.set_positions(positions)

        # Update velocities
        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # new_forces = atoms.get_forces()
        epot, new_forces, stress, is_dft = OFT.energyandforces(atoms)
        # e_pot_cur, new_forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________
        new_cell_forces = lattice_derivative_stress(atoms.get_cell(complete=False), stress)

        atoms.set_velocities(velocities + 0.5 * dt * ((forces + new_forces) / masses))

        # Update velocities of the cell atoms
        cell_atoms.velocities = cell_velocities + 0.5 * dt * ((cell_forces + new_cell_forces) / cell_masses)

        forces = new_forces
        cell_forces = new_cell_forces


        sign = int(np.sign(epot_old - epot))
        if sign_old != sign:
            sign_old = sign
            n_change += 1

        epot_old = epot

        if verbose:
            e_kin = 0.5 * np.sum(masses * atoms.get_velocities() * atoms.get_velocities())
            e_kin += 0.5 * np.sum(cell_masses * cell_atoms.velocities * cell_atoms.velocities)
            md_msg = "MD STEP:  {:d}   e_pot: {:1.5f}  e_kin:  {:1.5f}   e_tot:  {:1.5f}".format(i, epot, e_kin,
                                                                                                 epot + e_kin)
            print(md_msg)
            filename = "MD{0:05d}.ascii".format(i)
            write(filename, atoms)

        if i > 10000:
            warning_msg = "Geometry did not converge in {:d} optimizations steps".format(i)
            warnings.warn(warning_msg, FutureWarning)
            break

        i += 1

    return None


def escape_trial(atoms, OFT, dt, T):
    """
    Escape loop to find a new minimum
    Input:
        atoms: ASE atoms object
            ase atoms object containing postitions and a calculator
        dt: float
            timestep for the MD
        T: float
            temperature in (K) to scale the velocities accordingly
    Return:
        ASE atoms object containing a new minimum
    """

    # BAZANT TEST
    # ___________________________________________________________________________________________________________
    # e_pot_curr = atoms.get_potential_energy()
    # e_pot_curr, forces, deralat, stress_tensor = energyandforces(atoms)
    e_pot_curr, new_forces, stress, is_dft = OFT.energyandforces(atoms)
    # ___________________________________________________________________________________________________________
    escape = 0.0
    beta_s = 1.1
    while escape < 1e-3:

        MaxwellBoltzmannDistribution(atoms, temperature_K=T)

        if True in atoms.pbc:
            mass = .75 * np.sum(atoms.get_masses()) / 10.
            cell_atoms = cell_atom(mass=mass, positions=atoms.get_cell())
            cell_atoms.set_velocities_boltzmann(temperature=T)
            #velocities, cell_velocities = vcs_soften(atoms, cell_atoms, 20)
            #atoms.set_velocities(velocities)
            #cell_atoms.velocities = cell_velocities
            vcsmd(atoms, OFT, cell_atoms, dt, verbose=True)
            reshape_cell2(atoms, 6)
            # reshape_cell(atoms,3)
            vcs_optimizer(atoms, verbose=False)
        else:
            velocities = soften(atoms, 20)
            atoms.set_velocities(velocities)
            md(atoms, dt, verbose=False)
            optimizer(atoms, verbose=False)

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # e_pot = atoms.get_potential_energy()
        e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________
        escape = abs(e_pot - e_pot_curr)
        T *= beta_s
        print("TEMPARATUR:   ", T, escape, e_pot_curr)

    return


def in_history(e_pot_cur, history, ):
    i = 0
    for s in history:
        e_pot = s[0]
        e_diff = abs(e_pot - e_pot_cur)
        if e_diff < 1e-3:
            i += 1
    return i + 1

def in_history_fp(fp1, history, epot ):
    i = 0
    for s in history:
        fp2 = s[5]
        fp_dist = fp_distance(fp1, fp2)
        ediff = abs(epot-s[0])
        #print(fp_dist, ediff)
        if fp_dist < 100.:
            i += 1
    return i + 1


def optimizer(atoms, OFT, initial_step_size=0.0001, nhist_max=10, alpha_min=1e-5, eps_subsp=1e-4,
              max_force_threshold=0.00002, verbose=False):
    """
    Variable cell shape optimizer from Gubler et. al 2022 arXiv2206.07339
    Function optimizes atom positions and lattice vecotors

    Input:
        atoms: ASE atoms object
            contains the unoptimized atoms object
        initial_step_size: float
            initial step size to start the geometry optimization
        nhist_max: int
            maximum structures which are stored in the history MORITZ?
        alpha_min: float
            minimal step size during the optimization
        eps_subsop: float
            MORITZ?
        max_force_threshold: float
            convergence criterion at which maximal force component convergence is reached.
        verbose: bool
            If true the optimization will print the iteration step, the energy and the maximal force component of the
            current iteration. Furthermore v_sim files are written in each iteration.

    Return:
        optimized structure in atoms object

    """

    nat = atoms.get_positions().shape[0]
    optim = free_or_fixed_cell_sqnm.free_sqnm(nat=nat, initial_step_size=initial_step_size, nhist_max=nhist_max,
                                              alpha_min=alpha_min, eps_subsp=eps_subsp)

    max_force_comp = 100
    i = 0
    while max_force_comp > max_force_threshold:
        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # energy = atoms.get_potential_energy()# and also deralat
        # forces = atoms.get_forces()
        # deralat = lattice_derivative(atoms)
        # energy, forces, deralat, stress_tensor = energyandforces(atoms)
        energy, forces, stress, is_dft = OFT.energyandforces(atoms)
        deralat = lattice_derivative_stress(atoms.get_cell(complete=False),stress)

        # ___________________________________________________________________________________________________________


        i += 1
        max_force_comp = np.max(forces)
        positions = atoms.get_positions()

        new_positions = optim.optimizer_step(positions.T, energy, forces.T)
        atoms.set_positions(new_positions.T)

        if verbose:
            opt_msg = "OPT Step: {:d}   energy: {:1.8f}  max_force_comp:  {:1.5e}".format(i, energy, max_force_comp)
            print(opt_msg)
            filename = "OPT{0:05d}.ascii".format(i)
            write(filename, atoms)

        if i > 10000:
            warning_msg = "Geometry did not converge in {:d} optimizations steps".format(i)
            warnings.warn(warning_msg, FutureWarning)
            break

    return None


def vcs_optimizer(atoms, OFT,initial_step_size=0.01, nhist_max=10, lattice_weight=2, alpha_min=1e-3, eps_subsop=1e-4,
                  max_force_threshold=0.0002, verbose=False):
    """
    Variable cell shape optimizer from Gubler et. al 2022 arXiv2206.07339
    Function optimizes atom positions and lattice vecotors

    Input:
        atoms: ASE atoms object
            contains the unoptimized atoms object
        initial_step_size: float
            initial step size to start the geometry optimization
        nhist_max: int
            maximum structures which are stored in the history MORITZ?
        lattice_weight: float
            lattice weight MORITZ?
        alpha_min: float
            minimal step size during the optimization
        eps_subsop: float
            MORITZ?
        max_force_threshold: float
            convergence criterion at which maximal force component convergence is reached.
        verbose: bool
            If true the optimization will print the iteration step, the energy and the maximal force component of the
            current iteration. Furthermore v_sim files are written in each iteration.

    Return:
        optimized structure in atoms object

    """

    # Get nessecairy parameters from atoms object
    nat = atoms.get_positions().shape[0]
    init_lat = atoms.get_cell().T
    # Initialize optimizer class
    optim = periodic_sqnm.periodic_sqnm(nat, init_lat, initial_step_size, nhist_max, lattice_weight, alpha_min,
                                        eps_subsop)

    # Initialize futher parameters
    max_force_comp = 100
    i = 0
    while max_force_comp > max_force_threshold:
        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        # energy = atoms.get_potential_energy()# and also deralat
        # forces = atoms.get_forces()
        # deralat = lattice_derivative(atoms)
        # energy, forces, deralat, stress_tensor = energyandforces(atoms)
        energy, forces, stress, is_dft = OFT.energyandforces(atoms)
        deralat = lattice_derivative_stress(atoms.get_cell(complete=False),stress)

        # ___________________________________________________________________________________________________________

        i += 1
        max_force_comp = np.max(forces)
        max_latderiv_comp = np.max(deralat)
        max_force_comp = np.maximum(max_force_comp, max_latderiv_comp)
        positions = atoms.get_positions()
        lattice = atoms.get_cell().T

        new_positions, new_lattice = optim.optimizer_step(positions.T, lattice, energy, forces.T, deralat)
        atoms.set_positions(new_positions.T)
        atoms.set_cell(new_lattice.T, scale_atoms=False, apply_constraint=False)

        if verbose:
            opt_msg = "OPT Step: {:d}   energy: {:1.8f}  max_force_comp:  {:1.5e}  gain ratio:   {:1.8f}".format(i, energy, max_force_comp, optim.optimizer.gainratio)
            print(opt_msg)
            filename = "OPT{0:05d}.ascii".format(i)
            write(filename, atoms)

        if i > 10000:
            warning_msg = "Geometry did not converge in {:d} optimizations steps".format(i)
            warnings.warn(warning_msg, FutureWarning)
            break

    return


def main():
    # Initial settings and paths
    #path = "/home/marco/NequipMH/data/38/"
    #filename = "acc-00000.xyz"
    # path = "/home/marco/NequipMH/data/"
    # filename = "LJC.extxyz"
    path = "./"
    filename = "Si_in3.extxyz"

    # model_path = " "
    dt = 0.01
    e_diff = 3.
    alpha_a = 0.95
    alpha_r = 1.05
    n_acc = 0
    n_min = 0
    history = []
    acc_min_list = []
    T = 2500
    beta_decrease = 1. / 1.03
    beta_increase = 1.03
    enhanced_feedback = False
    compare_energy = False

    # Read local minimum input file
    atoms = read(path + filename)
    atoms.set_positions(np.random.normal(atoms.get_positions(),0.01))
    #cell=np.array([[100,0,0],[0,100,0],[0,0,100]])
    #atoms.set_cell(cell)


    # set up GP hyperparameters
    kernels = ['twobody', 'threebody'] # use 2+3 body kernel
    parameters = {'cutoff_twobody': 5.0,
                  'cutoff_threebody': 3.5}
    pm = ParameterHelper(
        kernels = kernels,
        random = True,
        parameters=parameters
    )

    hm = pm.as_dict()
    hyps = hm['hyps']
    cut = hm['cutoffs']
    # print('hyps', hyps)

    gp_model = GaussianProcess(
        kernels = kernels,
        component = 'mc', # If you are using ASE, please set to "mc" no matter for single-component or multi-component
        hyps = hyps,
        cutoffs = cut,
        hyp_labels = ['sig2','ls2','sig3','ls3','noise'],
        opt_algorithm = 'L-BFGS-B',
        n_cpus = 7
    )


    pseudo_dir = '/home/marco/NequipMH/src/python/'
    pseudopotentials = {'Si': 'Si.UPF'}

    # PAY ATTENTION TO REMOVE PWSCF.SAVE IF RESTARTED
    input_data = {
       'system': {
           'ecutwfc': 60,
           'ecutrho': 280},
       'disk_io': 'low',
       'electrons': {
           'mixing_beta' : 0.4,
           'startingpot' : 'file',
       }
    }  # automatically put into 'control'

    dft_calc = Espresso(pseudo_dir=pseudo_dir, pseudopotentials=pseudopotentials,
                   tstress=True, tprnfor=True, kpts=(3, 3, 3), input_data=input_data)



    # dft_calc = LennardJones()
    # dft_calc.parameters.epsilon = 1.0
    # dft_calc.parameters.sigma = 1.0
    # dft_calc.parameters.rc = 6.0
    # atoms.calc = dft_calc
    # print(atoms.get_potential_energy())
    # quit()


    # atoms.calc = dft_calc
    # print(atoms.get_potential_energy())
    # quit()


    OFT = mh_otf(gp_model, dft_calc, read=False)

    vcs_optimizer(atoms, OFT,verbose=True)
    quit()
    atoms_cur = atoms.copy()

    # BAZANT TEST
    # ___________________________________________________________________________________________________________
    # e_pot_cur = atoms.get_potential_energy()
    #e_pot_cur, forces, deralat, stress_tensor = energyandforces(atoms)
    e_pot_cur, forces, stress, is_dft = OFT.energyandforces(atoms)
    # ___________________________________________________________________________________________________________
    fp = get_OMFP(atoms)
    history.append((e_pot_cur, 1, T, e_diff, "A", fp))

    for i in range(10000):
        escape_trial(atoms, OFT, dt, T)

        filename = "MIN" + str(n_min).zfill(5) + ".ascii"
        write(filename, atoms)

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        e_pot = atoms.get_potential_energy()
        #e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________

        n_min += 1
        log_msg = "LOG:  {:d}  Epot:  {:1.5f}   E_diff:  {:1.5f}    Temp:   {:1.5f} ".format(n_min, e_pot, e_diff, T)
        print(log_msg)

        if abs(e_pot_cur - e_pot) < e_diff:
            e_pot_cur = e_pot
            e_diff *= alpha_a
            n_acc += 1
            acc_min_list.append([atoms.copy(), e_pot_cur])
            acc_min_list = sorted(acc_min_list, key=lambda x: x[1])
            for k, acc_min in enumerate(acc_min_list):
                atom = acc_min[0]
                # print("DEBUGGY    ", acc_min[1])
                filename = "ACC" + str(k).zfill(5) + ".ascii"
                write(filename, atom)
            atoms_cur = atoms.copy()
            acc_rej = "A"
        else:
            e_diff *= alpha_r
            acc_rej = "R"
            # print(atoms.get_pote ntial_energy()-e_pot_cur)



        fp = get_OMFP(atoms)
        if compare_energy:
            n_visits = in_history(e_pot, history)
        else:
            n_visits = in_history_fp(fp, history, e_pot)


        if n_visits > 1:
            if enhanced_feedback:
                T = T * beta_increase * (1. + 4. * np.log(float(n_visits)))
            else:
                T *= beta_increase
        else:
            T *= beta_decrease

        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        e_pot = atoms.get_potential_energy()
        #e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________
        history.append((e_pot, n_visits, T, e_diff, acc_rej, fp))

        if i % 1 == 0:
            f = open("history.dat", "w")
            for s in history:
                history_msg = "{:1.5f}  {:d}  {:1.5f}  {:1.5f}  {:s} \n".format(s[0], s[1], s[2], s[3], s[4])
                f.write(history_msg)
            f.close()

        atoms = atoms_cur.copy()
        # BAZANT TEST
        # ___________________________________________________________________________________________________________
        e_pot = atoms.get_potential_energy()
        #e_pot, forces, deralat, stress_tensor = energyandforces(atoms)
        # ___________________________________________________________________________________________________________


if __name__ == '__main__':
    main()

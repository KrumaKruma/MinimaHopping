from minimahopping.minhop import Minimahopping
from ase.calculators.eam import EAM
from ase.cluster import Icosahedron
from mpi4py import MPI
from ase.cluster.wulff import wulff_construction

atoms = wulff_construction('Na', surfaces=[(1, 0, 0), (0, 1, 0),(0, 0, 1)], energies=[0.001, 0.001, 0.15],
                       size=13, # maximum number of atoms
                       structure='bcc', rounding='above')
calculator = EAM(potential="Na_v2.eam.fs")
atoms.calc = calculator
fnrm = 5e-3
minima_threshold = 1e-4
with Minimahopping(atoms, fmax=fnrm, fingerprint_threshold=minima_threshold, verbose_output=False, T0=2000, dt0=0.1, use_MPI=True) as mh:
    mh(totalsteps=100)




import numpy as np

from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.stress import full_3x3_to_voigt_6_stress
from numba import jit, cuda
#from sklearn.neighbors import NearestNeighbors
import scipy.spatial




class LennardJones(Calculator):
    

    implemented_properties = ['energy', 'forces']
    # implemented_properties += ['stress', 'stresses']  # bulk properties
    default_parameters = {
        'epsilon': 1.0,
        'sigma': 1.0,
        'rc': None,
        'ro': None,
        'smooth': False,
    }
    nolabel = True

    def __init__(self, **kwargs):


        Calculator.__init__(self, **kwargs)

        if self.parameters.rc is None:
            self.parameters.rc = 3 * self.parameters.sigma

        if self.parameters.ro is None:
            self.parameters.ro = 0.66 * self.parameters.rc

        self.nl = None

    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        e, f = LennardJones.force(atoms.get_positions())

        self.results['energy'] = e
        self.results['forces'] = f

    @staticmethod 
    @jit(nopython=True, parallel=False)
    def force(ats):
        nat = ats.shape[0]
        e = 0.0
        f = np.zeros(ats.shape)

        for i in range(nat):
            for j in range(i):
                dr = ats[i,:] - ats[j,:]
                dd = np.sum(dr**2)
                dd2 = 1.0 / dd
                dd6 = dd2 * dd2 * dd2
                dd12 = dd6 * dd6
                e += 4.0 * (dd12 - dd6)
                tt = 24.0 * dd2 * (2.0 * dd12 - dd6)
                t = dr * tt
                f[i,:] += t
                f[j,:] -= t
        return e, f

# not computing e for the moment
@cuda.jit
def force_cuda(ats, neis, nneis, f, eat):
    threadId = cuda.threadIdx.x
    blockId = cuda.blockIdx.x
    blockDim = cuda.blockDim.x
    i = blockDim * blockId + threadId
    nat = ats.shape[0]
    if i < nat:
        f[i,0] = 0.
        f[i,1] = 0.
        eat[i] = 0.
        #for j in range(i):
        #neis = tree.query_ball_point(ats[i,:], rc)
        for ij in range(nneis[i]): #neis[i,:nneis[i]]:
            j = neis[i,ij]
            if True: #i < j:
                drx = ats[i,0] - ats[j,0]
                dry = ats[i,1] - ats[j,1]
                dd = drx**2 + dry**2
                dd2 = 1.0 / dd
                dd6 = dd2 * dd2 * dd2
                dd12 = dd6 * dd6
                #e += 4.0 * (dd12 - dd6) #+ lj_rc
                eat[i] += 4.0 * (dd12 - dd6)
                tt = 24.0 * dd2 * (2.0 * dd12 - dd6)
                f[i,0] += tt * drx
                f[i,1] += tt * dry
                #f[j,:] -= t





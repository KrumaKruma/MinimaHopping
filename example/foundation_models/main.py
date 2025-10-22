from minimahopping.minhop import Minimahopping
from make_structure import get_random_packed
from tensorpotential.calculator import grace_fm
from mpi4py import MPI
import os
from ase.io import write
from sqnm.vcsqnm_for_ase import aseOptimizer


# Choose the composition to search.
COMPOSITION = f'(CaTiO3)4'

# If MPI is enabled, use N+1 SLURM processes but only N GPUS. 
USE_MPI = False


def main():
    
    init_structure = get_random_packed(COMPOSITION, scale=1.5)


    if USE_MPI:
        rank = MPI.COMM_WORLD.Get_rank()
        if rank == 0: # No GPU for the master process
            print(f'Rank {rank} is not using GPU')
            os.environ["CUDA_VISIBLE_DEVICES"] = ""
        else:
            print(f'Rank {rank} is using GPU {rank - 1}')
            os.environ["CUDA_VISIBLE_DEVICES"] = f"{rank - 1}"

    if (not USE_MPI) or (rank > 0):
        init_structure.calc = get_calc()
        opt = aseOptimizer(init_structure,
                           vc_relax=True, 
                           force_tol=1.e-4,
                           )
        opt.optimize()
        init_structure.wrap()

    write('initial_structure.extxyz', init_structure)
            

    with Minimahopping(init_structure, 
                       symprec=1e-5,
                       verbose_output=False, 
                       T0=500, 
                       dt0=0.1, 
                       use_MPI=USE_MPI, 
                       mdmin=3,
                       fixed_cell_simulation=False,
                       fmax=0.002,
                       energy_threshold=0.002,
                       fingerprint_threshold=0.005,
                       write_graph_output=False,
                       alpha_reject=1.05,
                       alpha_accept=0.7,
                       beta_increase=1.1,
                       beta_decrease=0.8,
                       enhanced_feedback=True,
                       collect_md_data=False,
                       ) as mh:
        mh(totalsteps=100000)

# Add your favouorite universal potential here
def get_calc():
    return grace_fm('GRACE-2L-OAM')



if __name__ == '__main__':
    main()

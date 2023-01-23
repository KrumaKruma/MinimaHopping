from mpi4py import MPI
import time
import sirius_ase.ase_simulation
import sirius_ase.siriusCalculator
import numpy as np
import json
import traceback
import sys
import ase.io
import os
from minimahopping.minhop import Minimahopping

numberOfSlaves = 3
slaveRankSize = 2
globalNumberOfProcesses = numberOfSlaves * slaveRankSize + 1
print(globalNumberOfProcesses)
comm_world = MPI.COMM_WORLD # this communicator contains the entire mpi world

rank = comm_world.Get_rank()
size = comm_world.Get_size()

if globalNumberOfProcesses != size:
    print('wrong number of mpi processes given to program. Expected number of processe, ', globalNumberOfProcesses)
    comm_world.Abort()
    quit()

# all processes with same color will be put in a group later

if rank == 0: # master group
    color = 0

if rank > 0: # first group with two processes 
    color = 1

if rank > 2: # second group with two processes
    color = 2

if rank > 4: # third group with two processes
    color = 3

# comm_world.barrier()

# create group communicators based on color

group_communicator = comm_world.Split(color, rank)
group_rank = group_communicator.Get_rank()
group_size = group_communicator.Get_size()

if rank == 0:
    print('group_rank, group_size, rank, size')

# comm_world.barrier()

# time.sleep(0.1)
msg = 'group rank: {:d}   group size: {:d}  rank: {:d}  size: {:d}'.format(group_rank, group_size, rank, size)
print(msg)


# Set up calculator for minimahopping
siriusJsonFileName = sys.argv[1]
structufileName = sys.argv[2]
if not os.path.exists(siriusJsonFileName):
    print('json file does not exist')
    quit()

atoms = ase.io.read(filename=structufileName, parallel=False, index = 0)

f = open(siriusJsonFileName)
jsonparams = json.load(f)
f.close()
try:
    pp_files = jsonparams["unit_cell"]["atom_files"]
    pw_cutoff = jsonparams['parameters']["pw_cutoff"]
    gk_cutoff = jsonparams['parameters']["gk_cutoff"]
    functionals = jsonparams['parameters']['xc_functionals']
    kpoints = jsonparams['parameters']['ngridk']
    kshift = jsonparams['parameters']["shiftk"]
    if "atom_types" in jsonparams["unit_cell"]:
        jsonparams["unit_cell"].pop("atom_types")
    jsonparams["unit_cell"].pop("atom_files")
except KeyError:
    print("required parameter was missing")
    traceback.print_exc()
    quit()

# If not master rank give group comunicator to sirius calculator
if rank != 0:
    # give the group communicator to the sirius calculator...
    atoms.calc = sirius_ase.siriusCalculator.SIRIUS(atoms, pp_files, functionals, kpoints, kshift, pw_cutoff, gk_cutoff, jsonparams, group_communicator)

# Start MPI Minimahopping
with Minimahopping(atoms, verbose_output=True, T0=500, dt=0.1, use_MPI=True, totalWorkers=numberOfSlaves) as mh:

    # mpi example that will run for two minutes and assumes that 1 process is used as the server and the rest (n-1) as clients.
    # mh = Minimahopping(atoms, verbose_output=False, 
    #    T0=2000, dt=0.1, use_MPI=True, fingerprint_threshold=5e-4, run_time='0-01:00:00', totalWorkers=numberOfSlaves)
    mh(totalsteps=5)




atoms.calc.close()



comm_world.Barrier()

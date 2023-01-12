import sirius_ase.ase_simulation
import sirius_ase.siriusCalculator
import numpy as np
import json
import traceback
import sys
import ase.io
import os
from minimahopping.minhop import Minimahopping


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

atoms.calc = sirius_ase.siriusCalculator.SIRIUS(atoms, pp_files, functionals, kpoints, kshift, pw_cutoff, gk_cutoff, jsonparams)

with Minimahopping(atoms, verbose_output=True, T0=2000, dt=0.1, use_MPI=False) as mh:

    # mpi example that will run for two minutes and assumes that 1 process is used as the server and the rest (n-1) as clients.
    # from mpi4py import MPI
    # totalWorkers = MPI.COMM_WORLD.Get_size() - 1
    # mh = Minimahopping(initial_configuration, verbose_output=False, 
    #     T0=2000, dt=0.1, use_MPI=True, fingerprint_threshold=5e-4, run_time='0-00:02:00', totalWorkers=totalWorkers)
    mh(totalsteps=5)




atoms.calc.close()




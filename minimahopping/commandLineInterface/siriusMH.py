from mpi4py import MPI
import numpy as np
import json
import traceback
import ase.io
from minimahopping.minhop import Minimahopping
from ase.build import bulk
import logging
import argparse
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mhJsonFilename", help="mh json filename")
    parser.add_argument("siriusJsonFilename", help="sirius json filename")
    parser.add_argument("structureFilename", help="structure filename")
    parser.add_argument("--numberOfSlaves", help="number of minima hopping slaves", type=int, default=1)
    parser.add_argument("--slaveRankSize", help="number of ranks per minima hopping slave", type=int, default=1)
    parser.add_argument("--distributedStartingStructure", help="""
                        If this argument is present, the input structure must contain a list of numberOfSlaves structures.
                        Each slave will get a different structure to start with.
                        """, action="store_true")
    parser.add_argument("--totalSteps", help="total number of steps", type=int, default=-1, required=False)
    parser.add_argument("--pressureGPa", help="pressure in GPa", type=float, default=0.0, required=False)
    parser.add_argument("--dispersion", help="json file containing the information about dispersion correction. Allowed are D3 and D4 dispersion corrections", type=str, default=None, required=False)
    args = parser.parse_args()

    numberOfSlaves = args.numberOfSlaves
    slaveRankSize = args.slaveRankSize
    siriusJsonFileName = args.siriusJsonFilename
    structufileName = args.structureFilename

    if numberOfSlaves < 1 or slaveRankSize < 1:
        raise ValueError('numberOfSlaves and slaveRankSize must be greater than 0')

    if numberOfSlaves == 1:
        globalNumberOfProcesses = slaveRankSize
        useSlave = False
    else:
        globalNumberOfProcesses = numberOfSlaves * slaveRankSize + 1
        useSlave = True
    comm_world = MPI.COMM_WORLD

    rank = comm_world.Get_rank()
    size = comm_world.Get_size()
    if useSlave and rank == 0:
        isMaster = True
    else:
        isMaster = False

    if globalNumberOfProcesses != size:
        raise ValueError('wrong number of mpi processes given to program. Expected number of processes: ', globalNumberOfProcesses)

    color = (rank + slaveRankSize - 1) // slaveRankSize
    print("Mpi rank: ", rank, " is assigned to color: ", color, flush=True)

    if useSlave:
        print("Creating group communicator", rank, flush=True)
        group_communicator = comm_world.Split(color, rank)
    else:
        group_communicator = comm_world
    group_rank = group_communicator.Get_rank()
    group_size = group_communicator.Get_size()

    # only import sirius after arguments are parsed and sanity checks are done.
    import sirius_ase.ase_simulation
    import sirius_ase.siriusCalculator

    if useSlave and args.distributedStartingStructure:
        print("Warning: I did not test this. It might not work.")
        atoms = ase.io.read(filename=structufileName, parallel=False, index = color - 1)
    else:
        atoms = ase.io.read(filename=structufileName, parallel=False, index = 0)

    f = open(siriusJsonFileName)
    jsonparams = json.load(f)
    f.close()
    # check sirius json file for required parameters
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

    f = open(args.mhJsonFilename)
    mhJson = json.load(f)
    mhJson["use_MPI"] = useSlave
    f.close()

    if args.dispersion:

        with open(args.dispersion) as f:
            allowed_dispersion = ["d3", "d4"]
            disp_dict = json.load(f)
            dispersion_functional = disp_dict["functional"]
            dispersion_method = disp_dict["method"]
            # convert dispersion_method string to lower case:
            if dispersion_method.lower() not in allowed_dispersion:
                raise ValueError("dispersion_method must be one of: ", allowed_dispersion)
            # check if dampening key exists:
            dispersion_damp = ""
            if "dampening" in disp_dict:
                dispersion_damp = disp_dict["dampening"]

        if dispersion_method.lower() == 'd3':
            from dftd3.ase import DFTD3
            if dispersion_damp != "":
                disp_calc = DFTD3(method=dispersion_functional, damping=dispersion_damp)
            else:
                disp_calc = DFTD3(method=dispersion_functional)
        elif dispersion_method.lower() == 'd4':
            from dftd4.ase import DFTD4
            if dispersion_damp != "":
                disp_calc = DFTD4(method=dispersion_functional, damping=dispersion_damp)
            else:
                disp_calc = DFTD4(method=dispersion_functional)


    try:
        if not isMaster: # If not master rank give group comunicator to sirius calculator
            print("Creating calculator on rank", rank, flush=True)
            sirius_calc = sirius_ase.siriusCalculator.SIRIUS(atoms, pp_files, functionals, kpoints, kshift, pw_cutoff, gk_cutoff, jsonparams, communicator=group_communicator, pressure_giga_pascale=args.pressureGPa)
            print("Calculator created on rank", rank, flush=True)

            if args.dispersion:
                from ase.calculators.mixing import SumCalculator
                atoms.calc = SumCalculator([disp_calc, sirius_calc])
            else:
                atoms.calc = sirius_calc

        with Minimahopping(atoms, **mhJson) as mh:
            if args.totalSteps > 0:
                mh(args.totalSteps)
            else:
                mh()
    finally:
        # make sure that not a slave slave tries to close itselve.
        if useSlave and not isMaster and group_rank == 0:
            print("Closing calculator on rank", rank, flush=True)
            atoms.calc.close()
        elif not useSlave and rank == 0:
            print("Closing calculator on rank", rank, flush=True)
            atoms.calc.close()

if __name__ == "__main__":
    main()

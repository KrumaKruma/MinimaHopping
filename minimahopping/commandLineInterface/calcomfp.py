import argparse
from ase.io import read
from minimahopping.mh import parameters
from minimahopping.mh import minimum
import minimahopping.logging.logger as mhlogging
import logging
import os
import numpy as np

def main():

    mhlogging.setupLogger(logging.INFO)

    defaultParams = parameters.minimaHoppingParameters()

    parser = argparse.ArgumentParser(description ='Read a list of input files and calculate the overlap matrix fingerprint and write it to a file.')

    parser.add_argument('-i', '--inputfile1', dest ='file',
                action ='store', help ='input filename of file.', required=True)
    # parse indeices
    parser.add_argument('--index', dest='index', action='store', type=str,
        help='index of first structure in input file', default=0, required= False)
    
    parser.add_argument('-o --output', dest='output', action='store', type=str,
        help='output directory', default='./', required= False)

    # parse additional OMFP parameters
    parser.add_argument('--n_S_orbitals', dest='n_S_orbitals', action='store', type=int,
        help='number of s orbitals for constructing the OMFP', default=1, required= False)
    parser.add_argument('--n_P_orbitals', dest='n_P_orbitals', action='store', type=int,
        help='number of p orbitals for constructing the OMFP', default=1, required= False)
    parser.add_argument('--width_cutoff', dest='width_cutoff', action='store', type=float,
        help='cutoff for the OMFP', default=defaultParams.width_cutoff, required= False)
    parser.add_argument(
        "--exclude",
        dest='exclude',
        action='store', 
        nargs="*",
        type=str,
        default=defaultParams.exclude,
        help='List of elements to exclude in the OMFP.',
        required= False
    )
    
    args = parser.parse_args()

    atoms = read(args.file, index=args.index)
    
    if '.ascii' in args.file:
        atoms.pbc = [True, True, True]

    if not type(atoms) == list:
        atoms = [atoms]
    
    i = 0

    fname_in = args.file
    basename = os.path.basename(fname_in)
    basename = os.path.splitext(basename)[0]

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    for at in atoms:
        m1 = minimum.Minimum(at, 0.0, args.n_S_orbitals, args.n_P_orbitals, args.width_cutoff,
        T=0, ediff=0, exclude=args.exclude)

        maxNatInEnv = m1.maxNatInEnv

        fp =np.zeros((len(at), maxNatInEnv))
        for j in range(len(at)):
            fp[j, :len(m1.fp[j])] = m1.fp[j]

        with open(args.output + '/' + basename + '_omfp_' + str(i) + '.txt', 'w') as f:
            f.write("# fingerprint indes, atom1, atom2, ... \n")
            for ii in range(maxNatInEnv):
                f.write(str(ii) + ' ' + ' '.join([str(fp[j, ii]) for j in range(len(at))]) + '\n')
        i += 1


if __name__ =='__main__':
    main()
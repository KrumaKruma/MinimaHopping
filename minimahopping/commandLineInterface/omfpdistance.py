import argparse
from ase.io import read
from minimahopping.mh import parameters
from minimahopping.mh import minimum
import logging

def main():

    logging.root.setLevel(logging.INFO)

    defaultParams = parameters.minimaHoppingParameters()

    parser = argparse.ArgumentParser(description ='Read two input files and calculate the overlap matrix fingerprint distance between them.')

    parser.add_argument('-i1', '--inputfile1', dest ='file1',
                action ='store', help ='input filename of first file.', required=True)
    parser.add_argument('-i2', '--inputfile2', dest ='file2',
                action ='store', help ='input filename of second file.', required=True)
    # parse indeices
    parser.add_argument('--index1', dest='index1', action='store', type=int,
        help='index of first structure in input file', default=0, required= False)
    parser.add_argument('--index2', dest='index2', action='store', type=int,
        help='index of second structure in input file', default=0, required= False)

    # parse additional OMFP parameters
    parser.add_argument('--n_S_orbitals', dest='n_S_orbitals', action='store', type=int,
        help='number of s orbitals for constructing the OMFP', default=defaultParams.n_S_orbitals, required= False)
    parser.add_argument('--n_P_orbitals', dest='n_P_orbitals', action='store', type=int,
        help='number of p orbitals for constructing the OMFP', default=defaultParams.n_P_orbitals, required= False)
    parser.add_argument('--width_cutoff', dest='width_cutoff', action='store', type=float,
        help='cutoff for the OMFP', default=defaultParams.width_cutoff, required= False)
    parser.add_argument('--maxnatsphere', dest='maxnatsphere', action='store', type=int,
        help='Truncation lentgh of the OMFP length', default=defaultParams.maxnatsphere, required= False)
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

    atom1 = read(args.file1, index=args.index1)
    atom2 = read(args.file2, index=args.index2)
    
    if '.ascii' in args.file1:
        atom1.pbc = [True, True, True]
    if '.ascii' in args.file2:
        atom2.pbc = [True, True, True]

    m1 = minimum.Minimum(atom1, 0.0, args.n_S_orbitals, args.n_P_orbitals, args.width_cutoff, args.maxnatsphere,
        T=0, ediff=0, exclude=args.exclude, calculate_fingerprint=True)
    m2 = minimum.Minimum(atom2, 0.0, args.n_S_orbitals, args.n_P_orbitals, args.width_cutoff, args.maxnatsphere,
        T=0, ediff=0, exclude=args.exclude, calculate_fingerprint=True)

    fingerprintDistance = m1.fingerprint_distance(m2)
    print(fingerprintDistance)


if __name__ =='__main__':
    main()
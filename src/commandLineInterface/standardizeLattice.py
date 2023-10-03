import argparse
from ase.io import read, write
import spglib
from ase.symbols import Symbols
from ase.atoms import Atoms

def main():

    parser = argparse.ArgumentParser(description ='Read an input files and standardizes periodic cell.')

    parser.add_argument('-i', '--inputfile', dest ='infile',
                action ='store', help ='input filename of first file.', required=True)
    parser.add_argument('-o', '--outputfile', dest ='outfile',
                action ='store', help ='output filename', required=True)

    parser.add_argument('--index', dest='index', action='store', type=int,
        help='index of structure in input file. Default is 0', default=0, required= False)
    parser.add_argument('--spglib_tol', dest='spglib_tol', action='store', type=float,
        help='Tolerance of spglib. Default is 1e-5', default=1e-5, required= False)
    parser.add_argument('--primitive', help='Find primitive cell', action='store_true')
    parser.add_argument('--spacegroup', help='Find and print spacegroup', action='store_true')

    args = parser.parse_args()

    findPrimitive = args.primitive
    findSpacegroup = args.spacegroup

    atoms = read(args.infile, index=args.index)

    atoms.pbc = [True, True, True]


    lattice, scaled_positions, numbers = spglib.standardize_cell(atoms, to_primitive=findPrimitive, no_idealize=False, symprec=args.spglib_tol, angle_tolerance=-1.0)


    if findPrimitive:
        symbs = Symbols(numbers)
        atoms = Atoms(symbs, scaled_positions=scaled_positions, cell=lattice)
    else:
        atoms.set_cell(lattice)
        atoms.set_scaled_positions(scaled_positions)

    atoms.pbc = [True, True, True]

    write(args.outfile, atoms)

    if findSpacegroup:
        print(spglib.get_spacegroup(atoms, symprec=args.spglib_tol, angle_tolerance=-1.0, symbol_type=0))


if __name__ =='__main__':
    main()
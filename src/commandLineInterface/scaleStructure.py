import argparse
from ase.io import read, write

def main():

    parser = argparse.ArgumentParser(description ='Read an input file and scale their coordinates and lattice vectors.')

    parser.add_argument('-i', '--inputfile', dest ='file',
                action ='store', help ='input filename of first file.', required=True)
    parser.add_argument('-o', '--outputfile', dest ='outfile', default=None,
                action ='store', help ='output filename', required=False)
    parser.add_argument('-s', '--scale', dest='scale', action='store', type=float,
        help='scale, default is 1', default=1, required= False)

    parser.add_argument('--index', dest='index', action='store', type=int,
        help='index of structure in input file. Default is 0', default=0, required= False)

    args = parser.parse_args()

    if args.outfile == None:
        args.outfile = args.file

    atoms = read(args.file, index = args.index)

    atoms.set_positions(atoms.get_positions() * args.scale)

    try:
        atoms.set_cell(atoms.get_cell() * args.scale)
    except:
        pass

    write(args.outfile, atoms)


    
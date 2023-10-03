import argparse
from ase.io import read, write

def main():
    parser = argparse.ArgumentParser(description ='Read a number of input files, sort them by energy and then write them back into an output file')
    parser.add_argument('-i', '--inputfile', dest ='inputFilename',
                action ='store', help ='input filename.', required=True)
    parser.add_argument('-o', '--outputfile', dest ='outputFilename',
                action ='store', help ='output filename.', required=True)

    parser.add_argument('--index', dest='index', action='store', type=str,
        help='index of structure in input file. Default is :', default=':', required= False)
    
    args = parser.parse_args()

    atomList = read(args.inputFilename, index = args.index)
    atomList.sort(key= lambda x: x.info['energy'])

    write(args.outputFilename, atomList)


if __name__ =='__main__':
    main()
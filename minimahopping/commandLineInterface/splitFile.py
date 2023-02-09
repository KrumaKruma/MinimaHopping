import argparse
from ase.io import read, write

def main():
    parser = argparse.ArgumentParser(description ='Read a long input files and split them up into many small files')
    parser.add_argument('-i', '--inputfile', dest ='inputFilename',
                action ='store', help ='input filename.', required=True)
    parser.add_argument('-m', '--max', dest='max', action='store', help='Maximal number of files',
        required=False, type=int, default=1000)
    
    args = parser.parse_args()

    atomList = read(args.inputFilename, index = ':')

    for i in range(args.max):
        fname = f"{i:07d}.ascii"
        write(fname, atomList[i])


if __name__ =='__main__':
    main()
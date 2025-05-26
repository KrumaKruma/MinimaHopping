import argparse
from ase.io import read, write, formats
import numpy as np

def main():
    parser = argparse.ArgumentParser(description ='Read a long input files and split them up into many small files')
    parser.add_argument('-i', '--inputfile', dest ='inputFilename',
                action ='store', help ='input filename.', required=True)
    parser.add_argument('-m', '--max', dest='max', action='store', help='Maximal number of files',
        required=False, type=int, default=-1)
    parser.add_argument('-f', '--format', dest='format', action='store', help='Output file format',
        required=False, default='.extxyz')
    
    args = parser.parse_args()

    # Remove the leading dot from the format string
    file_format = args.format.lstrip('.')

    # Check if the specified format is supported
    if file_format not in formats.ioformats.keys():
        raise ValueError(f"Unsupported file format: {args.format}. Please choose a format supported by ASE.")


    atomList = read(args.inputFilename, index = ':')

    mymax = args.max
    if mymax < 0:
        mymax = len(atomList)
    else:
        mymax = min(mymax, len(atomList))
   
    for i in range(mymax):
        fname = f"{i:07d}.{file_format}"
        write(fname, atomList[i])


if __name__ =='__main__':
    main()
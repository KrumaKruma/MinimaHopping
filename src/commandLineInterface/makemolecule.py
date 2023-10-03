import numpy as np
from ase.io import read, write
from sklearn.metrics import pairwise_distances
#from ase import units
import minimahopping.mh.periodictable as periodictable
import argparse


def get_rcovs(elements):
    rcovs = []
    for element in elements:
        rcovs.append(periodictable.getRcov_n(element))
    rcovs = np.array(rcovs)
    return rcovs


def get_minimal_pairwise_distances(atoms):
    # Initializations
    nat = len(atoms)
    lattice = atoms.get_cell() 
    positions = atoms.get_positions() 
    minimal_distances = np.ones((nat,nat)) * np.inf

    for ix in range(-1,2,1):
        for iy in range(-1,2,1):
            for iz in range(-1,2,1):
                positions_w = positions - ix*lattice[0,:]  - iy*lattice[1,:] - iz*lattice[2,:]
                distances = pairwise_distances(positions, positions_w)
                minimal_distances = np.minimum(minimal_distances, distances)

    assert np.inf not in minimal_distances, "error in minimal distances. Found infinite distance"

    return minimal_distances


def write_log(molecule_index, molecule_size):
    log_msg = "Molecule {:d} has {:d} atoms".format(molecule_index, molecule_size)
    print(log_msg)


def write_structures(atoms, molecule_atoms, molecule_index):
    filename = "molecule" + str(molecule_index).zfill(3) + ".extxyz"
    write(filename, atoms[molecule_atoms])


def get_molecules(atoms, distances, rcovs, factor_cov, verbose = True):
    nat = len(atoms)
    number_of_molecules = 0
    molecule_index = np.zeros((nat), dtype = int)
    belongs_to = np.zeros(nat, dtype=bool)

    for iat in range(nat):
        if not belongs_to[iat]:
            belongs_to[iat] = True
            molecule_atoms = []
            molecule_atoms.append(iat)
            number_of_molecules += 1
            molecule_index[iat] = number_of_molecules
            molecule_size = 1
            for kat in range(nat):
                if belongs_to[kat] and molecule_index[kat] == number_of_molecules:
                    for jat in range(nat):
                        cutoff_distance = factor_cov*(rcovs[kat]+rcovs[jat])
                        if not belongs_to[jat] and distances[kat,jat] < cutoff_distance:
                            molecule_size = molecule_size + 1
                            belongs_to[jat] = True
                            molecule_index[jat] = number_of_molecules
                            molecule_atoms.append(jat)
            if verbose:
                write_log(number_of_molecules, molecule_size)
                write_structures(atoms, molecule_atoms, number_of_molecules)
    
    return belongs_to


def main():

    parser = argparse.ArgumentParser(description ='Splitting molecular crystals into different molecules')
    parser.add_argument('-i', '--inputfile', dest ='inputFilename',
                action ='store', help ='input filename.', required=True)
    parser.add_argument('--factor_cov', dest='factor_cov', 
                action='store', type=float, help='multiplication factor for detecting covalent bonds (factor_cov * (rcov1 + rcov2))', default=1.1, required= False)
    parser.add_argument('-o', '--outputfile', dest ='outputFilename',
                action ='store', help ='output filename (only the header of the file).', required=False)
    parser.add_argument('--index', dest='index', 
                action='store', type=int, help='index of structure in input file', default=0, required= False)
    
    args = parser.parse_args()

    atoms = read(args.inputFilename, index=args.index)

    # Get the pairwise distances of the atoms in atomic units
    distances = get_minimal_pairwise_distances(atoms)

    # Get the covalent radii of the elements
    elements = atoms.get_atomic_numbers()
    rcovs = get_rcovs(elements)

    # Initializations
    belongs_to = get_molecules(atoms, distances, rcovs, args.factor_cov)






if __name__ == '__main__':
    main()

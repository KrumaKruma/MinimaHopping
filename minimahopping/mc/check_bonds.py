import numpy as np
from scipy.spatial import distance_matrix

def get_mindist(atoms, threshold_distance):
    positions = atoms.get_positions()
    distances = distance_matrix(positions, positions)
    distances = distances + np.eye(len(atoms), len(atoms))*10
    min_distance = np.min(distances)

    if min_distance < threshold_distance:
        is_chrashed = True
    else:
        is_chrashed = False
    return is_chrashed
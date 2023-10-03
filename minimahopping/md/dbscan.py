from ase.io import read, write
from sklearn.cluster import DBSCAN
import minimahopping.mh.periodictable as periodictable
import numpy as np


def get_com(positions, mass):
    total_mass = np.sum(mass)
    mass_3d = np.vstack([mass] * 3).T
    weighted_positions = positions*mass_3d
    com = np.sum(weighted_positions, axis=0)
    com /= total_mass
    return com


def dbscan(eps, positions,):
    db = DBSCAN(eps=eps, min_samples=1).fit(positions)
    labels = db.labels_
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    assert n_noise == 0, "Some atoms in DBSCAN were recognized as noise"
    return labels, n_clusters

def get_eps(elements):
    rcovs = []
    for element in elements:
        rcovs.append(periodictable.getRcov_n(element))
    rcovs = np.array(rcovs)
    rcov_mean = np.mean(rcovs)
    eps = 2.*rcov_mean
    return eps


def one_cluster(positions, elements):
    eps = get_eps(elements)
    labels, n_clusters = dbscan(eps, positions)
    if n_clusters > 1:
        is_one_cluster = False
    else:
        is_one_cluster = True
    return is_one_cluster

    

def adjust_velocities(positions, velocities,elements, masses):
    com = get_com(positions, masses)
    eps = get_eps(elements)
    mass_3d = np.vstack([masses] * 3).T
    _e_kin = 0.5 * np.sum(mass_3d * velocities * velocities)
    _v_average = np.sqrt((2.*_e_kin)/(np.mean(masses) * velocities.shape[0]))
    shifts = np.zeros(positions.shape)
    labels, n_clusters = dbscan(eps, positions)
    if n_clusters > 1:
        for i in range(n_clusters):
            indices = np.where(labels == i)
            cluster_pos = positions[indices,:][0]
            cluster_mass = masses[indices]
            com_cluster = get_com(cluster_pos, cluster_mass)
            shift = (com - com_cluster)/np.linalg.norm(com - com_cluster)
            indices = np.where(labels == i)
            shifts[indices,:] = shift * _v_average

    return shifts





# def test():
#     atoms = read("acc.extxyz")
#     positions = atoms.get_positions()
#     elements = atoms.get_atomic_numbers()
#     masses = atoms.get_masses()
#     com = get_com(positions, masses)
#     rcovs = []
#     for element in elements:
#         rcovs.append(periodictable.getRcov_n(element))
#     rcovs = np.array(rcovs)
#     rcov_mean = np.mean(rcovs)
#     eps = 1.*rcov_mean
#     db = DBSCAN(eps=eps, min_samples=1).fit(positions)
#     labels = db.labels_
#     n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
#     n_noise = list(labels).count(-1)
#     assert n_noise == 0, "Some atoms in DBSCAN were recognized as noise"
#     if n_clusters > 1:
#         while n_clusters > 1:
#             for i in range(n_clusters):
#                 indices = np.where(labels == i)
#                 cluster_pos = positions[indices,:][0]
#                 cluster_mass = masses[indices]
#                 com_cluster = get_com(cluster_pos, cluster_mass)
#                 shift = (com - com_cluster)/np.linalg.norm(com - com_cluster)
#                 indices = np.where(labels == i)
#                 positions[indices,:] += 0.1 * shift
            
#             atoms.set_positions(positions)
#             write("toto.extxyz", atoms, append=True)
#             db = DBSCAN(eps=eps, min_samples=1).fit(positions)
#             labels = db.labels_
#             n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
#             n_noise = list(labels).count(-1)
#             assert n_noise == 0, "Some atoms in DBSCAN were recognized as noise"




# if __name__ == '__main__':
#     main()
import numpy as np
import scipy as sp
import pandas as pd
import scipy.spatial.distance as dist
from scipy import sparse
from scipy.spatial import cKDTree as KDTree


def RaDec2XYZ(ra, dec):
    """
    From (ra,dec) -> unit vector on the sphere
    """

    rar = np.radians(ra)
    decr = np.radians(dec)

    x = np.cos(rar) * np.cos(decr)
    y = np.sin(rar) * np.cos(decr)
    z = np.sin(decr)
    vec = np.array([x, y, z]).T

    return vec


def get_vec_distances(ra, dec, comoving_dist):
    '''
    Converts RA,DEC, and Comoving Distances into a single vector in 3D space
    '''
    vec_unit = RaDec2XYZ(ra, dec)
    vec_dist = (vec_unit.T * comoving_dist).T
    return vec_dist


def getPairs(data_frame, max_sep=20, results_loc='PAIRS_sparse_dist.npz', save_frame=False):
    '''
    Takes a data frame with RA, DEC, and COMOVING distances, and pairs
    up those vectors with the closest vector under a given maximum separation
    and store the results in an array where:
    results[0] = First Item
    results[1] = Second Item
    results[2] = Distance between First Item and Second Item
    '''
    assert(data_frame['RA'])
    assert(data_frame['DEC'])
    assert(data_frame['COMOVING'])

    vec_distances = []
    for index in range(len(data_frame)):
        vec_distances.append(get_vec_distances(
            data_frame['RA'][index], data_frame['DEC'][index], data_frame['COMOVING'][index]))

    data_frame['x_vec'] = pd.Series(np.transpose(vec_distances)[0])
    data_frame['y_vec'] = pd.Series(np.transpose(vec_distances)[1])
    data_frame['z_vec'] = pd.Series(np.transpose(vec_distances)[2])

    vec = [data_frame['x_vec'][0], data_frame['y_vec']
           [0], data_frame['z_vec'][0]]
    vec_frame = np.vstack(
        (data_frame['x_vec'], data_frame['y_vec'], data_frame['z_vec']))
    vec_frame = vec_frame.T

    frame_tree = KDTree(vec_frame)
    dist_matr = frame_tree.sparse_distance_matrix(
        frame_tree, max_distance=max_sep, p=2.0, output_type='ndarray')

    dist_u = dist_matr[dist_matr['i'] < dist_matr['j']]
    result = sparse.coo_matrix(
        (dist_u['v'], (dist_u['i'], dist_u['j'])), (len(dist_u), len(dist_u)))

    if save_frame == True:
        data_frame.to_csv('DES_REDMAGIC_Manipulated.csv')

    sp.sparse.save_npz(results_loc, result)

    pairs_matrix = np.vstack((dist_u['i'], dist_u['j'], dist_u['v']))

    return(pairs_matrix)

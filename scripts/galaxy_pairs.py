import numpy as np
import scipy as sp
import pandas as pd
import scipy.spatial.distance as dist
from scipy import sparse
from scipy.spatial import cKDTree as KDTree
import sptpol_software.observation as obs
from sptpol_software import *


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


def getPairs(data_frame, max_sep=20, results_loc='PAIRS_sparse_dist.npz', save_frame=False, output='DES_REDMAGIC_Manipulated.csv'):
    '''
    Takes a data frame with RA, DEC, and COMOVING distances, and pairs
    up those vectors with the closest vector under a given maximum separation
    and store the results in an array where:
    results[0] = First Item
    results[1] = Second Item
    results[2] = Distance between First Item and Second Item
    '''
    if 'RA' not in data_frame.columns:
        return("Error: RA not in Data Frame")
    if 'DEC' not in data_frame.columns:
        return("Error: DEC not in Data Frame")
    if 'COMOVING' not in data_frame.columns:
        return("Error: COMOVINGb not in Data Frame")

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
        data_frame.to_csv(output)

    sp.sparse.save_npz(results_loc, result)

    pairs_matrix = np.vstack((dist_u['i'], dist_u['j'], dist_u['v']))

    return(pairs_matrix)


def euclideanDistance(instance1, instance2, length):
    '''
    Returns the Euclidean Distance between two vectors of a given length
    '''
    distance = 0
    for x in range(length):
        distance += pow((instance1[x] - instance2[x]), 2)
    return np.sqrt(distance)


def get_subarray(array, centre, sqr_radius):
    '''
    Gets Sub Array with with an input centre and half width from an input Array
    '''
    x_cen = centre[0]
    y_cen = centre[1]
    sl_x = slice(x_cen - sqr_radius, x_cen + sqr_radius)
    sl_y = slice(y_cen - sqr_radius, y_cen + sqr_radius)

    return(array[sl_x, sl_y])


def get_midpoint(ra_dec_1, ra_dec_2):
    '''
    Find the midpoint between two points in array space
    '''
    pt1 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_1, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
    pt2 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_2, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))

    X1 = pt1[0]
    X2 = pt2[0]

    Y1 = pt1[1]
    Y2 = pt2[1]

    return(((X1 + X2) / 2, (Y1 + Y2) / 2))


def get_rotn_angle(ra_dec_1, ra_dec_2):
    '''
    Return angle needed to rotate array based on dot-product of ra_dec vectors
    '''
    pt1 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_1, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
    pt2 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_2, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))

    X1 = pt1[0]
    X2 = pt2[0]

    Y1 = pt1[1]
    Y2 = pt2[1]

    m = (Y2 - Y1) / (X2 - X1)

    return(np.arctan(m))


def stack(y_map, galaxy_catalogue, pairs):
    '''
    Take input Y-map, galaxy catalogue, and list of pairs, and stacks them on top of each other
    '''
    for index, row in pairs.iterrows():
        galaxy_1 = row['galaxy_index_1']
        galaxy_2 = row['galaxy_index_2']

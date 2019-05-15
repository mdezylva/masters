import numpy as np
import pandas as pd
import scipy as sp
import scipy.spatial.distance as dist
from scipy import sparse
from scipy.spatial import cKDTree as KDTree
import pdb
import astropy.units as u 
from mpdaf.obj import Image,WCS

import sptpol_software
import sptpol_software.observation as obs
from sptpol_software import *

import matplotlib.pyplot as plt


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

def getPairs(data_frame, max_sep=20, query_type = 1,results_loc='PAIRS_sparse_dist.npz', save_frame=False, output='DES_REDMAGIC_Manipulated.csv'):
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

    # vec = [data_frame['x_vec'][0], data_frame['y_vec']
    #    [0], data_frame['z_vec'][0]]
    vec_frame = np.vstack(
        (data_frame['x_vec'], data_frame['y_vec'], data_frame['z_vec']))
    vec_frame = vec_frame.T

    frame_tree = KDTree(vec_frame)
    if query_type == 1:
        dist_matr = frame_tree.sparse_distance_matrix(
            frame_tree, max_distance=max_sep, p=2.0, output_type='ndarray')
    elif query_type == 2:
            dist_matr = frame_tree.query_pairs(r=max_sep, p=2.0, output_type='ndarray')
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

def get_subarray(array, centre, sqr_radius,max_size = 120):
    '''
    Gets Sub Array with with an input centre in array space and half width from an input Array
    '''
    x_cen = centre[0]
    y_cen = centre[1]
    if len(array)< x_cen + sqr_radius:
        padded_array = np.pad(array,min(abs(np.shape(array)[0]-(x_cen + sqr_radius)),np.shape(array)[0]-max_size/2),mode="constant")   

        return(padded_array)
    else:    
        sl_x = slice(x_cen - sqr_radius, x_cen + sqr_radius)
        sl_y = slice(y_cen - sqr_radius, y_cen + sqr_radius)

        return(array[sl_x, sl_y])

def get_midpoint(ra_dec_1, ra_dec_2, debug = False):
    '''
    Find the midpoint between two points in array space
    '''
    pt1 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_1, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
    pt2 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_2, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))

    X1 = float(pt1[0][0][0])
    X2 = float(pt2[0][0][0])
    
    Y1 = float(pt1[0][1][0])
    Y2 = float(pt2[0][1][0])

    if debug:
        print("Distance Between Pairs = ", str(np.sqrt((X2-X1)**2 + (Y2 - Y1)**2)))
        pdb.set_trace()
    return((abs(X1 + X2) / 2, (abs(Y1 + Y2) / 2)))

def get_rotn_angle(ra_dec_1, ra_dec_2, debug = False):
    '''
    Return angle needed to rotate array based on gradient of ra_dec vectors
    '''
    pt1 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_1, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
    
    pt2 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_2, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))

    if debug:
        print("Point 1 = " + str(pt1))
        print("Point 2 = " + str(pt2))

    if (not pt1[1][0][0]) or (not pt1[1][1][0]) or (not pt2[1][0][0]) or (not pt2[1][1][0]):
        raise ValueError("Point not inside bounds of map")

    X1 = float(pt1[0][0][0])
    X2 = float(pt2[0][0][0])

    Y1 = float(pt1[0][1][0])
    Y2 = float(pt2[0][1][0])
    if debug:
        print("X1 = " + str(X1))
        print("X2 = " + str(X2))
        print("Y1 = " + str(Y1))
        print("Y2 = " + str(Y2))
    if (X2-X1) == 0:
        # print(X2-X1)
        # print(ra_dec_1)
        # print(ra_dec_2)
        # pdb.set_trace()
        return(90.)
    
    m = (Y2 - Y1) / (X2 - X1)
    if debug:
        print(m)
    return(np.degrees(np.arctan(m)))
 
def cut_out_pair(pair, y_map, galaxy_catalogue, sqr_radius ,debug = False):
    '''
    Takes an input pair and a Compton Y-Map, and extract the pair as a sub map
    '''
    first_point = pair[0]
    second_point = pair[1]

    ra_1 = galaxy_catalogue.loc[first_point]['RA']
    dec_1 = galaxy_catalogue.loc[first_point]['DEC']

    ra_2 = galaxy_catalogue.loc[second_point]['RA']
    dec_2 = galaxy_catalogue.loc[second_point]['DEC']

    point_1 = (ra_1, dec_1)
   # print(point_1)
    point_2 = (ra_2, dec_2)
   # print(point_2)

    midpoint = get_midpoint(point_1, point_2)
   # print(midpoint)
    midpoint = np.array(midpoint).astype(int)

    return(get_subarray(y_map, midpoint, sqr_radius))

def extract_ra_dec(galaxy_index,galaxy_catalogue):
    '''
    Takes input galaxy index from get_pairs, and returns corresponding RA and DEC
    '''
    ra = galaxy_catalogue.loc[galaxy_index]['RA']
    dec = galaxy_catalogue.loc[galaxy_index]['DEC']
    return((ra,dec))
    
def calc_array_scale_factor(ra_dec_1,ra_dec_2, scaled_coords):
    """
    Calculates the factor by which to scale the array
    """
    pt1 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_1, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
    pt2 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_2, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))

    X1 = float(pt1[0][0][0])
    X2 = float(pt2[0][0][0])
    
    Y1 = float(pt1[0][1][0])
    Y2 = float(pt2[0][1][0])

    sep = np.sqrt((X2-X1)**2 + (Y2 - Y1)**2)

    scale_fac = 100.0/sep
    
    return(scale_fac)

def rescale_array(array,ra_dec_1,ra_dec_2,sqr_radius,scaled_width=100.):
    '''
    Rescales Array 
    '''
    pt1 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_1, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
    pt2 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_2, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))

    X1 = float(pt1[0][0][0])
    X2 = float(pt2[0][0][0])
    
    Y1 = float(pt1[0][1][0])
    Y2 = float(pt2[0][1][0])

    sep = np.sqrt((X2-X1)**2 + (Y2 - Y1)**2)
    scale_fact = scaled_width/sep
    if sep == 0:
        return(0)
    rescaled_array = sp.ndimage.zoom(array,scale_fact)

    centre = [len(rescaled_array)/2,len(rescaled_array)/2]

    new_array = get_subarray(array=rescaled_array,centre = centre,sqr_radius = min(len(rescaled_array)/2,sqr_radius))

    return(new_array)

def stack_pairs(y_map, galaxy_catalogue, pairs, size_of_cutout=100, debug = False):
    '''
    Take input Y-map, galaxy catalogue, and list of pairs, and stacks them on top of each other
    returning a stacked array
    '''

    output = np.ndarray(shape = (size_of_cutout*2, size_of_cutout*2))
    if debug:
        print("Output Array Shape = " + str(np.shape(output)))

    for index, row in pairs.iterrows():
        galaxy_1 = row['galaxy_index_1']
        galaxy_2 = row['galaxy_index_2']
        pair = [galaxy_1, galaxy_2]

        cut_array = cut_out_pair(pair, y_map, galaxy_catalogue,30,debug=True)

        if debug:
            print("Cut Array Shape = " + str(np.shape(cut_array)))

        gal_1_coords = extract_ra_dec(galaxy_1, galaxy_catalogue)
        gal_2_coords = extract_ra_dec(galaxy_2, galaxy_catalogue)

        angle = get_rotn_angle(gal_1_coords, gal_2_coords)
        rot_array = sp.ndimage.rotate(cut_array, 90-angle, reshape=True)
        prev_cell = float()
        if debug:
            print("Rotated Array Shape = " + str(np.shape(rot_array)))
        scaled_array = rescale_array(rot_array,gal_1_coords,gal_2_coords,len(output)/2)
        if debug:
            print("Scaled Array Shape = " + str(np.shape(scaled_array)))
            if np.shape(scaled_array) != np.shape(output) and np.shape(scaled_array) != ():
                pdb.set_trace()
        output += scaled_array

        print("Added pair " + str(index))

    return(output)

def stack_pairs_V2(y_map, galaxy_catalogue, pairs, size_of_cutout=70, debug = False,save= False):
    '''
    Take input Y-map, galaxy catalogue, and list of pairs, and stacks them on top of each other returning a stacked array
    '''
    output = np.zeros(shape = (120,120))
    num_rejected = 0 
    if debug:
        print("Initial Output Shape : " + str(np.shape(output)))
    for index, row in pairs.iterrows():
        galaxy_1 = row['galaxy_index_1']
        galaxy_2 = row['galaxy_index_2']
        pair = [galaxy_1, galaxy_2]

        ra_1 = galaxy_catalogue.loc[galaxy_1]['RA']
        dec_1 = galaxy_catalogue.loc[galaxy_1]['DEC']

        ra_2 = galaxy_catalogue.loc[galaxy_2]['RA']
        dec_2 = galaxy_catalogue.loc[galaxy_2]['DEC']

        pt1 = sptpol_software.observation.sky.ang2Pix(
        (ra_1,dec_1), [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
        pt2 = sptpol_software.observation.sky.ang2Pix(
        (ra_2,dec_2), [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))

        X1 = float(pt1[0][0][0])
        X2 = float(pt2[0][0][0])
    
        Y1 = float(pt1[0][1][0])
        Y2 = float(pt2[0][1][0])

    

    
        midpoint = (int(abs(X1 + X2) / 2), int((abs(Y1 + Y2) / 2)))
        if debug:
            print("Midpoint = " +str(midpoint))
        cut_array = get_subarray(y_map, midpoint, 30)
        if debug:
            print("x1 = " + str(X1))
            # print("-x1 = " + str(len(cut_array)-X1))
            print("y1 = " + str(Y1))
            # print("-y1 = " + str(len(cut_array)-Y1))
            print("x2 = " + str(X2))
            # print("-x2 = " + str(len(cut_array)-X2))
            print("y2 = " + str(Y2))
            # print("-y2 = " + str(len(cut_array)-Y2))
            pdb.set_trace()
        if (float(X2)-float(X1) == 0):
            angle = 90
        else:
            angle = np.degrees(np.arctan((float(Y2)-float(Y1))/(float(X2)-float(X1))))

        if debug:
            print("Angle = " + str(angle) + ' or ' + str(90-angle))

        rot_array =  sp.ndimage.rotate(cut_array, 90-angle, reshape=True)

        sep = np.sqrt((X2-X1)**2 + (Y2 - Y1)**2)
        if sep < 2:
            num_rejected += 1
            continue
        scale_fac = 80.0/sep

        if debug:
            print("Separation = " + str(sep))
            print("Scale Factor = " + str(scale_fac))
            pdb.set_trace()
        
        rescaled_array = sp.ndimage.zoom(rot_array,scale_fac)

        centre = [len(rescaled_array)/2,len(rescaled_array)/2]

        

        re_cut_array = get_subarray(rescaled_array,centre,60)

        if debug:
            print("Output Shape: " + str(np.shape(output)))
            print("Re Cut Array:" + str(np.shape(re_cut_array)))
            pdb.set_trace()

        if np.shape(output) != np.shape(re_cut_array):
            num_rejected += 1
            continue

        if abs(np.mean(re_cut_array)) > 1e-4:
            num_rejected += 1
            continue

        output = np.add(output,re_cut_array)

        print("Added pair " + str(index))

        if index%1000 == 0 and save:
            plt.imshow(output)
            # plt.title("Output")
            # plt.show()
            # plt.imshow(re_cut_array)
            # plt.show()
            filename = "output_" + str(index) + ".png"
            plt.savefig(filename)
            # np.savetxt(filename,output,delimiter = ',')
            # pdb.set_trace()
    print("Number of Cut Pairs = " + str(num_rejected))
    return(output)

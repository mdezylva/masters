import numpy as np
import pandas as pd
import scipy as sp
import scipy.spatial.distance as dist
from scipy import sparse
from scipy.spatial import cKDTree as KDTree
import pdb
import astropy.units as u 
from mpdaf.obj import Image,WCS
from PIL import Image
import sptpol_software
import sptpol_software.observation as obs
from sptpol_software import *
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import healpy as hp


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

def normalize(num, lower=0, upper=360, b=False):
    """Normalize number to range [lower, upper) or [lower, upper].
    Parameters
    ----------
    num : float
        The number to be normalized.
    lower : int
        Lower limit of range. Default is 0.
    upper : int
        Upper limit of range. Default is 360.
    b : bool
        Type of normalization. Default is False. See notes.
    Returns
    -------
    n : float
        A number in the range [lower, upper) or [lower, upper].
    Raises
    ------
    ValueError
      If lower >= upper.
    Notes
    -----
    If the keyword `b == False`, then the normalization is done in the
    following way. Consider the numbers to be arranged in a circle,
    with the lower and upper ends sitting on top of each other. Moving
    past one limit, takes the number into the beginning of the other
    end. For example, if range is [0 - 360), then 361 becomes 1 and 360
    becomes 0. Negative numbers move from higher to lower numbers. So,
    -1 normalized to [0 - 360) becomes 359.
    If the keyword `b == True`, then the given number is considered to
    "bounce" between the two limits. So, -91 normalized to [-90, 90],
    becomes -89, instead of 89. In this case the range is [lower,
    upper]. This code is based on the function `fmt_delta` of `TPM`.
    Range must be symmetric about 0 or lower == 0.
    Examples
    --------
    >>> normalize(-270,-180,180)
    90.0
    >>> import math
    >>> math.degrees(normalize(-2*math.pi,-math.pi,math.pi))
    0.0
    >>> normalize(-180, -180, 180)
    -180.0
    >>> normalize(180, -180, 180)
    -180.0
    >>> normalize(180, -180, 180, b=True)
    180.0
    >>> normalize(181,-180,180)
    -179.0
    >>> normalize(181, -180, 180, b=True)
    179.0
    >>> normalize(-180,0,360)
    180.0
    >>> normalize(36,0,24)
    12.0
    >>> normalize(368.5,-180,180)
    8.5
    >>> normalize(-100, -90, 90)
    80.0
    >>> normalize(-100, -90, 90, b=True)
    -80.0
    >>> normalize(100, -90, 90, b=True)
    80.0
    >>> normalize(181, -90, 90, b=True)
    -1.0
    >>> normalize(270, -90, 90, b=True)
    -90.0
    >>> normalize(271, -90, 90, b=True)
    -89.0
    """
    from math import floor, ceil
    # abs(num + upper) and abs(num - lower) are needed, instead of
    # abs(num), since the lower and upper limits need not be 0. We need
    # to add half size of the range, so that the final result is lower +
    # <value> or upper - <value>, respectively.
    res = num
    if not b:
        if lower >= upper:
            raise ValueError("Invalid lower and upper limits: (%s, %s)" %
                             (lower, upper))

        res = num
        if num > upper or num == lower:
            num = lower + abs(num + upper) % (abs(lower) + abs(upper))
        if num < lower or num == upper:
            num = upper - abs(num - lower) % (abs(lower) + abs(upper))

        res = lower if num == upper else num
    else:
        total_length = abs(lower) + abs(upper)
        if num < -total_length:
            num += ceil(num / (-2 * total_length)) * 2 * total_length
        if num > total_length:
            num -= floor(num / (2 * total_length)) * 2 * total_length
        if num > upper:
            num = total_length - num
        if num < lower:
            num = -total_length - num

        res = num

    res *= 1.0  # Make all numbers float, to be consistent

    return res

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

def getUnPairs(data_frame, max_sep=20):
    '''
    Takes a data frame with RA, DEC, and COMOVING distances, and pairs
    up those vectors with the closest vector under a given maximum separation
    and store the results in an array where:
    results[0] = First Item
    results[1] = Second Item
    results[2] = Distance between First Item and Second Item
    '''
    if 'NORM_RA' not in data_frame.columns:
        return("Error: Normalised (-180,180] RA not in Data Frame")
    if 'DEC' not in data_frame.columns:
        return("Error: DEC not in Data Frame")

    np.asarray([df['NORM_RA'].as_matrix().T,df['DEC'].as_matrix().T]).T

    frame_tree = scipy.spatial.cKDTree(sky_locs, leafsize=100)

    for item in frame_tree:
        Result = YourTreeName.query(item, k=1, distance_upper_bound=6)

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
        if len(array) < 0:
            return(np.zeros(shape=(max_size,max_size) ))

        padded_array = np.pad(array,(x_cen+sqr_radius)-len(array),mode="constant")   
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
    output = np.zeros(shape = (169,169))
    num_rejected = 0 

    sep_cut = 0
    shape_cut = 0
    mean_cut = 0
    pairs = pairs.reset_index()

    cuts = []
    rotations  = []
    rescales = []
    recuts = []
 
    if debug:
        print("Initial Output Shape : " + str(np.shape(output)))
    for index, row in pairs.iterrows():
        galaxy_1 = row['galaxy_index_1'].astype('int32')
        galaxy_2 = row['galaxy_index_2'].astype('int32')
        pair = [galaxy_1, galaxy_2]

        ra_1 = galaxy_catalogue.loc[galaxy_1]['RA']
        dec_1 = galaxy_catalogue.loc[galaxy_1]['DEC']

        ra_2 = galaxy_catalogue.loc[galaxy_2]['RA']
        dec_2 = galaxy_catalogue.loc[galaxy_2]['DEC']

        pt1 = sptpol_software.observation.sky.ang2Pix(
        (ra_1,dec_1), [0, -57.5], reso_arcmin=0.25, map_pixel_shape=np.array([6000, 11280]))
        pt2 = sptpol_software.observation.sky.ang2Pix(
        (ra_2,dec_2), [0, -57.5], reso_arcmin=0.25, map_pixel_shape=np.array([6000, 11280]))

        X1 = float(pt1[0][0][0])
        X2 = float(pt2[0][0][0])
    
        Y1 = float(pt1[0][1][0])
        Y2 = float(pt2[0][1][0])

	midpoint = (int(abs(X1 + X2) / 2), int((abs(Y1 + Y2) / 2)))
	#points = SkyCoord([(str(ra_1) + " " + str(dec_1)), (str(ra_2) + " " + str(dec_2))], unit=(u.deg, u.deg))        
	#midpoint = SkyCoord(points.data.mean(), representation='unitspherical', frame=points)
	#midpoint_coord = [midpoint.ra.degree,midpoint.dec.degree]
	
	#midpoint = [np.mean([ra_1,ra_2]),np.mean([dec_1,dec_2])]
		
	if debug:
            print("Midpoint = " +str(midpoint))
        cut_array = get_subarray(y_map, midpoint, 200)

	if debug:
            print("x1 = " + str(X1))
            # print("-x1 = " + str(len(cut_array)-X1))
            print("y1 = " + str(Y1))
            # print("-y1 = " + str(len(cut_array)-Y1))
            print("x2 = " + str(X2))
            # print("-x2 = " + str(len(cut_array)-X2))
            print("y2 = " + str(Y2))
            # print("-y2 = " + str(len(cut_array)-Y2))
            #pdb.set_trace()
        if (float(X2)-float(X1) == 0):
            angle = 90
        else:
            angle = np.degrees(np.arctan((float(Y2)-float(Y1))/(float(X2)-float(X1))))

        if debug:
            print("Angle = " + str(angle) + ' or ' + str(90-angle))

        rot_array =  sp.ndimage.rotate(cut_array,90-angle, reshape=True)
       

	sep = np.sqrt((X2-X1)**2 + (Y2 - Y1)**2)
        #sep_cut = 0
	#ang_sep = points[0].separation(points[1]).arcmin
	
	#if sep < 2:
        #    sep_cut += 1
        #    continue
        scale_fac = 80.0/sep

        
        
	if debug:
            print("Separation = " + str(sep))
            print("Scale Factor = " + str(scale_fac))
            pdb.set_trace()
        
        rescaled_array = sp.ndimage.zoom(rot_array,scale_fac)

        centre = [len(rescaled_array)/2,len(rescaled_array)/2]

        

        re_cut_array = get_subarray(rescaled_array,centre,85)


        if debug:
            print("Output Shape: " + str(np.shape(output)))
            print("Re Cut Array:" + str(np.shape(re_cut_array)))
            pdb.set_trace()

	#shape_cut = 0
        if np.shape(output) != np.shape(re_cut_array):

	    re_cut_array = re_cut_array[:-1,:-1]
            
	    #print("Pair No: " + str(index))
	    #print("Output Shape: ")
	    #print(np.shape(output))
	    #print("Re Cut Array Shape: ")
	    #print(np.shape(re_cut_array))
	    #pdb.set_trace()
	    if np.shape(output) != np.shape(re_cut_array):
	    	shape_cut += 1
            	continue

	#mean_cut = 0
        if abs(np.mean(re_cut_array)) > 1e-5:
            mean_cut += 1
            continue

        flipped = np.fliplr(re_cut_array)
        av = np.add(flipped,re_cut_array)/2
        #flipped2 = np.flipud(re_cut_array)
	output = np.add(output,av)
        #output = np.add(output,flipped2)

        #if index%10000 == 0:
        #    print("Added pair " + str(index))
	#plt.switch_backend('Agg')
	#print("Added pair " + str(index))
	if index < 100:
	    
	    cuts.append(cut_array)
	    rescales.append(rescaled_array)
	    recuts.append(re_cut_array)
	    rotations.append(rot_array)
#
	    plt.imshow(output)
	    filename = "fh_"+str(index)+"_pairs"
	    plt.colorbar()
	    plt.title(filename)
	    plt.savefig(filename+'.png')
	    plt.close()
	    centre = len(output)/2+1
	    central_line = output[centre,:]
	    plt.plot(central_line)
	    plt.title(filename)
	    plt.savefig(filename+"_central_slice.png")
        if index%1000 == 0 and save:
            #plt.switch_backend('Agg')
	    plt.imshow(output/index)
            print("Added pair " + str(index))
	    filename = 'output_' + str(index) + '_pairs'
            plt.title("Output")
	    plt.colorbar()
            # plt.imshow(re_cut_array)
            # plt.show()
            # filename = "output_" + str(index) + ".png"
            plt.savefig(filename+".png")
	    plt.close()
	    centre = len(output)/2+1
	    central_line = output[centre,:]/index
	    plt.plot(central_line)
	    plt.title("Central Line: "+str(index))
	    plt.savefig(filename+"_cen_slice.png")
            np.savetxt(filename+".txt",output,delimiter = ',')
            # pdb.set_trace()
    print("Number of Cut Pairs (Sep) = " + str(sep_cut))
    print("Number of Cut Pairs (Shape) = " + str(shape_cut))
    print("Number of Cut Pairs (Mean) = " + str(mean_cut))
    dic = {}
    dic['cuts'] = np.asarray(cuts)
    dic['rotated'] = np.asarray(rotations)
    dic['sep'] = np.asarray(rescales)
    dic['recuts'] = np.asarray(recuts)
    return(output,dic)

def stack_pairs_hp(y_map, galaxy_catalogue, pairs, size_of_cutout=70, debug = False,save= False):
    '''
    Take input Y-map in HEALPIX format, galaxy catalogue, and list of pairs, and stacks them on top of each other returning a stacked array
    '''
    output = np.zeros(shape = (169,169))
    num_rejected = 0 
    #Keep track of how many pairs are cut for what reason
    sep_cut = 0
    shape_cut = 0
    mean_cut = 0
    #Ensure the index of the DF isn't out of order due to cuts made 
    pairs = pairs.reset_index()
 
    if debug:
        print("Initial Output Shape : " + str(np.shape(output)))
    # Iterate over all indexes and rows
    for index, row in pairs.iterrows():
        galaxy_1 = row['galaxy_index_1'].astype('int32')
        galaxy_2 = row['galaxy_index_2'].astype('int32')
        pair = [galaxy_1, galaxy_2]
	
	# RAs and DECs are stored in the Row Now
        ra_1 = row['RA_1']
        dec_1 = row['DEC_1']

        ra_2 = row['RA_2']
        dec_2 = row['DEC_2']

	#Operations to find the midpoint coordinate. Might need to move this into the Pair Creation Scripts
	#Currently doesn't work due to Software Versioning on CLOUD!

	points = SkyCoord([(str(ra_1) + " " + str(dec_1)), (str(ra_2) + " " + str(dec_2))], unit=(u.deg, u.deg))        
	#midpoint = SkyCoord(points.data.mean(), representation='unitspherical', frame=points)
	#midpoint_coord = [midpoint.ra.degree,midpoint.dec.degree]
	
	#Naive approach	
	midpoint = [np.mean([ra_1,ra_2]),np.mean([dec_1,dec_2])]
	# Work out the rotation angle before cutting out from HEALPIX array	
	if debug:
            print("RA 1 = " + str(ra_1))
            # print("-x1 = " + str(len(cut_array)-X1))
            print("DEC 1 = " + str(dec_1))
            # print("-y1 = " + str(len(cut_array)-Y1))
            print("RA 2 = " + str(ra_2))
            # print("-x2 = " + str(len(cut_array)-X2))
            print("DEC 2 = " + str(dec_2))
            # print("-y2 = " + str(len(cut_array)-Y2))
            #pdb.set_trace()
        if (float(ra_2)-float(ra_1) == 0):
            angle = 90
        else:
            angle = np.degrees(np.arctan((float(dec_2)-float(dec_1))/(float(ra_2)-float(ra_1))))

        if debug:
            print("Angle = " + str(angle) + ' or ' + str(90-angle))

        #rot_array =  sp.ndimage.rotate(cut_array,90-angle, reshape=True)

	#Append Angle to Midpoint Rotation variable
	midpoint.append(angle)

	if debug:
            print("Midpoint = " +str(midpoint))
	    pdb.set_trace()
	
	#Cut out array from HEALPIX Array, including midpoint and rotation angle
        cut_array_masked = hp.gnomview(y_map, rot=midpoint,  xsize=100, reso=1, return_projected_map=True)
	cut_array = cut_array_masked.data
	plt.close()	

	if debug:
	    pdb.set_trace()
	
	#Calculate the Angular Separation of the two points
	ang_sep = points[0].separation(points[1]).arcmin
	
	#Rescale the HALOs so they are 60 arcmins apart
        scale_fac = 60.0/ang_sep
        
	if debug:
            print("Separation = " + str(ang_sep))
            print("Scale Factor = " + str(scale_fac))
            pdb.set_trace()
        
	#Rescale Array
        rescaled_array = sp.ndimage.zoom(cut_array,scale_fac)

        centre = [len(rescaled_array)/2,len(rescaled_array)/2]

        #Resclice for adding arrays
        re_cut_array = get_subarray(rescaled_array,centre,85)

        if debug:
            print("Output Shape: " + str(np.shape(output)))
            print("Re Cut Array:" + str(np.shape(re_cut_array)))

        if np.shape(output) != np.shape(re_cut_array):

	    re_cut_array = re_cut_array[:-1,:-1]
            
	    if np.shape(output) != np.shape(re_cut_array):
	    	shape_cut += 1
            	continue

        if abs(np.mean(re_cut_array)) > 1e-5:
            mean_cut += 1
            continue

        output = np.add(output,re_cut_array)
        flipped = np.fliplr(re_cut_array)
        flipped2 = np.flipud(re_cut_array)
	output = np.add(output,flipped)
        output = np.add(output,flipped2)

	plt.switch_backend('Agg')
	
	#Processing operations
        if index%1000 == 0 and save:
	    plt.imshow(output)
            print("Added pair " + str(index))
	    filename = 'output_' + str(index) + '_pairs'
            plt.title("Output")
            plt.savefig(filename+".png")
	    plt.close()
	    centre = len(output)/2+1
	    central_line = output[centre,:]
	    plt.plot(central_line)
	    plt.title("Central Line: "+str(index))
	    plt.savefig(filename+"_cen_slice.png")
            np.savetxt(filename+".txt",output,delimiter = ',')
            # pdb.set_trace()
    print("Number of Cut Pairs (Sep) = " + str(sep_cut))
    print("Number of Cut Pairs (Shape) = " + str(shape_cut))
    print("Number of Cut Pairs (Mean) = " + str(mean_cut))

    return(output)
def null_stack_rand_rot(y_map, galaxy_catalogue, pairs, size_of_cutout=70, debug = False,save= False):
    '''
    Take input Y-map, galaxy catalogue, and list of pairs, and stacks them on top of each other returning a stacked array
    '''
    output = np.zeros(shape = (169,169))
    num_rejected = 0 

    sep_cut = 0
    shape_cut = 0
    mean_cut = 0
    pairs = pairs.reset_index()

    cuts = []
    rotations  = []
    rescales = []
    recuts = []
 
    if debug:
        print("Initial Output Shape : " + str(np.shape(output)))
    for index, row in pairs.iterrows():
        galaxy_1 = row['galaxy_index_1'].astype('int32')
        galaxy_2 = row['galaxy_index_2'].astype('int32')
        pair = [galaxy_1, galaxy_2]

        ra_1 = galaxy_catalogue.loc[galaxy_1]['RA']
        dec_1 = galaxy_catalogue.loc[galaxy_1]['DEC']

        ra_2 = galaxy_catalogue.loc[galaxy_2]['RA']
        dec_2 = galaxy_catalogue.loc[galaxy_2]['DEC']

        pt1 = sptpol_software.observation.sky.ang2Pix(
        (ra_1,dec_1), [0, -57.5], reso_arcmin=0.25, map_pixel_shape=np.array([6000, 11280]))
        pt2 = sptpol_software.observation.sky.ang2Pix(
        (ra_2,dec_2), [0, -57.5], reso_arcmin=0.25, map_pixel_shape=np.array([6000, 11280]))

        X1 = float(pt1[0][0][0])
        X2 = float(pt2[0][0][0])
    
        Y1 = float(pt1[0][1][0])
        Y2 = float(pt2[0][1][0])

	midpoint = (int(abs(X1 + X2) / 2), int((abs(Y1 + Y2) / 2)))
	#points = SkyCoord([(str(ra_1) + " " + str(dec_1)), (str(ra_2) + " " + str(dec_2))], unit=(u.deg, u.deg))        
	#midpoint = SkyCoord(points.data.mean(), representation='unitspherical', frame=points)
	#midpoint_coord = [midpoint.ra.degree,midpoint.dec.degree]
	
	#midpoint = [np.mean([ra_1,ra_2]),np.mean([dec_1,dec_2])]
		
	if debug:
            print("Midpoint = " +str(midpoint))
        cut_array = get_subarray(y_map, midpoint, 200)

	if debug:
            print("x1 = " + str(X1))
            # print("-x1 = " + str(len(cut_array)-X1))
            print("y1 = " + str(Y1))
            # print("-y1 = " + str(len(cut_array)-Y1))
            print("x2 = " + str(X2))
            # print("-x2 = " + str(len(cut_array)-X2))
            print("y2 = " + str(Y2))
            # print("-y2 = " + str(len(cut_array)-Y2))
            #pdb.set_trace()

        angle = np.random.rand()*90.0
        if debug:
            print("Angle = " + str(angle) + ' or ' + str(90-angle))

        rot_array =  sp.ndimage.rotate(cut_array,90-angle, reshape=True)
       

	sep = np.sqrt((X2-X1)**2 + (Y2 - Y1)**2)
        #sep_cut = 0
	#ang_sep = points[0].separation(points[1]).arcmin
	
	#if sep < 2:
        #    sep_cut += 1
        #    continue
        scale_fac = 80.0/sep

        
        
	if debug:
            print("Separation = " + str(sep))
            print("Scale Factor = " + str(scale_fac))
            pdb.set_trace()
        
        rescaled_array = sp.ndimage.zoom(rot_array,scale_fac)

        centre = [len(rescaled_array)/2,len(rescaled_array)/2]

        

        re_cut_array = get_subarray(rescaled_array,centre,85)


        if debug:
            print("Output Shape: " + str(np.shape(output)))
            print("Re Cut Array:" + str(np.shape(re_cut_array)))
            pdb.set_trace()

	#shape_cut = 0
        if np.shape(output) != np.shape(re_cut_array):

	    re_cut_array = re_cut_array[:-1,:-1]
            
	    #print("Pair No: " + str(index))
	    #print("Output Shape: ")
	    #print(np.shape(output))
	    #print("Re Cut Array Shape: ")
	    #print(np.shape(re_cut_array))
	    #pdb.set_trace()
	    if np.shape(output) != np.shape(re_cut_array):
	    	shape_cut += 1
            	continue

	#mean_cut = 0
        if abs(np.mean(re_cut_array)) > 1e-5:
            mean_cut += 1
            continue

        flipped = np.fliplr(re_cut_array)
        av = np.add(flipped,re_cut_array)/2
        #flipped2 = np.flipud(re_cut_array)
	output = np.add(output,av)
        #output = np.add(output,flipped2)

        #if index%10000 == 0:
        #    print("Added pair " + str(index))
	#plt.switch_backend('Agg')
	#print("Added pair " + str(index))
	if index < 100:
	    
	    cuts.append(cut_array)
	    rescales.append(rescaled_array)
	    recuts.append(re_cut_array)
	    rotations.append(rot_array)
#
	    plt.imshow(output)
	    filename = "fh_"+str(index)+"_pairs"
	    plt.colorbar()
	    plt.title(filename)
	    plt.savefig(filename+'.png')
	    plt.close()
	    centre = len(output)/2+1
	    central_line = output[centre,:]
	    plt.plot(central_line)
	    plt.title(filename)
	    plt.savefig(filename+"_central_slice.png")
        if index%1000 == 0 and save:
            #plt.switch_backend('Agg')
	    plt.imshow(output/index)
            print("Added pair " + str(index))
	    filename = 'output_' + str(index) + '_pairs'
            plt.title("Output")
	    plt.colorbar()
            # plt.imshow(re_cut_array)
            # plt.show()
            # filename = "output_" + str(index) + ".png"
            plt.savefig(filename+".png")
	    plt.close()
	    centre = len(output)/2+1
	    central_line = output[centre,:]/index
	    plt.plot(central_line)
	    plt.title("Central Line: "+str(index))
	    plt.savefig(filename+"_cen_slice.png")
            np.savetxt(filename+".txt",output,delimiter = ',')
            # pdb.set_trace()
    print("Number of Cut Pairs (Sep) = " + str(sep_cut))
    print("Number of Cut Pairs (Shape) = " + str(shape_cut))
    print("Number of Cut Pairs (Mean) = " + str(mean_cut))
    dic = {}
    dic['cuts'] = np.asarray(cuts)
    dic['rotated'] = np.asarray(rotations)
    dic['sep'] = np.asarray(rescales)
    dic['recuts'] = np.asarray(recuts)
    return(output,dic)
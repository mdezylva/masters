#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'Notebooks'))
	print(os.getcwd())
except:
	pass

#%%
import os 
#data_location = input("Enter Path Location of data")
os.chdir("/home/mitchell/Documents/masters/masters/scripts/")
import galaxy_pairs
import cmb
os.chdir("/home/mitchell/Documents/masters/masters/data")
cwd = os.getcwd()
print(cwd)


#%%
import numpy as np
import scipy as sp
import astropy as ap
import glob
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.constants as const
from astropy import constants as ap_const
import scipy.ndimage 
import PIL
import sptpol_software as sps
from astropy.io import fits
from astropy.cosmology import Planck15
from astropy import units as u
import scipy.spatial.distance as dist
print(Planck15)
cosmo= Planck15


#%%
from sptpol_software.util.tools import stat
from sptpol_software.observation import *
import sptpol_software.observation as obs
import sptpol_software as sps
import sptpol_software.observation.sky
from sptpol_software.util import files


#%%
T_cmb = 2.725
freqs_ghz = [93.2000,147.700]
beam_norm_correction = [1./0.99673, 1./0.99470, 1./1. ] 
cal_factors = [0.9097,0.7765] # no 220 is why last is zero. 
#pol_cal_factors = pol_cal_factors_800 * beam_norm_correction


#%%
glob.glob('*.h5')


#%%
map_150ghz = files.read("ra0dec-57p5_sum5000_150ghz.h5")
map_90ghz = files.read("ra0dec-57p5_sum5000_090ghz.h5")


#%%
# map_150ghz = pol_cal_factors_800[0]*map_150ghz
# map_90ghz = pol_cal_factors_800[1]*map_90ghz
# map_150ghz.mapper_arguments


#%%
# map_150ghz.getSubmap([5,5], center_offset=[0,0], units='degree').drawImage(bw=False,vmax=0)
# map_90ghz.getSubmap([10,10], center_offset=[0,0], units='degree').drawImage(bw=False)


#%%
# diff_T = map_150ghz-map_90ghz
# diff_T.getSubmap([5,5], center_offset=[1,1], units='degree').drawImage(bw=False)


#%%
# def convert_ghz_to_y(freq_power):
#     x = freq_power/56.85  # x = h v / k_B T_CMB
#     return((x/np.tanh(x/2)) - 4)
print(cmb.convert_ghz_to_y(90))
print(cmb.convert_ghz_to_y(150))


#%%
freq_range = np.arange(20,1000,1)
x = freq_range/56.85
y_param = [(ex/np.tanh(ex/2) - 4) for ex in x]
plt.semilogx(freq_range,y_param,label='Y')
I_v = [(ex**4 * np.exp(ex))/(np.exp(ex)-1)**2 * (ex/np.tanh(ex/2) - 4) for ex in x]
plt.semilogx(freq_range,I_v,label='Intensity')
plt.axhline(0, color='black')
plt.title('Frequency Dependance of SZ Effect')
plt.legend()


#%%
# map_150ghz_array = map_150ghz.getTOnly().map
# print(type(map_150ghz))


#%%
# freq_scaling = convert_ghz_to_y(150)-convert_ghz_to_y(90)
# print(freq_scaling)
# freq_scaling_fac = 1.0/freq_scaling
# print(freq_scaling_fac)


#%%
#y_map = sptpol_software.observation.sky.Map
# y_map = freq_scaling_fac*(map_150ghz-map_90ghz)
#/(convert_ghz_to_y(150)-convert_ghz_to_y(90))


#%%
# y_map.getSubmap([5,5], center_offset=[0,0], units='degree').drawImage(bw=True)
# print("Functionally a factor of 2 conversion, look at scale")
# diff_T.getSubmap([2,2], center_offset=[0,0], units='degree').drawImage(bw=False)
y_map_array = cmb.get_y_map([map_150ghz,map_90ghz],cal_factors,freqs_ghz)


#%%
# y_map.drawImage(bw=False,vmax=-0.00001,vmin=0.00001)


#%%
# np.shape(y_map)


#%%
# y_map.writeToHDF5('y_map.h5',overwrite=True,use_compression=False)
# y_map_array = y_map.getTOnly().map
print(y_map_array)
# np.save('y_map', y_map_array)


#%%
plt.imshow(y_map_array,vmax=0.00001,vmin=-0.00001)


#%%
# red_mapper_cat = fits.open("DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits")
from astropy.table import Table
dat = Table.read('DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits', format='fits')
df = dat.to_pandas()


#%%
# red_mapper_cat.info()
df


#%%
# com_dists = cosmo.comoving_distance(df['ZREDMAGIC']).to_value(u.Mpc)
df['COMOVING'] = pd.Series(cosmo.comoving_distance(df['ZREDMAGIC']).to_value(u.Mpc))
df['COMOVING_E'] = pd.Series(cosmo.comoving_distance(df['ZREDMAGIC_E']).to_value(u.Mpc))


#%%
# help(galaxy_pairs.getPairs)


#%%
# print(df)
pairs = galaxy_pairs.getPairs(df,20)
# First index is first item in pair
# Second index in second item in pair
# Third Index is distance apart from each other


#%%
pairs_df = pd.DataFrame(pairs.T,columns = ['galaxy_index_1','galaxy_index_2','Sep'])


#%%
# print(pairs)
pairs_df['galaxy_index_1'] = pairs_df.galaxy_index_1.astype(int)
pairs_df['galaxy_index_2'] = pairs_df.galaxy_index_2.astype(int)


#%%
pairs_df


#%%
# import random
# index = random.randint(1,len(df))
# pos1 = [df['RA'][pairs_df['First_Loc'][index]],df['DEC'][pairs_df['First_Loc'][index]]]
# print(index)
# print(pos1)
# # pos2 = [df['RA'][pairs_df['Second_Loc'][0]],df['DEC'][pairs_df['Second_Loc'][0]]]
# # print(pos2)


#%%
# df.loc[index]


#%%
# # df.loc[pairs_df['Second_Loc'][0]]
# def SubtractRaDec(radec_1, radec_2):
# #     assert((len(radec_1) == 2)) #,"First enetered RaDec does not have two components")
# #     assert((len(radec_2) == 2)) #,"Second enetered RaDec does not have two components")

#     ra_1 = radec_1[0]
#     dec_1 = radec_1[1]
#     ra_2 = radec_2[0]
#     dec_2 = radec_2[1]
    
#     ra_diff = 180 - abs(abs(ra_1 - ra_2) - 180)
    
#     dec_diff = 180 - abs(abs(dec_1 - dec_2) - 180)
    
#     return([ra_diff,dec_diff])


#%%
# # rel_pos_1 = [0-pos1[0],-57.5-pos1[1]]
# rel_pos_1 = SubtractRaDec([0,-57.5],pos1)
# print(rel_pos_1)
# # print(rel_pos_1[0])
# # [15.089213,5.669675000000012]


#%%
# y_map.getSubmap([1,1], center_offset=[rel_pos_1[0],rel_pos_1[1]], units='degree').drawImage(bw=False)


#%%

    


#%%
# print(np.max(df['RA']))
# print(np.min(df['RA']))
# print(np.max(df['DEC']))
# print(np.min(df['DEC']))


#%%
# np.savetxt('pairs.csv',pairs.T,delimiter=',')


#%%
# vec1 = galaxy_pairs.get_vec_distances(df['RA'][int(pairs[0][0])],df['DEC'][int(pairs[0][0])],df['COMOVING'][int(pairs[0][0])])
# vec2 = galaxy_pairs.get_vec_distances(df['RA'][int(pairs[1][0])],df['DEC'][int(pairs[1][0])],df['COMOVING'][int(pairs[1][0])])
# print(galaxy_pairs.euclideanDistance(vec1,vec2,3))


#%%
# df.to_csv('DES_REDMAGIC_Manipulatedo.csv')


#%%
# vec_unit_test = RaDec2XYZ(df['RA'],df['DEC'])
# print(vec_unit_test)
# print(np.shape(vec_unit_test))
# vec_unit_test2 = vec_unit_test[0]
# print(vec_unit_test2)


#%%
# help(sptpol_software.observation.sky.ang2Pix)


#%%
# sptpol_software.observation.sky.ang2Pix(rel_pos_1,[0,-57.5],reso_arcmin=0.25,map_pixel_shape=[1320, 2520])


#%%
# def get_subarray(array,centre, sqr_radius):
#     ''' 
#     Gets Sub Array with with an input centre and half width from an input Array
#     '''
#     x_cen = centre[0]
#     y_cen = centre[1]
#     sl_x = slice(x_cen-sqr_radius,x_cen+sqr_radius)
#     sl_y = slice(y_cen-sqr_radius,y_cen+sqr_radius)
    
#     return(array[sl_x,sl_y])


#%%
plt.imshow(galaxy_pairs.get_subarray(y_map_array,[1275,500],50),vmax=0.00001,vmin=-0.0001)


#%%
top_left_corner = np.array([330,52])
top_right_corner = np.array([330,2468])
bottom_right_corner = np.array([1275,2000])
bottom_left_corner = np.array([1275,2520-2000])


#%%
edges = np.vstack((top_left_corner,top_right_corner,bottom_left_corner,bottom_right_corner))


#%%
edges


#%%



#%%
# temp1 = sptpol_software.observation.sky.pix2Ang(top_left_corner,np.array([0,-57.5]),reso_arcmin=1,map_pixel_shape=np.array([1320, 2520]))
# print(temp1[0]%360,temp1[1]%360)
# print(temp1)


#%%
sptpol_software.observation.sky.ang2Pix((0, -57.5),[0,-57.5],reso_arcmin=1,map_pixel_shape=np.array([1320,2520]))


#%%
temp2 = sptpol_software.observation.sky.pix2Ang(top_right_corner,np.array([0,-57.5]),reso_arcmin=1,map_pixel_shape=np.array([1320, 2520]))


#%%
temp3 = sptpol_software.observation.sky.pix2Ang(bottom_left_corner,np.array([0,-57.5]),reso_arcmin=1,map_pixel_shape=np.array([1320, 2520]))


#%%
temp4 = sptpol_software.observation.sky.pix2Ang(bottom_right_corner,np.array([0,-57.5]),reso_arcmin=1,map_pixel_shape=np.array([1320, 2520]))


#%%
edges_ang = np.array([temp1,temp2,temp3,temp4])
edges_ang = edges_ang.astype(int)
ra_range = np.array([edges_ang[0][0],edges_ang[1][0]])
dec_range = np.array([edges_ang[1][1],edges_ang[2][1]])
print(ra_range)
print(dec_range)


#%%
# from ..constants import DTOR, RTOD
# def pixel_2_angle(pixel_coords,ra_dec_centre,reso_arcmin,map_pixel_shape):
    
    
#     pixel_coords = (y_coord, x_coord) = pixel_coords[0].astype(float), pixel_coords[1].astype(float)
#     n_pixels = map_pixel_shape.astype(float)
    
#     y_coord = (y_coord + 0.5 - 0.5 * n_pixels[0]) * reso_arcmin / 60
#     x_coord = (x_coord + 0.5 - 0.5 * n_pixels[1]) * reso_arcmin / 60
    
    


#%%
# print(max(df['RA']))
# print(min(df['RA']))
# print(max(df['DEC']))
# print(min(df['DEC']))


#%%
new_df = df[(df.DEC < dec_range[0]) & (df.DEC > dec_range[1]) & (df.RA > ra_range[0]) & (df.RA < ra_range[1])]
new_df = new_df.reset_index(drop=True)


#%%
cut_pairs = galaxy_pairs.getPairs(new_df,20)
cut_pairs_df = pd.DataFrame(cut_pairs.T,columns = ['galaxy_index_1','galaxy_index_2','Sep'])
cut_pairs_df['galaxy_index_1'] = cut_pairs_df.galaxy_index_1.astype(int)
cut_pairs_df['galaxy_index_2'] = cut_pairs_df.galaxy_index_2.astype(int)
cut_pairs_df


#%%
help(galaxy_pairs.get_subarray)


#%%
def get_midpoint(ra_dec_1, ra_dec_2):
    '''
    Find the midpoint between two points in array space
    '''
    pt1 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_1, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
    pt2 = sptpol_software.observation.sky.ang2Pix(
        ra_dec_2, [0, -57.5], reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))

    X1 = float(pt1[0][0])
    X2 = float(pt2[0][0])

    Y1 = float(pt1[0][1])
    Y2 = float(pt2[0][1])

    return(((X1 + X2) / 2., (Y1 + Y2) / 2.))


#%%
def cut_out_pair(pair,y_map,galaxy_catalogue):
    '''
    Takes an input pair and a Compton Y-Map, and extract the pair as a sub map
    '''
    first_point = pair[0]
    second_point = pair[1]
    
    ra_1 = galaxy_catalogue.loc[first_point]['RA']
    dec_1 = galaxy_catalogue.loc[first_point]['DEC']
    
    ra_2 = galaxy_catalogue.loc[second_point]['RA']
    dec_2 = galaxy_catalogue.loc[second_point]['DEC']
    
    point_1 = (ra_1,dec_1)
    point_2 = (ra_2,dec_2)
    
    midpoint = get_midpoint(point_1,point_2)
    print(np.array(midpoint).astype(int))
    midpoint = np.array(midpoint).astype(int)
    return(galaxy_pairs.get_subarray(y_map,midpoint,50))
    
    


#%%
test_pair = [cut_pairs_df.loc[1000]['galaxy_index_1'],cut_pairs_df.loc[1000]['galaxy_index_2']]


#%%
test_pair


#%%
def stack_pairs(y_map, galaxy_catalogue, pairs):
    '''
    Take input Y-map, galaxy catalogue, and list of pairs, and stacks them on top of each other
    returning a stacked array
    '''
    size_of_cutout = 60
    output = np.ndarray([int(size_of_cutout / 2.), int(size_of_cutout / 2.)])
    for index, row in pairs.iterrows():
        galaxy_1 = row['galaxy_index_1']
        galaxy_2 = row['galaxy_index_2']
        pair = [galaxy_1, galaxy_2]

        cut_array = galaxy_pairs.cut_out_pair(pair, y_map,
                                 galaxy_catalogue, int(size_of_cutout / 2.))

        gal_1_coords = galaxy_pairs.extract_ra_dec(galaxy_1,galaxy_catalogue)
        gal_2_coords = galaxy_pairs.extract_ra_dec(galaxy_2,galaxy_catalogue)

        angle = galaxy_pairs.get_rotn_angle(gal_1_coords,gal_2_coords)
        rot_array = sp.ndimage.rotate(cut_array, angle, reshape=False)
        output += rot_array

    return(output)


#%%
gal_coords = galaxy_pairs.extract_ra_dec(cut_pairs_df.loc[1]['galaxy_index_1'],df)


#%%
galaxy_pairs.extract_ra_dec(gal_coords,df)


#%%
help(galaxy_pairs.extract_ra_dec)


#%%
output = stack_pairs(y_map_array,df,cut_pairs_df)


#%%
plt.imshow(galaxy_pairs.cut_out_pair(test_pair,y_map_array,df),vmax=0.000001,vmin=-0.000001)


#%%




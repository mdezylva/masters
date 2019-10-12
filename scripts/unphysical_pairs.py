from astropy.table import Table
import scipy.spatial.distance as dist
from astropy import units as u
from astropy.cosmology import Planck15
from astropy.io import fits
from astropy.coordinates import SkyCoord
import sptpol_software as sps
from sptpol_software.util.tools import stat
from sptpol_software.observation import *
import sptpol_software.observation as obs
import sptpol_software as sps
import sptpol_software.observation.sky
from sptpol_software.util import files
import PIL
import pickle
from astropy import constants as ap_const
import scipy.constants as const
from scipy.integrate import quad
import matplotlib.pyplot as plt
import pandas as pd
import glob
import astropy as ap
import scipy as sp
import numpy as np
import os
import time 
import pdb
os.chdir("/home/mitchell/Documents/masters/masters/scripts/")
import cmb
import galaxy_pairs
#data_location = input("Enter Path Location of data")
os.chdir("/home/mitchell/Documents/masters/masters/data")
cwd = os.getcwd()

import scipy.ndimage
cosmo = Planck15

# Correction Factors for sptPohl
T_cmb = 2.725
freqs_ghz = [93.2000, 147.700]
beam_norm_correction = [1./0.99673, 1./0.99470, 1./1.]
cal_factors = [0.9097, 0.7765]  # no 220 is why last is zero.
#pol_cal_factors = pol_cal_factors_800 * beam_norm_correction

max_trsv = 14/cosmo.h
min_trsv = 6/cosmo.h
min_los = 40/cosmo.h
max_los = 200/cosmo.h
# Read in Galaxy Catalogue as a dataframe 
dat = Table.read('DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits', format='fits')
df = dat.to_pandas()

# PreProcess DataFrame
print("Computing Comoving Distances from Redshifts in Galaxy Catalogue ...")
start = time.time()
df['COMOVING'] = pd.Series(
    cosmo.comoving_distance(df['ZREDMAGIC']).to_value(u.Mpc))
df['COMOVING_E'] = pd.Series(
    cosmo.comoving_distance(df['ZREDMAGIC_E']).to_value(u.Mpc))
#print("Pickling Data Frame with Comoving Coords")

#df.__module__ = "pandas"
#df.to_pickle("data_frame.pickle")

#uncut_dict = df.to_dict()
#pickle_out = open("data_frame.pickle","wb")
#pickle.dump(uncut_dict, pickle_out)
#pickle_out.close()
end = time.time()

print("Done in " + str(end-start) + " seconds")
# Establish Cuts 
start = end

top_left_corner = np.array([330, 52])
top_right_corner = np.array([330, 2468])
bottom_right_corner = np.array([1275, 2000])
bottom_left_corner = np.array([1275, 2520-2000])

edges = np.vstack((top_left_corner, top_right_corner,
                   bottom_left_corner, bottom_right_corner))

# Convert Edges into Angles 
temp1 = sptpol_software.observation.sky.pix2Ang(top_left_corner, np.array(
    [0, -57.5]), reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
temp2 = sptpol_software.observation.sky.pix2Ang(top_right_corner, np.array(
    [0, -57.5]), reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
temp3 = sptpol_software.observation.sky.pix2Ang(bottom_left_corner, np.array(
    [0, -57.5]), reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
temp4 = sptpol_software.observation.sky.pix2Ang(bottom_right_corner, np.array(
    [0, -57.5]), reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))


edges_ang = np.array([temp1, temp2, temp3, temp4])
edges_ang = edges_ang.astype(int)
ra_range = np.array([edges_ang[0][0], edges_ang[1][0]])
dec_range = np.array([edges_ang[1][1], edges_ang[2][1]])
# print(ra_range)
# print(dec_range)


cut_df = df[((df.DEC < dec_range[0]-0.04) & (df.DEC > dec_range[1]+0.04))
            & ((df.RA > ra_range[0] +0.04 ) & (df.RA < ra_range[1] -0.04))]
cut_df = cut_df.reset_index(drop=True)

print("Saving Cut Galaxy Catalogue to Pickle")
# cut_df.to_pickle('unphysical_cat.pkl')

# cut_dict = cut_df.to_dict()
# pickle_out = open("cut_df.pickle","wb")
# pickle.dump(cut_dict, pickle_out)
# pickle_out.close()

# half_cat = cut_df.take(np.arange(len(cut_df/10)))
gal_1 = []
gal_2 = []
gal_1_z = []
gal_2_z = []
gal_1_ra = []
gal_1_dec = []
gal_2_ra = []
gal_2_dec = []
comov_trv_list = []
del_theta = []
sep_los = []
sep_trv = []


print("Constructing Pairs")
cut_pairs = galaxy_pairs.getPairs(cut_df, max_sep=100,query_type=1)
# for index1, row in cut_df.iterrows():
#     for index2 in range(index1,len(cut_df)):
#         coord1 = SkyCoord(ra= cut_df['RA'].iloc[index1]*u.degree,dec = cut_df['DEC'].iloc[index1]*u.degree)
#         coord2 SkyCoord(ra= cut_df['RA'].iloc[index2]*u.degree,dec = cut_df['DEC'].iloc[index2]*u.degree)
        
#         # Quickest cut to make is on Angular Separation
#         # If they are more than 20 arcmins in separation, skip the pair calculation
#         ang_sep = coord_1.separation(coord_2).arcmin

#         if ang_sep > 20:
#             continue

#         # Make cut based on Transverse Separation 
#         z1 = cut_df['ZREDMAGIC'].iloc[index1]
#         z2 = cut_df['ZREDMAGIC'].iloc[index2]
#         avg_z = np.mean((z1,z2))

#         comov_trv = cosmo.comoving_transverse_distance(avg_z)

#         trv_sep = ang_sep*comov_trv

#         if (trv_sep > max_trsv) or (trv_sep < min_trsv):
#             continue 

#         # Work out line of sight separation now
#         gal_1_comov = cut_df['COMOVING'].iloc[index1]
#         gal_2_comov = cut_df['COMOVING'].iloc[index2]
#         comov_sep = pd.Series(abs(gal_1_comov - gal_2_comov))   

#         # Cut if line of sight separation is too large
#         if (comov_sep < min_los) or (comov_sep > max_los):
#             continue  

#         gal_1.append(index1)
#         gal_2.append(index2)
#         gal_1_z.append(z1)
#         gal_2_z.append(z2)
#         gal_1_ra.append(cut_df['RA'].iloc[index1])
#         gal_1_dec.append(cut_df['DEC'].iloc[index1])
#         gal_2_ra.append(cut_df['RA'].iloc[index2])
#         gal_2_dec.append(cut_df['DEC'].iloc[index2])
#         comov_trv_list.append(comov_trv)
#         del_theta.append(ang_sep)
#         sep_los.append(comov_sep)
#         sep_trv.append(trv_sep)


# data = list(zip(gal_1, gal_2, gal_1_z, gal_2_z, gal_1_ra, gal_1_dec, gal_2_ra, gal_2_dec, comov_trv_list, del_theta, sep_los, sep_trv))
# cut_pairs_df = pd.DataFrame(data, columns = ['galaxy_index_1', 'galaxy_index_2', 'Z1', 'Z2', 'RA_1', 'DEC_1', 'RA_2', 'DEC_2', 'COMOV_TRV', 'DEL_THETA', 'SEP_LOS', 'SEP_TRV'])
# #pdb.set_trace()
cut_pairs_df = pd.DataFrame(
    cut_pairs.T, columns=['galaxy_index_1', 'galaxy_index_2', 'Sep'])
cut_pairs_df['galaxy_index_1'] = cut_pairs_df.galaxy_index_1.astype(int)
cut_pairs_df['galaxy_index_2'] = cut_pairs_df.galaxy_index_2.astype(int)

print(len(cut_pairs_df))


print("Saving Cut Galaxy Catalogue to Pickle")
cut_pairs_df.__module__ = 'pandas'
cut_pairs_df.to_pickle('unphysical_cat_4.pkl')


# cut_dict = cut_pairs_df.to_dict()
# pickle.dump(cut_dict, pickle_out)
# pickle_out = open("cut_pairs.pickle","wb")
# pickle_out.close()
print("Computing Comoving Separations...")
gal_1_comov = cut_df['COMOVING'].iloc[cut_pairs_df['galaxy_index_1']]
gal_2_comov = cut_df['COMOVING'].iloc[cut_pairs_df['galaxy_index_2']]
comov_sep = pd.Series(abs(gal_1_comov.values - gal_2_comov.values))

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Calculating Comoving Separation...")
cut_pairs_df['SEP_LOS'] = comov_sep

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

gal_1_z = cut_df['ZREDMAGIC'].iloc[cut_pairs_df['galaxy_index_1']]
gal_2_z = cut_df['ZREDMAGIC'].iloc[cut_pairs_df['galaxy_index_2']]

print("Appending Redshifts of Pair to Dataframe...")
cut_pairs_df['Z1'] = pd.Series(gal_1_z.values)
cut_pairs_df['Z2'] = pd.Series(gal_2_z.values)

end = time.time()
print("Done in " + str(end-start) + " seconds")
# start = end

print("Computing Average Redshift...")
avg_redshift = pd.Series(np.mean([gal_1_z,gal_2_z],axis=0))

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Computing Transverse Comoving Distance at this redshift...")
transverse_comoving = cosmo.comoving_transverse_distance(avg_redshift).to_value(u.Mpc)

cut_pairs_df['COMOV_TRV'] = transverse_comoving

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Appending (RA,DEC) to Dataframe...")
cut_pairs_df['RA_1'] = pd.Series(cut_df['RA'].iloc[cut_pairs_df['galaxy_index_1']].values)
cut_pairs_df['RA_2'] = pd.Series(cut_df['RA'].iloc[cut_pairs_df['galaxy_index_2']].values)
cut_pairs_df['DEC_1'] = pd.Series(cut_df['DEC'].iloc[cut_pairs_df['galaxy_index_1']].values)
cut_pairs_df['DEC_2'] = pd.Series(cut_df['DEC'].iloc[cut_pairs_df['galaxy_index_2']].values)

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Storing as Sky Coordinates...")
coord_1 = SkyCoord(ra= cut_pairs_df['RA_1']*u.degree,dec = cut_pairs_df['DEC_1']*u.degree)
coord_2 = SkyCoord(ra= cut_pairs_df['RA_2']*u.degree,dec = cut_pairs_df['DEC_2']*u.degree)

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

t_sep = np.empty((len(coord_1)))

print("Computing Angular Separations...")
for index in range(len(coord_1)):
    t_sep[index] = coord_1[index].separation(coord_2[index]).radian

cut_pairs_df['DEL_THETA'] = pd.Series(t_sep)

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Computing Transverse Separation...")
cut_pairs_df['SEP_TRV'] = cut_pairs_df['DEL_THETA'] * cut_pairs_df['COMOV_TRV']

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Done!")

print("Cutting Pairs based on separation criteria... ")
max_trsv = 14/cosmo.h
min_trsv = 6/cosmo.h 
max_los = 40/cosmo.h

print("Max Line of Sight Separation: " + str(max_los) + " Mpc" )
print("Transverse Separation Range: " + str(max_trsv) + ' to ' + str(min_trsv) + " Mpc")


# output_df = cut_pairs_df[ (cut_pairs_df.SEP_TRV < max_trsv ) & (cut_pairs_df.SEP_TRV > max_trsv) & (cut_pairs_df.SEP_LOS < max_los)]
cut_pairs_df.reset_index()
cut_pairs_df.head()




print(' Cutting pairs based on separation conditions')
cut_pairs_df = cut_pairs_df[ (cut_pairs_df.SEP_TRV < max_trsv) & (cut_pairs_df.SEP_TRV > min_trsv) & (cut_pairs_df.SEP_LOS > max_los)]

cut_pairs_df = cut_pairs_df.reset_index()

print("Pickling List of Pairs...")
cut_pairs_df.to_pickle('unphysical_pairs_4.pkl')

cut_pairs_df.head()

print("Done! ")
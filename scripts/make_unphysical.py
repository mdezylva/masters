from astropy.table import Table
import scipy.spatial.distance as dist
from astropy import units as u
from astropy.cosmology import Planck13
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
os.chdir("/home/mdz/scripts")
import cmb
import galaxy_pairs
#data_location = input("Enter Path Location of data")
os.chdir("/data53/cr/mitchell")
cwd = os.getcwd()

import scipy.ndimage
cosmo = Planck13

# Correction Factors for sptPohl
T_cmb = 2.725
freqs_ghz = [93.2000, 147.700]
beam_norm_correction = [1./0.99673, 1./0.99470, 1./1.]
cal_factors = [0.9097, 0.7765]  # no 220 is why last is zero.
#pol_cal_factors = pol_cal_factors_800 * beam_norm_correction

# # Read in signal maps
# map_150ghz = files.read("ra0dec-57p5_sum5000_150ghz.h5")
# map_90ghz = files.read("ra0dec-57p5_sum5000_090ghz.h5")

# # Produce Y Map from Signal Maps
# y_map_array = cmb.get_y_map([map_150ghz, map_90ghz], cal_factors, freqs_ghz)

# plt.imshow(y_map_array, vmax=0.00001, vmin=-0.00001)

# Read in Galaxy Catalogue as a dataframe 
#dat = Table.read('DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits', format='fits')
#df = dat.to_pandas()

df = pd.read_pickle('DES_cat.pkl')

# PreProcess DataFrame
print("Computing Comoving Distances from Redshifts in Galaxy Catalogue ...")
start = time.time()
df['COMOVING'] = pd.Series(
    cosmo.comoving_distance(df['ZREDMAGIC']))
df['COMOVING_E'] = pd.Series(
    cosmo.comoving_distance(df['ZREDMAGIC_E']))
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

norm_ra = []
for i in range(len(df)):
    norm_ra.append(galaxy_pairs.normalize(df['RA'][i],-180,180))
df['NORM_RA'] = pd.Series(norm_ra)


### SPT POL
ra_range = np.array([-32,32])
dec_range = np.array([-52,-67])



cut_df = df[((df.DEC < dec_range[0]) & (df.DEC > dec_range[1]))
            & ((df.NORM_RA > ra_range[0]) & (df.NORM_RA < ra_range[1]))]
cut_df = cut_df.reset_index(drop=True)

#print("Saving Cut Galaxy Catalogue to Pickle")
#cut_df.to_pickle('DES_cat.pkl')

# cut_dict = cut_df.to_dict()
# pickle_out = open("cut_df.pickle","wb")
# pickle.dump(cut_dict, pickle_out)
# pickle_out.close()

sky_locs = np.asarray([cut_df['NORM_RA'].as_matrix().T,cut_df['DEC'].as_matrix().T]).T

frame_tree = scipy.spatial.cKDTree(sky_locs, leafsize=100)

result_set = frame_tree.query_pairs(r=0.5)

result_list = list(result_set)

result_array = np.asarray(result_list)

cut_pairs_df = pd.DataFrame(data=result_array,columns = ['galaxy_index_1','galaxy_index_2'])

#cut_pairs = galaxy_pairs.getUnPairs(cut_df, max_sep=0.5)
#pdb.set_trace()
#cut_pairs_df = pd.DataFrame(
#    cut_pairs.T, columns=['galaxy_index_1', 'galaxy_index_2', 'Sep'])
cut_pairs_df['galaxy_index_1'] = cut_pairs_df.galaxy_index_1.astype(int)
cut_pairs_df['galaxy_index_2'] = cut_pairs_df.galaxy_index_2.astype(int)

print("Num of Pairs Produced on Footprint = " + str(len(cut_pairs_df)))

# print("Saving Cut Galaxy Catalogue to Pickle")
#cut_pairs_df.__module__ = 'pandas'
# cut_pairs_df.to_pickle('cut_catalogue.pkl')


# cut_dict = cut_pairs_df.to_dict()
# pickle_out = open("cut_pairs.pickle","wb")
# pickle.dump(cut_dict, pickle_out)
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


min_los = 100/cosmo.h
max_los = 200/cosmo.h


cut_pairs_df = cut_pairs_df[(cut_pairs_df.SEP_LOS < max_los) & (cut_pairs_df.SEP_LOS > min_los)]

cut_pairs_df = cut_pairs_df.reset_index()

print("Num of Pairs after LOS Cuts = " + str(len(cut_pairs_df)))

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
start = end

print("Computing Average Redshift...")
avg_redshift = pd.Series(np.mean([gal_1_z,gal_2_z],axis=0))

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Computing Transverse Comoving Distance at this redshift...")
transverse_comoving = cosmo.comoving_transverse_distance(avg_redshift)

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
print("Max Line of Sight Separation: " + str(max_los) + " Mpc" )
print("Transverse Separation Range: " + str(max_trsv) + ' to ' + str(min_trsv) + " Mpc")


# output_df = cut_pairs_df[ (cut_pairs_df.SEP_TRV < max_trsv ) & (cut_pairs_df.SEP_TRV > max_trsv) & (cut_pairs_df.SEP_LOS < max_los)]
cut_pairs_df.reset_index()
cut_pairs_df.head()


max_trsv = 14/cosmo.h
min_trsv = 6/cosmo.h

print(' Cutting pairs based on separation conditions')
cut_pairs_df = cut_pairs_df[ (cut_pairs_df.SEP_TRV < max_trsv) & (cut_pairs_df.SEP_TRV > min_trsv)
	& (cut_pairs_df.SEP_LOS < max_los) & (cut_pairs_df.SEP_LOS > min_los)]

cut_pairs_df = cut_pairs_df.reset_index()


print("Num of Pairs after TRV Cuts = " + str(len(cut_pairs_df)))
print("Pickling List of Pairs...")
cut_pairs_df.to_pickle('unphysical_pairs_pol_2.pkl')

cut_pairs_df.head()

print("Done! ")
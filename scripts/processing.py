from astropy.table import Table
import scipy.spatial.distance as dist
from astropy import units as u
from astropy.cosmology import Planck15
from astropy.io import fits
import sptpol_software as sps
from sptpol_software.util.tools import stat
from sptpol_software.observation import *
import sptpol_software.observation as obs
import sptpol_software as sps
import sptpol_software.observation.sky
from sptpol_software.util import files
import PIL
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

# Read in signal maps
map_150ghz = files.read("ra0dec-57p5_sum5000_150ghz.h5")
map_90ghz = files.read("ra0dec-57p5_sum5000_090ghz.h5")

# Produce Y Map from Signal Maps
y_map_array = cmb.get_y_map([map_150ghz, map_90ghz], cal_factors, freqs_ghz)

plt.imshow(y_map_array, vmax=0.00001, vmin=-0.00001)

# Read in Galaxy Catalogue as a dataframe 
dat = Table.read('DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits', format='fits')
df = dat.to_pandas()

# PreProcess DataFrame
df['COMOVING'] = pd.Series(
    cosmo.comoving_distance(df['ZREDMAGIC']).to_value(u.Mpc))
df['COMOVING_E'] = pd.Series(
    cosmo.comoving_distance(df['ZREDMAGIC_E']).to_value(u.Mpc))


# Establish Cuts 

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


cut_df = df[((df.DEC < dec_range[0]) & (df.DEC > dec_range[1]))
            & ((df.RA > ra_range[0]) & (df.RA < ra_range[1]))]
cut_df = cut_df.reset_index(drop=True)


cut_pairs = galaxy_pairs.getPairs(cut_df, max_sep=14,query_type=1)
cut_pairs_df = pd.DataFrame(
    cut_pairs.T, columns=['galaxy_index_1', 'galaxy_index_2', 'Sep'])
cut_pairs_df['galaxy_index_1'] = cut_pairs_df.galaxy_index_1.astype(int)
cut_pairs_df['galaxy_index_2'] = cut_pairs_df.galaxy_index_2.astype(int)


output = galaxy_pairs.stack_pairs_V2(y_map_array, cut_df, cut_pairs_df,debug=False)

print(output)
plt.imshow(output)#,vmax=1,vmin=-1)
plt.show()
from astropy.table import Table
import scipy.spatial.distance as dist
from astropy import units as u
from astropy.cosmology import Planck13
from astropy.io import fits
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
from astropy.coordinates import SkyCoord
import astropy as ap
import scipy as sp
import numpy as np
import healpy as hp
import os
import pdb
import gzip
os.chdir("/home/mdz/scripts/")
import cmb
import galaxy_pairs
#data_location = input("Enter Path Location of data")
os.chdir("/data53/cr/mitchell/")
cwd = os.getcwd()
print(cwd)
import scipy.ndimage
cosmo = Planck13

# Correction Factors for sptPohl
T_cmb = 2.725
freqs_ghz = [93.2000, 147.700]
beam_norm_correction = [1./0.99673, 1./0.99470, 1./1.]
cal_factors = [0.9097, 0.7765]  # no 220 is why last is zero.
#pol_cal_factors = pol_cal_factors_800 * beam_norm_correction
print(glob.glob('*.h5'))
print(glob.glob('*.fits'))
# Read in signal maps
print("... Loading Map(s) ...")
#map_150ghz = files.read("ra0dec-57p5_sum5000_150ghz.h5")
f1 = gzip.open('map_90.pkl.gz','rb')
map_90 = pickle.load(f1)
map_90ghz = map_90[0]['T'].map
#map_90ghz = files.read("ra0dec-57p5_sum5000_90ghz.h5")
f2 = gzip.open('map_150.pkl.gz','rb')
map_150 = pickle.load(f2)
map_150ghz = map_150[0] 

#y_map_array = hp.read_map('Planck_SPT_100_353_ymapdec18_ymap_minimum_variance1_cib4_nside_8192.fitsmap.fits')


# Produce Y Map from Signal Maps
print("... Creating y-map ...")
#y_map_array = cmb.get_y_map([map_150ghz, map_90ghz], cal_factors, freqs_ghz)
map_90ghz = map_90ghz*cal_factors[0]
map_150ghz = map_150ghz*cal_factors[1]
freq_scaling = 1.0 /(cmb.convert_ghz_to_y(150) - cmb.convert_ghz_to_y(90))
diff_map = map_150ghz - map_90ghz

y_map_array = freq_scaling*diff_map

#plt.imshow(y_map_array, vmax=0.00001, vmin=-0.00001)
print("... Opening galaxy cataloge ...")

# Read in Galaxy Catalogue as a dataframe 
#dat = Table.read('DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits', format='fits')
#df = dat.to_pandas()
#with fits.open('DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits') as data:
#    df = pd.DataFrame(data[1].data)
print(glob.glob("*.pkl"))

#df = pd.read_pickle('data_frame.pickle')

# PreProcess DataFrame
#print(" ... Converting to Comoving Mpc ...")
#df['COMOVING'] = pd.Series(cosmo.comoving_distance(df['ZREDMAGIC']).to(u.Mpc).value)
#df['COMOVING_E'] = pd.Series(cosmo.comoving_distance(df['ZREDMAGIC_E']).to(u.Mpc).value)


# Establish Cuts 
#print("... cutting galaxies outside sptPol area ...")
#top_left_corner = np.array([330, 52])
#top_right_corner = np.array([330, 2468])
#bottom_right_corner = np.array([1275, 2000])
#bottom_left_corner = np.array([1275, 2520-2000])

#edges = np.vstack((top_left_corner, top_right_corner,
#                   bottom_left_corner, bottom_right_corner))

# Convert Edges into Angles 
#temp1 = sptpol_software.observation.sky.pix2Ang(top_left_corner, np.array(
#    [0, -57.5]), reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
#temp2 = sptpol_software.observation.sky.pix2Ang(top_right_corner, np.array(
#    [0, -57.5]), reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
#temp3 = sptpol_software.observation.sky.pix2Ang(bottom_left_corner, np.array(
#    [0, -57.5]), reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))
#temp4 = sptpol_software.observation.sky.pix2Ang(bottom_right_corner, np.array(
#    [0, -57.5]), reso_arcmin=1, map_pixel_shape=np.array([1320, 2520]))


#edges_ang = np.array([temp1, temp2, temp3, temp4])
#edges_ang = edges_ang.astype(int)
#ra_range = np.array([edges_ang[0][0], edges_ang[1][0]])
#dec_range = np.array([edges_ang[1][1], edges_ang[2][1]])
# print(ra_range)
# print(dec_range)

#print("... Cutting dataframe ...")
#cut_df = df[((df.DEC < dec_range[0]-0.04) & (df.DEC > dec_range[1]+0.04))
#            & ((df.RA > ra_range[0] +0.04 ) & (df.RA < ra_range[1] -0.04))]
#cut_df = cut_df.reset_index(drop=True)
print("Loading in cut dataframe")
#pickle_in = open("cut_df.pickle")
#cut_df_dict = pickle.load(pickle_in)
#cut_df = pd.DataFrame.from_dict(cut_df_dict)

cut_df = pd.read_pickle('cut_catalogue.pkl')

#cut_pairs = galaxy_pairs.getPairs(cut_df, max_sep=20,query_type=1)
#cut_pairs_df = pd.DataFrame(
#    cut_pairs.T, columns=['galaxy_index_1', 'galaxy_index_2', 'Sep'])
#cut_pairs_df['galaxy_index_1'] = cut_pairs_df.galaxy_index_1.astype(int)
#cut_pairs_df['galaxy_index_2'] = cut_pairs_df.galaxy_index_2.astype(int)

print("Number of Galaxies in Catalogue = " + str(len(cut_df.index)))

print("Loading in cut data_frame pickle")
#pickle_in = open("cut_pairs.pickle","rb")
#cut_pairs_dict = pickle.load(pickle_in)
#print("..done..")
#print("..converting to data frame..")
#cut_pairs_df = pd.DataFrame.from_dict(cut_pairs_dict)

max_trsv = 14/cosmo.h
min_trsv = 6/cosmo.h 
max_los = 10/cosmo.h

data_frame = pd.read_pickle('galaxy_pairs.pkl')

cut_pairs_df = data_frame[ (data_frame.SEP_TRV < max_trsv) & (data_frame.SEP_TRV > min_trsv) & (data_frame.SEP_LOS < max_los)]

pdb.set_trace()

print("Number of Pairs = " +str(len(cut_pairs_df.index)))
print("... Stacking Galaxies ...")
plt.switch_backend('Agg')
output = galaxy_pairs.stack_pairs_V2(y_map_array, galaxy_catalogue=cut_df, pairs=cut_pairs_df,debug=False,save=True)
sp.misc.imsave("/home/mdz/output.png",output)
print(output)
#plt.imshow(output)#,vmax=1,vmin=-1)
#plt.show()

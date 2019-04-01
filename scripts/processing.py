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
print(cwd)

import scipy.ndimage
print(Planck15)
cosmo = Planck15

# Correction Factors for sptPol
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

import numpy as np
import os
import random
import matplotlib.pyplot as plt
from IPython import embed
import cPickle as pickle
import gzip
import healpy as hp
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from IPython import embed
# plt.switch_backend('Agg')
from tqdm import tqdm

RA_SPTSZ_min  = 0
RA_SPTSZ_max  = 165
DEC_SPTSZ_min = -65
DEC_SPTSZ_max = -40

RA_SPTpol_min  = 30
RA_SPTpol_max  = 100
DEC_SPTpol_min = -65
DEC_SPTpol_max = -50

def fn_load_halo(fname, Mmin=1e14, Mmax=1e15, zmin=0., zmax=1., nobj=None, like_SPTSZ=0):
	"""
	Returns a RA, DEC, Z with redmapper catalog info.
	"""
	cat = fits.open(fname)[1].data
	cat = cat[(cat.REDSHIFT >= zmin) & (cat.REDSHIFT <= zmax)]
	cat = cat[(cat.M200 >= Mmin) & (cat.M200 <= Mmax)]
	if like_SPTSZ:
		cat = cat[(cat.RA >= RA_SPTSZ_min) & (cat.RA <= RA_SPTSZ_max)]
		cat = cat[(cat.DEC >= DEC_SPTSZ_min) & (cat.DEC <= DEC_SPTSZ_max)]

	if nobj is not None:
		cat = cat[np.random.randint(0,len(cat),size=nobj)]

	ra    = cat.RA
	dec   = cat.DEC
	zs    = cat.REDSHIFT
	v_los = cat.VLOS
	M200  = cat.M200

	return ra, dec, zs, v_los, M200

# Params 'n' stuff
flenders_cat    = '/home/mitchell/Documents/masters/masters/data/simulated/catalog_fullsky_zsm1.fits'
tsz_map_file    = 'data/simulated/tsz150_R13.fits'
# cutout_dic_name = 'kSZ_tSZ_SPTSZ_footprint_zsm1.pkl.gz'
nside           = 8192
boxsize         = 50   # arcmin
reso            = 0.25 #arcmin

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("...loading Flender's halo catalogue...")
ra, dec, z, v_los, M200 = fn_load_halo(flenders_cat, like_SPTSZ=1, nobj=None)
print("...# CLUSTERS = %d"%len(ra))
print("...done...")


print("...loading tSZ map...")
tsz_map = hp.read_map(tsz_map_file)
tsz_map_smooth = hp.smoothing(tsz_map, fwhm=np.radians(1.2/60.))
print("...done...")

cuts = {}
cuts['cutouts_tSZ'] = []
cuts['RA']      = []
cuts['DEC']     = []
cuts['Z']       = []
cuts['VLOS']    = []
cuts['M200']    = []
cuts['SPT-SZ']  = []
cuts['SPTpol']  = []

# coord = SkyCoord(ra=ra, dec=dec, unit='deg').transform_to('galactic')
# l = coord.l.value
# b = coord.b.value

# Loop over sources
print("...starting loop over objects...")
for i in tqdm(range(ra.size)):
	cut_tSZ = hp.gnomview(tsz_map, rot=[ra[i],dec[i]], xsize=boxsize, reso=reso, return_projected_map=True)
	plt.close()
	cuts['cutouts_tSZ'].append(cut_tSZ)
	cuts['RA'].append(ra[i])
	cuts['DEC'].append(dec[i])
	cuts['Z'].append(z[i])
	cuts['VLOS'].append(v_los[i])
	cuts['M200'].append(M200[i])

	if ra[i] > RA_SPTSZ_min and ra[i] < RA_SPTSZ_max and dec[i] > DEC_SPTSZ_min and dec[i] < DEC_SPTSZ_max:
		cuts['SPT-SZ'].append(True)
	else:
		cuts['SPT-SZ'].append(False)

	if ra[i] > RA_SPTpol_min and ra[i] < RA_SPTpol_max and dec[i] > DEC_SPTpol_min and dec[i] < DEC_SPTpol_max:
		cuts['SPTpol'].append(True)
	else:
		cuts['SPTpol'].append(False)

cuts['cutouts_kSZ'] = np.asarray(cuts['cutouts_kSZ'])
cuts['cutouts_tSZ'] = np.asarray(cuts['cutouts_tSZ'])
cuts['RA']  = np.asarray(cuts['RA'])
cuts['DEC'] = np.asarray(cuts['DEC'])
cuts['Z'] = np.asarray(cuts['Z'])
cuts['VLOS'] = np.asarray(cuts['VLOS'])
cuts['M200'] = np.asarray(cuts['M200'])
cuts['SPT-SZ']  = np.asarray(cuts['SPT-SZ'])
cuts['SPTpol']  = np.asarray(cuts['SPTpol'])
print("...done...")

#embed()

pickle.dump(cuts, gzip.open(cutout_dic_name, 'w'), protocol=2)
embed()

import numpy as np
import pandas as pd
import pickle
import os 
import glob 
import time
import astropy as ap

import galaxy_pairs

from astropy.cosmology import Planck13
from astropy import units as u
from astropy.coordinates import SkyCoord




os.chdir('/home/mitchell/Documents/masters/masters/data/simulated')
cwd = os.getcwd()
cosmo = Planck13

glob.glob('*')

cat = pd.read_pickle('flender_cat.pkl')

print("Cutting Catalogue so that it behaves like SPTpol..")

start = time.time()

RA_SPTpol_min  = 30
RA_SPTpol_max  = 100
DEC_SPTpol_min = -65
DEC_SPTpol_max = -50

cat = cat[(cat.RA >= RA_SPTpol_min) & (cat.RA <= RA_SPTpol_max)]
cat = cat[(cat.DEC >= DEC_SPTpol_min) & (cat.DEC <= DEC_SPTpol_max)]
cat = cat.reset_index(drop=True)


end = time.time()
print("Completed in: "+ str(end-start) + " seconds")

print("Computing Comoving Distances from Redshifts in Galaxy Catalogue ...")
start = time.time()

cat['COMOVING'] = pd.Series(
	cosmo.comoving_distance(cat['REDSHIFT']))

end = time.time()
print("Completed in: "+ str(end-start) + " seconds")


print("Producing Galaxy Pairs...")

start = time.time()
pairs = galaxy_pairs.getPairs(cat, max_sep=30,query_type=1)
#pdb.set_trace()
pairs_df = pd.DataFrame(
    pairs.T, columns=['galaxy_index_1', 'galaxy_index_2', 'Sep'])
pairs_df['galaxy_index_1'] = pairs_df.galaxy_index_1.astype(int)
pairs_df['galaxy_index_2'] = pairs_df.galaxy_index_2.astype(int)

end = time.time()
print("Completed in: "+ str(end-start) + " seconds")


## Change this Code
print("Computing Comoving Separations...")
gal_1_comov = cat['COMOVING'].iloc[pairs_df['galaxy_index_1']]
gal_2_comov = cat['COMOVING'].iloc[pairs_df['galaxy_index_2']]
comov_sep = pd.Series(abs(gal_1_comov.values - gal_2_comov.values))

end = time.time()
print("Completed in: " + str(end-start) + " seconds")
start = end

print("Calculating Comoving Separation...")
pairs_df['SEP_LOS'] = comov_sep

end = time.time()
print("Completed in: " + str(end-start) + " seconds")
start = end

gal_1_z = cat['REDSHIFT'].iloc[pairs_df['galaxy_index_1']]
gal_2_z = cat['REDSHIFT'].iloc[pairs_df['galaxy_index_2']]

print("Appending Redshifts of Pair to Dataframe...")
pairs_df['Z1'] = pd.Series(gal_1_z.values)
pairs_df['Z2'] = pd.Series(gal_2_z.values)

end = time.time()
print("Completed in: " + str(end-start) + " seconds")
start = end

print("Computing Average Redshift...")
avg_redshift = pd.Series(np.mean([gal_1_z,gal_2_z],axis=0))

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Computing Transverse Comoving Distance at this redshift...")
transverse_comoving = cosmo.comoving_transverse_distance(avg_redshift).to_value(u.Mpc)

pairs_df['COMOV_TRV'] = transverse_comoving

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Appending (RA,DEC) to Dataframe...")
pairs_df['RA_1'] = pd.Series(cat['RA'].iloc[pairs_df['galaxy_index_1']].values)
pairs_df['RA_2'] = pd.Series(cat['RA'].iloc[pairs_df['galaxy_index_2']].values)
pairs_df['DEC_1'] = pd.Series(cat['DEC'].iloc[pairs_df['galaxy_index_1']].values)
pairs_df['DEC_2'] = pd.Series(cat['DEC'].iloc[pairs_df['galaxy_index_2']].values)

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Storing as Sky Coordinates...")
coord_1 = SkyCoord(ra= pairs_df['RA_1']*u.degree,dec = pairs_df['DEC_1']*u.degree)
coord_2 = SkyCoord(ra= pairs_df['RA_2']*u.degree,dec = pairs_df['DEC_2']*u.degree)

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

t_sep = np.empty((len(coord_1)))

print("Computing Angular Separations...")
for index in range(len(coord_1)):
    t_sep[index] = coord_1[index].separation(coord_2[index]).radian

pairs_df['DEL_THETA'] = pd.Series(t_sep)

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Computing Transverse Separation...")
pairs_df['SEP_TRV'] = pairs_df['DEL_THETA'] * pairs_df['COMOV_TRV']

end = time.time()
print("Done in " + str(end-start) + " seconds")
start = end

print("Done!")

print("Cutting Pairs based on separation criteria... ")
max_trsv = 14/cosmo.h
min_trsv = 6/cosmo.h 
max_los = 10/cosmo.h

print("Max Line of Sight Separation: " + str(max_los) + " Mpc" )
print("Transverse Separation Range: " + str(max_trsv) + ' to ' + str(min_trsv) + " Mpc")

print(' Cutting pairs based on separation conditions')
pairs_df = pairs_df[ (pairs_df.SEP_TRV < max_trsv) & (pairs_df.SEP_TRV > min_trsv) & (pairs_df.SEP_LOS < max_los)]

pairs_df = pairs_df.reset_index()

print("Pickling List of Pairs...")
pairs_df.to_pickle('flender_pairs.pkl')


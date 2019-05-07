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
from mpdaf.obj import Image,WCS
import sptpol_software as sps
from astropy.io import fits
from astropy.cosmology import Planck15
from astropy import units as u
import scipy.spatial.distance as dist
print(Planck15)
cosmo= Planck15
cwd = os.getcwd()
print(cwd)
from sptpol_software.util.tools import stat
from sptpol_software.observation import *
import sptpol_software.observation as obs
import sptpol_software as sps
import sptpol_software.observation.sky
from sptpol_software.util import files


#%%
size_of_array = int(input("Enter the size of the array: "))
array_1= np.zeros(shape = (size_of_array,size_of_array))
array_2 = np.zeros(shape = (size_of_array,size_of_array))


#%%
source_loc_1 = np.array([np.random.randint(0,len(array_1)-10),np.random.randint(0,len(array_1)-10)])
source_loc_2 = np.array([np.random.randint(0,len(array_2)-10),np.random.randint(0,len(array_2)-10)])
print(source_loc_1)
print(source_loc_2)


#%%
x1 = source_loc_1[0]
print("x1 = " + str(x1))
print("-x1 = " + str(len(array_1)-x1))
y1 = source_loc_1[1]
print("y1 = " + str(y1))
print("-y1 = " + str(len(array_1)-y1))

x2 = source_loc_2[0]
print("x2 = " + str(x2))
print("-x2 = " + str(len(array_1)-x2))

y2 = source_loc_2[1]
print("y2 = " + str(y2))
print("-y2 = " + str(len(array_1)-y2))

# source = makeGaussian(size_of_source,size_of_source)
# array[x1:x1+source.shape[0], y1:y1+source.shape[1]] = source
# array[-x1-source.shape[0]:-x1, -y1-source.shape[1]:-y1] = source


#%%
# plt.imshow(array)


#%%
# source_x1 = x1+(source.shape[0]/2)
# source_y1 = y1+(source.shape[1]/2)
# source_x2 = -x1-(source.shape[0]/2)
# source_y2 = -y1+(source.shape[1]/2)

angle_1 = np.degrees(np.arctan((float(len(array_1)-y1)-float(y1))/(float(len(array_1)-x1)-float(x1))))
angle_2 =  np.degrees(np.arctan((float(len(array_2)-y2)-float(y2))/(float(len(array_2)-x2)-float(x2))))


#%%
print(angle_1)
print(angle_2)


#%%



#%%



#%%



#%%
array_1[source_loc_1[0]][source_loc_1[1]] = 1
array_1[-source_loc_1[0]][-source_loc_1[1]] = 1
gauss_array_1 = sp.ndimage.gaussian_filter(array_1,2)
array_2[source_loc_2[0]][source_loc_2[1]] = 1
array_2[-source_loc_2[0]][-source_loc_2[1]] = 1
gauss_array_2 = sp.ndimage.gaussian_filter(array_2,2)


#%%
plt.imshow(gauss_array_1)
plt.show()
plt.imshow(gauss_array_2)
plt.show()


#%%
rot_array_1 =  sp.ndimage.rotate(gauss_array_1, 90-angle_1, reshape=True)
rot_array_2 =  sp.ndimage.rotate(gauss_array_2, 90-angle_2, reshape=True)


#%%
plt.imshow(rot_array_1)
plt.show()
plt.imshow(rot_array_2)
plt.show()


#%%
sep_1 = np.sqrt(((len(array_1)-x1)-x1)**2 + ((len(array_1)-y1) - y1)**2)
print(sep_1)
scale_fac_1 = 100.0/sep_1
sep_2 = np.sqrt(((len(array_2)-x2)-x2)**2 + ((len(array_2)-y2) - y2)**2)
print(sep_2)
scale_fac_2 = 100.0/sep_2


#%%
print(scale_fac_1)
print(scale_fac_2)


#%%
rescaled_array_1 = sp.ndimage.zoom(rot_array_1,scale_fac_1)
plt.imshow(rescaled_array_1)
plt.show()
test = np.pad(rescaled_array_1,30,mode="constant")
plt.imshow(test)
plt.show()


#%%
rescaled_array_2 = sp.ndimage.zoom(rot_array_2,scale_fac_2)
plt.imshow(rescaled_array_2)


#%%
centre_1 = [len(rescaled_array_1)/2,len(rescaled_array_1)/2]
print(centre_1)
centre_2 = [len(rescaled_array_2)/2,len(rescaled_array_2)/2]
print(centre_2)


#%%
np.shape(rescaled_array_1)


#%%
np.shape(rescaled_array_2)


#%%
abs(len(rescaled_array_2)-(len(rescaled_array_2)/2 + 70))


#%%



#%%
output_1 = galaxy_pairs.get_subarray(rescaled_array_1,centre_1,min(len(rescaled_array_1),70))
output_2 = galaxy_pairs.get_subarray(rescaled_array_2,centre_2,min(len(rescaled_array_2),70))


#%%
plt.imshow(output_1)
plt.show()
plt.imshow(output_2)
plt.show()


#%%
final_output = output_1+output_2
plt.imshow(final_output)


#%%
np.shape(output_1)


#%%
np.shape(output_2)


#%%




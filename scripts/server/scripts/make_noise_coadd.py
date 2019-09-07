import glob
import os
import pdb

from sptpol_software.util import files
import sptpol_software.observation as obs
from sptpol_software.observation import *
import sptpol_software.analysis.maps as maps
import numpy as np


def removeTWeight(init_map,in_place=True):
    if in_place:
        new_map = init_map
    else:
        new_map = init_map.copy()

    # If the map is unweighted, leave it be
    if not getattr(new_map,'weighted_map',False):
	return(new_map)

    ind = new_map.weight >0
    new_map.map[ind] /= new_map.weight[ind]
    
    new_map.weighted_map = False
    new_map.setMapAttr('weighted_map',False)

    return(new_map)


def make_noise():
    
    
    dir='/data56/lmocanu/lens500d/ra0hdec-57p5/IDF_V3/03032016'
    #os.chdir('/data56/lmocanu/lens500d/ra0hdec-57p5/IDF_V3/03032016/left')
    #cwd = os.getcwd()
    #print cwd

    output_dir = '/home/mdz/data/coadds/noise'

    base = 'ra0dec-57p5'
#    sum_fname = output_dir+base+'_sum'
#    dif_fname = output_dir+base+'_dif'
    fstrs = ['_090ghz','_150ghz']
    nfiles = len(fstrs)

    j = 0 
    realisation_size = 100
    first_list = sorted(glob.glob(dir+'/[left,right]*/*_090ghz.h5'))
    second_list = sorted(glob.glob(dir+'/[left,right]*/*_150ghz.h5'))

	#slist = glob.glob('*'+fstrs[i]+'.h5')
        #right_list = glob.glob('*'+fstrs[i]+'.h5' )  
	#dmap = 0
#	print flist
    num_realisations = 20
    for mapcnt in range(len(first_list)):
	print '%s of %s: %s' %(mapcnt, len(first_list), first_list[mapcnt])
	tmp_90 = files.read(first_list[mapcnt])
	print '%s of %s: %s' %(mapcnt, len(second_list), second_list[mapcnt])
	tmp_150 = files.read(second_list[mapcnt])
        if (j==0):
	    j = 1
	    noise_map = (tmp_150 - tmp_90)
	    noise_map.weight = (tmp_150.weight-tmp_90.weight)
	elif (((mapcnt%realisation_size) == 0) and (mapcnt != 0)):
            print 'removing weight'
            noise_map.removeWeight()
	    #pdb.set_trace()
            print 'writing realisation %s to hf5' %(mapcnt/realisation_size) 
     	    noise_map.writeToHDF5(output_dir+"/noise_real_"+str(mapcnt/realisation_size)+'.h5',overwrite=True,use_compression=False)
            j=0
	    continue
	elif (mapcnt%8000 == 0):
	    break
        else:
            noise_map += (tmp_150-tmp_90)
            noise_map.weight += (tmp_150.weight-tmp_90.weight)
         #   dmap = dmap + tmp
         #   dmap2=dmap.copy()
         #   dmap2=removeTWeight(dmap2)
         #   dmap20=dmap10.copy()
         #   dmap20=removeTWeight(dmap20)
         #   dmap3=removeTWeight(tmp)
         #   rmses[j]=np.std(dmap2.map[2600:3000,1400:1600])
         #   print j,rmses[j],np.std(dmap20[2600:3000,1400:1600]),np.std(dmap3[2600:3000,1400:1600]),rmses[0]/np.sqrt(j+1)
            j+=1
#        pdb.set_trace()
 
       # smap = removeTWeight(smap)

#        dmap = removeTWeight(dmap10)
#        dmap.removeWeight()
        #pdb.set_trace()

 #   try:
#	os.chdir(output_dir) 
        #os.mkdir(output_dir)
#    except OSError: 
        # If output directory already exists, do nothing
        # Output directory not found
#	print 'Output Directory Not Found'
#	pass
     
#    dmap10.writeToHDF5(dif_ofile+fstrs[i]+'.h5',overwrite=True,use_compression=False)
    


if __name__=="__main__":
     # If run from commandline, automatically execute
     make_noise()

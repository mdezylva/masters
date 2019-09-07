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


def make_coadd():
    
    
    dir='/data53/cr/ra0hdec-57p5/clusterlensing_maps_20170910/bundle_maps'
    #'/data56/lmocanu/lens500d/ra0hdec-57p5/IDF_V3/03032016'
    #os.chdir('/data56/lmocanu/lens500d/ra0hdec-57p5/IDF_V3/03032016/left')
    #cwd = os.getcwd()
    #print cwd

    output_dir = '/data53/cr/mitchell/'

    base = 'ra0dec-57p5'
    sum_fname = output_dir+base+'_sum'
    dif_fname = output_dir+base+'_dif'
    fstrs = ['_90ghz','_150ghz']
    nfiles = len(fstrs)

    j = 0 
    map_cutoff = 5000
    for i in range(nfiles):
	print fstrs[i]
       	flist = sorted(glob.glob(dir+'/*'+fstrs[i]+'_*.h5'))
	#slist = glob.glob('*'+fstrs[i]+'.h5')
        #right_list = glob.glob('*'+fstrs[i]+'.h5' )  
	#dmap = 0
	print flist
 
        for mapcnt,fname in enumerate(flist):
	    print '%s of %s: %s' %(mapcnt, len(flist), fname)
	    tmp = files.read(fname)
            if (mapcnt == 0):
		sum_map = tmp
		sum_map.weight = tmp.weight
	    elif (mapcnt == map_cutoff):
		break
            else:
                sum_map += tmp
                sum_map.weight += tmp.weight
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
        print 'removing weight'
 
       # smap = removeTWeight(smap)

#        dmap = removeTWeight(dmap10)
        sum_map.removeWeight()
#        dmap.removeWeight()
#        pdb.set_trace()

 #   try:
#	os.chdir(output_dir) 
        #os.mkdir(output_dir)
#    except OSError: 
        # If output directory already exists, do nothing
        # Output directory not found
#	print 'Output Directory Not Found'
#	pass
     
        print 'writing to hf5'
     	sum_map.writeToHDF5(sum_fname+str(map_cutoff)+fstrs[i]+'.h5',overwrite=True,use_compression=False)
#    dmap10.writeToHDF5(dif_ofile+fstrs[i]+'.h5',overwrite=True,use_compression=False)
    


if __name__=="__main__":
     # If run from commandline, automatically execute
     make_coadd()

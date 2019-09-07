import glob
import os
import pdb

from sptpol_software.util import files
import numpy as np

def removeTWeight(imap,in_place=True):
    if in_place:
        new_map = imap
    else:
        new_map = imap.copy()

        # If this map is already unweighted, we have nothing to do.  
    if not getattr(new_map, 'weighted_map', False):
        return new_map

    ind = new_map.weight > 0
    new_map.map[ind] /= new_map.weight[ind]
    
    new_map.weighted_map = False
    new_map.setMapAttr('weighted_map', False)

    return new_map


def coadd_szpol_pol():
    dir='/data56/lmocanu/lens500d/ra0hdec-57p5/IDF_V3/03032016'
    odir=dir+'coadds/'
    try: 
        os.mkdir(odir)
    except OSError: 
        #do nothing
        pass
    base='ra0hdec-57.5'
    sum_ofile=odir+base+'_sum_'
    dif_ofile=odir+base+'_dif_'
    fstrs=['150ghz','90ghz']
    nf=2

    for i in range(nf):

        slist = glob.glob(dir+'bundle_maps/'+base+'_'+fstrs[i]+'*.h5')
        dlist = glob.glob(dir+'bundle_lrjacks/'+base+'_'+fstrs[i]+'*.h5')


        
        smap = 0
        dmap = 0
        j=0
        for fname in slist:
            print fname
            #tmp = files.read(fname)
            #smap = smap + tmp
            #print j,np.std(smap.map[2600:3000,1400:1600])
            #print j,np.std(tmp[2600:3000,1400:1600])
            #print j,smap.weight[2800,1500],smap.weight[2800,1500]/(j+1)
            
            j=j+1

        N=len(slist)
        rmses=np.zeros(N)
        j=0
        dmap=0

        for fname in dlist:
            print fname
            tmp = files.read(fname)
            if (j == 0):
                dmap10=tmp.copy()
            else:
                dmap10.map=dmap10.map+tmp.map
                dmap10.weight=dmap10.weight+tmp.weight
            dmap = dmap + tmp
            dmap2=dmap.copy()
            dmap2=removeTWeight(dmap2)
            dmap20=dmap10.copy()
            dmap20=removeTWeight(dmap20)
            dmap3=removeTWeight(tmp)
            rmses[j]=np.std(dmap2.map[2600:3000,1400:1600])
            print j,rmses[j],np.std(dmap20[2600:3000,1400:1600]),np.std(dmap3[2600:3000,1400:1600]),rmses[0]/np.sqrt(j+1)
            j=j+1
        pdb.set_trace()
        print 'removing weight'
 
        smap = removeTWeight(smap)

        dmap = removeTWeight(dmap)
#        smap.removeWeight()
#        dmap.removeWeight()
        pdb.set_trace()

        print 'writing to hf5'
        smap.writeToHDF5(sum_ofile+fstrs[i]+'.h5',overwrite=True,use_compression=False)
        dmap.writeToHDF5(dif_ofile+fstrs[i]+'.h5',overwrite=True,use_compression=False)
        
			
if __name__=="__main__":
    coadd_szpol_pol()

import sptpol_software.util.files as files
import sptpol_software.analysis.maps as maps
import pickle
import gzip
import numpy as np
import glob
import sys
import os
import time

#fdpath = '/data57/jtsayre/ra0hdec-57p5/500dBB_1p5arcmin/left'
fdpath = '/data57/jtsayre/ra0hdec-57p5/500dBB_1p5arcmin/'
maplist = sorted(glob.glob('%s/[left,right]*/*150ghz.h5' % (fdpath)))
# print len(maplist), sys.exit()
for mapcnt, currmapname in enumerate(maplist):
    print '%s of %s: %s' % (mapcnt, len(maplist), currmapname)
    currmap = files.read(currmapname)
    if mapcnt == 0:
        BUNDLES = currmap
    else:
        BUNDLES += currmap

BUNDLES.removeWeight()
print 'done'

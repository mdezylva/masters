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
import healpy as hp
import os
os.chdir("/home/mitchell/Documents/masters/masters/scripts/")
import cmb
import galaxy_pairs
#data_location = input("Enter Path Location of data")
os.chdir("/home/mitchell/Documents/masters/masters/data/simulated")
cwd = os.getcwd()

import scipy.ndimage
cosmo = Planck15


tsz_150ghz = files.read("tsz150_R13.fits")


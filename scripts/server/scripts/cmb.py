from sptpol_software.util import files
import sptpol_software.observation.sky
import sptpol_software.observation as obs
from sptpol_software.observation import *
from sptpol_software.util.tools import stat
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
import sptpol_software as sps
from astropy.io import fits
from astropy.cosmology import Planck13
from astropy import units as u
import scipy.spatial.distance as dist
print(Planck13)
cosmo = Planck13


def convert_ghz_to_y(freq_power):
    '''
    Converts frequency power to Compton Y parameter
    '''
    x = freq_power / 56.85  # x = h v / k_B T_CMB
    return((x / np.tanh(x / 2)) - 4)


def get_y_map(maps, cal_factors, frequencies, save=False):
    '''
    Takes an array of maps, and calibration factors and converts it to a Compton y-map
    '''
    assert(len(maps) == len(cal_factors))
    assert(len(frequencies) == len(maps))

    for index in range(len(maps)):
        maps[index] = maps[index] * cal_factors[index]

    difference_map = maps[0] - maps[1]

    freq_scaling = 1.0 / \
        (convert_ghz_to_y(frequencies[0]) - convert_ghz_to_y(frequencies[1]))

    y_map = freq_scaling * difference_map

    #y_map_array = y_map.getTOnly().map

    if save:
        np.save('y_map', y_map_array)

    return(y_map_array)

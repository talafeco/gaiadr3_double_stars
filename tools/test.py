import os
import numpy as np
import math
import sys
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from astropy.table import Column, MaskedColumn
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
import warnings
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta

pmval = float(sys.argv[1]) / 100

def calcPmCategory(pmfact):
    print(pmfact)
    pmCommon = ()
    if pmfact >= 0.8:
        pmCommon = 'CPM'
    elif 0.4 <= pmfact < 0.8:
        pmCommon = 'SPM'
    elif pmfact < 0.4:
        pmCommon = 'DPM'
    return pmCommon

print(calcPmCategory(pmval))
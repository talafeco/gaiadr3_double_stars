#! /usr/bin/python3

# Usage: phot_gaia_data <image_file>

import os
import sys
import numpy as np
import math
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, FK5, ICRS
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta


# Configuration for the ATIK camera

dummyObservationDate = "2022-01-01T12:00:00"
file_name = sys.argv[1]

fitsFileDate = ''
fitsHeader = fits.open(file_name)[0].header
key_to_lookup = 'DATE-OBS'
if key_to_lookup in fitsHeader:
    fitsFileDate = fitsHeader['DATE-OBS']
else:
    fitsFileDate = np.nan

mywcs = WCS(fitsHeader)

photo_left_upper = mywcs.pixel_to_world(0, 0, 1)
photo_left_lower = mywcs.pixel_to_world(fitsHeader['NAXIS2'], 0, 1)
photo_right_upper = mywcs.pixel_to_world(0, fitsHeader['NAXIS1'], 1)
photo_right_lower = mywcs.pixel_to_world(fitsHeader['NAXIS2'], fitsHeader['NAXIS1'], 1)

print(photo_left_upper, photo_left_lower, photo_right_upper, photo_right_lower)
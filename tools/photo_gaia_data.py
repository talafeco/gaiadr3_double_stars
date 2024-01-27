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
from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
Gaia.ROW_LIMIT = -1 # To return an unlimited number of rows

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

photo_left_upper = SkyCoord.from_pixel(0, 0, mywcs, origin=0, mode='all')
photo_left_lower = SkyCoord.from_pixel(fitsHeader['NAXIS2'], 0, mywcs, origin=0, mode='all')
photo_right_upper = SkyCoord.from_pixel(0, fitsHeader['NAXIS1'], mywcs, origin=0, mode='all')
photo_right_lower = SkyCoord.from_pixel(fitsHeader['NAXIS2'], fitsHeader['NAXIS1'], mywcs, origin=0, mode='all')
photo_center = SkyCoord(fitsHeader['RA'] * u.degree, fitsHeader['DEC'] * u.degree)

photo_width = photo_left_upper.separation(photo_right_upper)
photo_height = photo_left_upper.separation(photo_left_lower)

print(photo_center.to_string('hmsdms'), photo_width, photo_height)

gaia_photo_catalog = Gaia.query_object_async(photo_center, photo_width / 2, photo_height / 2)

print(gaia_photo_catalog)
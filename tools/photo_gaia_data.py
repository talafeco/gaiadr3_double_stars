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

file_name = sys.argv[1]

fitsHeader = fits.open(file_name)[0].header
mywcs = WCS(fitsHeader)

photo_left_upper = SkyCoord.from_pixel(0, 0, mywcs, origin=0, mode='all')
#photo_left_lower = SkyCoord.from_pixel(fitsHeader['NAXIS2'], 0, mywcs, origin=0, mode='all')
#photo_right_upper = SkyCoord.from_pixel(0, fitsHeader['NAXIS1'], mywcs, origin=0, mode='all')
photo_right_lower = SkyCoord.from_pixel(fitsHeader['NAXIS2'], fitsHeader['NAXIS1'], mywcs, origin=0, mode='all')
photo_center = SkyCoord(fitsHeader['RA'] * u.degree, fitsHeader['DEC'] * u.degree)

#photo_width = photo_left_upper.separation(photo_right_upper)
#photo_height = photo_left_upper.separation(photo_left_lower)
#photo_diagonal = photo_left_upper.separation(photo_right_lower)
photo_radius = photo_left_upper.separation(photo_right_lower) / 2

print('Center of photo (hours / decimal degree): ', photo_center.to_string('hmsdms'), '/', photo_center.to_string('decimal'),
      #'\nCenter of photo (decimal degree): ', photo_center.to_string('decimal'),
      #'\nWidth of photo: ', photo_width,
      #'\nHeight of photo: ', photo_height,
      '\nRadius of photo: ', photo_radius)

# gaia_photo_catalog = Gaia.query_object_async(photo_center, photo_width / 2, photo_height / 2)
gaia_photo_catalog = Gaia.cone_search_async(photo_center, radius=u.Quantity(photo_radius))
r = gaia_photo_catalog.get_results()

print(r)
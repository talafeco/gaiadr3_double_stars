#! /usr/bin/python3

# Usage: gdr3_comp <Image_folder>

import sys
import os
import numpy as np
import math
from astropy.io import fits
from astropy.table import QTable, Table, Column
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky, search_around_sky, position_angle, angular_separation
from astropy.table import Table, vstack, hstack
import astropy.units as u
from astropy.time import Time, TimeDelta

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
Gaia.ROW_LIMIT = -1 # To return an unlimited number of rows

# Define the directory containing the FITS files
workingDirectory = sys.argv[1]

# Define the magnitude limit of the image used during Gaia DR3 queries
gaia_dr3_magnitude_limit = 18
dao_sigma = 3.0
dao_fwhm = 8.0
dao_threshold = 12.0
possible_distance = 30000.0 # AU
search_cone = 0.001 # Decimal degree


directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and (f.endswith('.new') or f.endswith('.fit') or f.endswith('.fits'))]
print(files)


def extract_sources(fits_file):
    # Get file base data
    hdu = fits.open(fits_file)
    mywcs = WCS(hdu[0].header)
    file_header = hdu[0].header
    
    # Set observation date and time
    fitsFileDate = ''
    key_to_lookup_a = 'DATE-OBS'
    key_to_lookup_b = 'DATE'
    if key_to_lookup_a in file_header:
        fitsFileDate = file_header['DATE-OBS']
    elif key_to_lookup_b in file_header:
        fitsFileDate = file_header['DATE']
    else:
        fitsFileDate = np.nan
        
    # Estimate the background and background noise
    data = hdu[0].data
    mean, median, std = sigma_clipped_stats(data, sigma = dao_sigma)  

    daofind = DAOStarFinder(fwhm=dao_fwhm, threshold=dao_threshold*std)  
    sources = daofind(data - median)
    ra2, dec2 = mywcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 1)
    
    sources.add_column(ra2, name='ra_deg') 
    sources.add_column(dec2, name='dec_deg')
    sources.add_column(fitsFileDate, name='image_date')
    sources.add_column(fits_file, name='file')
    
    return sources

def create_master_catalog():
    fits_file = workingDirectory + '/' + files[0]
    hdu = fits.open(fits_file)
    mywcs = WCS(hdu[0].header)
    file_header = hdu[0].header
    photo_left_upper = SkyCoord.from_pixel(0, 0, mywcs, origin=0, mode='all')
    photo_right_lower = SkyCoord.from_pixel(file_header['NAXIS2'], file_header['NAXIS1'], mywcs, origin=0, mode='all')
    photo_center = SkyCoord(file_header['CRVAL1'] * u.degree, file_header['CRVAL2'] * u.degree)
    photo_radius = photo_left_upper.separation(photo_right_lower) / 2
    print('Center of photo: ', photo_center.to_string('hmsdms'), '/', photo_center.to_string('decimal'),
      '\nRadius of photo: ', photo_radius)
    gaia_photo_catalog = Gaia.cone_search_async(photo_center, radius=u.Quantity(photo_radius) * 1.1)
    gaia_stars = gaia_photo_catalog.get_results()
    
    return gaia_stars

master_catalog = create_master_catalog()

file_counter = 0

for fitsFile in files:
    # 1. Read the list of sources extracted from an image (fits) file
    file_counter = file_counter + 1
    print('\n\n### Processing file', file_counter, 'out of', len(files),': ', fitsFile, '###')
    sources = extract_sources(workingDirectory + '/' + fitsFile)

    sources_catalog = SkyCoord(ra=sources['ra_deg']*u.degree, dec=sources['dec_deg']*u.degree, frame='fk5')
    idxw, idxs, wsd2d, wsd3d = search_around_sky(master_catalog, sources_catalog, search_cone*u.deg)
    composit_catalog = hstack([master_catalog[idxw], sources[idxs]])
    print('### Composit Catalog No. ', file_counter, '\n', composit_catalog)


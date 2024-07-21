#! /usr/bin/python3

# Usage: gdr3_comp <Image_folder>

'''
To achieve the task, you'll need to combine several steps using astropy, numpy, and astroquery. Below is a detailed Python script to perform the required tasks:

Read multiple FITS files from a directory.
Extract sources and their magnitudes from the images.
Estimate background and compare source magnitudes to the background.
Define source coordinates using the WCS information from the FITS images.
Query Gaia DR3 for sources in the same area.
Summarize information in a QTable.
Calculate the background limiting magnitude.
Compare measured magnitudes to Gaia DR3 magnitudes and collect any major differences.

Explanation:
Reading FITS Files: The script iterates through all .fits files in the specified directory.
Extracting Sources: Uses DAOStarFinder to find sources in each image and calculate their WCS coordinates.
Querying Gaia DR3: Uses astroquery to fetch Gaia DR3 data for the region covered by each image.
Calculating Background Limiting Magnitude: Based on the background standard deviation and the number of sources.
Comparing Magnitudes: Calculates the difference between measured magnitudes and Gaia magnitudes, identifying significant discrepancies.
Storing Results: Aggregates all sources and major differences into QTables and saves them as FITS files.
'''

import sys
import os
import numpy as np
from astropy.io import fits
from astropy.table import QTable, Table
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from astroquery.gaia import Gaia
import astropy.units as u

# Define the directory containing the FITS files
fits_dir = sys.argv[1]

# Initialize lists to store results
all_sources = []
major_differences = []

# Function to extract sources from a FITS image
def extract_sources(fits_file):
    with fits.open(fits_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)
        
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
        sources = daofind(data - median)
        
        if sources is None:
            return None, None, None
        
        sources['ra'], sources['dec'] = wcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 0)
        return sources, wcs, data

# Function to query Gaia DR3 for sources in the same area
def query_gaia_sources(wcs, data_shape):
    ny, nx = data_shape
    corners = np.array([[0, 0], [0, ny], [nx, 0], [nx, ny]])
    ra_dec_corners = wcs.all_pix2world(corners, 0)
    
    ra_min, dec_min = np.min(ra_dec_corners, axis=0)
    ra_max, dec_max = np.max(ra_dec_corners, axis=0)
    
    query = f"""
    SELECT * FROM gaiadr3.gaia_source
    WHERE ra BETWEEN {ra_min} AND {ra_max}
    AND dec BETWEEN {dec_min} AND {dec_max}
    """
    job = Gaia.launch_job(query)
    gaia_sources = job.get_results()
    return gaia_sources

# Function to calculate background limiting magnitude
def calculate_limiting_magnitude(data, std, wcs, num_sources):
    bkg_mag = -2.5 * np.log10(std / np.sqrt(num_sources))
    return bkg_mag

# Read FITS files and process each one
for fits_file in os.listdir(fits_dir):
    if fits_file.endswith('.new'):
        fits_path = os.path.join(fits_dir, fits_file)
        
        sources, wcs, data = extract_sources(fits_path)
        if sources is None:
            continue
        
        gaia_sources = query_gaia_sources(wcs, data.shape)
        sources['gaia_id'] = ''
        sources['gaia_g_mag'] = np.nan
        
        for i, source in enumerate(sources):
            ra, dec = source['ra'], source['dec']
            idx = np.argmin(np.sqrt((gaia_sources['ra'] - ra)**2 + (gaia_sources['dec'] - dec)**2))
            
            gaia_source = gaia_sources[idx]
            print('Gaia Source ID: ', gaia_source['SOURCE_ID'])
            
            sources['gaia_id'][i] = gaia_source['SOURCE_ID']
            sources['gaia_g_mag'][i] = gaia_source['phot_g_mean_mag']
        
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        bkg_limiting_mag = calculate_limiting_magnitude(data, std, wcs, len(sources))
        sources['bkg_limiting_mag'] = bkg_limiting_mag
        
        sources['measured_mag'] = -2.5 * np.log10(sources['flux'] / bkg_limiting_mag)
        sources['mag_difference'] = sources['measured_mag'] - sources['gaia_g_mag']
        
        major_diff = sources[np.abs(sources['mag_difference']) > 1]
        major_differences.append(major_diff)
        
        all_sources.append(sources)

# Convert results to QTables
all_sources_table = QTable(np.hstack(all_sources))
major_differences_table = QTable(np.hstack(major_differences))

# Save results to files
# all_sources_table.write('all_sources_summary.fits', format='fits', overwrite=True)
all_sources_table.write('all_sources_summary.txt',  format='ascii', overwrite=True, delimiter=',')

# major_differences_table.write('major_differences_summary.fits', format='fits', overwrite=True)
major_differences_table.write('major_differences_summary.txt', format='ascii', overwrite=True, delimiter=',')
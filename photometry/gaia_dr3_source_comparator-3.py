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

# Feladatok:
# Hozzáadni a Gaia koordináták és a mérési koordináták különbségét

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
from astropy.coordinates import match_coordinates_sky
import astropy.units as u
from astropy.time import Time, TimeDelta

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
Gaia.ROW_LIMIT = -1 # To return an unlimited number of rows

# Define the directory containing the FITS files
fits_dir = sys.argv[1]

# Define the magnitude limit of the image used during Gaia DR3 queries
gaia_dr3_magnitude_limit = 18

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
        # daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
        daofind = DAOStarFinder(fwhm=8.0, threshold=12.0*std)
        sources = daofind(data - median)
        
        if sources is None:
            return None, None, None
        
        sources['ra'], sources['dec'] = wcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 0)

        photo_left_upper = SkyCoord.from_pixel(0, 0, wcs, origin=0, mode='all')
        photo_right_lower = SkyCoord.from_pixel(header['NAXIS2'], header['NAXIS1'], wcs, origin=0, mode='all')
        photo_center = SkyCoord(header['CRVAL1'] * u.degree, header['CRVAL2'] * u.degree)
        photo_radius = photo_left_upper.separation(photo_right_lower) / 2

        

        return sources, wcs, data, photo_center, photo_radius

# Function to query Gaia DR3 for sources in the same area
'''def query_gaia_sources(wcs, data_shape, limit):
    ny, nx = data_shape
    corners = np.array([[0, 0], [0, ny], [nx, 0], [nx, ny]])
    ra_dec_corners = wcs.all_pix2world(corners, 0)
    
    ra_min, dec_min = np.min(ra_dec_corners, axis=0)
    ra_max, dec_max = np.max(ra_dec_corners, axis=0)

    gaia_mag_limit = limit
    
    query = f"""
    SELECT * FROM gaiadr3.gaia_source
    WHERE ra BETWEEN {ra_min} AND {ra_max}
    AND dec BETWEEN {dec_min} AND {dec_max}
    """
    job = Gaia.launch_job(query)
    gaia_sources = job.get_results()
    return gaia_sources'''

# Function to query Gaia DR3 for sources in the same area
def query_gaia_sources(photo_center, photo_radius):
    print('Center of photo: ', photo_center.to_string('hmsdms'), '/', photo_center.to_string('decimal'),
      '\nRadius of photo: ', photo_radius)
    gaia_photo_catalog = Gaia.cone_search_async(photo_center, radius=u.Quantity(photo_radius))
    gaia_sources = gaia_photo_catalog.get_results()
    return gaia_sources


# Function to calculate background limiting magnitude
def calculate_limiting_magnitude(measured_mags, dr3_mags):
    # bkg_mag = -2.5 * np.log10(math.fabs(std) / np.sqrt(num_sources))
    bkg_mag = np.mean(dr3_mags - measured_mags)
    return bkg_mag

# Function to calculate Star positions based on Gaia DR3 coordinates and proper motion
def calcCurrentDR3Coord(date, star_ra, star_dec, star_pr_ra, star_pr_dec):
    date_time = Time(date, format='jyear')
    star = SkyCoord(ra = star_ra * u.degree,
                dec = star_dec * u.degree,
                pm_ra_cosdec = star_pr_ra * u.mas/u.yr,
                pm_dec = star_pr_dec * u.mas/u.yr,
                frame = 'icrs',
                obstime=Time('2016-01-01 00:00:00.0'))
    starUpdCoord = star.apply_space_motion(new_obstime=date_time)
    return starUpdCoord

# Get Gaia DR3 catalog data based on the first image WCS coordinates
fits_master = os.listdir(fits_dir)[0]
fits_master_file = os.path.join(fits_dir, fits_master)
master_sources, master_wcs, master_data, master_photo_center, master_radius = extract_sources(fits_master)
# gaia_sources = query_gaia_sources(master_wcs, master_data.shape, gaia_dr3_magnitude_limit)

gaia_sources = query_gaia_sources(master_photo_center, master_radius * 1.1)
print('### gaia_sources: ', len(gaia_sources))
gaia_catalog = SkyCoord(gaia_sources['ra']*u.degree, gaia_sources['dec']*u.degree, frame='fk5')

# Read FITS files and process each one
for fits_file in os.listdir(fits_dir):
    if fits_file.endswith('.new'):
        fits_path = os.path.join(fits_dir, fits_file)
        sources, wcs, data, photo_center, radius = extract_sources(fits_path)
        if sources is None:
            continue

        # print('### SOURCES TALBE ###\n', sources)

        idx, d2d, d3d = sources.match_to_catalog_sky(gaia_catalog)

        sources.add_column(np.empty, name='gaia_designation')
        sources.add_column(np.empty, name='gaia_source_id')
        sources.add_column(np.empty, name='object_id')
        sources['gaia_bp_mag'] = 0.0
        sources['gaia_g_mag'] = 0.0
        sources['gaia_rp_mag'] = 0.0
        sources['gaia_ra'] = 0.0
        sources['gaia_dec'] = 0.0
        sources.add_column(fits_file, name='file_name')
        
        # for i, source in enumerate(sources):
            # ra, dec = source['ra'], source['dec']
            # idx = np.argmin(np.sqrt((gaia_sources['ra'] - ra)**2 + (gaia_sources['dec'] - dec)**2))
            # print('### IDX 1. ###\n', idx)
            # gaia_source = gaia_sources[idx]
            # sources['gaia_id'][i] = gaia_source['DESIGNATION']
            # sources['gaia_g_mag'][i] = gaia_source['phot_g_mean_mag']
        
        '''
        for i, source in enumerate(sources):
            source_ra, source_dec = wcs.all_pix2world([[source ['xcentroid'], source ['ycentroid']]], 0)[0]   
            source_catalog_coords = SkyCoord(ra=source_ra*u.degree, dec=source_dec*u.degree, frame='fk5')  
            
            idx, d2d, d3d = match_coordinates_sky(source_catalog_coords, gaia_catalog)

            gaia_source = gaia_sources[idx]
            '''
            # sep_source = SkyCoord(gaia_source['ra']*u.degree, gaia_source['dec']*u.degree, frame='fk5')
            # sep_target = SkyCoord(ra*u.degree, dec*u.degree, frame='fk5')
            # separation_final = sep_source.separation(sep_target)*u.arcsecond
            
            
            
            
        print('gaia_star: ', gaia_source['DESIGNATION'], '\n', idx, sep_target)
        #print('gaia_star: ', gaia_source['DESIGNATION'], '\n', idx, d2d)
        sources['gaia_designation'] = gaia_source['DESIGNATION']
        sources['gaia_source_id'] = str(gaia_source['SOURCE_ID'])
        sources['object_id'] = str(gaia_source['SOURCE_ID']) + str(source['file_name'])
        sources['gaia_g_mag'] = gaia_source['phot_g_mean_mag']
        sources['gaia_bp_mag'] = gaia_source['phot_bp_mean_mag']
        sources['gaia_rp_mag'] = gaia_source['phot_rp_mean_mag']
        sources['gaia_ra'] = gaia_source['ra']
        sources['gaia_dec'] = gaia_source['dec']

            # sources['gaia_bp_rp'][i] = gaia_source['bp_rp']
            # sources['gaia_bp_g'][i] = gaia_source['bp_g']
            # sources['gaia_g_rp'][i] = gaia_source['g_rp']


        
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        bkg_limiting_mag = calculate_limiting_magnitude(sources['mag'], sources['gaia_g_mag'])
        sources['bkg_limiting_mag'] = bkg_limiting_mag
        sources['measured_mag'] = bkg_limiting_mag + sources['mag']
        sources['mag_difference'] = sources['gaia_g_mag'] - sources['measured_mag']
        
        major_diff = sources[np.abs(sources['mag_difference']) > 1]
        major_differences.append(major_diff)
        
        all_sources.append(sources)

# Convert results to QTables
all_sources_table = QTable(np.hstack(all_sources))
major_differences_table = QTable(np.hstack(major_differences))

# Create table for the aggregated data
designation, sourceid, flux, mag, ra_deg, de_deg, dif_arcsec, bp_mag, g_mag, rp_mag, background_mag, measured_mag, measured_mag_err, mag_difference, mag_difference_err = np.array([], dtype=str), np.array([], dtype=str), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64), np.array([], dtype=np.float64)

aggregated_table = Table([designation, sourceid, flux, mag, ra_deg, de_deg, dif_arcsec, bp_mag, g_mag, rp_mag, background_mag, measured_mag, measured_mag_err, mag_difference, mag_difference_err], names=('Gaia DR3 Designation', 'Gaia DR3 Source ID', 'Flux', 'Magnitude', 'RA', 'DEC', 'Position Deviation', 'Gaia DR3 Bp Magnitude', 'Gaia DR3 G Magnitude', 'Gaia DR3 Rp Magnitude', 'Background Limiting Magnitude', 'Measured Magnitude', 'Std. Error of Measured Magnitude', 'Magnitude Difference', 'Std. Error of Magnitude Difference'))

# Average measurements
reportTable_by_object = all_sources_table.group_by('gaia_source_id')

# print('reportTable_by_object: ', reportTable_by_object.info)

count = 1
for star in reportTable_by_object.groups:
    # print(len(star['gaia_designation']))
    # print(len(os.listdir(fits_dir)) / 2)

    if len(star['gaia_designation']) >= len(os.listdir(fits_dir)) / 2:
        
        '''print('\n### Star index:', count, '###')
        count = count + 1
        print('\nStar DR3 Designation: ', star['gaia_designation'][0])
        print('\nGaia G Magnitude: ', star['gaia_g_mag'].groups.aggregate(np.mean))
        print('\nMeasured Magnitude: ', star['measured_mag'].groups.aggregate(np.mean))
        print('\nMagnitude Difference: ', star['mag_difference'].groups.aggregate(np.mean))'''

        dr3_current_coord = SkyCoord(ra = star[0]['gaia_ra'] * u.degree, dec = star[0]['gaia_dec'] * u.degree)
        measured_coord = SkyCoord(ra = star['ra'].groups.aggregate(np.mean) * u.degree, dec = star['dec'].groups.aggregate(np.mean) * u.degree,)
        coord_error = dr3_current_coord.separation(measured_coord).arcsecond

        print('Star: ', star)

        aggregated_table.add_row([star[0]['gaia_designation'], star[0]['gaia_source_id'], star['flux'].groups.aggregate(np.mean), star['mag'].groups.aggregate(np.mean), star['ra'].groups.aggregate(np.mean), star['dec'].groups.aggregate(np.mean), coord_error, star[0]['gaia_bp_mag'], star[0]['gaia_g_mag'], star[0]['gaia_rp_mag'], star['bkg_limiting_mag'].groups.aggregate(np.mean), star['measured_mag'].groups.aggregate(np.mean), star['measured_mag'].groups.aggregate(np.std), star['mag_difference'].groups.aggregate(np.mean), star['mag_difference'].groups.aggregate(np.std)])

        # print('\nStar DR3 Designation: ', star['gaia_designation'])
        # print('\nStar DR3 Designation: ', star['gaia_designation'])


# Save results to files
# all_sources_table.write('all_sources_summary.fits', format='fits', overwrite=True)
all_sources_table.write('all_sources_summary.csv',  format='ascii', overwrite=True, delimiter=',')

# major_differences_table.write('major_differences_summary.fits', format='fits', overwrite=True)
major_differences_table.write('major_differences_summary.csv', format='ascii', overwrite=True, delimiter=',')

aggregated_table.write('aggregated_table.csv', format='ascii', overwrite=True, delimiter=',')


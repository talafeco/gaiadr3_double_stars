#! /usr/bin/python3

# Usage: gdr3-comp <Image_folder>

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
from astropy.utils.masked import Masked
import matplotlib.pyplot as plt
import matplotlib as mpl


Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default
Gaia.ROW_LIMIT = -1 # To return an unlimited number of rows

# Define the directory containing the FITS files
workingDirectory = sys.argv[1]

# Define the magnitude limit of the image used during Gaia DR3 queries
gaia_dr3_magnitude_limit = 20
dao_sigma = 3.0
dao_fwhm = 8.0
dao_threshold = 12.0
search_cone = 0.00028 # Decimal degree
image_limit = 2000

directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and (f.endswith('.new') or f.endswith('.fit') or f.endswith('.fits'))]
print(files)

file_counter = 0


def extract_sources(fits_file):
    # Get file base data
    hdu = fits.open(fits_file)
    mywcs = WCS(hdu[0].header)
    file_header = hdu[0].header
    # image = hdu[0].data
    
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
    
    return sources, hdu, mywcs

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

# Function to calculate background limiting magnitude
def calculate_limiting_magnitude(measured_mags, dr3_mags):
    # bkg_mag = -2.5 * np.log10(math.fabs(std) / np.sqrt(num_sources))
    bkg_mag = np.mean(dr3_mags.data) - np.mean(measured_mags)

    return bkg_mag

def normalize_colors(colors):
    color_min = colors.min()
    color_max = colors.max()
    
    # Normalize to the range -1 to 1
    normalized_range = [(2 * (num - color_min) / (color_max - color_min)) - 1 for num in colors]
    
    return normalized_range



def plot_sources(wcs, hdu, source_list, composite_list, filename):
    ax = plt.subplot(projection=wcs, label='overlays')
    ax.imshow(hdu[0].data, vmax=image_limit, vmin=0, cmap='Greys')
    ax.coords.grid(True, color='white', ls='solid')
    ax.coords[0].set_axislabel('Right Ascension (J2000)')
    ax.coords[1].set_axislabel('Declination (J2000)')
    ax.scatter(source_list['xcentroid'], source_list['ycentroid'], s=20, edgecolor="white", color='none', linewidths=1)
    # Calculate color map of identified sources based on the magnitude difference
    normalized_color_scale = normalize_colors(composite_list['mag_difference'])
    ax.scatter(composite_list['xcentroid'], composite_list['ycentroid'], s=20, facecolor='none', c=normalized_color_scale, cmap='jet', alpha=0.4, linewidths=1)
    im = ax.scatter(composite_list['xcentroid'], composite_list['ycentroid'], s=20, facecolor='none', c=normalized_color_scale, cmap='jet', alpha=0.4, linewidths=1)
    plt.colorbar(im)
    plt.title(filename)
    plt.savefig(str(workingDirectory + '/' + filename + '.jpg').replace(' ', ''), dpi=300.0, bbox_inches='tight', pad_inches=0.2)
    plt.close()

def plot_sum(wcs, hdu, sum_list, filename):
    ax = plt.subplot(projection=wcs, label='overlays')
    ax.imshow(hdu[0].data, vmax=image_limit, vmin=0, cmap='Greys')
    ax.coords.grid(True, color='white', ls='solid')
    ax.coords[0].set_axislabel('Right Ascension (J2000)')
    ax.coords[1].set_axislabel('Declination (J2000)')
    # Calculate color map of identified sources based on the magnitude difference
    normalized_color_scale = normalize_colors(sum_list['mag_difference'])
    ax.scatter(sum_list['ra_deg'], sum_list['dec_deg'], s=20, facecolor='none', c=normalized_color_scale, cmap='jet', alpha=0.4, linewidths=1, transform=ax.get_transform('fk5'))
    im = ax.scatter(sum_list['ra_deg'], sum_list['dec_deg'], s=20, facecolor='none', c=normalized_color_scale, cmap='jet', alpha=0.4, linewidths=1, transform=ax.get_transform('fk5'))
    plt.colorbar(im)
    plt.title(filename)
    plt.savefig(filename, dpi=300.0, bbox_inches='tight', pad_inches=0.2)
    plt.close()

gaia_catalog = create_master_catalog()

measured_objects = Table()

#measurements_table_dtype = {'designation': np.str_, 'ra': np.float64, 'dec': np.float64, 'parallax': np.float64, 'parallax_error': np.float64, 'pm': np.float64, 'pmra': np.float64, 'pmdec': np.float64, 'ruwe': np.float64, 'phot_g_mean_mag': np.float64, 'phot_bp_mean_mag': np.float64, 'phot_rp_mean_mag': np.float64, 'bp_rp': np.float64, 'bp_g': np.float64, 'g_rp': np.float64, 'radial_velocity': np.float64, 'radial_velocity_error': np.float64, 'phot_variable_flag': np.str_, 'teff_gspphot': np.float64, 'distance_gspphot': np.float64, 'mag': np.float64, 'ra_deg': np.float64, 'dec_deg': np.float64, 'separation': np.float64, 'measured_mag': np.float64, 'mag_difference': np.float64}
measurements_table_dtype = (np.str_, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.str_, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, np.int32)


measurements_table = QTable(dtype=measurements_table_dtype, names=('designation', 'ra', 'dec', 'parallax', 'parallax_error', 'pm', 'pmra', 'pmdec', 'ruwe', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', 'phot_variable_flag', 'teff_gspphot', 'distance_gspphot', 'mag', 'ra_deg', 'dec_deg', 'separation', 'measured_mag', 'measured_mag_err', 'mag_difference', 'mag_difference_err', 'number_of_measurements'), meta={'name': 'object table'})

for fitsFile in files:
    # 1. Read the list of sources extracted from an image (fits) file
    file_counter = file_counter + 1
    print('\n\n### Processing file', file_counter, 'out of', len(files),': ', fitsFile, '###')
    sources = extract_sources(workingDirectory + '/' + fitsFile)
    print('Number of sources: ', len(sources[0]))
    master_catalog = SkyCoord(ra=gaia_catalog['ra'], dec=gaia_catalog['dec'], frame='icrs')
    sources_catalog = SkyCoord(ra=sources[0]['ra_deg']*u.degree, dec=sources[0]['dec_deg']*u.degree, frame='fk5')
    idxw, idxs, wsd2d, wsd3d = search_around_sky(master_catalog, sources_catalog, search_cone*u.deg)
    composit_catalog = hstack([gaia_catalog[idxw], sources[0][idxs]])
    
    composit_catalog.add_column(wsd2d.arcsecond, name='separation', index=-1)

    bkg_limiting_mag = calculate_limiting_magnitude(composit_catalog['mag'], composit_catalog['phot_g_mean_mag'])
    composit_catalog['bkg_limiting_mag'] = bkg_limiting_mag
    composit_catalog.add_column(bkg_limiting_mag + composit_catalog['mag'], name='measured_mag')
    composit_catalog.add_column(composit_catalog['phot_g_mean_mag'].data - composit_catalog['measured_mag'], name='mag_difference')

    sources[0].write('sources-' + str(file_counter) + '.csv', format='ascii', overwrite=True, delimiter=',')
    composit_catalog.write('composit_catalog-' + str(file_counter) + '.csv', format='ascii', overwrite=True, delimiter=',')
    plot_sources(sources[2], sources[1], sources[0], composit_catalog, 'source_plot_' + str(file_counter))
    measured_objects = vstack([measured_objects, composit_catalog])

grouped_objects = measured_objects.group_by(['DESIGNATION'])

# Write code to aggregate values of each stars

count = 1
for measured_object in grouped_objects.groups:
    measurements_table.add_row([
        measured_object[0]['DESIGNATION'],
        measured_object[0]['ra'],
        measured_object[0]['dec'],
        measured_object[0]['parallax'],
        measured_object[0]['parallax_error'],
        measured_object[0]['pm'],
        measured_object[0]['pmra'],
        measured_object[0]['pmdec'],
        measured_object[0]['ruwe'],
        measured_object[0]['phot_g_mean_mag'],
        measured_object[0]['phot_bp_mean_mag'],
        measured_object[0]['phot_rp_mean_mag'],
        measured_object[0]['bp_rp'],
        measured_object[0]['bp_g'],
        measured_object[0]['g_rp'],
        measured_object[0]['radial_velocity'],
        measured_object[0]['radial_velocity_error'],
        measured_object[0]['phot_variable_flag'],
        measured_object[0]['teff_gspphot'],
        measured_object[0]['distance_gspphot'],
        float(measured_object['mag'].mean()),
        float(measured_object['ra_deg'].mean()),
        float(measured_object['dec_deg'].mean()),
        float(measured_object['separation'].mean()),
        float(measured_object['measured_mag'].mean()),
        float(measured_object['measured_mag'].std()),
        float(measured_object['mag_difference'].mean()),
        float(measured_object['mag_difference'].std()),
        int(len(measured_object['mag']))
    ])

    count = count + 1

sum_plot = extract_sources(workingDirectory + '/' + files[0])
plot_sum(sum_plot[2], sum_plot[1], measurements_table, 'sum_plot.jpg')

measurements_table.write('measurements_table' + '.csv', format='ascii', overwrite=True, delimiter=',')
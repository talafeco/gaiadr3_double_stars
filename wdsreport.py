#!/usr/bin/python3
# WDS Report tool to measure double stars on astronomical images based on Gaia DR3 data
# Version: 1.0
# Usage: wdsmeasure -d [IMAGE FOLDER] -[OPTIONS]
# -d <image_folder>             Define the image folder path
# -o    --orbit_calculations    Calculate the historical orbit of the double star
# -g    --gaia_measurements     Calculate the double star attributes based on the Gaia DR3 data release
# -O    --offline               Use the offline DR3 star data instead of the online TAP service query to calculate the double star attributes

# Tasks
# Calculate the double star's precise coordinates for current epoch
# Get all double stars from wds, which should be found on the images
# Search double stars based on magnitude difference, separation and position angle, if not found based on the coordinates

''' Errors:  fizikai jellemzők egyeznek az Excellel.
szeparáció átszámítás: wdsreport 37 171 AU - OK
Excel 38 708 AU ( képlet: WTD Dist*Gaia Sep)
szökési sebesség:            wdsreport  0,561 - OK
Excel   0,55 (ezzel együtt lehet élni!)
sajátmozgás (rPM)           wdsreport   CPM
xcel    SPM (0,31) ez azért gond, mert nem CPM!
Harshaw faktor                   wdsreport 0,352
Excel  0,4457 végülis ez is az is „??” eredményez
Historic velocity                 wdsreport max orbit 0,5571
Excel max orbit 0,43
wdsreport obs orbit 23,6803
Excel obs orbit 3,1109 optikai így is, úgy is, de az eltérés nagy, lehet a képlet hibás? '''

# 2024-07-27: updated only data[0] is counted due to rgb images
# check color fits file, then slice and measure independently: https://astronomy.stackexchange.com/questions/32222/converting-an-rgb-image-to-fits-astropy
# Check, which double stars should be on the image?
# it there should be a double star, but not found, search for similar separation + similar mag difference doubles
# design the measurement and calculation functions independently, process the images and write measurements, if there is no response from gaia in 10s or the star not found
# Upodate gaia search, to get a minimal distance in arcsecs

import pprint
import argparse
import os
import sys
import numpy as np
import numpy.ma as ma
import math
import sys
import datetime
import configparser
from time import gmtime, strftime
from dsmodules import dscalculation
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.coordinates import match_coordinates_sky, search_around_sky, position_angle, angular_separation
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from astropy.table import Column, MaskedColumn
from astropy.utils.masked import MaskedNDArray
from photutils.detection import DAOStarFinder
from astropy.wcs import WCS
import astropy.units as u
import warnings
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta
warnings.filterwarnings("ignore")
from astroquery.gaia import Gaia
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="A program with command line options.")

# Add optional argument --gaia-measurements
parser.add_argument(
    '-G', '--gaia_measurements',
    action='store_true',  # This makes it a flag (True when present, False when absent)
    help="Do double star gaia calculations"
)

# Add optional argument --orbit_calculations
parser.add_argument(
    '-O', '--offline',
    action='store_true',  # This makes it a flag (True when present, False when absent)
    help="Do double star historic orbit calculations"
)

# Add optional argument --orbit_calculations
parser.add_argument(
    '-o', '--orbit_calculations',
    action='store_true',  # This makes it a flag (True when present, False when absent)
    help="Do double star historic orbit calculations"
)

# Add another example option (e.g., input file)
parser.add_argument(
    '-d', '--directory',
    type=str,
    help="Path of the directory containing the fits files"
)

parser.add_argument(
    '-f', '--find_doubles_on_image',
    action='store_true',  # This makes it a flag (True when present, False when absent)
    help="Search for double stars in the WDS database, which should be on the image no matter, if they are too faint, close, etc."
)

parser.add_argument(
    '-y', '--generate_wds_list',
    action='store_true',  # This makes it a flag (True when present, False when absent)
    help="Write the list of double stars in the WDS database, which should be on the image into a file."
)

args = parser.parse_args()

# Set working directory to read Double Star images
# workingDirectory = sys.argv[1]
workingDirectory = args.directory

# Read configuration
# Create a ConfigParser object
config = configparser.ConfigParser()

# Define the configuration file
config_file = 'C:\\Users\\gerge\\Documents\\Github\\gaiadr3_double_stars\\config.ini'

# Read the configuration file
config.read(config_file)

# Configuration of the camera
dao_sigma = float(config['source_detection']['dao_sigma'])
dao_fwhm = float(config['source_detection']['dao_fwhm'])
dao_threshold = float(config['source_detection']['dao_threshold'])

# Configurations for calculations
possible_distance = float(config['calculations']['possible_distance'])
search_cone = float(config['calculations']['search_cone'])
dummyObservationDate = config['calculations']['dummyObservationDate']
gaia_dr3_epoch = float(config['calculations']['gaia_dr3_epoch'])
auToParsec = float(config['calculations']['auToParsec'])

# Gravitational constant is convenient if measure distances in parsecs (pc), velocities in kilometres per second (km/s) and masses in solar units M
gravConst = float(config['calculations']['gravConst'] )
image_limit = float(config['imageplot']['image_limit'])

# Constant to calculate star luminosity and mass
sun_luminosity = 3.0128 * (10 ** 28)
sun_absolute_luminosity = 3.828 * (10 ** 26)

# Constant variables
hipparcos_file = Table.read(config['data']['hipparcos_file'], format='ascii')

# Insert the downloaded wds file path here
wds_file = config['data']['wds_file']

# Library of map segments
segment_lib = config['data']['segment_lib']

# Create WDS table
wds_converters = {  '2000 Coord': np.str_,
                    'Discov': np.str_,
                    'Comp': np.str_,
                    'Date (first)': np.str_,
                    'Date (last)': np.str_,
                    'Obs': np.str_,
                    'PA_f': np.float64,
                    'PA_l': np.float64,
                    'Sep_f': np.float64,
                    'Sep_l': np.float64,
                    'Mag_A': np.str_,
                    'Mag_B': np.str_,
                    'Spectral A/B': np.str_,
                    'PM_A_ra': np.str_,
                    'PM_A_dec': np.str_,
                    'PM_B_ra': np.str_,
                    'PM_B_dec': np.str_,
                    'D. Number': np.str_,
                    'Notes': np.str_,
                    'Coord (RA)': np.str_,
                    'Coord (DEC)': np.str_
                    }
wds_data = Table.read(wds_file,
                      names=('2000 Coord', 'Discov', 'Comp', 'Date (first)', 'Date (last)', 'Obs',
                             'PA_f', 'PA_l', 'Sep_f', 'Sep_l', 'Mag_A',
                             'Mag_B', 'Spectral A/B', 'PM_A_ra', 'PM_A_dec',
                             'PM_B_ra', 'PM_B_dec', 'D. Number', 'Notes', 'Coord (RA)',
                             'Coord (DEC)'
                        ),
                      converters=wds_converters,
                      format='ascii.fixed_width',
                      header_start=2, data_start=5,
                      col_starts=(0, 10, 17, 23, 28, 33, 38, 42, 46, 52, 58, 64, 70, 80, 84, 89, 93, 98, 107, 112, 121),
                      col_ends=(9, 16, 21, 26, 31, 36, 40, 44, 50, 56, 61, 68, 78, 83, 87, 92, 96, 105, 110, 120, 129),
                      )

###################################################################################################################################

print('\n### Reading WDS database ###')

dsfilename = np.array([], dtype=str)
dswds_identifier = np.array([], dtype=str)
dsdiscovr = np.array([], dtype=str)
dscomp = np.array([], dtype=str)
dstheta = np.array([], dtype=np.float64)
dsrho = np.array([], dtype=np.float64)
dsmag_pri = np.array([], dtype=np.float64)
dsmag_sec = np.array([], dtype=np.float64)
dsmag_diff = np.array([], dtype=np.float32)
dsspectra = np.array([], dtype=str)
dspm_a_ra = np.array([], dtype=np.float64)
dspm_a_dec = np.array([], dtype=np.float64)
dspm_b_ra = np.array([], dtype=np.float64)
dspm_b_dec = np.array([], dtype=np.float64)
dsra = np.array([], dtype=str)
dsdec = np.array([], dtype=str)
dsdegra = np.array([], dtype=np.float64)
dsdegdec = np.array([], dtype=np.float64)
dsobjectid = np.array([], dtype=str)
dspaactual = np.array([], dtype=np.float64)
dssepactual = np.array([], dtype=np.float64)
dsmagdiff = np.array([], dtype=np.float64)
dsTable = QTable([dswds_identifier, dsdiscovr, dscomp, dstheta, dsrho, dsmag_pri, dsmag_sec, dsmag_diff, dsspectra, dspm_a_ra, dspm_a_dec, dspm_b_ra, dspm_b_dec, dsra, dsdec, dsdegra, dsdegdec, dsobjectid, dspaactual, dssepactual, dsmagdiff], names=('wds_identifier', 'discovr', 'comp', 'theta', 'rho', 'mag_pri', 'mag_sec', 'dsmag_diff', 'spectra', 'pm_a_ra', 'pm_a_dec', 'pm_b_ra', 'pm_b_dec', 'ra_hms', 'dec_dms', 'ra_deg', 'dec_deg', 'object_id', 'dspaactual', 'dssepactual', 'dsmagdiff'), meta={'name': 'ds table'})

# Define Report table
reportw_identifier = np.array([], dtype=str)
reportdate = np.array([], dtype=str)
reporttheta = np.array([], dtype=np.float64)
reportthetaerr = np.array([], dtype=np.float64)
reportrho = np.array([], dtype=np.float64)
reportrhoerr = np.array([], dtype=np.float64)
reportmag_pri = np.array([], dtype=np.float64)
reportmag_prierr = np.array([], dtype=np.float64)
reportmag_sec = np.array([], dtype=np.float64)
reportmag_secerr = np.array([], dtype=np.float64)
reportfilter =  np.array([], dtype=str)
reportfilterfwhm = np.array([], dtype=str)
reporttelescopeap = np.array([], dtype=str)
reportnights = np.array([], dtype=str)
reportrefcode = np.array([], dtype=str)
reporttech = np.array([], dtype=str)
reportcat = np.array([], dtype=str)
preccoord = np.array([], dtype=str)
reportTable = QTable([reportw_identifier, reportdate, reporttheta, reportthetaerr, reportrho, reportrhoerr, reportmag_pri, reportmag_prierr, reportmag_sec, reportmag_secerr, reportfilter, reportfilterfwhm, reporttelescopeap, reportnights, reportrefcode, reporttech, reportcat, preccoord], names=('wds_identifier', 'date_of_obs', 'mean_theta', 'mean_theta_err', 'mean_rho', 'mean_rho_err', 'mag_pri', 'mag_pri_err', 'mag_sec', 'mag_sec_err', 'filter', 'filter_fwhm', 'telescope_ap', 'nights_of_obs', 'reference_code', 'tech_code', 'catalog_code', 'precise_cordinates'), meta={'name': 'report table'})

print('\n### Creating filelist ###')
directoryContent = os.listdir(workingDirectory)
print('Working directory: ', workingDirectory)

# Add fit, fits file extensions too
files = [f for f in directoryContent if os.path.isfile(workingDirectory+'/'+f) and (f.endswith('.new') or f.endswith('.fit') or f.endswith('.fits'))]
print('Files:', files)

# Create the WDS catalog table
print('Creating wds catalog - ' , datetime.datetime.now())
wdsTable = dscalculation.create_wds_table(wds_data)
print('WDS angle formatting done, creating catalog... - ' , datetime.datetime.now())
wds_catalog = SkyCoord(ra=wdsTable['Coord (RA) hms'], dec=Angle(wdsTable['Coord (DEC) dms']), unit='hour, degree', frame="icrs") #
print('WDS Catalog created successfully - ' , datetime.datetime.now())
sources_ds = Table()

file_counter = 0

### Run source detection, collect star data to Qtable
print('\n### Running source detection ###\n' , datetime.datetime.now())
for fitsFile in files:
    # 1. Read the list of sources extracted from an image (fits) file
    file_counter = file_counter + 1
    print('\n\n### Processing file', file_counter, 'out of', len(files),': ', fitsFile, '###')
    fitsFileName = workingDirectory + '/' + fitsFile
    hdu, image_wcs, fits_header, fits_data, fits_file_date = dscalculation.get_fits_data(fitsFileName)

    if args.find_doubles_on_image:
        image_center, image_radius = dscalculation.calculate_photo_center(image_wcs, fits_header)
        doubles_on_image_table = dscalculation.catalog_search_in_image(image_wcs, fits_header, image_center, image_radius, wds_catalog, wdsTable)
        print(doubles_on_image_table)
        print(doubles_on_image_table.info)
        #wds_doubles = SkyCoord.to_pixel(doubles_on_image_table['Coord (RA) hms'], doubles_on_image_table['Coord (DEC) dms'])
        if args.generate_wds_list:
            doubles_on_image_table.write(workingDirectory + fitsFileName[:-4] + '_wds_pairs.csv', format='ascii', overwrite=True, delimiter=',')

    sources_ds = dscalculation.get_sources_from_image(sources_ds, wds_catalog, fits_data, fits_header, fitsFileName, fits_file_date, image_wcs, wdsTable, dao_sigma, dao_fwhm, dao_threshold, search_cone)
    image_plane = dscalculation.define_image_plane(image_wcs, fits_header)

    #dscalculation.plot_image_with_frame(hdu, wds_doubles, image_plane, 'testing.txt', image_limit, workingDirectory)

upd_sources_ds = sources_ds[sources_ds['rho_measured'] != 0]
upd_sources_ds_by_object = upd_sources_ds.group_by(['2000 Coord', 'Discov', 'Comp'])

#dscalculation.exit_if_no_doubles_found(upd_sources_ds_by_object)

print('### Updated sources DS table grouped by WDS Identifier, Discoverer and Components ###')
print(upd_sources_ds_by_object)
upd_sources_ds_by_object.write(workingDirectory + '/double_stars.csv', format='ascii', overwrite=True, delimiter=',')
objectMean = upd_sources_ds_by_object.groups.aggregate(np.mean)

count = 1
for ds in upd_sources_ds_by_object.groups:
    print('\n#---------------------------------------------------------------------------------------------------------------------#')
    print('\n### Group index:', count, '###')
    count = count + 1
    # Search component in the offline Gaia DR3 database
    print('DS: ', ds)
    firstFitsImageFileName = ds['file'][0]
    wds_double_star = dscalculation.wds_measurement(ds)
    dscalculation.imagePlot(firstFitsImageFileName, workingDirectory, wds_double_star.pairObjectId, wds_double_star.starActualRa1, wds_double_star.starActualDec1, wds_double_star.starActualRa2, wds_double_star.starActualDec2, image_limit)
    dscalculation.write_wds_report(ds, wds_double_star, workingDirectory)
    gaiaAStar, gaiaBStar, searchKey = 0, 0, 0

    if args.gaia_measurements and args.offline:
        gaiaAStar, gaiaBStar = dscalculation.get_gaia_dr3_data_offline(ds, segment_lib)
        searchKey = 'designation'
    elif args.gaia_measurements and args.offline is False:
        gaiaAStar, gaiaBStar = dscalculation.get_gaia_dr3_data(ds, search_cone)
        searchKey = 'DESIGNATION'
        
    if gaiaAStar and gaiaBStar:
        gaia_ds = dscalculation.gaia_calculations(gaiaAStar, gaiaBStar, ds, searchKey)
        gaiaData = gaia_ds.gaiaData
        dscalculation.hrdPlot(wds_double_star.pairObjectId, workingDirectory, gaia_ds.pairAbsMag1, gaia_ds.pairAbsMag2, gaia_ds.pairBVIndexA, gaia_ds.pairBVIndexB, hipparcos_file)
        dscalculation.write_gaia_report(ds, wds_double_star, gaia_ds, workingDirectory)

        if args.orbit_calculations:
            historic_orbit_calculation = dscalculation.calculate_historical_orbit(gaia_ds, wds_double_star, ds)
            dscalculation.write_historic_orbit_report(wds_double_star, historic_orbit_calculation, workingDirectory)

reportTable.write(workingDirectory + '/double_stars_wds_format.txt', format='ascii', overwrite=True, delimiter=',')
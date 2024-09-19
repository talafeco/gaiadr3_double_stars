#!/usr/bin/python3
# WDS Report tool to measure double stars on astronomical images based on Gaia DR3 data
# Version: 1.0


# Tasks
# Calculate the double star's precise coordinates for current epoch
# Get all double stars from wds, which should be found on the images
# Search double stars based on magnitude difference, separation and position angle, if not found based on the coordinates


# 2024-07-27: updated only data[0] is counted due to rgb images
# check color fits file, then slice and measure independently: https://astronomy.stackexchange.com/questions/32222/converting-an-rgb-image-to-fits-astropy
# Check, which double stars should be on the image?
# it there should be a double star, but not found, search for similar separation + similar mag difference doubles
# design the measurement and calculation functions independently, process the images and write measurements, if there is no response from gaia in 10s or the star not found
# Upodate gaia search, to get a minimal distance in arcsecs

# 40-50 ivmásodperccel beljebb mérjen - ok
# Listázza ki fileba a képen szereplő kettősöket - ok
# Limitálni a magnitúdókat, amiket feldob
# plotolni a wds-ben szereplő párokat a képre

import os
import sys
import numpy as np
import sys
import datetime
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.coordinates import match_coordinates_sky, search_around_sky
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from photutils.detection import DAOStarFinder
from astropy.wcs import WCS
import astropy.units as u
import warnings
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta
warnings.filterwarnings("ignore")
from regions import RegularPolygonPixelRegion, RectangleSkyRegion, RectanglePixelRegion
from dsmodules import dscalculation
import configparser

# Read configuration
# Create a ConfigParser object
config = configparser.ConfigParser()

# Read the configuration file
config.read('config.ini')

dao_sigma = float(config['source_detection']['dao_sigma'])
dao_fwhm = float(config['source_detection']['dao_fwhm'])
dao_threshold = float(config['source_detection']['dao_threshold'])

# Configurations for calculations
possible_distance = float(config['calculations']['possible_distance'])
search_cone = float(config['calculations']['search_cone'])
dummyObservationDate = config['calculations']['dummyObservationDate']

# Gravitational constant is convenient if measure distances in parsecs (pc), velocities in kilometres per second (km/s) and masses in solar units M
image_limit = float(config['imageplot']['image_limit'])

# Define frame percent
frame_percentage = float(config['frame']['frame_percentage'])

# Define the limiting magnitue of the components
limiting_magnitude_primary = 13
limiting_magnitude_secondary = 13

# Insert the downloaded wds file path here
wds_file = config['data']['wds_file']

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

# Set working directory to read Double Star images
workingDirectory = sys.argv[1]

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

    sources_ds = dscalculation.get_sources_from_image(sources_ds, wds_catalog, fits_data, fits_header, fitsFileName, fits_file_date, image_wcs, wdsTable, dao_sigma, dao_fwhm, dao_threshold)

upd_sources_ds = sources_ds[sources_ds['rho_measured'] != 0]

print('upd_sources_ds info: ', upd_sources_ds.info)

upd_sources_ds_by_object = upd_sources_ds.group_by(['2000 Coord', 'Discov', 'Comp'])

print('### Updated sources DS table grouped by WDS Identifier, Discoverer and Components ###')

upd_sources_ds_by_object.write(workingDirectory + '/double_stars.csv', format='ascii', overwrite=True, delimiter=',')

dscalculation.exit_if_no_doubles_found(upd_sources_ds_by_object)

objectMean = upd_sources_ds_by_object.groups.aggregate(np.mean)

count = 1
for ds in upd_sources_ds_by_object.groups:
    print('\n#---------------------------------------------------------------------------------------------------------------------#')
    print('\n### Group index:', count, '###')
    # print(ds)
    count = count + 1

    double_star = dscalculation.wds_measurement(ds)

    reportName = (workingDirectory + '/' + double_star.pairObjectId + '.txt').replace(' ', '')
    reportFile = open(reportName, "a")

    firstFitsImageFileName = ds['file'][0]
    dscalculation.imagePlot(firstFitsImageFileName, workingDirectory, double_star.pairObjectId, double_star.starActualRa1, double_star.starActualDec1, double_star.starActualRa2, double_star.starActualDec2)

    # Print temp data
    print('\n### COMPONENTS ###')
    print('WDS Identifier:', ds[0]['2000 Coord'], ds[0]['Discov'], ds[0]['Comp'])
    print('Magnitude Pri: ' + str(ds[0]['Mag_A']))
    print('Magnitude Sec: ' + str(ds[0]['Mag_B']))
    print('PA last: ', str(ds[0]['PA_l']))
    print('Sep last: ',  str(ds[0]['Sep_l']))
    print('Date of observation (human readable): ', Time(ds['image_date'].data).mean())
    print('Date of observation (Julian date): ' + double_star.dateOfObservation)
    print('Precise coordinates (J2000): ' + double_star.preciseCoord)
    print('\nTheta measurements\n') # , ds['dspaactual']
    print('Mean:', double_star.pairMeanTheta)
    print('Error:', double_star.pairMeanThetaErr)
    print('\nRho measurements\n') # , ds['dssepactual']
    print('Mean:', double_star.pairMeanRho)
    print('Error:', double_star.pairMeanRhoErr)
    print('\nMagnitude difference measurements\n') # , ds['dsmagdiff']
    print('Mean:', double_star.pairMagDiff)
    print('Error:', double_star.pairMagDiffErr)
    
    # Write results to file
    reportTable.add_row([ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp']), double_star.dateOfObservation, double_star.pairMeanTheta, double_star.pairMeanThetaErr, double_star.pairMeanRho, double_star.pairMeanRhoErr, np.nan, np.nan, double_star.pairMagDiff, double_star.pairMagDiffErr, 'Filter wawelenght', 'filter FWHM', '0.2', '1', 'TLB_2024', 'C', '7', double_star.preciseCoord])
    reportFile.write('### WDS Data ###')
    reportFile.write('\nWDS Identifier: ' + ds[0]['2000 Coord'])
    reportFile.write('\nDiscoverer and components: ' + str(ds[0]['Discov']) + ' ' + str(ds[0]['Comp']))
    reportFile.write('\nMagnitude Pri: ' + str(ds[0]['Mag_A']))
    reportFile.write('\nMagnitude Sec: ' + str(ds[0]['Mag_B']))
    reportFile.write('\nPA last: ' + str(ds[0]['PA_l']))
    reportFile.write('\nSep last: ' +  str(ds[0]['Sep_l']))
    reportFile.write('\n\n### Measurements ###')
    reportFile.write('\nDate of observation (human readable): ' + str(Time(ds['image_date'].data).mean()))
    reportFile.write('\nDate of observation: ' + double_star.dateOfObservation)
    reportFile.write('\nPrecise coordinates (J2000): ' + double_star.preciseCoord)
    reportFile.write('\n\nPosition angle:')
    reportFile.write('\nTheta measurements' + str(ds['theta_measured'].degree))
    reportFile.write('\nMean: ' + str(dscalculation.roundNumber(double_star.pairMeanTheta)))
    reportFile.write('\nError: ' + str(dscalculation.roundNumber(double_star.pairMeanThetaErr)))
    reportFile.write('\n\nSeparation:')
    reportFile.write('\nRho measurements\n' + str(ds['rho_measured'].arcsec))
    reportFile.write('\nMean: ' + str(dscalculation.roundNumber(double_star.pairMeanRho)))
    reportFile.write('\nError: ' + str(dscalculation.roundNumber(double_star.pairMeanRhoErr)))
    reportFile.write('\n\nMagnitude measurements\n'  + str(ds['mag_diff']))
    reportFile.write('\nMean: ' + str(dscalculation.roundNumber(double_star.pairMagDiff)))
    reportFile.write('\nError: ' + str(dscalculation.roundNumber(double_star.pairMagDiffErr)))

    
    reportFile.write('\n\n### WDS form:\n')
    wdsform = str(ds[0]['2000 Coord']) + ',' + double_star.dateOfObservation + ',' +  str(dscalculation.roundNumber(double_star.pairMeanTheta)) + ',' +  str(dscalculation.roundNumber(double_star.pairMeanThetaErr)) + ',' +  str(dscalculation.roundNumber(double_star.pairMeanRho)) + ',' +  str(dscalculation.roundNumber(double_star.pairMeanRhoErr)) + ',' +  'nan' + ',' +  'nan' + ',' +  str(dscalculation.roundNumber(double_star.pairMagDiff)) + ',' +  str(dscalculation.roundNumber(double_star.pairMagDiffErr)) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TLB_2023' + ',' +  'C' + ',' + '7'+ ',' + str(dscalculation.getPreciseCoord(double_star.starActualRa1, double_star.starActualDec1, fits_file_date))
    reportFile.write(str(wdsform))

reportTable.write(workingDirectory + '/double_stars_wds_format.txt', format='ascii', overwrite=True, delimiter=',')
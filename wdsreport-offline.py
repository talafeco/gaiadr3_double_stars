#!/usr/bin/python3
# WDS Report tool to measure double stars on astronomical images based on Gaia DR3 data
# Version: 1.0
# Usage: wdsreport <image_folder>

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

# Set working directory to read Double Star images
workingDirectory = sys.argv[1]

# Read configuration
# Create a ConfigParser object
config = configparser.ConfigParser()

# Read the configuration file
config.read('config.ini')

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

    sources_ds = dscalculation.get_sources_from_image(sources_ds, wds_catalog, fits_data, fits_header, fitsFileName, fits_file_date, image_wcs, wdsTable, dao_sigma, dao_fwhm, dao_threshold)

upd_sources_ds = sources_ds[sources_ds['rho_measured'] != 0]
upd_sources_ds_by_object = upd_sources_ds.group_by(['2000 Coord', 'Discov', 'Comp'])

print(upd_sources_ds_by_object.info)

print('### Updated sources DS table grouped by WDS Identifier, Discoverer and Components ###')
print(upd_sources_ds_by_object)

upd_sources_ds_by_object.write(workingDirectory + '/double_stars.csv', format='ascii', overwrite=True, delimiter=',')


objectMean = upd_sources_ds_by_object.groups.aggregate(np.mean)

count = 1
for ds in upd_sources_ds_by_object.groups:
    print('\n#---------------------------------------------------------------------------------------------------------------------#')
    print('\n### Group index:', count, '###')
    # print(ds)
    count = count + 1
    #pairObjectId = ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp'])

    # Search component in the offline Gaia DR3 database
    print('DS: ', ds)
    gaiaAStar, gaiaBStar = dscalculation.get_gaia_dr3_data_offline(ds, segment_lib)

    print('GaiaStar: ', gaiaAStar['designation'], gaiaBStar['designation'])

    if gaiaAStar and gaiaBStar:
        print('Gaia A star: ' + str(gaiaAStar['designation']))
        print('Gaia B star: ' + str(gaiaBStar['designation']))

        wds_double_star = dscalculation.wds_measurement(ds)
        gaia_ds = dscalculation.gaia_calculations(gaiaAStar, gaiaBStar, ds)
        historic_orbit_calculation = dscalculation.calculate_historical_orbit(gaia_ds, wds_double_star, ds)

        print('### Historic Orbit Calculation ###')
        print(pprint.pprint(vars(historic_orbit_calculation)))

        reportName = (workingDirectory + '/' + wds_double_star.pairObjectId + '.txt').replace(' ', '')
        reportFile = open(reportName, "a")

        gaiaData = gaia_ds.gaiaData
        
        dscalculation.hrdPlot(wds_double_star.pairObjectId, workingDirectory, gaia_ds.pairAbsMag1, gaia_ds.pairAbsMag2, gaia_ds.pairBVIndexA, gaia_ds.pairBVIndexB)
        
        print('firstFitsImageFileName =', ds['file'][0])
        firstFitsImageFileName = ds['file'][0]
        print(str(firstFitsImageFileName), str(wds_double_star.pairObjectId), str(gaia_ds.pairRaA), str(gaia_ds.pairDecA), str(gaia_ds.pairRaB), str(gaia_ds.pairDecB))
        dscalculation.imagePlot(firstFitsImageFileName, workingDirectory, wds_double_star.pairObjectId, gaia_ds.pairRaA, gaia_ds.pairDecA, gaia_ds.pairRaB, gaia_ds.pairDecB)
    
        # Print temp data
        print('\n### COMPONENTS ###')
        print('\nWDS Identifier:', ds[0]['2000 Coord'], ds[0]['Discov'], ds[0]['Comp'])
        print('\nDate of observation: ' + wds_double_star.dateOfObservation)
        print('\nPrecise coordinates (J2000): ' + wds_double_star.preciseCoord)
        print('\nComponent A:', gaia_ds.pairDesignationA)
        print('Component B:', gaia_ds.pairDesignationB)
        print('\nCalculated coordinates')
        print('\nComponent A DR3 2016:', gaia_ds.pairRaA, gaia_ds.pairDecA)
        print('Component A DR3 on date:', gaia_ds.pairACurrentCoord.ra.degree, gaia_ds.pairACurrentCoord.dec.degree)
        print('Component A measured:', wds_double_star.pairAMeasuredCoord.ra.degree, wds_double_star.pairAMeasuredCoord.dec.degree)
        print('Component A error (on date - measured):', gaia_ds.pairACoordErr.arcsecond)
        print('\nComponent B DR3 2016:', gaia_ds.pairRaB, gaia_ds.pairDecB)
        print('Component B DR3 on date:', gaia_ds.pairBCurrentCoord.ra.degree, gaia_ds.pairBCurrentCoord.dec.degree)
        print('Component B measured:', wds_double_star.pairBMeasuredCoord.ra.degree, wds_double_star.pairBMeasuredCoord.dec.degree)
        print('Component B error (on date - measured):', gaia_ds.pairBCoordErr.arcsecond)
        print('2016 Calculated Position angle / Separation: ', SkyCoord(ra=gaia_ds.pairRaA*u.degree, dec=gaia_ds.pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=gaia_ds.pairRaB*u.degree, dec=gaia_ds.pairDecB*u.degree, frame='icrs')).degree, SkyCoord(ra=gaia_ds.pairRaA*u.degree, dec=gaia_ds.pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=gaia_ds.pairRaB*u.degree, dec=gaia_ds.pairDecB*u.degree, frame='icrs')).arcsecond)
        #print('Current Calculated Position angle / Separation: ', SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).degree, SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).arcsecond)
        print('Current Calculated Position angle / Separation: ', gaia_ds.pairACurrentCoord.position_angle(gaia_ds.pairBCurrentCoord).degree, gaia_ds.pairACurrentCoord.separation(gaia_ds.pairBCurrentCoord).arcsecond)
        print('\nTheta measurements\n') # , ds['dspaactual']
        print('Mean:', wds_double_star.pairMeanTheta)
        print('Error:', wds_double_star.pairMeanThetaErr)
        print('\nRho measurements\n') # , ds['dssepactual']
        print('Mean:', wds_double_star.pairMeanRho)
        print('Error:', wds_double_star.pairMeanRhoErr)
        print('\nMagnitude measurements\n') # , ds['dsmagdiff']
        print('Mean:', wds_double_star.pairMagDiff)
        print('Error:', wds_double_star.pairMagDiffErr)
        print('\n\nParallax factor:', gaia_ds.pairParallaxFactor, '%')
        print('Proper motion factor:', gaia_ds.pairPmFactor * 100, '%')
        print('Proper motion category:', gaia_ds.pairPmCommon)
        print('Absolute magnitude A:', gaia_ds.pairAbsMag1)
        print('Absolute magnitude B:', gaia_ds.pairAbsMag2)
        print('Luminosity A:', gaia_ds.pairLum1)
        print('Luminosity B:', gaia_ds.pairLum2)
        print('Luminosity Alternate A:', gaia_ds.pairAltLum1)
        print('Luminosity Alternate B:', gaia_ds.pairAltLum2)
        print('Mass A:', gaia_ds.pairMass1)
        print('Mass B:', gaia_ds.pairMass2)
        print('BV index A:', gaia_ds.pairBVIndexA, 'B:', gaia_ds.pairBVIndexB)
        print('Radial velocity of the stars', 'A:', gaia_ds.pairRadVelA, 'km/s (Err:', gaia_ds.pairRadVelErrA, 'km/s)', 'B:', gaia_ds.pairRadVelB, 'km/s (Err:', gaia_ds.pairRadVelErrB, 'km/s)')
        print('Radial velocity ratio A:', gaia_ds.pairRadVelRatioA, '%')
        print('Radial velocity ratio B:', gaia_ds.pairRadVelRatioB, '%')
        print('Separation:', gaia_ds.pairSepPar, 'parsec,', gaia_ds.pairSepPar * 206265, 'AU')
        print('Pair Escape velocity:', gaia_ds.pairEscapeVelocity, 'km/s')
        print('Pair Relative velocity:', gaia_ds.pairRelativeVelocity, 'km/s')
        '''print('### Pair historical orbit calculations ###')
        print('Historic criterion: ', gaia_ds.pair_orbit[0])
        print('Max orbit velolicy: ', gaia_ds.pair_orbit[1])
        print('Observed velocity: ', gaia_ds.pair_orbit[2])'''
        print('Pair Harshaw factor:', gaia_ds.pairHarshawFactor)
        print('Pair Harshaw physicality:', gaia_ds.pairHarshawPhysicality)
        print('Pair binarity:', gaia_ds.pairBinarity)
        print('Analysis finished: ' , datetime.datetime.now())
        
        # Write results to file
        reportTable.add_row([ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp']), wds_double_star.dateOfObservation, wds_double_star.pairMeanTheta, wds_double_star.pairMeanThetaErr, wds_double_star.pairMeanRho, wds_double_star.pairMeanRhoErr, np.nan, np.nan, wds_double_star.pairMagDiff, wds_double_star.pairMagDiffErr, 'Filter wawelenght', 'filter FWHM', '0.2', '1', 'TLB_2024', 'C', '7', wds_double_star.preciseCoord])
        reportFile.write('### WDS Data ###')
        reportFile.write('\nWDS Identifier: ' + ds[0]['2000 Coord'])
        reportFile.write('\nDiscoverer and components: ' + str(ds[0]['Discov']) + ' ' + str(ds[0]['Comp']))
        reportFile.write('\nMagnitude Pri: ' + str(ds[0]['Mag_A']))
        reportFile.write('\nMagnitude Sec: ' + str(ds[0]['Mag_B']))
        reportFile.write('\nPA last: ' + str(ds[0]['PA_l']))
        reportFile.write('\nSep last: ' +  str(ds[0]['Sep_l']))
        reportFile.write('\n\n### Gaia DR3 Data ###')
        reportFile.write('\nMain star: ' + gaia_ds.pairDesignationA)
        reportFile.write('\nCompanion: ' + gaia_ds.pairDesignationB)
        reportFile.write('\nPair G magnitudes A: ' + str(dscalculation.roundNumber(gaia_ds.pairMagA)) + ' B: ' + str(dscalculation.roundNumber(gaia_ds.pairMagB)))
        reportFile.write('\nPosition angle: ' + str(dscalculation.roundNumber(gaia_ds.pairDR3Theta)))
        reportFile.write('\nSeparation: ' + str(dscalculation.roundNumber(gaia_ds.pairDR3Rho)))
        reportFile.write('\nMagnitude difference: ' + str(dscalculation.roundNumber(gaia_ds.pairGMagDiff)))
        reportFile.write('\nPrecise coordinates (J2000): ' + wds_double_star.preciseCoord)
        reportFile.write('\nDate of observation: ' + wds_double_star.dateOfObservation)
        reportFile.write('\n\nCalculated coordinates')
        reportFile.write('\nComponent A DR3 2016: ' + str(gaia_ds.pairRaA) + ' ' + str(gaia_ds.pairDecA))
        reportFile.write('\nComponent A DR3 on date: ' + str(gaia_ds.pairACurrentCoord.ra.degree) + ' ' + str(gaia_ds.pairACurrentCoord.dec.degree))
        reportFile.write('\nComponent A measured: ' + str(wds_double_star.pairAMeasuredCoord.ra.degree) + ' ' + str(wds_double_star.pairAMeasuredCoord.dec.degree))
        reportFile.write('\nComponent A error (on date - measured): ' + str(gaia_ds.pairACoordErr.arcsecond))
        reportFile.write('\nComponent B DR3 2016: ' + str(gaia_ds.pairRaB) + ' ' + str(gaia_ds.pairDecB))
        reportFile.write('\nComponent B DR3 on date: ' + str(gaia_ds.pairBCurrentCoord.ra.degree) + ' ' + str(gaia_ds.pairBCurrentCoord.dec.degree))
        reportFile.write('\nComponent B measured: ' + str(wds_double_star.pairBMeasuredCoord.ra.degree) + ' ' + str(wds_double_star.pairBMeasuredCoord.dec.degree))
        reportFile.write('\nComponent B error (on date - measured): ' + str(gaia_ds.pairBCoordErr.arcsecond))
        reportFile.write('\n\n2016 Calculated Position angle / Separation: '  + str(SkyCoord(ra=gaia_ds.pairRaA*u.degree, dec=gaia_ds.pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=gaia_ds.pairRaB*u.degree, dec=gaia_ds.pairDecB*u.degree, frame='icrs')).degree) + ' ' + str(SkyCoord(ra=gaia_ds.pairRaA*u.degree, dec=gaia_ds.pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=gaia_ds.pairRaB*u.degree, dec=gaia_ds.pairDecB*u.degree, frame='icrs')).arcsecond))
        reportFile.write('\nCurrent Calculated Position angle / Separation: ' + str(gaia_ds.pairACurrentCoord.position_angle(gaia_ds.pairBCurrentCoord).degree) + ' ' + str(gaia_ds.pairACurrentCoord.separation(gaia_ds.pairBCurrentCoord).arcsecond))
        reportFile.write('\n\n### Measurements ###')
        reportFile.write('\nPosition angle:')
        reportFile.write('\nTheta measurements' + str(ds['theta_measured'].degree))
        reportFile.write('\nMean: ' + str(dscalculation.roundNumber(wds_double_star.pairMeanTheta)))
        reportFile.write('\nError: ' + str(dscalculation.roundNumber(wds_double_star.pairMeanThetaErr)))
        reportFile.write('\n\nSeparation:')
        reportFile.write('\nRho measurements\n' + str(ds['rho_measured'].arcsec))
        reportFile.write('\nMean: ' + str(dscalculation.roundNumber(wds_double_star.pairMeanRho)))
        reportFile.write('\nError: ' + str(dscalculation.roundNumber(wds_double_star.pairMeanRhoErr)))
        reportFile.write('\n\nMagnitude measurements\n'  + str(ds['mag_diff']))
        reportFile.write('\nMean: ' + str(dscalculation.roundNumber(wds_double_star.pairMagDiff)))
        reportFile.write('\nError: ' + str(dscalculation.roundNumber(wds_double_star.pairMagDiffErr)))
        reportFile.write('\n\n### Calculated attributes ###')
        reportFile.write('\nSeparation (Measured): ' + str(dscalculation.roundNumber(wds_double_star.pairMeanRho)))
        reportFile.write('\nPosition angle (Measured): ' + str(dscalculation.roundNumber(wds_double_star.pairMeanTheta)))
        reportFile.write('\nMagnitude difference (Measured): ' + str(dscalculation.roundNumber(wds_double_star.pairMagDiff)) + ' (Err: ' + str(dscalculation.roundNumber(wds_double_star.pairMagDiffErr)) + ')')
        reportFile.write('\nAbsolute magnitude A: ' + str(dscalculation.roundNumber(gaia_ds.pairAbsMag1)))
        reportFile.write('\nAbsolute magnitude B: ' + str(dscalculation.roundNumber(gaia_ds.pairAbsMag2)))
        reportFile.write('\nLuminosity A: ' + str(dscalculation.roundNumber(gaia_ds.pairLum1)))
        reportFile.write('\nLuminosity B: ' + str(dscalculation.roundNumber(gaia_ds.pairLum2)))
        reportFile.write('\nRad A: ' + str(dscalculation.roundNumber(gaia_ds.pairRad1)))
        reportFile.write('\nRad B: ' + str(dscalculation.roundNumber(gaia_ds.pairRad2)))
        reportFile.write('\nMass A: ' + str(dscalculation.roundNumber(gaia_ds.pairMass1)))
        reportFile.write('\nMass B: ' + str(dscalculation.roundNumber(gaia_ds.pairMass2)))
        reportFile.write('\nBV index A: ' + str(dscalculation.roundNumber(gaia_ds.pairBVIndexA)) + ' B: ' + str(dscalculation.roundNumber(gaia_ds.pairBVIndexB)))
        reportFile.write('\nRadial velocity of the stars ' + 'A:' + str(dscalculation.roundNumber(gaia_ds.pairRadVelA)) + 'km/s (Err:' + str(dscalculation.roundNumber(gaia_ds.pairRadVelErrA)) + 'km/s)' + ' B:' + str(dscalculation.roundNumber(gaia_ds.pairRadVelB)) + 'km/s (Err:' + str(dscalculation.roundNumber(gaia_ds.pairRadVelErrB)) + 'km/s)')
        reportFile.write('\nRadial velocity ratio A: ' + str(dscalculation.roundNumber(gaia_ds.pairRadVelRatioA)) + ' %')
        reportFile.write('\nRadial velocity ratio B: ' + str(dscalculation.roundNumber(gaia_ds.pairRadVelRatioB)) + ' %')
        reportFile.write('\nSeparation: ' + str(dscalculation.roundNumber(gaia_ds.pairDistance[2] * auToParsec)) + ' parsec, ' + str(dscalculation.roundNumber((gaia_ds.pairDistance[2]))) + ' AU')
        reportFile.write('\nPair Escape velocity: ' + str(dscalculation.roundNumber(gaia_ds.pairEscapeVelocity)) + ' km/s')
        reportFile.write('\nPair Relative velocity: ' + str(dscalculation.roundNumber(gaia_ds.pairRelativeVelocity)) + ' km/s')
        reportFile.write('\n\n### Analysis ###')
        reportFile.write('\nParallax factor: ' + str(dscalculation.roundNumber(gaia_ds.pairParallaxFactor)) + ' %')
        reportFile.write('\nProper motion factor: ' + str(dscalculation.roundNumber(gaia_ds.pairPmFactor) * 100) + ' %')
        reportFile.write('\nProper motion category: '+ str(gaia_ds.pairPmCommon))
        reportFile.write('\nPair Harshaw factor: ' + str(dscalculation.roundNumber(gaia_ds.pairHarshawFactor)))
        reportFile.write('\nPair Harshaw physicality: ' + str(gaia_ds.pairHarshawPhysicality))
        reportFile.write('\nPair binarity: ' + str(gaia_ds.pairBinarity))
        
        # new function - orbit calculation
        '''reportFile.write('\n\n### Pair historical orbit calculations ###')
        reportFile.write('\nHistoric criterion: ' + str(gaia_ds.pair_orbit[0]))
        reportFile.write('\nDelta theta: ' + str(gaia_ds.pair_orbit[4]))
        reportFile.write('\nDelta rho: ' + str(gaia_ds.pair_orbit[5]))
        reportFile.write('\nDelta time: ' + str(gaia_ds.pair_orbit[6]))
        reportFile.write('\nMax orbit velolicy: ' + str(gaia_ds.pair_orbit[1]))
        reportFile.write('\nObserved velocity: ' + str(gaia_ds.pair_orbit[2]))
        reportFile.write('\nInput data variables: ' + str(gaia_ds.pair_orbit[3]))'''
        
        reportFile.write('\n\n### WDS form:\n')
        wdsform = str(ds[0]['2000 Coord']) + ',' + wds_double_star.dateOfObservation + ',' +  str(dscalculation.roundNumber(wds_double_star.pairMeanTheta)) + ',' +  str(dscalculation.roundNumber(wds_double_star.pairMeanThetaErr)) + ',' +  str(dscalculation.roundNumber(wds_double_star.pairMeanRho)) + ',' +  str(dscalculation.roundNumber(wds_double_star.pairMeanRhoErr)) + ',' +  'nan' + ',' +  'nan' + ',' +  str(dscalculation.roundNumber(wds_double_star.pairMagDiff)) + ',' +  str(dscalculation.roundNumber(wds_double_star.pairMagDiffErr)) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TLB_2023' + ',' +  'C' + ',' + '7'+ ',' + str(dscalculation.getPreciseCoord(gaia_ds.pairRaA, gaia_ds.pairDecA, fits_file_date))
        reportFile.write(str(wdsform))
        reportFile.write('\n\n### Gaia data:\n')
        reportFile.write(str(gaiaData))

reportTable.write(workingDirectory + '/double_stars_wds_format.txt', format='ascii', overwrite=True, delimiter=',')
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

wdsTable = hstack([wds_data, dscalculation.calculate_wds_ra_hourangle(wds_data['Coord (RA)'])])
wdsTable.rename_column('col0', 'Coord (RA) hms')
wdsTable = hstack([wdsTable, dscalculation.calculate_wds_dec_hourangle(wds_data['Coord (DEC)'])])
wdsTable.rename_column('col0', 'Coord (DEC) dms')
wdsTable = hstack([wdsTable, dscalculation.create_unique_id(wds_data['2000 Coord'], wds_data['Discov'])])
wdsTable.rename_column('col0', 'Unique ID')
wdsTable = dscalculation.delete_invalid_lines_wds(wdsTable)

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
    hdu = fits.open(fitsFileName)
    mywcs = WCS(hdu[0].header, naxis=2)
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
    # data[0]
    data = hdu[0].data
    mean, median, std = sigma_clipped_stats(data, sigma = dao_sigma)  

    daofind = DAOStarFinder(fwhm=dao_fwhm, threshold=dao_threshold*std)  
    sources = daofind(data - median)
    ra2, dec2 = mywcs.all_pix2world(sources['xcentroid'], sources['ycentroid'], 1)
    sources.add_column(ra2, name='ra_deg') 
    sources.add_column(dec2, name='dec_deg')
    sources.add_column(fitsFileDate, name='image_date')
    sources.add_column(fitsFileName, name='file')
    photo_center, photo_radius = dscalculation.calculate_photo_center(mywcs, file_header)
    sources_catalog = SkyCoord(ra=sources['ra_deg']*u.degree, dec=sources['dec_deg']*u.degree, frame='fk5')
    idxw, idxs, wsd2d, wsd3d = search_around_sky(wds_catalog, sources_catalog, search_cone*u.deg)
    composit_catalog = hstack([wdsTable[idxw]['2000 Coord', 'Discov', 'Comp', 'Date (first)', 'PA_f', 'PA_l', 'Sep_f', 'Sep_l', 'Mag_A', 'Mag_B'], sources[idxs]['id', 'mag', 'ra_deg', 'dec_deg']])
    companion_catalog = SkyCoord(ra=composit_catalog['ra_deg'] * u.degree, dec=composit_catalog['dec_deg'] * u.degree).directional_offset_by(composit_catalog['PA_l']*u.degree, composit_catalog['Sep_l']*u.arcsec)
    idxs2, d2ds2, d3ds2 = match_coordinates_sky(companion_catalog, sources_catalog)
    composit_catalog2 = hstack([composit_catalog, sources[idxs2]]) #['id', 'mag', 'ra_deg', 'dec_deg']

    sources_pa = SkyCoord(ra=composit_catalog2['ra_deg_1']*u.degree, dec=composit_catalog2['dec_deg_1']*u.degree).position_angle(SkyCoord(ra=composit_catalog2['ra_deg_2']*u.degree, dec=composit_catalog2['dec_deg_2']*u.degree)).to(u.deg)
    sources_sep = SkyCoord(ra=composit_catalog2['ra_deg_1']*u.degree, dec=composit_catalog2['dec_deg_1']*u.degree).separation(SkyCoord(ra=composit_catalog2['ra_deg_2']*u.degree, dec=composit_catalog2['dec_deg_2']*u.degree)).to(u.arcsec)
    sources_mag_diff = composit_catalog2['mag_2'] - composit_catalog2['mag_1']
    
    composit_catalog2.add_column(sources_pa, name='theta_measured')
    composit_catalog2.add_column(sources_sep, name='rho_measured')
    composit_catalog2.add_column(sources_mag_diff, name='mag_diff')
    sources_ds = vstack([sources_ds, composit_catalog2])

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
    pairObjectId = ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp'])

    # Search component in the offline Gaia DR3 database
    print('DS: ', ds)
    gaiaAStar, gaiaBStar = dscalculation.get_gaia_dr3_data_offline(ds, segment_lib)

    print('GiaiaStar: ', gaiaAStar['designation'], gaiaBStar['designation'])

    if gaiaAStar and gaiaBStar:
        print('Gaia A star: ' + str(gaiaAStar['designation']))
        print('Gaia B star: ' + str(gaiaBStar['designation']))
        pairDistanceMinA = dscalculation.calcDistanceMin(float(gaiaAStar['parallax']), float(gaiaAStar['parallax_error']))
        pairDistanceMinB = dscalculation.calcDistanceMin(float(gaiaBStar['parallax']), float(gaiaBStar['parallax_error']))
        
        # Calculate physical binarity

        # Creating empty arrays for Star related calculations
        #Set input data
        starRa1 = float(gaiaAStar['ra'])
        starDec1 = float(gaiaAStar['dec'])
        starRa2 = float(gaiaBStar['ra'])
        starDec2 = float(gaiaBStar['dec'])
        starCoord1 = SkyCoord(starRa1, starDec1, unit="deg")
        starCoord2 = SkyCoord(starRa2, starDec2, unit="deg")
        starParallax1 = float(gaiaAStar['parallax'])
        starParallaxError1 = float(gaiaAStar['parallax_error'])
                    
        # Calculate the widest possible separation for StarA
        possSep1 = possible_distance / dscalculation.calcDistanceMax(starParallax1, starParallaxError1)
        # rhoStar = rhoCalc(starRa1, starDec1, starRa2, starDec2)
        rhoStar = starCoord1.separation(starCoord2).arcsecond
        if possSep1 > rhoStar:
            #starId1 = gaiaAStar['solution_id']
            starName1 = gaiaAStar['designation']
            #starId2 = gaiaBStar['solution_id']
            starName2 = gaiaBStar['designation']
            starParallax2 = float(gaiaBStar['parallax'])
            starParallaxError2 = float(gaiaBStar['parallax_error'])
            starActualRa1 = float(ds['ra_deg_1'].mean())
            starActualDec1 = float(ds['dec_deg_1'].mean())
            starActualRa2 = float(ds['ra_deg_2'].mean())
            starActualDec2 = float(ds['dec_deg_2'].mean())
            starActualCoord1 = SkyCoord(starActualRa1, starActualDec1, unit="deg")
            starActualCoord2 = SkyCoord(starActualRa2, starActualDec2, unit="deg")

            rhoActual = starActualCoord1.separation(starActualCoord2).arcsecond
            starDistance1 = dscalculation.calcDistance(starParallax1)
            starDistanceMax1 = dscalculation.calcDistanceMax(starParallax1, starParallaxError1)
            starDistanceMin1 = dscalculation.calcDistanceMin(starParallax1, starParallaxError1)
            starDistanceRange1 = starDistanceMax1 - starDistanceMin1
            starDistance2 = dscalculation.calcDistance(starParallax2)
            starDistanceMax2 = dscalculation.calcDistanceMax(starParallax2, starParallaxError2)
            starDistanceMin2 = dscalculation.calcDistanceMin(starParallax2, starParallaxError2)
            starDistanceRange2 = starDistanceMax2 - starDistanceMin2

            # Check if stars shares a common distance range

            distanceCommon = ()
            if starDistanceMin1 < starDistanceMin2 < starDistanceMax1 or starDistanceMin2 < starDistanceMin1 < starDistanceMax2:
                distanceCommon = 'overlapping'
            else:
                distanceCommon = 'no'
        
        # Calculate attributes
        pairParallaxFactor, pairPmFactor, pairPmFactor, pairPmCommon, pairAbsMag1, pairAbsMag2, pairLum1, pairLum2, pairRad1, pairRad2, pairDR3Theta, pairDR3Rho, pairMass1, pairMass2, pairBVIndexA, pairBVIndexB, pairSepPar, pairEscapeVelocity, pairRelativeVelocity, pairHarshawFactor, pairHarshawPhysicality, pairBinarity = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 

        pairParallaxFactor = (dscalculation.calcParallaxFactor(gaiaAStar['parallax'], gaiaBStar['parallax'])) * 100
        pairPmFactor = (dscalculation.calcPmFactor(gaiaAStar['pmra'], gaiaAStar['pmdec'], gaiaBStar['pmra'], gaiaBStar['pmdec']))
        pairPmCommon = dscalculation.calcPmCategory(pairPmFactor * 100)
        pairAbsMag1 = dscalculation.calcAbsMag(float(gaiaAStar['phot_g_mean_mag']), float(gaiaAStar['parallax'])) # Calculate Absolute magnitude
        pairAbsMag2 = dscalculation.calcAbsMag(float(gaiaBStar['phot_g_mean_mag']), float(gaiaBStar['parallax'])) # Calculate Absolute magnitude
        pairLum1 = dscalculation.calcLuminosity(pairAbsMag1)
        pairLum2 = dscalculation.calcLuminosity(pairAbsMag2)
        pairAltLum1 = dscalculation.calcLuminosityAlternate(pairAbsMag1)
        pairAltLum2 = dscalculation.calcLuminosityAlternate(pairAbsMag2)
        pairRad1 = dscalculation.calcRadius(pairLum1, gaiaAStar['teff_gspphot'])
        pairRad2 = dscalculation.calcRadius(pairLum2, gaiaBStar['teff_gspphot'])
        # pairDR3Theta = thetaCalc(dscalculation.deltaRa(gaiaAStar['ra'], gaiaBStar['ra'], gaiaBStar['dec']), dscalculation.deltaDec(gaiaBStar['dec'], gaiaAStar['dec'])) + addThetaValue
        pairDR3Theta = starCoord1.position_angle(starCoord2).degree
        # pairDR3Rho = rhoCalc(gaiaAStar['ra'], gaiaAStar['dec'], gaiaBStar['ra'], gaiaBStar['dec'])
        pairDR3Rho = rhoStar
        pairMass1 = dscalculation.calcMass(pairLum1)
        pairMass2 = dscalculation.calcMass(pairLum2)
        pairBVIndexA = float(gaiaAStar['phot_bp_mean_mag']) - float(gaiaAStar['phot_rp_mean_mag'])
        pairBVIndexB = float(gaiaBStar['phot_bp_mean_mag']) - float(gaiaBStar['phot_rp_mean_mag'])
        pairSepPar2 = dscalculation.sepCalc(pairDistanceMinA, pairDistanceMinB, rhoStar) # Separation of the pairs in parsecs
        pairDistance = dscalculation.calc_average_distance(float(gaiaAStar['parallax']), float(gaiaAStar['parallax_error']), float(gaiaBStar['parallax']), float(gaiaBStar['parallax_error']), pairDR3Rho)
        pairSepPar = pairDistance[2] * auToParsec
        print('pairSepPar: ', pairSepPar)
        print('pairSepPar2: ', pairSepPar2)
        print('pairDistance: ', pairDistance[1])

        pairEscapeVelocity = dscalculation.calcEscapevelocity(pairMass1, pairMass2, pairSepPar, gravConst)
        pairRelativeVelocity = dscalculation.calcRelativeVelocity(float(gaiaAStar['pmra']), float(gaiaAStar['pmdec']), float(gaiaBStar['pmra']), float(gaiaBStar['pmdec']), gaiaAStar['radial_velocity'], gaiaBStar['radial_velocity'], pairDistanceMinA, pairDistanceMinB)
        pairHarshawFactor = dscalculation.calcHarshaw((pairParallaxFactor) / 100, (pairPmFactor))
        pairHarshawPhysicality = dscalculation.calcHarshawPhysicality(pairHarshawFactor * 100)
        pairBinarity = dscalculation.calcBinarity(pairRelativeVelocity, pairEscapeVelocity)
        
        # Calculate values for each pair based on the groups
        pairDesignationA = gaiaAStar['designation']
        pairDesignationB = gaiaBStar['designation']
        pairRaA = gaiaAStar['ra']
        pairDecA = gaiaAStar['dec']
        pairRaB = gaiaBStar['ra']
        pairDecB = gaiaBStar['dec']
        pairMagA = gaiaAStar['phot_g_mean_mag']
        pairMagB = gaiaBStar['phot_g_mean_mag']
        pairMeanTheta = float(ds['theta_measured'].degree.mean())
        pairMeanThetaErr = float(ds['theta_measured'].degree.std())
        pairMeanRho = float(ds['rho_measured'].arcsec.mean())
        pairMeanRhoErr = float(ds['rho_measured'].arcsec.std())
        pairMagnitudeA = ds['mag_1']
        pairMagnitudeB = ds['mag_2']
        pairGMagDiff = float(gaiaBStar['phot_g_mean_mag']) - float(gaiaAStar['phot_g_mean_mag'])
        pairMagDiff = float((ds['mag_diff']).mean())
        pairMagDiffErr = (ds['mag_diff']).std()
        pairRadVelA = dscalculation.convertStringToNan(gaiaAStar['radial_velocity'])
        pairRadVelErrA = dscalculation.convertStringToNan(gaiaAStar['radial_velocity_error'])
        pairRadVelB = dscalculation.convertStringToNan(gaiaBStar['radial_velocity'])
        pairRadVelErrB = dscalculation.convertStringToNan(gaiaBStar['radial_velocity_error'])
        pairRadVelRatioA = math.fabs(float(dscalculation.convertStringToNan(gaiaAStar['radial_velocity_error'])) / float(dscalculation.convertStringToNan(gaiaAStar['radial_velocity']))) * 100
        pairRadVelRatioB = math.fabs(float(dscalculation.convertStringToNan(gaiaBStar['radial_velocity_error'])) / float(dscalculation.convertStringToNan(gaiaBStar['radial_velocity']))) * 100
        pairDesA = str(gaiaAStar['designation'])
        pairDesB = str(gaiaBStar['designation'])
        
        # dateOfObservation = getUTC(fitsFileDate)
        print('dateOfObservation: ', Time(ds['image_date'].data))
        print('dateOfObservationMean: ', Time(ds['image_date'].data).mean())
        #dateOfObservation = getUTC(Time(ds['image_date']).mean())
        dateOfObservation = dscalculation.getUTC(Time(ds['image_date'].data).mean())
        
        pairACurrentCoord = dscalculation.calcCurrentDR3Coord(dateOfObservation, pairRaA, pairDecA, float(gaiaAStar['pmra']), float(gaiaAStar['pmdec']))
        pairBCurrentCoord = dscalculation.calcCurrentDR3Coord(dateOfObservation, pairRaB, pairDecB, float(gaiaBStar['pmra']), float(gaiaBStar['pmdec']))
        pairAMeasuredCoord = SkyCoord(ra=ds['ra_deg_1'].groups.aggregate(np.mean) * u.deg, dec=ds['dec_deg_1'].groups.aggregate(np.mean) * u.deg)
        pairBMeasuredCoord = SkyCoord(ra=ds['ra_deg_2'].groups.aggregate(np.mean) * u.deg, dec=ds['dec_deg_2'].groups.aggregate(np.mean) * u.deg)
        pairACoordErr = pairACurrentCoord.separation(pairAMeasuredCoord)
        pairBCoordErr = pairBCurrentCoord.separation(pairBMeasuredCoord)
        # Caculate the common distance from Earth
        
        if (pairMass1 is not None and pairMass2 and pairSepPar is not None and pairDistance[1] is not None and pairACurrentCoord.separation(pairBCurrentCoord).arcsecond is not None and ds[0]['Sep_f'] is not None and pairMeanRho is not None and ds[0]['PA_f'] is not None and pairMeanTheta is not None and ds[0]['Date (first)'] is not None and dateOfObservation):
            # OLD function to calculate the historical orbit values based on the first measurements found in WDS
            # pair_orbit = calc_historic_orbit(pairMass1, pairMass2, pairSepPar, pairDistance[1], pairACurrentCoord.separation(pairBCurrentCoord).arcsecond, ds[0]['Sep_f'], pairMeanRho, ds[0]['PA_f'], pairMeanTheta, ds[0]['Date (first)'], dateOfObservation)
            
            # NEW function to calculate the historical orbit values based on the calculated PA and SEP from Gaia DR3 on epoch 2016
            pair_orbit = dscalculation.calc_historic_orbit(pairMass1, pairMass2, pairSepPar, pairDistance[1], pairACurrentCoord.separation(pairBCurrentCoord).arcsecond, SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).arcsecond, pairMeanRho, SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).degree, pairMeanTheta, gaia_dr3_epoch, dateOfObservation)
        else:
            pair_orbit = ['Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.', 'Cannot be determined, missing data.']
        
        preciseCoord = str(dscalculation.getPreciseCoord(pairRaA, pairDecA, fitsFileDate))
        reportName = (workingDirectory + '/' + pairObjectId + '.txt').replace(' ', '')
        reportFile = open(reportName, "a")
        gaiaData = str(ds[0]['2000 Coord']) + ',' + str(ds[0]['Discov']) + ',' + str(gaiaAStar['pmra']) + ',' + str(gaiaAStar['pmdec']) + ',' + str(gaiaBStar['pmra']) + ',' + str(gaiaBStar['pmdec']) + ',' + str(gaiaAStar['parallax']) + ',' + str(gaiaBStar['parallax']) + ',' + str(dscalculation.calcDistance(gaiaAStar['parallax'])) + ',' + str(dscalculation.calcDistance(gaiaBStar['parallax'])) + ',' + str(gaiaAStar['radial_velocity']) + ',' + str(gaiaBStar['radial_velocity']) + ',' + 'pairRad1' + ',' + 'pairRad2' + ',' + str(pairLum1) + ',' + str(pairLum2) + ',' + str(gaiaAStar['teff_gspphot']) + ',' + str(gaiaBStar['teff_gspphot']) + ',' + str(gaiaAStar['phot_g_mean_mag']) + ',' + str(gaiaBStar['phot_g_mean_mag']) + ',' + str(gaiaAStar['phot_bp_mean_mag']) + ',' + str(gaiaBStar['phot_bp_mean_mag']) + ',' + str(gaiaAStar['phot_rp_mean_mag']) + ',' + str(gaiaBStar['phot_rp_mean_mag']) + ',' + str(pairDR3Theta) + ',' + str(pairDR3Rho) + ',' + str(gaiaAStar['ra']) + ',' + str(gaiaAStar['dec']) + ',' + str(gaiaBStar['ra']) + ',' + str(gaiaBStar['dec']) + ',' + str(gaiaAStar['parallax_error']) + ',' + str(gaiaBStar['parallax_error'])
        
        dscalculation.hrdPlot(pairObjectId, workingDirectory, pairAbsMag1, pairAbsMag2, pairBVIndexA, pairBVIndexB)
        
        print('firstFitsImageFileName =', ds['file'][0])
        firstFitsImageFileName = ds['file'][0]
        print(str(firstFitsImageFileName), str(pairObjectId), str(pairRaA), str(pairDecA), str(pairRaB), str(pairDecB))
        dscalculation.imagePlot(firstFitsImageFileName, workingDirectory, pairObjectId, pairRaA, pairDecA, pairRaB, pairDecB)
    
        # Print temp data
        print('\n### COMPONENTS ###')
        print('\nWDS Identifier:', ds[0]['2000 Coord'], ds[0]['Discov'], ds[0]['Comp'])
        print('\nDate of observation: ' + dateOfObservation)
        print('\nPrecise coordinates (J2000): ' + preciseCoord)
        print('\nComponent A:', pairDesignationA)
        print('Component B:', pairDesignationB)
        print('\nCalculated coordinates')
        print('\nComponent A DR3 2016:', pairRaA, pairDecA)
        print('Component A DR3 on date:', pairACurrentCoord.ra.degree, pairACurrentCoord.dec.degree)
        print('Component A measured:', pairAMeasuredCoord.ra.degree, pairAMeasuredCoord.dec.degree)
        print('Component A error (on date - measured):', pairACoordErr.arcsecond)
        print('\nComponent B DR3 2016:', pairRaB, pairDecB)
        print('Component B DR3 on date:', pairBCurrentCoord.ra.degree, pairBCurrentCoord.dec.degree)
        print('Component B measured:', pairBMeasuredCoord.ra.degree, pairBMeasuredCoord.dec.degree)
        print('Component B error (on date - measured):', pairBCoordErr.arcsecond)
        print('2016 Calculated Position angle / Separation: ', SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).degree, SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).arcsecond)
        #print('Current Calculated Position angle / Separation: ', SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).degree, SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).arcsecond)
        print('Current Calculated Position angle / Separation: ', pairACurrentCoord.position_angle(pairBCurrentCoord).degree, pairACurrentCoord.separation(pairBCurrentCoord).arcsecond)
        print('\nTheta measurements\n') # , ds['dspaactual']
        print('Mean:', pairMeanTheta)
        print('Error:', pairMeanThetaErr)
        print('\nRho measurements\n') # , ds['dssepactual']
        print('Mean:', pairMeanRho)
        print('Error:', pairMeanRhoErr)
        print('\nMagnitude measurements\n') # , ds['dsmagdiff']
        print('Mean:', pairMagDiff)
        print('Error:', pairMagDiffErr)
        print('\n\nParallax factor:', pairParallaxFactor, '%')
        print('Proper motion factor:', pairPmFactor * 100, '%')
        print('Proper motion category:', pairPmCommon)
        print('Absolute magnitude A:', pairAbsMag1)
        print('Absolute magnitude B:', pairAbsMag2)
        print('Luminosity A:', pairLum1)
        print('Luminosity B:', pairLum2)
        print('Luminosity Alternate A:', pairAltLum1)
        print('Luminosity Alternate B:', pairAltLum2)
        print('Mass A:', pairMass1)
        print('Mass B:', pairMass2)
        print('BV index A:', pairBVIndexA, 'B:', pairBVIndexB)
        print('Radial velocity of the stars', 'A:', pairRadVelA, 'km/s (Err:', pairRadVelErrA, 'km/s)', 'B:', pairRadVelB, 'km/s (Err:', pairRadVelErrB, 'km/s)')
        print('Radial velocity ratio A:', pairRadVelRatioA, '%')
        print('Radial velocity ratio B:', pairRadVelRatioB, '%')
        print('Separation:', pairSepPar, 'parsec,', pairSepPar * 206265, 'AU')
        print('Pair Escape velocity:', pairEscapeVelocity, 'km/s')
        print('Pair Relative velocity:', pairRelativeVelocity, 'km/s')
        print('### Pair historical orbit calculations ###')
        print('Historic criterion: ', pair_orbit[0])
        print('Max orbit velolicy: ', pair_orbit[1])
        print('Observed velocity: ', pair_orbit[2])
        print('Pair Harshaw factor:', pairHarshawFactor)
        print('Pair Harshaw physicality:', pairHarshawPhysicality)
        print('Pair binarity:', pairBinarity)
        print('Analysis finished: ' , datetime.datetime.now())
        
        # Write results to file
        reportTable.add_row([ds[0]['2000 Coord'] + ds[0]['Discov'] + str(ds[0]['Comp']), dateOfObservation, pairMeanTheta, pairMeanThetaErr, pairMeanRho, pairMeanRhoErr, np.nan, np.nan, pairMagDiff, pairMagDiffErr, 'Filter wawelenght', 'filter FWHM', '0.2', '1', 'TLB_2024', 'C', '7', preciseCoord])
        reportFile.write('### WDS Data ###')
        reportFile.write('\nWDS Identifier: ' + ds[0]['2000 Coord'])
        reportFile.write('\nDiscoverer and components: ' + str(ds[0]['Discov']) + ' ' + str(ds[0]['Comp']))
        reportFile.write('\nMagnitude Pri: ' + str(ds[0]['Mag_A']))
        reportFile.write('\nMagnitude Sec: ' + str(ds[0]['Mag_B']))
        reportFile.write('\nPA last: ' + str(ds[0]['PA_l']))
        reportFile.write('\nSep last: ' +  str(ds[0]['Sep_l']))
        reportFile.write('\n\n### Gaia DR3 Data ###')
        reportFile.write('\nMain star: ' + pairDesA)
        reportFile.write('\nCompanion: ' + pairDesB)
        reportFile.write('\nPair G magnitudes A: ' + str(dscalculation.roundNumber(pairMagA)) + ' B: ' + str(dscalculation.roundNumber(pairMagB)))
        reportFile.write('\nPosition angle: ' + str(dscalculation.roundNumber(pairDR3Theta)))
        reportFile.write('\nSeparation: ' + str(dscalculation.roundNumber(pairDR3Rho)))
        reportFile.write('\nMagnitude difference: ' + str(dscalculation.roundNumber(pairGMagDiff)))
        reportFile.write('\nPrecise coordinates (J2000): ' + preciseCoord)
        reportFile.write('\nDate of observation: ' + dateOfObservation)
        reportFile.write('\n\nCalculated coordinates')
        reportFile.write('\nComponent A DR3 2016: ' + str(pairRaA) + ' ' + str(pairDecA))
        reportFile.write('\nComponent A DR3 on date: ' + str(pairACurrentCoord.ra.degree) + ' ' + str(pairACurrentCoord.dec.degree))
        reportFile.write('\nComponent A measured: ' + str(pairAMeasuredCoord.ra.degree) + ' ' + str(pairAMeasuredCoord.dec.degree))
        reportFile.write('\nComponent A error (on date - measured): ' + str(pairACoordErr.arcsecond))
        reportFile.write('\nComponent B DR3 2016: ' + str(pairRaB) + ' ' + str(pairDecB))
        reportFile.write('\nComponent B DR3 on date: ' + str(pairBCurrentCoord.ra.degree) + ' ' + str(pairBCurrentCoord.dec.degree))
        reportFile.write('\nComponent B measured: ' + str(pairBMeasuredCoord.ra.degree) + ' ' + str(pairBMeasuredCoord.dec.degree))
        reportFile.write('\nComponent B error (on date - measured): ' + str(pairBCoordErr.arcsecond))
        reportFile.write('\n\n2016 Calculated Position angle / Separation: '  + str(SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').position_angle(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).degree) + ' ' + str(SkyCoord(ra=pairRaA*u.degree, dec=pairDecA*u.degree, frame='icrs').separation(SkyCoord(ra=pairRaB*u.degree, dec=pairDecB*u.degree, frame='icrs')).arcsecond))
        reportFile.write('\nCurrent Calculated Position angle / Separation: ' + str(pairACurrentCoord.position_angle(pairBCurrentCoord).degree) + ' ' + str(pairACurrentCoord.separation(pairBCurrentCoord).arcsecond))
        reportFile.write('\n\n### Measurements ###')
        reportFile.write('\nPosition angle:')
        reportFile.write('\nTheta measurements' + str(ds['theta_measured'].degree))
        reportFile.write('\nMean: ' + str(dscalculation.roundNumber(pairMeanTheta)))
        reportFile.write('\nError: ' + str(dscalculation.roundNumber(pairMeanThetaErr)))
        reportFile.write('\n\nSeparation:')
        reportFile.write('\nRho measurements\n' + str(ds['rho_measured'].arcsec))
        reportFile.write('\nMean: ' + str(dscalculation.roundNumber(pairMeanRho)))
        reportFile.write('\nError: ' + str(dscalculation.roundNumber(pairMeanRhoErr)))
        reportFile.write('\n\nMagnitude measurements\n'  + str(ds['mag_diff']))
        reportFile.write('\nMean: ' + str(dscalculation.roundNumber(pairMagDiff)))
        reportFile.write('\nError: ' + str(dscalculation.roundNumber(pairMagDiffErr)))
        reportFile.write('\n\n### Calculated attributes ###')
        reportFile.write('\nSeparation (Measured): ' + str(dscalculation.roundNumber(pairMeanRho)))
        reportFile.write('\nPosition angle (Measured): ' + str(dscalculation.roundNumber(pairMeanTheta)))
        reportFile.write('\nMagnitude difference (Measured): ' + str(dscalculation.roundNumber(pairMagDiff)) + ' (Err: ' + str(dscalculation.roundNumber(pairMagDiffErr)) + ')')
        reportFile.write('\nAbsolute magnitude A: ' + str(dscalculation.roundNumber(pairAbsMag1)))
        reportFile.write('\nAbsolute magnitude B: ' + str(dscalculation.roundNumber(pairAbsMag2)))
        reportFile.write('\nLuminosity A: ' + str(dscalculation.roundNumber(pairLum1)))
        reportFile.write('\nLuminosity B: ' + str(dscalculation.roundNumber(pairLum2)))
        reportFile.write('\nRad A: ' + str(dscalculation.roundNumber(pairRad1)))
        reportFile.write('\nRad B: ' + str(dscalculation.roundNumber(pairRad2)))
        reportFile.write('\nMass A: ' + str(dscalculation.roundNumber(pairMass1)))
        reportFile.write('\nMass B: ' + str(dscalculation.roundNumber(pairMass2)))
        reportFile.write('\nBV index A: ' + str(dscalculation.roundNumber(pairBVIndexA)) + ' B: ' + str(dscalculation.roundNumber(pairBVIndexB)))
        reportFile.write('\nRadial velocity of the stars ' + 'A:' + str(dscalculation.roundNumber(pairRadVelA)) + 'km/s (Err:' + str(dscalculation.roundNumber(pairRadVelErrA)) + 'km/s)' + ' B:' + str(dscalculation.roundNumber(pairRadVelB)) + 'km/s (Err:' + str(dscalculation.roundNumber(pairRadVelErrB)) + 'km/s)')
        reportFile.write('\nRadial velocity ratio A: ' + str(dscalculation.roundNumber(pairRadVelRatioA)) + ' %')
        reportFile.write('\nRadial velocity ratio B: ' + str(dscalculation.roundNumber(pairRadVelRatioB)) + ' %')
        reportFile.write('\nSeparation: ' + str(dscalculation.roundNumber(pairDistance[2] * auToParsec)) + ' parsec, ' + str(dscalculation.roundNumber((pairDistance[2]))) + ' AU')
        reportFile.write('\nPair Escape velocity: ' + str(dscalculation.roundNumber(pairEscapeVelocity)) + ' km/s')
        reportFile.write('\nPair Relative velocity: ' + str(dscalculation.roundNumber(pairRelativeVelocity)) + ' km/s')
        reportFile.write('\n\n### Analysis ###')
        reportFile.write('\nParallax factor: ' + str(dscalculation.roundNumber(pairParallaxFactor)) + ' %')
        reportFile.write('\nProper motion factor: ' + str(dscalculation.roundNumber(pairPmFactor) * 100) + ' %')
        reportFile.write('\nProper motion category: '+ str(pairPmCommon))
        reportFile.write('\nPair Harshaw factor: ' + str(dscalculation.roundNumber(pairHarshawFactor)))
        reportFile.write('\nPair Harshaw physicality: ' + str(pairHarshawPhysicality))
        reportFile.write('\nPair binarity: ' + str(pairBinarity))
        
        # new function - orbit calculation
        reportFile.write('\n\n### Pair historical orbit calculations ###')
        reportFile.write('\nHistoric criterion: ' + str(pair_orbit[0]))
        reportFile.write('\nDelta theta: ' + str(pair_orbit[4]))
        reportFile.write('\nDelta rho: ' + str(pair_orbit[5]))
        reportFile.write('\nDelta time: ' + str(pair_orbit[6]))
        reportFile.write('\nMax orbit velolicy: ' + str(pair_orbit[1]))
        reportFile.write('\nObserved velocity: ' + str(pair_orbit[2]))
        reportFile.write('\nInput data variables: ' + str(pair_orbit[3]))
        
        reportFile.write('\n\n### WDS form:\n')
        wdsform = str(ds[0]['2000 Coord']) + ',' + dateOfObservation + ',' +  str(dscalculation.roundNumber(pairMeanTheta)) + ',' +  str(dscalculation.roundNumber(pairMeanThetaErr)) + ',' +  str(dscalculation.roundNumber(pairMeanRho)) + ',' +  str(dscalculation.roundNumber(pairMeanRhoErr)) + ',' +  'nan' + ',' +  'nan' + ',' +  str(dscalculation.roundNumber(pairMagDiff)) + ',' +  str(dscalculation.roundNumber(pairMagDiffErr)) + ',' + 'Filter wawelenght' + ',' + 'filter FWHM' + ',' + '0.2' + ',' + '1' + ',' + 'TLB_2023' + ',' +  'C' + ',' + '7'+ ',' + str(dscalculation.getPreciseCoord(pairRaA, pairDecA, fitsFileDate))
        reportFile.write(str(wdsform))
        reportFile.write('\n\n### Gaia data:\n')
        reportFile.write(str(gaiaData))

reportTable.write(workingDirectory + '/double_stars_wds_format.txt', format='ascii', overwrite=True, delimiter=',')